#!/usr/bin/env python3
"""Numerical regression tests for the MaVaN + PNGB CLASS fork.

The tests deliberately exercise both the exact momentum hierarchy and
CLASS's optional late-time ncdm fluid closure.  They also verify the
beta -> 0 limit against the uncoupled scalar-field code path.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "mavan_validation"


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def load(name: str) -> np.ndarray:
    data = np.loadtxt(OUT / name)
    if data.ndim != 2 or not np.all(np.isfinite(data)):
        raise AssertionError(f"{name} is not a finite two-dimensional table")
    return data


def log_interpolated_max_relative(path_a: str, path_b: str) -> float:
    a = load(path_a)
    b = load(path_b)
    if np.any(a[:, 0] <= 0) or np.any(b[:, 0] <= 0):
        raise AssertionError("P(k) comparison requires positive k")
    if np.any(a[:, 1] <= 0) or np.any(b[:, 1] <= 0):
        raise AssertionError("P(k) comparison requires positive power")

    lower = max(a[:, 0].min(), b[:, 0].min())
    upper = min(a[:, 0].max(), b[:, 0].max())
    k = np.geomspace(lower, upper, 1000)
    power_a = np.exp(np.interp(np.log(k), np.log(a[:, 0]), np.log(a[:, 1])))
    power_b = np.exp(np.interp(np.log(k), np.log(b[:, 0]), np.log(b[:, 1])))
    return float(np.max(np.abs(power_a / power_b - 1.0)))


def check_background(name: str, beta: float, f_pngb: float) -> dict[str, float]:
    # Column layout is fixed by background_output_titles() for a one-species
    # MaVaN run.  Indices below are zero-based.
    table = load(name)
    if table.shape[1] != 32:
        raise AssertionError(f"{name} has {table.shape[1]} columns instead of 32")

    z = table[:, 0]
    loga = -np.log1p(z)
    conformal_hubble = table[:, 3] / (1.0 + z)

    rho_nu = table[:, 11]
    p_nu = table[:, 12]
    trace_nu = table[:, 13]
    mass_eV = table[:, 14]
    rho_phi = table[:, 17]
    p_phi = table[:, 18]
    phi = table[:, 20]
    phi_prime = table[:, 21]
    potential = table[:, 22]
    potential_prime = table[:, 23]
    potential_second = table[:, 24]
    omega_r = table[:, 28]
    omega_m = table[:, 29]

    expected_mass = mass_eV[-1] * np.exp(beta * (phi - phi[-1]))
    mass_law_error = float(np.max(np.abs(mass_eV / expected_mass - 1.0)))
    target_error_eV = float(abs(mass_eV[-1] - 0.06))

    v0_samples = potential / (1.0 + np.cos(phi / f_pngb))
    v0 = float(np.median(v0_samples))
    potential_error = float((v0_samples.max() - v0_samples.min()) / v0)
    potential_prime_error = float(
        np.max(np.abs(potential_prime + v0 / f_pngb * np.sin(phi / f_pngb)))
        / (v0 / f_pngb)
    )
    potential_second_error = float(
        np.max(np.abs(potential_second + v0 / f_pngb**2 * np.cos(phi / f_pngb)))
        / (v0 / f_pngb**2)
    )

    dphi_dloga = phi_prime / conformal_hubble
    drho_nu_dloga = np.gradient(rho_nu, loga, edge_order=2)
    drho_phi_dloga = np.gradient(rho_phi, loga, edge_order=2)
    source = beta * dphi_dloga * trace_nu
    residual_nu = drho_nu_dloga + 3.0 * (rho_nu + p_nu) - source
    residual_phi = drho_phi_dloga + 3.0 * (rho_phi + p_phi) + source
    residual_sum = (
        np.gradient(rho_nu + rho_phi, loga, edge_order=2)
        + 3.0 * (rho_nu + p_nu + rho_phi + p_phi)
    )

    # Very early on, rho_phi is almost exactly constant and a relative
    # finite-difference residual is dominated by printed-table precision.
    # The late mask covers the entire dynamically relevant scalar evolution.
    late = (z < 1.0e5) & (z > 1.0e-4)

    def p99_relative(residual: np.ndarray, scale: np.ndarray) -> float:
        ratio = np.abs(residual[late]) / np.maximum(np.abs(scale[late]), 1.0e-300)
        return float(np.quantile(ratio, 0.99))

    nu_conservation = p99_relative(residual_nu, rho_nu)
    phi_conservation = p99_relative(residual_phi, rho_phi)
    total_conservation = p99_relative(residual_sum, rho_nu + rho_phi)

    # At high redshift rho-3p is far below rho and direct subtraction is
    # intentionally inaccurate.  Compare only where subtraction is resolved.
    resolved_trace = trace_nu > 1.0e-8 * rho_nu
    trace_error = float(
        np.max(np.abs((rho_nu[resolved_trace] - 3.0 * p_nu[resolved_trace])
                      / trace_nu[resolved_trace] - 1.0))
    )

    checks = {
        "mass target [eV]": target_error_eV,
        "exponential mass law": mass_law_error,
        "PNGB potential": potential_error,
        "PNGB first derivative": potential_prime_error,
        "PNGB second derivative": potential_second_error,
        "neutrino continuity p99": nu_conservation,
        "scalar continuity p99": phi_conservation,
        "coupled-sum continuity p99": total_conservation,
        "stable trace cross-check": trace_error,
    }

    limits = {
        "mass target [eV]": 1.0e-6,
        "exponential mass law": 1.0e-10,
        "PNGB potential": 1.0e-9,
        "PNGB first derivative": 1.0e-9,
        "PNGB second derivative": 1.0e-9,
        "neutrino continuity p99": 2.0e-4,
        "scalar continuity p99": 2.0e-5,
        "coupled-sum continuity p99": 2.0e-4,
        "stable trace cross-check": 1.0e-3,
    }
    for label, value in checks.items():
        if value > limits[label]:
            raise AssertionError(
                f"{name}: {label}={value:.6e} exceeds {limits[label]:.6e}"
            )

    if not (0.0 < omega_r[-1] < 1.0e-3):
        raise AssertionError(f"{name}: unphysical Omega_r(0)={omega_r[-1]}")
    if not (0.2 < omega_m[-1] < 0.4):
        raise AssertionError(f"{name}: unphysical Omega_m(0)={omega_m[-1]}")

    print(f"\n{name}")
    for label, value in checks.items():
        print(f"  {label:31s} {value:.6e}")
    print(f"  {'Omega_r(0)':31s} {omega_r[-1]:.6e}")
    print(f"  {'Omega_m(0)':31s} {omega_m[-1]:.6e}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-build", action="store_true", help="reuse the existing binary")
    parser.add_argument("--no-run", action="store_true", help="only inspect existing outputs")
    args = parser.parse_args()

    if not args.no_build:
        run(["make", "class", "-j2"])

    if not args.no_run:
        # CLASS protects existing outputs by altering the requested root.  A
        # clean, validation-only directory guarantees that every table read
        # below was produced by this invocation (and is not a stale partial
        # file left by an interrupted earlier run).
        if OUT.exists():
            shutil.rmtree(OUT)
        OUT.mkdir(parents=True, exist_ok=True)
        for ini in (
            "mavan_pngb_exact.ini",
            "mavan_pngb_fluid.ini",
            "mavan_pngb_beta0.ini",
            "pngb_uncoupled.ini",
            "mavan_pngb_strong.ini",
            "mavan_pngb_strong_fluid.ini",
            "mavan_pngb_prior_edge.ini",
            "mavan_pngb_prior_edge_fluid.ini",
            "mavan_pngb_negative.ini",
            "mavan_pngb_cmb.ini",
            "mavan_pngb_theta.ini",
        ):
            run(["./class", f"examples/{ini}"])
    else:
        OUT.mkdir(parents=True, exist_ok=True)

    check_background("exact_background.dat", beta=1.0, f_pngb=1.0)
    check_background("strong_background.dat", beta=4.0, f_pngb=0.6)
    check_background("prior_edge_background.dat", beta=6.0, f_pngb=0.3)
    check_background("negative_background.dat", beta=-1.0, f_pngb=1.0)

    fluid_difference = log_interpolated_max_relative("exact_pk.dat", "fluid_pk.dat")
    strong_fluid_difference = log_interpolated_max_relative(
        "strong_pk.dat", "strong_fluid_pk.dat"
    )
    edge_fluid_difference = log_interpolated_max_relative(
        "prior_edge_pk.dat", "prior_edge_fluid_pk.dat"
    )
    beta0_difference = log_interpolated_max_relative("beta0_pk.dat", "uncoupled_pk.dat")
    cmb_cl = load("cmb_cl.dat")
    cmb_lensed_cl = load("cmb_cl_lensed.dat")
    theta_pk = load("theta_pk.dat")
    if np.any(cmb_cl[:, 0] < 2) or np.any(cmb_lensed_cl[:, 0] < 2):
        raise AssertionError("CMB output contains invalid multipoles")
    print(f"\n  exact vs fluid max relative P(k)  {fluid_difference:.6e}")
    print(f"  strong exact vs fluid relative P(k) {strong_fluid_difference:.6e}")
    print(f"  prior-edge exact vs fluid P(k)    {edge_fluid_difference:.6e}")
    print(f"  beta=0 vs uncoupled relative P(k) {beta0_difference:.6e}")
    print(f"  finite unlensed/lensed CMB rows    {len(cmb_cl)}/{len(cmb_lensed_cl)}")
    print(f"  acoustic-scale shooting P(k) rows  {len(theta_pk)}")

    if fluid_difference > 2.0e-3:
        raise AssertionError("late-time ncdm fluid closure differs from exact P(k) by >0.2%")
    if strong_fluid_difference > 2.0e-3:
        raise AssertionError("strong-coupling fluid closure differs from exact P(k) by >0.2%")
    if edge_fluid_difference > 1.0e-2:
        raise AssertionError("prior-edge fluid closure differs from exact P(k) by >1%")
    if beta0_difference > 1.0e-5:
        raise AssertionError("beta=0 code path does not recover the uncoupled P(k)")

    print("\nAll MaVaN + PNGB validation checks passed.")


if __name__ == "__main__":
    main()
