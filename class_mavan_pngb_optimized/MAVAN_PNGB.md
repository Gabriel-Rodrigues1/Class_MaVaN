# MaVaN + PNGB model in CLASS v3.3.4

This fork implements one massive neutrino species coupled to a canonical
scalar field through an exponential, field-dependent mass, together with a
PNGB potential. Both the homogeneous evolution and scalar perturbations are
included. The exact momentum-dependent ncdm hierarchy is the reference
implementation; CLASS's late-time ncdm fluid approximation is also supported.

The implementation follows the corrected Boltzmann hierarchy in the erratum
to Brookfield et al., *Cosmology of Mass-Varying Neutrinos Driven by
Quintessence: Theory and Observations*,
[arXiv:astro-ph/0512367](https://arxiv.org/abs/astro-ph/0512367), and the
background, stress-energy, and fluid equations used by da Fonseca, Barreiro,
and Nunes, *Relaxing cosmological constraints on current neutrino masses*,
[arXiv:2311.01803](https://arxiv.org/abs/2311.01803).

## Model and conventions

The public mass parameter is the physical mass today:

\[
  m_\nu(\phi)=m_{\nu,0}\exp\!\left[\beta(\phi-\phi_0)\right],
  \qquad \beta=\frac{d\ln m_\nu}{d\phi}.
\]

Here \(\phi\) and the PNGB scale \(f\) are measured in reduced-Planck units.
The standard CLASS input `m_ncdm` is \(m_{\nu,0}\), the physical mass today;
it is never varied by the shooting algorithm. CLASS stores the corresponding
dimensionless present mass \(M_{\rm ncdm}=m_{\nu,0}/T_{\nu,0}\) and evaluates
the mass history with an internal reference field \(\phi_0\).

Because \(\phi_0\) is an endpoint value needed during a forward integration,
CLASS solves the consistency equation

\[
 F_{\phi_0}=\phi(a=1;\phi_0)-\phi_0=0.
\]

When `Omega_scf = -1`, this condition and spatial closure are solved together
for \(\phi_0\) and \(\log_{10}V_0\) with one multidimensional Newton solve.
There is no neutrino-mass shooting and no nested root finder inside the
background integrator.

The scalar potential is

\[
  V(\phi)=V_0\left[1+\cos\left(\frac{\phi}{f}\right)\right],
\]

with

\[
 V_{,\phi}=-\frac{V_0}{f}\sin(\phi/f),\qquad
 V_{,\phi\phi}=-\frac{V_0}{f^2}\cos(\phi/f).
\]

`scf_parameters` is interpreted as

```text
f, log10(V0/Mpc^-2), phi_ini, phi_prime_ini
```

where CLASS stores \(V\) in \(M_{\rm Pl}^2/{\rm Mpc}^2\) and
`phi_prime_ini` is a conformal-time derivative. The second entry is normally
only a shooting guess: use `Omega_scf = -1` and `scf_tuning_index = 1`.

For \(\epsilon=(q^2+a^2m_\nu^2)^{1/2}\), the background neutrino moments are
evaluated from the usual phase-space integrals. In particular, the coupling
source is computed directly as

\[
 \rho_\nu-3p_\nu = \frac{C}{a^4}\int dq\,
 q^2\frac{a^2m_\nu^2}{\epsilon}f_0(q).
\]

This direct integral is essential in the relativistic regime. Subtracting
the separately computed \(\rho_\nu\) and \(3p_\nu\) causes catastrophic
cancellation and can make the background integrator fail.

With \({\cal H}=a'/a\), the exchange equations are

\[
 \rho_\nu'+3{\cal H}(\rho_\nu+p_\nu)
 =\beta\phi'(\rho_\nu-3p_\nu),
\]

\[
 \phi''+2{\cal H}\phi'+a^2V_{,\phi}
 =-a^2\beta(\rho_\nu^{\rm phys}-3p_\nu^{\rm phys})/M_{\rm Pl}^2.
\]

CLASS densities are \(\rho_{\rm phys}/(3M_{\rm Pl}^2)\), which is why the
implemented Klein-Gordon source contains an explicit factor of three.
The two exchange terms cancel in the total continuity equation.

## Perturbations

The implementation is currently restricted to synchronous gauge. In the
exact hierarchy, every occurrence of the comoving energy uses the local
effective mass. The corrected spatial mass-gradient term is added to the
dipole; it is not omitted as in the original, pre-erratum version of
Brookfield et al.

The stress-energy perturbations include the explicit response of
\(\epsilon\) to \(\delta\phi\):

\[
\delta\rho_\nu = \frac{C}{a^4}\int dq\,q^2 f_0
\left[\epsilon\Psi_0+
\beta\delta\phi\frac{a^2m_\nu^2}{\epsilon}\right],
\]

\[
\delta p_\nu = \frac{C}{3a^4}\int dq\,q^4\frac{f_0}{\epsilon}
\left[\Psi_0-
\beta\delta\phi\frac{a^2m_\nu^2}{\epsilon^2}\right].
\]

The perturbed Klein-Gordon equation uses the complete, cancellation-free
\(\delta\rho_\nu-3\delta p_\nu\). Since that trace already contains the
mass response proportional to \(\beta\delta\phi\), it automatically includes
the \(\beta^2\) contribution for a constant exponential coupling.

When the ncdm fluid approximation is active, its continuity and Euler
equations contain the coupling, \(\delta\phi'\), fifth-force, and modified
friction terms. The shear equation remains CLASS's closure ansatz. Use the
exact hierarchy for reference calculations:

```text
ncdm_fluid_approximation = 3
```

For production chains, the faster CLASS closure is

```text
ncdm_fluid_approximation = 2
```

The included tests find maximum exact-versus-fluid differences of
\(6.5\times10^{-4}\), \(6.9\times10^{-4}\), and \(7.9\times10^{-4}\) in
linear \(P(k)\) for \((\beta,f)=(1,1),(4,0.6),(6,0.3)\), respectively, over
\(k\lesssim0.2\,h\,{\rm Mpc}^{-1}\). These are regression results, not a
guarantee over every combination in the prior: validate the approximation
again around the final posterior support.

## Performance design

The implementation avoids repeated exponentials and mutable shared state in
the perturbation loops. Each background row stores the local effective mass,
the trace \(\rho_\nu-3p_\nu\), and the pressure-mass derivative moment. The
Boltzmann hierarchy and fluid closure read these cached quantities directly.
This is thread-safe and keeps the extra work outside the inner momentum and
multipole loops.

The scalar closure and \(\phi_0\) endpoint conditions use one safeguarded
multidimensional solve. At the fiducial example this takes 13 inexpensive
trial background evaluations, followed by the final background calculation;
it does not nest one shooting loop inside another.
The trust-region safeguard is enabled only for this MaVaN joint solve, so it
does not alter ordinary CLASS shooting.

For production CMB, lensing, BAO, and supernova chains, start with:

```text
ncdm_fluid_approximation = 2
nonlinear_min_k_max = 5
```

On the validation machine, three warm one-core Python-wrapper runs had median
times of 1.52 s for the supplied CMB+lensing+matter-power calculation with
`nonlinear_min_k_max = 5`, versus 2.10 s when it was forced to 20. Standalone
minimal linear-\(P(k)\) runs had median times of 0.50 s with the fluid closure
and 0.59 s with the exact hierarchy. These numbers are machine-dependent; use
the exact hierarchy around representative posterior points as a validation
calculation, not for the main chain.

## Minimal CLASS input

```ini
N_ncdm = 1
N_ur = 2.0308
deg_ncdm = 1
T_ncdm = 0.71611

mavan_ncdm = yes
beta_ncdm = 1.0
m_ncdm = 0.06

Omega_Lambda = 0
Omega_fld = 0
Omega_scf = -1
attractor_ic_scf = no
scf_tuning_index = 1
scf_parameters = 1.0, -7.0, 0.4, 0.0

gauge = synchronous
ncdm_fluid_approximation = 2
```

`m_ncdm` is the mass of one particle represented by the coupled distribution.
With `deg_ncdm = 1`, one massive eigenstate of 0.06 eV plus two effectively
massless states is represented by `N_ur = 2.0308`; with the default ncdm
temperature this gives \(N_{\rm eff}\simeq3.044\). Do not use `deg_ncdm = 3`
unless the intended model really has three degenerate coupled states, in which
case the total present mass is three times the per-particle mass.

Do not pass `Omega_ncdm`, `omega_ncdm`, or the internal
`mavan_phi_today`. The density is inferred from the present physical mass.

## Build and validate

Build the standalone CLASS executable with

```bash
make clean
make class -j4
```

Use the explicit `class` target. A plain parallel `make` can try to build the
native executable and Python extension concurrently on some systems.

For Cobaya, install the Python wrapper inside the same environment that runs
Cobaya:

```bash
CC=gcc CXX=g++ python -m pip install . --no-build-isolation
```

This environment must already contain NumPy, Cython, setuptools, and wheel.
Setting `CC` and `CXX` explicitly avoids inheriting a nonexistent compiler
from a Python installation built with Clang.

Run all supplied numerical checks with

```bash
python validation/test_mavan_pngb.py
```

The suite checks:

- the exponential mass law and the 0.06 eV boundary condition;
- both PNGB derivatives;
- separate neutrino and scalar continuity equations and conservation of
  their sum;
- the stable trace integral;
- physical present-day `Omega_m` and `Omega_r` bookkeeping;
- the \(\beta=0\) limit against an uncoupled PNGB model;
- exact-hierarchy versus fluid-closure linear matter power at moderate,
  strong, and prior-edge couplings;
- finite CMB, lensing, and matter-power output through the full pipeline;
- the optional three-target `theta_s_100` shooting path.

The MaVaN background output adds `rho_nu-3p_nu` and `m_nu_eff[eV]`, making the
mass history and coupling source directly auditable.

## Cobaya notes

An end-to-end example is in `examples/mavan_pngb_cobaya.yaml`. Sampling
`H0` directly is computationally preferable here. Scalar closure and the
endpoint condition for \(\phi_0\) already form a two-dimensional background
shooting problem; sampling `theta_s_100` adds a third target to every
likelihood call.

For a fixed present mass, keep `m_ncdm: 0.06` under CLASS `extra_args`. To
sample it, move `m_ncdm` to the Cobaya `params` block. The internal reference
field `mavan_phi_today` is generated by CLASS and must never be exposed to the
sampler.

A rectangular prior in `(beta_ncdm, f_pngb)` can include points for which the
chosen `phi_ini` and PNGB potential have no solution satisfying the requested
present scalar closure. At such points CLASS intentionally returns a shooting
failure; with Cobaya `stop_at_error: false` this is an ordinary rejected
proposal. Do not add a slower nested mass/root solve to hide these points.
After a pilot run, trimming demonstrably no-solution regions from the prior can
improve sampling efficiency without changing the viable model.

## Current scope

- one coupled ncdm species (`N_ncdm = 1`), with arbitrary multiplicity;
- constant exponential coupling \(\beta\);
- canonical PNGB scalar field;
- scalar modes in synchronous gauge;
- exact hierarchy plus an optional late-time fluid closure.

Tensor/vector MaVaN sources, a field-dependent \(\beta(\phi)\), multiple
independently coupled massive species, and Newtonian-gauge evolution are not
part of this implementation.
