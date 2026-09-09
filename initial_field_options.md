**Allowing the pNGB initial field to vary**

The existing implementation reads the initial field from the third entry of
`scf_parameters`. Changing this entry does not require a C-source change.
Keep the fourth entry zero, so the field initially has zero conformal-time
velocity. The present-day neutrino mass remains 0.06 eV, and the existing
internal solve for the field reference and potential amplitude is retained.
There is no additional shooting unknown when the initial field is sampled.

For this potential,

\[
\vartheta_{\rm ini}=\phi_{\rm ini}/f,\qquad
\left.\frac{V_{,\phi}}{V}\right|_{\rm ini}
=-\frac{1}{f}\tan\left(\frac{\phi_{\rm ini}}{2f}\right).
\]

Fixing `phi_ini=0.4` imposes `theta_ini=0.4/f`: over the supplied f prior,
the initial angle runs from approximately 0.2 to 1.33 radians. Freeing the
initial position tests this restriction. Near the hilltop, the slope is
approximately `-phi_ini/(2*f**2)`, which explains one source of correlation
between the initial field and f. The neutrino coupling adds further
dependence on the field history. A wider model can fit the data better while
also giving wider reconstruction bands and more difficult MCMC mixing.

**Suggested sequence**

1. Make a short constrained profile in the raw initial field, with a few
   values around 0.4, such as 0.2, 0.4, and 0.6. Extend it if the best fit
   occurs at a scan boundary or the fit changes substantially. Include
   negative-field starts/points when assessing the full signed-angle model
   in the sampling YAML. These values are an exploratory grid, not a
   recommended final prior or a complete coverage of possible solutions.
2. Reoptimize all the other cosmological and likelihood nuisance parameters
   at each value. Check the minimizer's success and agreement between starts.
   Record CLASS shooting failures separately: they do not establish that a
   physical solution is absent.
3. For Bayesian reconstruction bands that include initial-condition
   uncertainty, sample the extra parameter jointly and marginalize. A
   profile optimum does not provide that marginalization. Profiling inside
   every MCMC proposal also changes the statistical target and requires
   repeated model evaluations; it is not a simple speed optimization.

Profiling can diagnose the dependence on initial conditions without running
a full posterior chain at every fixed value. It still requires optimization
and can be costly. Sampling adds a slow cosmological parameter and its
correlations; it is not expected automatically to improve acceptance or
convergence. Use the exploratory results to improve starting points and
proposal geometry, and state the scientific prior independently.

**Sampling variant: `mavan_pngb_sample_initial_angle.yaml`**

This samples `theta_ini`, with a uniform prior on `[-pi, pi]`, and passes
`phi_ini=f_pngb*theta_ini` as the third CLASS input. The physical initial
field is also recorded as a derived chain column. The supplied f, beta,
and other cosmological priors, theory settings, and likelihoods are retained.

The signed interval covers one period of the pNGB potential. If the analysis
deliberately restricts the initial field to the positive branch, set the
minimum to zero and report that restriction. With the supplied asymmetric
beta prior, excluding negative initial fields cannot generally be justified
just by the even symmetry of the potential: reflection of the coupled
equations also reverses beta.

A uniform angle is a transparent prior choice, but is not equivalent to a
uniform field value independent of f. At fixed f it induces

\[
p(\phi_{\rm ini}\mid f)=\frac{1}{2\pi f},\qquad
|\phi_{\rm ini}|\leq\pi f.
\]

The `ref` and `proposal` widths are initial numerical choices, not results
of a tuning study. The reference is near the current solution for f near
one. Also initialize independent chains on the negative-angle side and
check whether they reach the same explored posterior; convergence among
chains initialized in one region cannot exclude a missed mode.

Run with the supplied new output prefix:

```bash
cobaya-run mavan_pngb_sample_initial_angle.yaml
```

Do not resume the old fixed-field chains as though they had sampled this
extra dimension. Once this new run exists, it can be resumed with its own
unchanged parameterization and output prefix.

The dynamic input mapping uses Cobaya's documented `drop`/`value` mechanism.
The recorded `phi_ini` is derived directly from the sampled inputs; it is
not sent to CLASS as an additional standalone input. See
[Cobaya parameterization](https://cobaya.readthedocs.io/en/latest/params_prior.html).

**Profiling template: `mavan_pngb_profile_phi_ini.yaml`**

This fixes the raw `phi_ini` to 0.4 and switches to BOBYQA minimization with
four initial attempts. To run another point, change the third entry in
`scf_parameters` and the output prefix. Keeping f free at fixed phi_ini is
different from fixing theta_ini, because theta_ini=phi_ini/f.

```bash
cobaya-run mavan_pngb_profile_phi_ini.yaml
```

The template intentionally has `ignore_prior: false`, preserving the
original nuisance/calibration prior terms. It therefore computes a
penalized, or constrained, profile. The explicit cosmological priors in
the attached input are uniform and contribute only constant densities
within their bounds; any nonuniform likelihood-component priors remain
in the objective.

For a pure likelihood profile Cobaya uses `ignore_prior: true`. That option
also drops nuisance prior densities, so external calibration information
represented as priors must be handled consistently if it is to remain part
of the profiled fit. See [Cobaya minimization](https://cobaya.readthedocs.io/en/latest/sampler_minimize.html).

Compare like objectives when constructing changes in fit quality. In
particular, Cobaya's GetDist-format header can show a posterior quantity
under a likelihood-style label; inspect the relevant likelihood and prior
contributions. Do not read a Bayesian error band directly from the sequence
of profile best-fit trajectories.

Both YAMLs use new output prefixes and `resume: false`. Their YAML syntax,
the transformation to `scf_parameters`, and retention of the supplied
theory and likelihood blocks have been checked. The custom Python CLASS
extension and cosmological likelihoods are not installed in this environment,
so no likelihood minimization or posterior sampling was performed here.
