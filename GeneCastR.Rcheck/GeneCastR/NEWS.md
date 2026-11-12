# GeneCastR 0.0.0.9000

First working version. The model is Wright-Fisher with genic selection and
Poisson flow out of the empty state, with an introduction time `t*`
marginalised analytically.

* `GCData` validates a survey panel onto the model's lattice.
* `GCFit$new(lead_in = )` sets where the window opens, in steps before the
  first survey, so a three-year lead-in at `M = 16` is `lead_in = 3 * 16`.
  This is the earliest arrival the model will entertain, and it is a prior
  choice rather than a property of the panel, which is why it sits on the fit
  rather than on `GCData`.
* `GCFit` runs blocked Metropolis-within-Gibbs against the **exact**
  likelihood -- no particle approximation, so this is not a pseudo-marginal
  scheme.
* `GCProjection` carries smoothed curves, projections, presence and
  exceedance probabilities, and the introduction times, as one row per deme
  and step.
* The founding count at `t*` is `Poisson(lambda)` conditioned on being
  non-zero -- the same law that governs every later re-entry -- so the first
  arrival is not a special case. `gc_control(founding = )` also offers
  `"one"` (exactly one copy) and `"rate"` (`t*` as the onset of the arrival
  rate rather than an arrival).
* `t*` and `T_first` are analytic. `T_est` comes from sampled paths and is
  the start of the excursion that produced the first detected count, so it
  is bounded by that detection and undefined where there was none. Their
  joint distribution, also free, has the probability that the detected
  lineage is the original arrival on its diagonal.
* `t*` carries a geometric prior by default, so arrival may postdate the
  record and `P(t* > T)` is learned rather than assumed zero. The hazard
  `theta` has an `Exponential(window)` prior, which makes
  `P(founded within the record)` uniform.
* `gc_control(neff = )` and `gc_control(selection = )` vary how much
  per-deme structure the model carries: shared or hierarchical effective
  sizes, and a shared or hierarchical growth advantage. **Both
  default to `"shared"`**, the simpler model, so a hierarchy is asked for
  rather than assumed.
* `sigma_u`, the spread of the per-deme growth advantages under
  `selection = "hierarchical"`, carries a half-`t` prior whose **mode is at
  zero** (`u_sigma_scale`, `u_sigma_df`; half-Cauchy by default), sampled by
  the Huang-Wand expansion so both conditionals stay conjugate. An
  `InvGamma` on a variance cannot reach zero and so cannot express complete
  sharing, which would leave the hierarchy held apart by the prior in exactly
  the comparison it is meant to inform. `sigma`, the spread of the effective
  sizes, still carries `InvGamma(sigma2_shape, sigma2_rate)`; use
  `neff = "shared"` to ask whether the demes differ at all.
* `gc_control()` names the settings it returns exactly as its own arguments,
  so a value can be read back under the name it was passed and each pair of
  bounds stays one vector.
* `gc_waic()` and `gc_waic_compare()`, computed at the deme level because
  the likelihood does not factorise over a deme's surveys.
* Nine plotting functions, each returning a ggplot.

Validated against the reference Python implementation: the likelihood and
the `t*` posterior agree to 1e-10, the smoother and exceedance to 6e-12.
