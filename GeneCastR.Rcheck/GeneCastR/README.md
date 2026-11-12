# GeneCastR

Introduction time and spread of genetic variants under Wright–Fisher dynamics.

Fits a Wright–Fisher model with genic selection and recurrent gene flow to
panels of allele-frequency surveys, with optional per-deme hierarchies over
the growth advantage and the effective size, and reports projected
frequencies, presence and exceedance probabilities, and posteriors over the
introduction time and the founding time of the lineage that was detected.

The latent allele count is held on a lattice, so the likelihood is evaluated
**exactly** by a forward–backward pass — no particle approximation — and the
introduction time is marginalised analytically rather than sampled.

## Status

Complete and checking clean (`R CMD check --as-cran`: 0 errors, 0 warnings).
See `../RCPP_PORT_SCOPE.md` for the specification, and `vignette("GeneCastR")`
for a worked example on the Ugandan panel.

Validated against the reference Python implementation: the likelihood and the
`t*` posterior agree to 1e-10, the smoother and exceedance to 6e-12, and the
sampler recovers the reference posterior for `u` and `nu` to within Monte-Carlo
error.

| phase | | |
|---|---|---|
| 0 | package skeleton, toolchain check, CI | done |
| 1 | `GCData` and its validation | done |
| 2 | kernel, backward table, likelihood, `t*` posterior | done |
| 3 | sampler, `GCFit` | done |
| 4 | smoother, projections, `GCProjection` | done |
| 5 | FFBS, `T_first`, `T_est`, joint | done |
| 6 | plots | done |
| 7 | OpenMP, caching | deferred |
| 8 | docs, vignette, `R CMD check` | done |

## Getting started

Needs R >= 4.1 and a C++ compiler. No Fortran toolchain is required --
the package links no BLAS or LAPACK, deliberately (see `src/Makevars`).

```r
install.packages(c("R6", "Rcpp", "roxygen2", "testthat", "pkgbuild",
                   "pkgload", "posterior", "ggplot2", "scales"))

Rcpp::compileAttributes(".")                      # regenerates RcppExports
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))   # then reads them
pkgload::load_all(".")
testthat::test_local(".")
```

Regenerating the stored fit (`data-raw/uganda_fit.R`) writes into the source
tree, while `gc_example_fit()` reads from the installed library -- reinstall
after running it.

Run those two in that order after changing any `// [[Rcpp::export]]`
signature: `compileAttributes()` rewrites `R/RcppExports.R`, and roxygen
reads the result. Installing without them leaves the `.Rd` describing the old
argument list, which `R CMD check` reports as a code/documentation mismatch.

**Do not benchmark a `load_all()` build.** `pkgbuild::compile_dll()`, which
`load_all()` calls, appends `-O0` after R's own `-O2`, so the compiled code
runs roughly six times slower than the installed package. Time things against
`R CMD INSTALL` instead:

```sh
R CMD INSTALL --preclean .
```

`man/` is empty in the repository because roxygen generates it; `document()`
is the first thing to run.

Check the compiled toolchain:

```r
gc_toolchain()
#> $openmp       TRUE / FALSE
#> $max_threads   n
#> $cxx_standard 202002
```

`openmp = FALSE` is expected on macOS with Apple clang alone, and means the
per-deme loop runs serially. That is correct but slower; phase 7 needs libomp
from <https://mac.r-project.org/openmp/> first.

## Speed

The 16-district Ugandan panel, one likelihood evaluation over all demes,
against the reference NumPy implementation on the same machine:

| | GeneCastR | NumPy |
|---|---|---|
| kernel construction | 0.55 ms | 5.0 ms |
| backward pass | 3.23 ms | 2.3 ms |
| **total, `M = 16`** | **3.8 ms** | 14.4 ms |
| **total, `M = 32`** | **19.4 ms** | 42.7 ms |

Two things account for it. The binomial rows are built by the ratio
recursion `pmf(j+1)/pmf(j) = ((N-j)/(j+1)) p/(1-p)`, one exponential per row
instead of one per entry; and the kernel is stored banded, which keeps about
half the entries and drops the rest, all of them below 1e-18 of their row's
maximum. The backward pass itself is still slower than NumPy's, which reaches
52 GFLOP/s through Accelerate's matrix coprocessor; banding is what closes
the gap.

## Fitting

```r
d   <- gc_example_data(M = 16)
fit <- GCFit$new(d, chains = 4, iter = 3000, lead_in = 6 * 16)
fit$summary()
fit$tstar()          # analytic, not sampled

pr <- fit$project(horizon = 3, thresholds = c(0.01, 0.05, 0.10),
                  n_paths = 300)
pr$curves           # one row per deme and step
pr$demes            # one row per deme
```

`lead_in` is where the window opens, counted in **steps** before the first
survey, so `6 * 16` is six years at `M = 16`. It is the earliest arrival the
model will entertain; the default is one unit of time. Re-anchoring rebuilds
the panel, so `fit$data` is the window that was fitted, and `fit` prints it.

`t*` carries a constant hazard by default, so `P(t* > T)` is positive and
learned: a deme can simply not have been reached yet. `fit$draws` reports
`theta` and, as a transform of it, `p_beyond_baseline` -- the chance a deme
with no data has not been reached. The per-deme posteriors, which use each
deme's own surveys, are a different quantity and live in
`pr$demes$p_tstar_beyond`.

`gc_control()` also varies how much per-deme structure the model carries --
`neff` for the effective sizes, `selection` for the growth advantage -- and
`gc_waic_compare()` scores the alternatives. On the Ugandan panel one `A`
for the whole panel wins: per-deme selection is penalised at 1.9 standard
errors even with a hierarchy sharing it, and at 2.4 without. Sixteen
effective sizes, by contrast, tie with one, because they are shrunk so hard
that they buy a single effective parameter between them.

`n_paths = 0` (the default) skips the path sampler. `T_first` is analytic
and costs nothing; `T_est` -- the start of the excursion that produced the
first detected count -- is the only expensive output. Under the default
founding law `T_first` is identically `t*`, so sampling it is a check rather
than a necessity: the two routes agree to 0.007 on the Ugandan panel.

`pr$curves` carries the median and the predictive band for the latent
frequency, the presence probability and its across-draw quantiles, an
exceedance column per threshold, and the introduction-time pmf. Two
distinctions it is easy to lose:

- `f_q*` are **predictive**: `theta`, the latent path and the Wright-Fisher
  noise are all integrated out. `p_present_q*` hold `theta` out and say how
  well the data pin a probability, not how uncertain the outcome is.
- The `f_*` columns are **unconditional**. The atom at zero carries
  `1 - p_present`, so a band resting on zero means absence is the likelier
  state rather than a frequency near zero.

A four-chain fit of the Ugandan panel runs at about 51 ms per sweep on one
core; `4 x 4000` takes four minutes on four. Cost grows as the square of the
lattice and the per-sweep time is lumpy, because the `N_eff` posterior has a
long right tail: the top 1% of draws take a quarter of the compute.

## Plots

Each returns a ggplot faceted by deme and adds nothing you cannot remove.

```r
gc_plot_frequency(pr)              # median, predictive bands, surveys
gc_plot_present(pr)                # P(present | X) and its across-draw band
gc_plot_exceed(pr)                 # P(f > x) per threshold
gc_plot_times(pr, "cdf")           # t* and T_est together
gc_plot_times(pr, "pmf")           # densities per unit time
gc_plot_joint(pr)                  # the two times paired
gc_plot_trace(fit)                 # chains, with Rhat and ESS in the labels

plot(pr)                           # frequency
plot(fit)                          # traces
```

`gc_plot_joint()` and `gc_plot_times(which = "test")` need
`project(n_paths = ...)` and say so if they do not have it;
`gc_plot_times()` falls back to `t*` alone.

## Example

```r
d <- gc_example_data(M = 16)
d
#> <GCData>
#>   demes    : 16 (Agago, Amolatar, Arua, ...)
#>   surveys  : 96 (42 non-zero)
#>   window   : 2015 to 2022   M = 16, T = 112 steps
#>   never positive: 2 deme(s) -- Jinja, Kanungu
```
