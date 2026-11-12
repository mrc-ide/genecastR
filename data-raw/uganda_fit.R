## The precomputed fit the vignette shows.
##
## A production fit takes about four minutes on four cores, which is too long
## to run while a package is checked, and a fit short enough to run there
## would show diagnostics nobody should quote.  So the real one is computed
## here, once, and shipped.
##
## Regenerate with:
##
##   Rscript data-raw/uganda_fit.R
##
## Reproducible: `cores` changes nothing, because each chain is seeded inside
## its own worker.
##
## This writes into the SOURCE tree.  `gc_example_fit()` reads through
## system.file(), which resolves to the installed library, so reinstall
## afterwards or you will keep testing the old object:
##
##   R CMD INSTALL --preclean .
library(GeneCastR)

## `lead_in` is in steps: 6 * 16 puts t0 six years before the first survey,
## in 2010.  The default would be one year, which asserts the variant cannot
## have arrived before 2015.  The `theta` prior self-scales to the window, so
## nothing else needs changing.
d <- gc_example_data(M = 16)
el <- system.time(
  fit <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, seed = 101,
                   cores = 4, lead_in = 6 * 16)
)
set.seed(17)
pr <- fit$project(horizon = 3, n_draws = 200, n_paths = 400,
                  thresholds = c(0.01, 0.05, 0.10))

print(fit)
cat(sprintf("fit wall-clock: %.0f s on 4 cores\n", el[["elapsed"]]))
saveRDS(list(fit = fit, pr = pr), "inst/extdata/uganda_fit.rds",
        compress = "xz")
