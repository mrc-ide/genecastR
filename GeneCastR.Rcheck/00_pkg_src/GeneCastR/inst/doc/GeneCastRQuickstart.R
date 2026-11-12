## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 9,
                      fig.height = 6.5, dpi = 96)
set.seed(1)

## -----------------------------------------------------------------------------
library(GeneCastR)

## -----------------------------------------------------------------------------
d <- gc_example_data(M = 16)
head(d$surveys, 4)

## ----eval = FALSE-------------------------------------------------------------
# d <- GCData$new(surveys, M = 16,
#                 cols = c(deme = "district", time = "year",
#                          x = "n_positive", n = "n_tested"))

## -----------------------------------------------------------------------------
d

## ----fit----------------------------------------------------------------------
stored <- gc_example_tstar_fit()$panel
fit <- stored$fit
pr <- stored$pr
fit

## ----eval = FALSE-------------------------------------------------------------
# fit <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, cores = 4,
#                  lead_in = 6 * 16,
#                  control = gc_control(tstar_min = 2012, tstar_mean = 2017))

## ----eval = FALSE-------------------------------------------------------------
# fit <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, cores = 4,
#                  lead_in = 6 * 16,
#                  control = gc_control(
#                    tstar_min  = c(Koboko = 2012, Arua = 2012, Kanungu = 2014),
#                    tstar_mean = c(Koboko = 2016, Arua = 2016, Kanungu = 2019)))

## ----fig.height = 5-----------------------------------------------------------
gc_plot_trace(fit)

## ----fig.height = 5-----------------------------------------------------------
gc_plot_hist(fit)

## ----fig.height = 3.6---------------------------------------------------------
gc_plot_hist(fit, "A", ncol = 1)

## -----------------------------------------------------------------------------
A <- as.numeric(posterior::as_draws_matrix(fit$draws)[, "A"])
round(stats::quantile(A, c(.025, .5, .975)), 3)
c(`P(A > 0)` = mean(A > 0), `P(A > 0.25)` = mean(A > 0.25))

## ----project------------------------------------------------------------------
pr

## ----eval = FALSE-------------------------------------------------------------
# pr <- fit$project(horizon = 3, n_paths = 400)

## -----------------------------------------------------------------------------
gc_plot_frequency(pr)

## -----------------------------------------------------------------------------
gc_plot_present(pr)

## -----------------------------------------------------------------------------
gc_plot_exceed(pr)

## -----------------------------------------------------------------------------
head(pr$curves[, c("deme", "time", "f_mean", "p_present", "p_exceed_5")], 4)

## -----------------------------------------------------------------------------
gc_plot_times(pr, "cdf")

## ----fig.height = 7.5---------------------------------------------------------
gc_plot_times(pr, "pmf")

## ----eval = FALSE-------------------------------------------------------------
# fit$summary()        # the table printed above
# fit$draws            # every posterior draw, as a posterior package object
# pr$curves            # one row per site and time step
# pr$demes             # one row per site

