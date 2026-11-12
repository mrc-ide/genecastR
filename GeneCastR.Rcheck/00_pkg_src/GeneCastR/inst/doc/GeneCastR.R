## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 9,
                      fig.height = 6.5, dpi = 96)
set.seed(1)

## -----------------------------------------------------------------------------
library(GeneCastR)

## -----------------------------------------------------------------------------
drift_per_year <- function(M, N_eff = 15.6, f = 0.3, N_lat = NULL) {
  nl <- if (is.null(N_lat)) gc_nlat(log(N_eff), M) else N_lat
  P <- gc_kernel(nl, 0)                       # neutral
  k <- round(f * nl); g <- (0:nl) / nl
  mean_f <- sum(P[k + 1L, ] * g)
  M * sum(P[k + 1L, ] * (g - mean_f)^2)       # variance per step, times M
}
c(M8 = drift_per_year(8), M16 = drift_per_year(16),
  M32 = drift_per_year(32), M128 = drift_per_year(128),
  target = 0.3 * 0.7 / 15.6)

## -----------------------------------------------------------------------------
vapply(c(8, 16, 32, 64), drift_per_year, numeric(1), N_lat = 250L)

## -----------------------------------------------------------------------------
d <- gc_example_data(M = 16)
d

## ----fit----------------------------------------------------------------------
stored <- gc_example_fit()
fit <- stored$fit
fit

## ----eval = FALSE-------------------------------------------------------------
# fit <- GCFit$new(d, chains = 4, iter = 4000, warmup = 2000, seed = 101,
#                  cores = 4, lead_in = 6 * 16, control = gc_control(selection = "shared"))

## ----fig.height = 5-----------------------------------------------------------
gc_plot_trace(fit)

## ----fig.height = 5-----------------------------------------------------------
gc_plot_hist(fit)

## ----fig.height = 6-----------------------------------------------------------
gc_plot_pairs(fit)

## ----fig.height = 3.6---------------------------------------------------------
gc_plot_hist(fit, "A", ncol = 1)

## -----------------------------------------------------------------------------
A <- as.numeric(posterior::as_draws_matrix(fit$draws)[, "A"])
round(stats::quantile(A, c(.025, .05, .25, .5, .75, .95, .975)), 3)

## -----------------------------------------------------------------------------
c(`P(A > 0)` = mean(A > 0), `P(A > 0.25)` = mean(A > 0.25))

## -----------------------------------------------------------------------------
fit$info[, c("chain", "accept_z", "accept_u", "accept_nu", "guard_rejects")]

## ----project------------------------------------------------------------------
pr <- stored$pr
pr

## ----eval = FALSE-------------------------------------------------------------
# pr <- fit$project(horizon = 3, n_draws = 200, n_paths = 400,
#                   thresholds = c(0.01, 0.05, 0.10))

## -----------------------------------------------------------------------------
names(pr$curves)

## -----------------------------------------------------------------------------
gc_plot_frequency(pr)

## -----------------------------------------------------------------------------
gc_plot_present(pr)

## -----------------------------------------------------------------------------
gc_plot_exceed(pr)

## -----------------------------------------------------------------------------
gc_plot_times(pr, "cdf")

## ----fig.height = 7.5---------------------------------------------------------
gc_plot_times(pr, "pmf")

## -----------------------------------------------------------------------------
J <- pr$joint
share <- vapply(seq_along(pr$demes$deme),
                function(k) sum(diag(J[, , k])), numeric(1)) /
         pmax(pr$demes$p_alive_T, 1e-12)
round(sort(stats::setNames(share, pr$demes$deme), decreasing = TRUE), 2)

## -----------------------------------------------------------------------------
ctl <- gc_control(selection = "hierarchical", mu_mean = log(30), tau = 0.5)
ctl$mu_mean

## ----eval = FALSE-------------------------------------------------------------
# fit <- GCFit$new(d, chains = 4, iter = 4000, cores = 4, lead_in = 6 * 16,
#                  control = ctl)

