pr <- function(...) .prior_list(gc_control(...), th_rate = 7)

small <- function() {
  ## three demes, short window: enough to exercise every block cheaply
  d <- data.frame(
    deme = rep(c("a","b","c"), each = 4L),
    time = rep(2016:2019, 3L),
    x    = c(0L, 1L, 4L, 6L,  0L, 0L, 1L, 2L,  0L, 0L, 0L, 0L),
    n    = rep(30L, 12L))
  GCData$new(d, M = 4L)
}

test_that("the sampler returns the documented shape", {
  d <- small()
  set.seed(2)
  r <- gc_mwg_run(d$deme_data(), M = d$M, T = d$T, n_iter = 60L, warm = 20L,
                  z0 = rep(log(15), 3L), u0 = 0.2, nu0 = 0.005,
                  step_z = 0.9, step_u = 0.25, step_nu = 0.6,
                  prior = pr(), target = 0.3)
  ## layout: D effective sizes, then u / mu / log_sigma2 / nu / theta, then
  ## D per-deme u, then mu_u, sigma_u, and the t* shape and mean offset.  The
  ## columns a given run does not use are NA rather than absent, so this
  ## matrix does not change width with the control settings.
  expect_equal(dim(r$draws), c(60L, 2L * 3L + 9L))
  expect_named(r$accept, c("z", "u", "nu", "theta", "shape", "mean"))
  expect_named(r$steps, c("z", "u", "nu", "theta", "shape", "mean"))
  expect_true(all(is.finite(r$draws[, seq_len(7L)])))
  expect_true(all(is.na(r$draws[, 8L])))      # theta: t* is uniform here
  expect_true(all(is.na(r$draws[, 12L:13L]))) # mu_u, sigma_u: A is shared
  expect_true(all(r$draws[, 7L] > 0))          # nu on the natural scale
  expect_true(all(r$accept >= 0 & r$accept <= 1))
  expect_true(all(r$steps > 0))
})

test_that("the guard is a pre-likelihood rejection and it bites", {
  d <- small()
  set.seed(3)
  wide <- gc_mwg_run(d$deme_data(), d$M, d$T, 60L, 20L, rep(log(15), 3L),
                     0.2, 0.005, 0.9, 0.25, 0.6, pr(), 0.3)
  ## a guard that admits almost nothing must reject almost everything
  set.seed(3)
  tight <- gc_mwg_run(d$deme_data(), d$M, d$T, 60L, 20L, rep(log(15), 3L),
                      0.2, 0.005, 0.9, 0.25, 0.6,
                      pr(neff_bounds = c(14, 16)), 0.3)
  expect_gt(tight$guard_rejects, wide$guard_rejects)
  ## and N_eff never leaves the bounds it was given
  expect_true(all(exp(tight$draws[, 1:3]) >= 14 - 1e-9))
  expect_true(all(exp(tight$draws[, 1:3]) <= 16 + 1e-9))
})

test_that("draws respect the hard bounds on u and log nu", {
  d <- small()
  set.seed(4)
  r <- gc_mwg_run(d$deme_data(), d$M, d$T, 80L, 30L, rep(log(15), 3L),
                  0.2, 0.005, 2.0, 2.0, 2.0,
                  pr(u_bounds = c(-0.5, 0.5), lognu_bounds = c(-6, -4)), 0.3)
  expect_true(all(r$draws[, 4L] >= -0.5 & r$draws[, 4L] <= 0.5))
  expect_true(all(log(r$draws[, 7L]) >= -6 & log(r$draws[, 7L]) <= -4))
})

test_that("adaptation happens only during warm-up", {
  d <- small()
  set.seed(5)
  none <- gc_mwg_run(d$deme_data(), d$M, d$T, 60L, 0L, rep(log(15), 3L),
                     0.2, 0.005, 0.9, 0.25, 0.6, pr(), 0.3)
  expect_equal(unname(none$steps)[1:3], c(0.9, 0.25, 0.6))
  set.seed(5)
  some <- gc_mwg_run(d$deme_data(), d$M, d$T, 60L, 55L, rep(log(15), 3L),
                     0.2, 0.005, 0.9, 0.25, 0.6, pr(), 0.3)
  expect_false(isTRUE(all.equal(unname(some$steps)[1:3], c(0.9, 0.25, 0.6))))
})

test_that("the chain is reproducible from the R seed", {
  d <- small()
  set.seed(7); a <- gc_mwg_run(d$deme_data(), d$M, d$T, 40L, 10L,
                               rep(log(15), 3L), 0.2, 0.005, 0.9, 0.25, 0.6,
                               pr(), 0.3)
  set.seed(7); b <- gc_mwg_run(d$deme_data(), d$M, d$T, 40L, 10L,
                               rep(log(15), 3L), 0.2, 0.005, 0.9, 0.25, 0.6,
                               pr(), 0.3)
  expect_identical(a$draws, b$draws)
  set.seed(8); cc <- gc_mwg_run(d$deme_data(), d$M, d$T, 40L, 10L,
                                rep(log(15), 3L), 0.2, 0.005, 0.9, 0.25, 0.6,
                                pr(), 0.3)
  expect_false(identical(a$draws, cc$draws))
})

test_that("bad arguments are rejected", {
  d <- small()
  expect_error(gc_mwg_run(d$deme_data(), d$M, d$T, 10L, 10L, rep(0, 3L),
                          0, 0.005, .9, .25, .6, pr(), .3), "0 <= warm")
  expect_error(gc_mwg_run(d$deme_data(), d$M, d$T, 10L, 2L, rep(0, 2L),
                          0, 0.005, .9, .25, .6, pr(), .3), "one entry per deme")
})

## ---- the R6 wrapper ----------------------------------------------------

test_that("GCFit assembles draws on the reported scales", {
  d <- small()
  ## `sigma` is a parameter only when the sizes are hierarchical, which is
  ## not the default, so this asks for it explicitly
  fit <- GCFit$new(d, chains = 2L, iter = 90L, warmup = 30L, seed = 3L,
                   control = gc_control(neff = "hierarchical",
                                        tstar = "geometric"), quiet = TRUE)
  expect_s3_class(fit$draws, "draws_array")
  expect_equal(posterior::ndraws(fit$draws), 120L)
  expect_equal(posterior::nchains(fit$draws), 2L)
  v <- posterior::variables(fit$draws)
  expect_true(all(c("u", "A", "nu", "mu", "sigma", "exp_mu",
                    "N_eff[a]", "N_eff[c]") %in% v))

  m <- posterior::as_draws_matrix(fit$draws)
  ## the derived scales are exactly the transformations they claim to be
  expect_equal(as.numeric(m[, "A"]), expm1(as.numeric(m[, "u"])))
  expect_equal(as.numeric(m[, "exp_mu"]), exp(as.numeric(m[, "mu"])))
  expect_equal(as.numeric(m[, "sigma"]),
               exp(as.numeric(m[, "log_sigma2"]) / 2))
  expect_true(all(m[, "N_eff[a]"] > 0))

  s <- fit$summary()
  expect_true(all(c("variable","median","q5","q95","rhat","ess_bulk")
                  %in% names(s)))
  ## the default t* prior is geometric, so its hazard and the implied
  ## baseline chance of not having been reached are reported too
  expect_setequal(s$variable, c("u", "A", "nu", "exp_mu", "sigma",
                                "theta", "p_beyond_baseline"))
  expect_setequal(GCFit$new(small(), chains = 1L, iter = 60L, warmup = 20L,
                            seed = 3L, quiet = TRUE,
                            control = gc_control(tstar = "uniform",
                                                 neff = "hierarchical")
                            )$summary()$variable,
                  c("u", "A", "nu", "exp_mu", "sigma"))
  expect_equal(nrow(fit$info), 2L)
  expect_true(is.numeric(attr(fit$info, "wall")))
})

test_that("GCFit$tstar() is a proper distribution per deme", {
  d <- small()
  fit <- GCFit$new(d, chains = 1L, iter = 60L, warmup = 20L, seed = 4L,
                   quiet = TRUE)
  ts <- fit$tstar(n_draws = 10L)
  expect_equal(nrow(ts), d$n_demes * d$T)
  ## under the default geometric prior the in-record mass falls short of
  ## one by exactly the chance the deme has not been reached
  agg <- as.numeric(tapply(ts$prob, ts$deme, sum)[d$demes])
  expect_equal(agg + unname(attr(ts, "p_beyond")[d$demes]),
               rep(1, d$n_demes), tolerance = 1e-10)
  expect_true(all(ts$prob >= 0))
  ## t* cannot postdate the first non-zero survey
  for (nm in names(d$cut)) {
    k <- d$cut[[nm]]
    if (!is.na(k) && k < d$T) {
      late <- ts$prob[ts$deme == nm & ts$step > k]
      expect_equal(sum(late), 0, info = nm)
    }
  }
  expect_equal(ts$time, d$t0 + ts$step / d$M)
})

test_that("GCFit rejects bad arguments", {
  d <- small()
  expect_error(GCFit$new(42), "must be a GCData")
  expect_error(GCFit$new(d, iter = 5L), "at least 10")
  expect_error(GCFit$new(d, iter = 50L, warmup = 50L), "\\[0, iter\\)")
})


test_that("a geometric t* adds a hazard block that actually moves", {
  d <- small()
  set.seed(31)
  r <- gc_mwg_run(d$deme_data(), d$M, d$T, 300L, 100L, rep(log(15), 3L),
                  0.2, 0.005, 0.9, 0.25, 0.6,
                  .prior_list(gc_control(), th_rate = d$span),
                  0.3, founding = 2L, tstar_p = -expm1(-0.3 / d$M))
  th <- r$draws[, 8L]
  expect_true(all(is.finite(th)))
  expect_true(all(th > 0))
  expect_gt(stats::sd(th), 0)                       # it is being sampled
  expect_gt(r$accept[["theta"]], 0)
  expect_true(all(log(th) > -8 & log(th) < 3))      # bounds respected
})

## --- lead_in ---------------------------------------------------------------

test_that("lead_in re-anchors the window in steps, not time units", {
  d <- small()                                   # M = 4, t0 = 2015 by default
  expect_equal(d$t0, 2015)

  f <- GCFit$new(d, chains = 1L, iter = 20L, warmup = 5L, lead_in = 3L * d$M,
                 quiet = TRUE)
  expect_equal(f$data$t0, 2013)                  # three time units, 12 steps
  expect_equal(f$data$span, 6)
  ## the panel the caller handed over is untouched
  expect_equal(d$t0, 2015)
})

test_that("lead_in = M reproduces the default window", {
  d <- small()
  f <- GCFit$new(d, chains = 1L, iter = 20L, warmup = 5L, lead_in = d$M,
                 quiet = TRUE)
  expect_equal(f$data$t0, d$t0)
  expect_equal(f$data$T, d$T)
})

test_that("lead_in and an explicit t0 give the same fit bit-for-bit", {
  dd <- data.frame(
    deme = rep(c("a","b","c"), each = 4L),
    time = rep(2016:2019, 3L),
    x    = c(0L, 1L, 4L, 6L,  0L, 0L, 1L, 2L,  0L, 0L, 0L, 0L),
    n    = rep(30L, 12L))
  by_t0   <- GCData$new(dd, M = 4L, t0 = 2013)
  by_lead <- GCData$new(dd, M = 4L)
  a <- GCFit$new(by_t0, chains = 2L, iter = 60L, warmup = 20L, seed = 5L,
                 quiet = TRUE)
  b <- GCFit$new(by_lead, chains = 2L, iter = 60L, warmup = 20L, seed = 5L,
                 lead_in = 12L, quiet = TRUE)
  expect_equal(a$data$deme_data(), b$data$deme_data())
  expect_equal(as.numeric(posterior::as_draws_matrix(a$draws)),
               as.numeric(posterior::as_draws_matrix(b$draws)))
})

test_that("lead_in must be a whole number of steps, at least one", {
  d <- small()
  expect_error(GCFit$new(d, lead_in = 2.5, quiet = TRUE), "whole number")
  expect_error(GCFit$new(d, lead_in = 0L, quiet = TRUE), "whole number")
  expect_error(GCFit$new(d, lead_in = c(4L, 8L), quiet = TRUE), "whole number")
  ## and the error says what a three-year lead-in would be at this M
  expect_error(GCFit$new(d, lead_in = 0L, quiet = TRUE), "lead_in = 12")
})

test_that("a moved lead-in is announced unless quiet", {
  d <- small()
  expect_message(
    GCFit$new(d, chains = 1L, iter = 20L, warmup = 5L, lead_in = 3L * d$M),
    "t0 moved from 2015 to 2013")
  expect_silent(
    GCFit$new(d, chains = 1L, iter = 20L, warmup = 5L, lead_in = d$M,
              quiet = TRUE))
})

test_that("gc_control names its settings as its arguments do", {
  cc <- gc_control(mu_mean = log(30), neff_bounds = c(4, 200))
  ## the trap this avoids: reading back a setting under the name you passed
  expect_equal(cc$mu_mean, log(30))
  expect_equal(cc$neff_bounds, c(4, 200))
  expect_true(all(names(cc) %in% names(formals(gc_control))))
  expect_setequal(names(cc), names(formals(gc_control)))
  ## and the compiled sampler still gets its own flat names
  pl <- .prior_list(cc, th_rate = 5)
  expect_equal(pl$mu0, log(30))
  expect_equal(c(pl$ne_lo, pl$ne_hi), c(4, 200))
  expect_equal(pl$th_rate, 5)
})

test_that("shared is the default for both A and the effective sizes", {
  cc <- gc_control()
  expect_identical(cc$selection, "shared")
  expect_identical(cc$neff, "shared")
})

test_that("a default fit reports one effective size and no sigma", {
  d <- small()
  f <- GCFit$new(d, chains = 2L, iter = 60L, warmup = 20L, quiet = TRUE)
  m <- posterior::as_draws_matrix(f$draws)
  ## sigma is always a column but is NA unless the sizes are hierarchical
  expect_true(all(is.na(m[, "sigma"])))
  expect_false("sigma" %in% f$summary()$variable)
  ## the panel-level size is named for the structure that was fitted: with
  ## one shared size there is no mean to take, so it is N_eff, not exp_mu
  v <- posterior::variables(f$draws)
  expect_true(all(c("N_eff", "log_N_eff") %in% v))
  expect_false(any(c("exp_mu", "mu") %in% v))
  expect_equal(as.numeric(m[, "N_eff"]), exp(as.numeric(m[, "log_N_eff"])))
  ## and every deme carries that one size
  ne <- m[, paste0("N_eff[", d$demes, "]")]
  expect_equal(apply(ne, 1, stats::sd), rep(0, nrow(m)), ignore_attr = TRUE)
  expect_equal(as.numeric(ne[, 1]), as.numeric(m[, "N_eff"]))
})

test_that("default plots skip variables the fitted structure does not have", {
  skip_if_not_installed("ggplot2")
  d <- small()
  f <- GCFit$new(d, chains = 2L, iter = 60L, warmup = 20L, quiet = TRUE)
  expect_s3_class(gc_plot_trace(f), "ggplot")
  expect_false("sigma" %in% levels(gc_plot_trace(f)$data$variable))
  expect_s3_class(gc_plot_hist(f), "ggplot")
  ## asking for nothing but an unused structure is an error, not an empty panel
  expect_error(gc_plot_trace(f, "sigma"), "not parameters of this fit")
})

## --- the half-t prior on sigma_u -------------------------------------------

test_that("sigma_u carries a half-Cauchy by default, not an inverse-gamma", {
  cc <- gc_control()
  expect_equal(cc$u_sigma_scale, 0.5)
  expect_equal(cc$u_sigma_df, 1)      # df = 1 is half-Cauchy
  expect_error(gc_control(u_sigma_scale = 0), "u_sigma_scale")
  expect_error(gc_control(u_sigma_df = -1), "u_sigma_df")
})

test_that("the u hierarchy can collapse toward a single shared advantage", {
  d <- small()
  f <- GCFit$new(d, chains = 2L, iter = 1500L, warmup = 500L, seed = 3L,
                 control = gc_control(selection = "hierarchical"),
                 quiet = TRUE)
  s <- as.numeric(posterior::as_draws_matrix(f$draws)[, "sigma_u"])
  expect_true(all(s > 0))
  ## An InvGamma(2, 1/2) on sigma_u^2 puts 5e-5 of its *prior* mass below
  ## 0.2, so the hierarchy could never collapse.  The half-t has its mode at
  ## zero, and the posterior duly reaches it.
  expect_lt(stats::quantile(s, 0.05), 0.2)
  expect_lt(min(s), 0.05)
})

test_that("mu_u shares the prior the shared advantage gets", {
  ## as sigma_u -> 0 the hierarchical model must become the shared one
  ## exactly, which needs mu_u ~ N(0, tau^2) and not a prior of its own
  d <- small()
  f <- GCFit$new(d, chains = 2L, iter = 400L, warmup = 100L, seed = 4L,
                 control = gc_control(selection = "hierarchical", tau = 0.3),
                 quiet = TRUE)
  m <- posterior::as_draws_matrix(f$draws)
  expect_true(all(abs(as.numeric(m[, "mu_u"])) < 4))
  ## the panel u is the average of the per-deme values
  ud <- m[, paste0("u[", d$demes, "]")]
  expect_equal(as.numeric(rowMeans(ud)), as.numeric(m[, "u"]))
})
