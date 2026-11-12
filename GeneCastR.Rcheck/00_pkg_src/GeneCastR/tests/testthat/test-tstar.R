## The prior on t* enters only as log-weights added to the conditional
## likelihoods, so these tests pin the weights and then check that routing a
## law through them changes nothing it should not.

small2 <- function(M = 4L) {
  GCData$new(data.frame(deme = rep(c("a", "b"), each = 4L),
                        time = rep(2016:2019, 2L),
                        x = c(0L, 1L, 4L, 6L, 0L, 0L, 0L, 0L),
                        n = rep(30L, 8L)), M = M, quiet = TRUE)
}

test_that("a negative binomial of size 1 is the geometric, exactly", {
  T <- 48L; M <- 4L
  p <- 1 - exp(-0.5 / M)                       # the hazard the package uses
  geo <- c((0:(T - 1L)) * log1p(-p) + log(p), T * log1p(-p))
  nb <- gc_tstar_negbin_logw(T, 1L, (1 - p) / p, 1)
  expect_equal(nb, geo, tolerance = 1e-12)
  ## and the whole likelihood, not just the weights
  d <- small2(M); dd <- d$deme_data()
  nl <- gc_nlat(rep(log(12), 2L), M); lam <- gc_lambda(0.02, rep(12, 2L))
  al <- gc_alpha(log(1.2), M)
  a <- gc_loglik_pt(nl, al, lam, d$T, dd, 2L, p)
  LW <- matrix(gc_tstar_negbin_logw(d$T, 1L, (1 - p) / p, 1), d$T + 1L, 2L)
  b <- gc_loglik_pt(nl, al, lam, d$T, dd, 2L, -1, LW)
  expect_equal(a, b)
})

test_that("the weights are a proper distribution, tail included", {
  for (tmin in c(1L, 5L, 30L)) for (mu in c(2, 40, 500)) for (r in c(0.3, 1, 7)) {
    lw <- gc_tstar_negbin_logw(60L, tmin, mu, r)
    expect_equal(sum(exp(lw)), 1, tolerance = 1e-12,
                 info = paste(tmin, mu, r))
    ## nothing before tmin
    if (tmin > 1L) expect_true(all(is.infinite(lw[seq_len(tmin - 1L)])))
    expect_true(is.finite(lw[tmin]))
  }
})

test_that("the tail is computed from pnbinom, not by subtraction", {
  ## with the mean far beyond the window almost all the mass is in the tail,
  ## where 1 - sum(pmf) would be catastrophic
  lw <- gc_tstar_negbin_logw(20L, 1L, 1e5, 1)
  expect_equal(exp(lw[21L]), pnbinom(19, size = 1, mu = 1e5,
                                     lower.tail = FALSE), tolerance = 1e-14)
  expect_gt(exp(lw[21L]), 0.999)
})

test_that("t_min shifts the law and the shape moves the tail", {
  tail <- function(tmin, mu, r)
    exp(gc_tstar_negbin_logw(192L, tmin, mu, r)[193L])
  ## later t_min, same offset mean -> more mass beyond the window
  expect_gt(tail(80L, 100, 1), tail(1L, 100, 1))
  ## the tail is not monotone in the shape at fixed mean, which is why the
  ## shape is not a reparameterisation of p_beyond
  t_small <- tail(1L, 64, 0.08); t_mid <- tail(1L, 64, 0.21)
  t_big <- tail(1L, 64, 5)
  expect_gt(t_mid, t_small)
  expect_gt(t_mid, t_big)
})

test_that("weights may differ by deme", {
  d <- small2(); dd <- d$deme_data()
  nl <- gc_nlat(rep(log(12), 2L), 4L); lam <- gc_lambda(0.02, rep(12, 2L))
  al <- gc_alpha(log(1.2), 4L)
  LW <- cbind(gc_tstar_negbin_logw(d$T, 1L, 20, 1),
              gc_tstar_negbin_logw(d$T, 9L, 20, 1))   # deme b cannot arrive early
  one <- matrix(gc_tstar_negbin_logw(d$T, 1L, 20, 1), d$T + 1L, 2L)
  a <- gc_loglik_pt(nl, al, lam, d$T, dd, 2L, -1, one)
  b <- gc_loglik_pt(nl, al, lam, d$T, dd, 2L, -1, LW)
  expect_equal(a[1], b[1])          # deme a saw the same prior
  expect_false(isTRUE(all.equal(a[2], b[2])))
  expect_error(gc_loglik_pt(nl, al, lam, d$T, dd, 2L, -1, LW[, 1, drop = FALSE]),
               "must be")
})

## --- the prior as GCFit resolves it ----------------------------------------

d3 <- function() GCData$new(data.frame(
  deme = rep(c("a", "b", "c"), each = 4L), time = rep(2016:2019, 3L),
  x = c(0L, 1L, 4L, 6L, 0L, 0L, 1L, 2L, 0L, 0L, 0L, 0L),
  n = rep(30L, 12L)), M = 4L, quiet = TRUE)

nbfit <- function(..., iter = 60L) {
  GCFit$new(d3(), chains = 1L, iter = iter, warmup = iter %/% 3L, seed = 4L,
            lead_in = 8L, control = gc_control(...), quiet = TRUE)
}

test_that("negbin is the default and its mean self-scales", {
  cc <- gc_control()
  expect_identical(cc$tstar, "negbin")
  expect_identical(cc$tstar_learn, "shape")
  expect_null(cc$tstar_min)
  expect_null(cc$tstar_mean)
  ## the default mean puts the median arrival at the end of the record when
  ## the shape is one, whatever the span -- so P(t* <= T) is a coin flip
  for (span in c(7, 12, 20)) {
    M <- 16L; T <- as.integer(M * span)
    lw <- gc_tstar_negbin_logw(T, 1L, M * span / log(2), 1)
    expect_equal(sum(exp(lw[seq_len(T)])), 0.5, tolerance = 0.01,
                 info = paste("span", span))
  }
})

test_that("t_min and the mean may be named by deme", {
  f <- nbfit(tstar_min = c(a = 2015, c = 2016), tstar_mean = 2018)
  ts <- GeneCastR:::.tstar_setup(f$control, f$data)
  ## t0 is 2016 - 8/4 = 2014; steps are M * (t_min - t0), floored at 1
  expect_equal(f$data$t0, 2014)
  expect_equal(unname(ts$tmin_step), c(4L, 1L, 8L))
  expect_false(ts$panel_wide)
  ## a per-deme prior makes the panel-level tail ill-defined, so it is not
  ## reported; the per-deme version lives on the projection
  expect_false("p_beyond_baseline" %in% posterior::variables(f$draws))
  expect_true("tstar_shape" %in% posterior::variables(f$draws))
})

test_that("the per-deme prior is validated against the data", {
  expect_error(nbfit(tstar_min = c(nope = 2015)), "not in the data")
  expect_error(nbfit(tstar_min = 2010), "before the window opens")
  ## deme a is first positive in 2017; a t_min after that contradicts it
  expect_error(nbfit(tstar_min = c(a = 2018)),
               "after the first non-zero survey")
  expect_error(nbfit(tstar_min = 2015, tstar_mean = 2014),
               "must be after")
  expect_error(nbfit(tstar_min = stats::setNames(2015, NULL)[c(1, 1)]),
               "single value or a vector named by deme")
  ## the shape is constrained to at least one
  expect_error(gc_control(tstar_shape = 0.5), "tstar_shape")
})

test_that("tstar_learn controls what moves, and none really freezes", {
  ## shape only
  a <- nbfit(tstar_mean = 2018)
  expect_true("tstar_shape" %in% posterior::variables(a$draws))
  expect_false("tstar_mean_offset" %in% posterior::variables(a$draws))
  ## both
  b <- nbfit(tstar_mean = 2018, tstar_learn = "both", iter = 120L)
  expect_true(all(c("tstar_shape", "tstar_mean_offset") %in%
                  posterior::variables(b$draws)))
  expect_true(all(as.numeric(posterior::as_draws_matrix(b$draws)[
    , "tstar_mean_offset"]) > 0))
  ## none: nothing about t* is a parameter, so nothing about it is reported.
  ## Setting the prior spread to zero instead would have pinned the shape to
  ## one rather than to the value asked for.
  c_ <- nbfit(tstar_mean = 2018, tstar_learn = "none", tstar_shape = 3)
  expect_false(any(c("tstar_shape", "tstar_mean_offset",
                     "p_beyond_baseline") %in%
                   posterior::variables(c_$draws)))
  ## and the frozen shape is the one requested: a fit at shape 3 must differ
  ## from one at shape 1
  e <- nbfit(tstar_mean = 2018, tstar_learn = "none", tstar_shape = 1)
  expect_false(isTRUE(all.equal(
    as.numeric(posterior::as_draws_matrix(c_$draws)[, "u"]),
    as.numeric(posterior::as_draws_matrix(e$draws)[, "u"]))))
})
