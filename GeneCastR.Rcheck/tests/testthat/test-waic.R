tiny <- function(control = gc_control(), iter = 150L, chains = 2L) {
  d <- GCData$new(data.frame(
    deme = rep(c("a", "b", "c"), each = 5L),
    time = rep(2016:2020, 3L),
    x    = c(0L, 1L, 4L, 6L, 8L,  0L, 0L, 0L, 1L, 2L,  0L, 0L, 0L, 0L, 0L),
    n    = rep(30L, 15L)), M = 4L)
  GCFit$new(d, chains = chains, iter = iter, warmup = iter %/% 3L, seed = 5L,
            control = control, quiet = TRUE)
}

## ---- shared effective size ---------------------------------------------

test_that("neff = 'shared' shares the demes into one size", {
  fit <- tiny(gc_control(neff = "shared"))
  m <- posterior::as_draws_matrix(fit$draws)
  ne <- as.matrix(m[, paste0("N_eff[", fit$data$demes, "]")])
  ## every draw has a single size for the whole panel
  expect_true(all(apply(ne, 1L, function(r) diff(range(r))) == 0))
  ## the between-deme spread is not a parameter and must not be invented
  expect_true(all(is.na(as.numeric(m[, "sigma"]))))
  expect_true(all(is.na(as.numeric(m[, "log_sigma2"]))))
  ## and it is dropped from the default summary rather than shown as NA
  expect_false("sigma" %in% fit$summary()$variable)
})

test_that("the hierarchy is still free to spread when asked for", {
  fit <- tiny(gc_control(neff = "hierarchical"))
  m <- posterior::as_draws_matrix(fit$draws)
  ne <- as.matrix(m[, paste0("N_eff[", fit$data$demes, "]")])
  expect_gt(mean(apply(ne, 1L, function(r) diff(range(r)))), 0)
  expect_true(all(is.finite(as.numeric(m[, "sigma"]))))
})

test_that("a shared fit still respects the N_eff guard", {
  fit <- tiny(gc_control(neff = "shared", neff_bounds = c(8, 12)))
  ne <- as.numeric(posterior::as_draws_matrix(fit$draws)[, "N_eff[a]"])
  expect_true(all(ne >= 8 - 1e-9 & ne <= 12 + 1e-9))
})

## ---- WAIC ---------------------------------------------------------------

test_that("gc_waic is computed at the deme level", {
  fit <- tiny()
  w <- gc_waic(fit, n_draws = 60L)
  ## one term per deme, not per survey: the likelihood marginalises each
  ## deme's whole path and does not factorise over its observations
  expect_equal(nrow(w$pointwise), fit$data$n_demes)
  expect_equal(w$pointwise$deme, fit$data$demes)
  expect_equal(w$waic, -2 * sum(w$pointwise$elpd), tolerance = 1e-12)
  expect_equal(w$elpd, w$lppd - w$p_waic, tolerance = 1e-12)
  expect_gt(w$p_waic, 0)
  expect_true(is.finite(w$se) && w$se > 0)
})

test_that("lppd is bounded above by the best single draw", {
  ## geometric here so the comparison likelihood can be rebuilt from theta
  ## alone; the negative binomial default is covered in test-tstar.R
  fit <- tiny(gc_control(tstar = "geometric"))
  w <- gc_waic(fit, n_draws = 60L)
  ## a posterior average of likelihoods cannot beat the best draw
  d <- fit$data
  m <- posterior::as_draws_matrix(fit$draws)
  idx <- unique(round(seq(1, nrow(m), length.out = 60L)))
  Z <- log(as.matrix(m[idx, seq_len(d$n_demes), drop = FALSE]))
  dimnames(Z) <- NULL
  best <- max(vapply(seq_along(idx), function(i)
    sum(gc_loglik_pt(gc_nlat(Z[i, ], d$M), gc_alpha(m[idx[i], "u"], d$M),
                     gc_lambda(as.numeric(m[idx[i], "nu"]), exp(Z[i, ])),
                     d$T, d$deme_data(), 2L,
                     -expm1(-as.numeric(m[idx[i], "theta"]) / d$M))),
    numeric(1L)))
  expect_lte(w$lppd, best + 1e-8)
})

test_that("the comparison orders by WAIC and pairs the standard error", {
  a <- tiny(gc_control(neff = "hierarchical"))
  b <- tiny(gc_control(neff = "shared"))
  cmp <- gc_waic_compare(hierarchical = a, shared = b, n_draws = 60L)
  expect_equal(nrow(cmp), 2L)
  expect_setequal(cmp$model, c("hierarchical", "shared"))
  expect_true(!is.unsorted(cmp$waic))            # best first
  expect_equal(cmp$d_waic[1L], 0)
  expect_equal(cmp$se_diff[1L], 0)               # a model against itself
  expect_gt(cmp$se_diff[2L], 0)
  expect_equal(cmp$d_waic[2L], cmp$waic[2L] - cmp$waic[1L], tolerance = 1e-12)
  ## the pairwise se is NOT the difference of the two marginal ones: it is
  ## computed from the paired per-deme contributions, which is the whole
  ## reason a non-best comparison needs its own
  expect_lt(cmp$se_diff[2L], max(cmp$se))
})

test_that("comparing a fit with itself gives exactly no difference", {
  fit <- tiny()
  cmp <- gc_waic_compare(one = fit, two = fit, n_draws = 40L)
  expect_equal(diff(cmp$waic), 0, tolerance = 1e-12)
  expect_equal(cmp$se_diff, c(0, 0))
})

test_that("WAIC works under either t* prior", {
  for (tp in c("geometric", "uniform")) {
    w <- gc_waic(tiny(gc_control(tstar = tp)), n_draws = 40L)
    expect_true(is.finite(w$waic), info = tp)
    expect_equal(nrow(w$pointwise), 3L, info = tp)
  }
})

## ---- per-deme selection -------------------------------------------------

test_that("selection = 'shared' keeps one growth advantage", {
  fit <- tiny()
  v <- posterior::variables(fit$draws)
  expect_false(any(grepl("^u\\[", v)))
  expect_true(all(c("u", "A") %in% v))
})

test_that("the hierarchy gives each deme its own, on the same scale", {
  for (sel in "hierarchical") {
    fit <- tiny(gc_control(selection = sel))
    m <- posterior::as_draws_matrix(fit$draws)
    ## strip the draws_matrix class: the numbers are what is being tested,
    ## not the attributes posterior carries along
    bare <- function(x) matrix(as.numeric(x), nrow = nrow(x))
    ud <- bare(as.matrix(m[, paste0("u[", fit$data$demes, "]")]))
    ad <- bare(as.matrix(m[, paste0("A[", fit$data$demes, "]")]))
    ## the same transform the shared model uses
    expect_equal(ad, expm1(ud), tolerance = 1e-12, info = sel)
    ## they genuinely differ between demes
    expect_gt(mean(apply(ud, 1L, function(r) diff(range(r)))), 0)
    ## and the reported scalar is the panel average
    expect_equal(as.numeric(m[, "u"]), unname(rowMeans(ud)),
                 tolerance = 1e-12, info = sel)
  }
})

test_that("only the hierarchical option carries hyperparameters", {
  h <- tiny(gc_control(selection = "hierarchical"))
  expect_true(all(c("mu_u", "sigma_u") %in% posterior::variables(h$draws)))
  m <- posterior::as_draws_matrix(h$draws)
  expect_true(all(as.numeric(m[, "sigma_u"]) > 0))
  expect_true(all(c("mu_u", "sigma_u") %in% h$summary()$variable))

  p <- tiny()
  expect_false(any(c("mu_u", "sigma_u") %in% posterior::variables(p$draws)))
})

test_that("WAIC handles every selection mode", {
  for (sel in c("shared", "hierarchical")) {
    w <- gc_waic(tiny(gc_control(selection = sel)), n_draws = 40L)
    expect_true(is.finite(w$waic), info = sel)
    expect_equal(nrow(w$pointwise), 3L, info = sel)
  }
})
