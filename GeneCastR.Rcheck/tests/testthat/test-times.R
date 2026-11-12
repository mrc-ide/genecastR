small <- function() {
  d <- data.frame(
    deme = rep(c("a", "b", "c"), each = 5L),
    time = rep(2016:2020, 3L),
    x    = c(0L, 1L, 4L, 6L, 8L,   0L, 0L, 0L, 1L, 2L,   0L, 0L, 0L, 0L, 0L),
    n    = rep(30L, 15L))
  GCData$new(d, M = 4L)
}

draws <- function(d, n = 1L) {
  list(Z = matrix(log(12), nrow = n, ncol = d$n_demes),
       U = rep(0.3, n), NU = rep(0.03, n))
}

test_that("the path sampler returns the documented shapes", {
  d <- small(); p <- draws(d)
  set.seed(1)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 200L)
  expect_equal(dim(tt$pmf_tfirst), c(d$n_demes, d$T))
  expect_equal(dim(tt$pmf_test), c(d$n_demes, d$T))
  expect_length(tt$p_alive, d$n_demes)
  expect_null(tt$joint)
  expect_true(all(tt$pmf_tfirst >= 0) && all(tt$pmf_test >= 0))
  ## p_alive is a sum of n_paths equal increments, so it lands within an
  ## ulp or two of 1 when every path survives
  expect_true(all(tt$p_alive >= 0 & tt$p_alive <= 1 + 1e-12))
})

test_that("T_first is a distribution; T_est exists only where detected", {
  d <- small(); p <- draws(d)
  set.seed(2)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 400L)
  ## every path is founded inside the window, so T_first always exists
  expect_equal(rowSums(tt$pmf_tfirst), rep(1, d$n_demes), tolerance = 1e-12)
  ## T_est is the start of the excursion that produced the first detected
  ## count, so it is a proper distribution where there was a detection and
  ## absent where there was none -- it does not condition on survival
  seen <- !is.na(d$cut)
  expect_equal(rowSums(tt$pmf_test)[seen], rep(1, sum(seen)),
               tolerance = 1e-12)
  expect_equal(rowSums(tt$pmf_test)[!seen], rep(0, sum(!seen)),
               tolerance = 1e-12)
  expect_true(all(tt$p_alive <= 1 + 1e-12))
})

test_that("T_est cannot postdate the detection that defines it", {
  d <- small(); p <- draws(d)
  set.seed(21)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 400L)
  for (i in seq_along(d$demes)) {
    k <- d$cut[[i]]
    if (!is.na(k) && k < d$T)
      expect_equal(sum(tt$pmf_test[i, (k + 1L):d$T]), 0, info = d$demes[i])
  }
})

test_that("T_est is at or after t*, path by path", {
  d <- small(); p <- draws(d)
  set.seed(22)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 400L,
                    joint = TRUE)
  for (k in seq_len(d$n_demes)) {
    m <- tt$joint[, , k]
    expect_equal(sum(m[lower.tri(m)]), 0, info = d$demes[k])
  }
})

test_that("more data means more certainty that the variant is still there", {
  d <- small(); p <- draws(d)
  set.seed(3)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 400L)
  ## deme a rises to 8/30, c is never seen
  expect_gt(tt$p_alive[1L], tt$p_alive[3L])
  expect_equal(tt$p_alive[1L], 1, tolerance = 0.02)
})

test_that("T_first reproduces the analytic t* posterior", {
  ## The two are the same quantity: the founding event deposits a copy at t*
  ## and nothing is present before it.  They are computed by completely
  ## different routes -- one a normalisation of the backward row, one from
  ## sampled paths -- so agreeing exercises the kernel, the filter and the
  ## path sampler against each other.
  d <- small(); p <- draws(d, 1L)
  set.seed(4)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 4000L)
  an <- gc_tstar_pt(gc_nlat(p$Z[1, ], d$M), gc_alpha(p$U[1], d$M),
                    gc_lambda(p$NU[1], exp(p$Z[1, ])), d$T,
                    d$deme_data())[, seq_len(d$T), drop = FALSE]
  for (i in seq_len(d$n_demes)) {
    worst <- max(abs(cumsum(tt$pmf_tfirst[i, ]) - cumsum(an[i, ])))
    ## 4000 paths: Monte-Carlo error of the sampled route is ~1/sqrt(n)
    expect_lt(worst, 0.03)
  }
})

test_that("t* cannot postdate the first non-zero survey, in the paths too", {
  d <- small(); p <- draws(d)
  set.seed(5)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 400L)
  for (i in seq_along(d$demes)) {
    k <- d$cut[[i]]
    if (!is.na(k) && k < d$T)
      expect_equal(sum(tt$pmf_tfirst[i, (k + 1L):d$T]), 0, info = d$demes[i])
  }
})

test_that("the sampler is reproducible from the R seed", {
  d <- small(); p <- draws(d)
  set.seed(6); a <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 100L)
  set.seed(6); b <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 100L)
  expect_identical(a$pmf_test, b$pmf_test)
  set.seed(7); cc <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 100L)
  expect_false(identical(a$pmf_test, cc$pmf_test))
})

## ---- the joint ---------------------------------------------------------

test_that("the joint has the right support and the right margins", {
  d <- small(); p <- draws(d)
  set.seed(8)
  tt <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 600L,
                    joint = TRUE)
  J <- tt$joint
  expect_equal(dim(J), c(d$T, d$T, d$n_demes))
  expect_true(all(J >= 0))
  for (k in seq_len(d$n_demes)) {
    m <- J[, , k]
    ## the surviving excursion cannot begin before the variant appears, so
    ## nothing sits below the diagonal -- structurally, not approximately
    expect_equal(sum(m[lower.tri(m)]), 0, info = d$demes[k])
    ## margins are the two one-dimensional pmfs
    expect_equal(colSums(m), tt$pmf_test[k, ], tolerance = 1e-12)
    expect_equal(sum(m), sum(tt$pmf_test[k, ]), tolerance = 1e-12)
    ## the T_first margin is the pmf restricted to paths that were detected,
    ## so it is dominated by the unrestricted one
    expect_true(all(rowSums(m) <= tt$pmf_tfirst[k, ] + 1e-12))
  }
})

test_that("asking for the joint does not change the margins", {
  d <- small(); p <- draws(d)
  set.seed(9); a <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 200L)
  set.seed(9); b <- gc_times_pt(p$Z, p$U, p$NU, d$M, d$T, d$deme_data(), 200L,
                                joint = TRUE)
  expect_identical(a$pmf_tfirst, b$pmf_tfirst)
  expect_identical(a$pmf_test, b$pmf_test)
  expect_identical(a$p_alive, b$p_alive)
})

test_that("a deme detected at once keeps its mass on the diagonal", {
  ## positive at the very first survey, strongly selected, large: there is
  ## no room for an earlier arrival to have come and gone, so the lineage
  ## that was detected is the lineage that arrived
  dd <- GCData$new(data.frame(deme = "a", time = 2016:2020,
                              x = c(10L, 15L, 20L, 24L, 26L),
                              n = rep(30L, 5L)), M = 4L)
  set.seed(10)
  tt <- gc_times_pt(matrix(log(40), 1L, 1L), 1.0, 0.001, dd$M, dd$T,
                    dd$deme_data(), 400L, joint = TRUE)
  m <- tt$joint[, , 1L]
  expect_equal(sum(m), 1, tolerance = 1e-12)
  expect_gt(sum(diag(m)) / sum(m), 0.95)
})

test_that("a late detection leaves room for turnover before it", {
  ## never seen until the last survey: earlier arrivals have had years to
  ## come and go, so T_est should sit well after t*
  dd <- GCData$new(data.frame(deme = "a", time = 2016:2024,
                              x = c(rep(0L, 8L), 6L), n = rep(30L, 9L)),
                   M = 4L)
  set.seed(11)
  tt <- gc_times_pt(matrix(log(25), 1L, 1L), 0.3, 0.05, dd$M, dd$T,
                    dd$deme_data(), 600L)
  md <- function(p) which(cumsum(p) / sum(p) >= 0.5)[1L]
  expect_gt(md(tt$pmf_test[1L, ]), md(tt$pmf_tfirst[1L, ]))
})

## ---- integration with GCProjection -------------------------------------

fitted <- function() {
  d <- suppressMessages(gc_example_data(M = 16))
  GCFit$new(d, chains = 1L, iter = 120L, warmup = 40L, seed = 12L,
            quiet = TRUE)
}

test_that("n_paths = 0 leaves the T_est columns absent, not wrong", {
  pr <- fitted()$project(horizon = 1, n_draws = 6L, n_paths = 0L)
  expect_true(all(is.na(pr$curves$pmf_test)))
  expect_true(all(is.na(pr$curves$cdf_test)))
  expect_true(all(is.na(pr$demes$p_alive_T)))
  ## and everything else is unaffected
  expect_false(any(is.na(pr$curves$p_present)))
  expect_false(any(is.na(pr$curves$pmf_tfirst)))
})

test_that("n_paths > 0 fills them consistently", {
  fit <- fitted()
  set.seed(11)
  pr <- fit$project(horizon = 1, n_draws = 6L, n_paths = 120L)
  cv <- pr$curves
  expect_false(any(is.na(cv$pmf_test)))
  ## zero beyond the window: t* is supported on 1..T, so T_est is too
  expect_true(all(cv$pmf_test[cv$phase == "projected"] == 0))
  ## a proper distribution where the deme was detected, absent where it was
  ## not; the cdf is its running sum either way
  agg <- as.numeric(tapply(cv$pmf_test, cv$deme, sum)[pr$demes$deme])
  seen <- !is.na(pr$demes$cut_time)
  expect_equal(agg[seen], rep(1, sum(seen)), tolerance = 1e-12)
  expect_equal(agg[!seen], rep(0, sum(!seen)), tolerance = 1e-12)
  for (nm in pr$demes$deme) {
    s <- cv[cv$deme == nm, ]
    expect_equal(s$cdf_test, cumsum(s$pmf_test), tolerance = 1e-12)
  }
  ## T_est is at or after t* on every path, so its CDF cannot lead.  The
  ## two sides are computed differently here -- cdf_tfirst is the analytic
  ## t* posterior, cdf_test comes from sampled paths -- so the comparison
  ## carries the sampler's Monte-Carlo error, about 1/sqrt(draws * paths).
  ## The exact, pathwise version of this statement is tested on the joint,
  ## where no mass may sit below the diagonal.
  mc <- 3 / sqrt(6 * 120)
  expect_true(all(cv$cdf_test <= cv$cdf_tfirst + mc))
})

test_that("the projection cache distinguishes n_paths", {
  fit <- fitted()
  a <- fit$project(horizon = 1, n_draws = 6L, n_paths = 0L)
  b <- fit$project(horizon = 1, n_draws = 6L, n_paths = 50L)
  expect_false(identical(a, b))
  expect_identical(b, fit$project(horizon = 1, n_draws = 6L, n_paths = 50L))
})
