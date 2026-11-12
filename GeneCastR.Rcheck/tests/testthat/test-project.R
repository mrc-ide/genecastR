fx <- function(f) utils::read.csv(testthat::test_path("fixtures", f),
                                  stringsAsFactors = FALSE)

setting <- function(d, nm = "median") {
  s <- fx("settings.csv")
  s <- s[s$setting == nm, ]
  s[match(d$demes, s$deme), ]
}

## ---- parity with the reference smoother --------------------------------

test_that("the smoother matches the reference on the real panel", {
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting(d)
  ## z, not a z rebuilt from the rounded N_lat: lambda = nu * N_eff is built
  ## from the continuous effective size, and the two differ by up to half a
  ## lattice point
  expect_equal(gc_nlat(s$z, 16), as.integer(s$Nlat))

  r <- gc_project_pt(matrix(s$z, nrow = 1L), s$u[1L], s$nu[1L],
                     M = 16L, T = 112L, H = 48L, demes = d$deme_data(),
                     thresholds = c(0.01, 0.05, 0.10), nbins = 250L)
  ref <- fx("smooth.csv")
  for (i in seq_along(d$demes)) {
    q <- ref[ref$deme == d$demes[i], ]
    q <- q[order(q$step), ]
    expect_equal(as.numeric(r$pres[i, ]), q$pres, tolerance = 1e-10,
                 info = d$demes[i])
    expect_equal(as.numeric(r$fmean[i, ]), q$fmean, tolerance = 1e-10,
                 info = d$demes[i])
  }
})

test_that("exceedance matches the reference, computed on the lattice", {
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting(d)
  thr <- c(0.01, 0.05, 0.10)
  r <- gc_project_pt(matrix(s$z, nrow = 1L), s$u[1L], s$nu[1L], 16L, 112L,
                     48L, d$deme_data(), thr, 250L)
  ref <- fx("exceed.csv")
  for (i in seq_along(d$demes)) for (j in seq_along(thr)) {
    q <- ref[ref$deme == d$demes[i] & abs(ref$threshold - thr[j]) < 1e-12, ]
    q <- q[order(q$step), ]
    expect_equal(as.numeric(r$exceed[j, , i]), q$prob, tolerance = 1e-10,
                 info = paste(d$demes[i], thr[j]))
  }
})

test_that("exceedance at zero is exactly the presence probability", {
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting(d)
  r <- gc_project_pt(matrix(s$z, nrow = 1L), s$u[1L], s$nu[1L], 16L, 112L,
                     12L, d$deme_data(), c(0, 0.05), 250L)
  expect_equal(as.numeric(r$exceed[1L, , ]), as.numeric(t(r$pres)),
               tolerance = 1e-12)
  ## and it decreases in the threshold
  expect_true(all(r$exceed[1L, , ] >= r$exceed[2L, , ] - 1e-15))
})

test_that("the smoother uses later surveys, so it is not a filter", {
  ## two demes identical up to step 32; one is then seen at 10/30, the other
  ## stays at zero.  A filter would give them the same state at step 16.
  mk <- function(late) {
    GCData$new(data.frame(deme = "a", time = 2016:2019,
                          x = c(0L, 0L, late, late), n = rep(30L, 4L)),
               M = 4L)
  }
  s <- function(dd) {
    gc_project_pt(matrix(log(12), nrow = 1L), 0.3, 0.02, 4L, dd$T, 0L,
                  dd$deme_data(), numeric(0), 50L)$pres[1L, 1L]
  }
  expect_gt(s(mk(10L)), s(mk(0L)))
})

test_that("the projection is the kernel applied to the smoothed state", {
  ## With T = 1 the introduction is certain at step 1 and deposits one copy,
  ## so the smoothed lattice distribution there is exactly a point mass at
  ## k = 1.  Every projected step is then a row of a power of the kernel,
  ## which can be built independently from the fixture-checked pieces.
  N <- 40L; alpha <- 0.02; lambda <- 0.3; H <- 12L
  dd <- list(list(step = 1L, x = 0L, n = 1L, cut = NA_integer_))
  r <- gc_project_pt(matrix(log(N), nrow = 1L), alpha * 0, 0, M = 1L, T = 1L,
                     H = H, demes = dd, thresholds = numeric(0), nbins = 50L)

  ## rebuild with the same parameters the C++ derives: alpha from u, lambda
  ## from nu * N_eff
  u <- 1 * log1p(alpha); nu <- lambda / N
  r <- gc_project_pt(matrix(log(N), nrow = 1L), u, nu, 1L, 1L, H, dd,
                     numeric(0), 50L)
  P <- gc_kernel(N, alpha)
  P[1L, ] <- gc_poisson_row(N, lambda)
  p <- c(0, 1, rep(0, N - 1L))
  for (h in seq_len(H)) {
    p <- as.numeric(p %*% P)
    expect_equal(r$pres[1L, 1L + h], sum(p[-1L]), tolerance = 1e-10,
                 info = paste("step", h))
    expect_equal(r$fmean[1L, 1L + h], sum(p[-1L] * (seq_len(N) / N)),
                 tolerance = 1e-10, info = paste("step", h))
  }
  expect_equal(r$pres[1L, 1L], 1, tolerance = 1e-12)
})

test_that("p_present steps at the join, and in the direction the model says", {
  ## Not a defect: the final survey constrains the smoother and not the
  ## projection, and immigration lifts whatever mass sits at zero.
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting(d)
  r <- gc_project_pt(matrix(s$z, nrow = 1L), s$u[1L], s$nu[1L], 16L, 112L,
                     4L, d$deme_data(), numeric(0), 250L)
  jump <- r$pres[, 113L] - r$pres[, 112L]
  pinned <- r$pres[, 112L] > 1 - 1e-6
  ## pinned by a positive final survey: can only fall
  expect_true(all(jump[pinned] <= 1e-9))
  ## carrying mass at zero: must rise, by at most the one-step immigration
  ## probability applied to that mass
  lam <- s$nu[1L] * exp(s$z)
  room <- (1 - r$pres[, 112L]) * (1 - exp(-lam))
  expect_true(all(jump[!pinned] > 0))
  expect_true(all(jump[!pinned] <= room[!pinned] + 1e-9))
})

## ---- the assembled tables ----------------------------------------------

fitted <- function() {
  d <- suppressMessages(gc_example_data(M = 16))
  GCFit$new(d, chains = 1L, iter = 120L, warmup = 40L, seed = 9L,
            quiet = TRUE)
}

test_that("curves has the documented schema", {
  pr <- fitted()$project(horizon = 1, n_draws = 8L,
                         thresholds = c(0.01, 0.05), probs = c(.025, .5, .975))
  nm <- names(pr$curves)
  expect_true(all(c("deme", "step", "time", "phase", "f_mean",
                    "f_q2.5", "f_q50", "f_q97.5",
                    "p_present", "p_present_q2.5", "p_present_q50",
                    "p_present_q97.5", "p_exceed_1", "p_exceed_5",
                    "pmf_tfirst", "cdf_tfirst", "pmf_test", "cdf_test")
                  %in% nm))
  expect_false(any(grepl(" ", nm)))
  expect_equal(nrow(pr$curves), 16L * (112L + 16L))
  expect_equal(levels(pr$curves$phase), c("smoothed", "projected"))
  expect_equal(sum(pr$curves$phase == "projected"), 16L * 16L)
  expect_equal(pr$curves$time, pr$settings$t0 + pr$curves$step / 16)
  expect_equal(pr$thresholds, c(0.01, 0.05))
})

test_that("the frequency quantiles are ordered and unconditional", {
  pr <- fitted()$project(horizon = 1, n_draws = 8L)
  cv <- pr$curves
  expect_true(all(cv$f_q2.5 <= cv$f_q25 + 1e-12))
  expect_true(all(cv$f_q25 <= cv$f_q50 + 1e-12))
  expect_true(all(cv$f_q50 <= cv$f_q75 + 1e-12))
  expect_true(all(cv$f_q75 <= cv$f_q97.5 + 1e-12))
  expect_true(all(cv$f_mean >= 0 & cv$f_mean <= 1))
  ## the atom at zero carries 1 - p_present, so any level below it is zero
  low <- cv$p_present < 0.975
  expect_true(all(cv$f_q2.5[low] == 0))
})

test_that("p_present and its quantiles are different objects", {
  pr <- fitted()$project(horizon = 1, n_draws = 20L)
  cv <- pr$curves
  expect_true(all(cv$p_present >= 0 & cv$p_present <= 1))
  expect_true(all(cv$p_present_q2.5 <= cv$p_present_q50 + 1e-12))
  expect_true(all(cv$p_present_q50 <= cv$p_present_q97.5 + 1e-12))
  ## the mean of P(present | X, theta) over draws is P(present | X), so it
  ## lies inside the across-draw interval
  expect_true(all(cv$p_present >= cv$p_present_q2.5 - 1e-9))
  expect_true(all(cv$p_present <= cv$p_present_q97.5 + 1e-9))
})

test_that("pmf_tfirst is the analytic t* posterior", {
  fit <- fitted()
  pr <- fit$project(horizon = 1, n_draws = 10L)
  ts <- fit$tstar(n_draws = 10L)
  cv <- pr$curves[pr$curves$phase == "smoothed", ]
  key <- paste(cv$deme, cv$step)
  expect_equal(cv$pmf_tfirst, ts$prob[match(key, paste(ts$deme, ts$step))],
               tolerance = 1e-12)
  expect_equal(sum(cv$pmf_tfirst), sum(ts$prob), tolerance = 1e-10)
  ## the in-record mass plus the probability of not having been reached
  ## sum to one; under a uniform t* the second term is zero
  agg <- as.numeric(tapply(pr$curves$pmf_tfirst, pr$curves$deme,
                           sum)[pr$demes$deme])
  expect_equal(agg + pr$demes$p_tstar_beyond, rep(1, 16L), tolerance = 1e-12)
  ## the fit's own t* method does the same accounting
  expect_equal(unname(attr(ts, "p_beyond")[pr$demes$deme]) +
                 as.numeric(tapply(ts$prob, ts$deme, sum)[pr$demes$deme]),
               rep(1, 16L), tolerance = 1e-10)
  expect_true(all(pr$curves$pmf_tfirst[pr$curves$phase == "projected"] == 0))
  ## the cdf is its running sum within each deme
  for (nm in unique(pr$curves$deme)) {
    s <- pr$curves[pr$curves$deme == nm, ]
    expect_equal(s$cdf_tfirst, cumsum(s$pmf_tfirst), tolerance = 1e-12)
  }
})

test_that("t* cannot postdate the first non-zero survey", {
  fit <- fitted()
  pr <- fit$project(horizon = 1, n_draws = 10L)
  cut <- fit$data$cut
  for (nm in names(cut)) {
    if (!is.na(cut[[nm]]) && cut[[nm]] < fit$data$T) {
      late <- pr$curves$pmf_tfirst[pr$curves$deme == nm &
                                   pr$curves$step > cut[[nm]]]
      expect_equal(sum(late), 0, info = nm)
    }
  }
})

test_that("T_est columns are absent rather than wrong", {
  pr <- fitted()$project(horizon = 1, n_draws = 8L)
  expect_true(all(is.na(pr$curves$pmf_test)))
  expect_true(all(is.na(pr$curves$cdf_test)))
  expect_true(all(is.na(pr$demes$p_alive_T)))
})

test_that("the per-deme table is at its own grain", {
  fit <- fitted()
  pr <- fit$project(horizon = 1, n_draws = 8L)
  expect_equal(nrow(pr$demes), 16L)
  expect_equal(pr$demes$deme, fit$data$demes)
  expect_true(all(c("n_surveys", "n_nonzero", "cut_time", "N_eff_median",
                    "p_present_final", "f_median_final") %in% names(pr$demes)))
  last <- pr$curves[pr$curves$step == max(pr$curves$step), ]
  expect_equal(pr$demes$p_present_final,
               last$p_present[match(pr$demes$deme, last$deme)])
  expect_equal(nrow(pr$surveys), fit$data$n_surveys)
})

test_that("project() memoises on its arguments", {
  fit <- fitted()
  a <- fit$project(horizon = 1, n_draws = 8L)
  b <- fit$project(horizon = 1, n_draws = 8L)
  expect_identical(a, b)
  cdiff <- fit$project(horizon = 2, n_draws = 8L)
  expect_false(identical(a, cdiff))
  expect_equal(nrow(cdiff$curves), 16L * (112L + 32L))
  fit$clear_cache()
  expect_false(identical(a, fit$project(horizon = 1, n_draws = 8L)))
})

test_that("a zero horizon is the smoother alone", {
  pr <- fitted()$project(horizon = 0, n_draws = 6L)
  expect_equal(nrow(pr$curves), 16L * 112L)
  expect_true(all(pr$curves$phase == "smoothed"))
})

test_that("bad projection arguments are rejected", {
  fit <- fitted()
  expect_error(fit$project(horizon = -1), "non-negative")
  expect_error(fit$project(thresholds = 1.5))
  expect_error(fit$project(probs = c(0, 0.5)))
})
