fx <- function(f) {
  p <- testthat::test_path("fixtures", f)
  utils::read.csv(p, stringsAsFactors = FALSE)
}

## ---- structural properties, no fixture needed --------------------------

test_that("the kernel is row-stochastic and absorbing at both ends", {
  for (N in c(2L, 5L, 40L, 201L)) {
    for (alpha in c(-0.5, 0, 0.0217, 1.5)) {
      P <- gc_kernel(N, alpha)
      expect_equal(dim(P), c(N + 1L, N + 1L))
      expect_equal(rowSums(P), rep(1, N + 1L), tolerance = 1e-13)
      expect_true(all(P >= 0))
      expect_equal(P[1, 1], 1)                 # k = 0 is absorbing
      expect_equal(P[N + 1L, N + 1L], 1)       # k = N is absorbing
      expect_equal(sum(P[1, -1]), 0)
      expect_equal(sum(P[N + 1L, -(N + 1L)]), 0)
    }
  }
})

test_that("selection moves the kernel's mean the way psi says", {
  N <- 100L
  k <- 20L
  for (alpha in c(-0.3, 0, 0.3)) {
    P <- gc_kernel(N, alpha)
    f <- k / N
    expect_equal(sum(P[k + 1L, ] * (0:N)) / N,
                 f * (1 + alpha) / (1 + alpha * f), tolerance = 1e-12)
  }
})

test_that("alpha is rejected at and below -1", {
  expect_error(gc_kernel(10L, -1), "exceed -1")
  expect_error(gc_kernel(10L, -2), "exceed -1")
  expect_error(gc_kernel(0L, 0), "at least 1")
})

test_that("the Poisson row sums to exactly one, tail folded onto N", {
  for (N in c(2L, 5L, 40L, 201L)) {
    for (lam in c(0, 1e-8, 0.27, 5, 1000)) {
      r <- gc_poisson_row(N, lam)
      expect_length(r, N + 1L)
      expect_true(all(r >= 0))
      ## exactly one: a leak of 1e-12 compounds over T steps
      expect_equal(sum(r), 1, tolerance = 1e-15)
    }
  }
  ## at a mean far beyond the lattice, essentially all mass is folded onto N
  r <- gc_poisson_row(40L, 1000)
  expect_gt(r[41L], 0.999)
  ## lambda = 0 cannot import
  expect_equal(gc_poisson_row(40L, 0), c(1, rep(0, 40)))
})

test_that("the emission is a proper likelihood over the lattice", {
  e <- gc_emission(20L, 0L, 10L)
  expect_equal(e[1L], 1)                        # k = 0 explains x = 0 exactly
  expect_equal(e[21L], 0)                       # k = N cannot give x = 0
  e <- gc_emission(20L, 10L, 10L)
  expect_equal(e[1L], 0)
  expect_equal(e[21L], 1)
  ## summing over x at fixed k is 1
  s <- rowSums(vapply(0:12, function(x) gc_emission(20L, x, 12L),
                      numeric(21L)))
  expect_equal(s, rep(1, 21L), tolerance = 1e-13)
  expect_error(gc_emission(20L, 5L, 4L), "0 <= x <= n")
})

## ---- parity with the reference implementation --------------------------

test_that("the kernel matches the reference entry by entry", {
  d <- fx("kernel.csv")
  for (key in unique(paste(d$N, d$alpha))) {
    s <- d[paste(d$N, d$alpha) == key, ]
    N <- s$N[1L]; alpha <- s$alpha[1L]
    P <- gc_kernel(as.integer(N), alpha)
    expect_equal(P[cbind(s$i + 1L, s$j + 1L)], s$value,
                 tolerance = 1e-12,
                 info = paste("N =", N, "alpha =", alpha))
  }
})

test_that("the Poisson row matches the reference", {
  d <- fx("poisson_row.csv")
  for (key in unique(paste(d$N, d$lambda))) {
    s <- d[paste(d$N, d$lambda) == key, ]
    r <- gc_poisson_row(as.integer(s$N[1L]), s$lambda[1L])
    expect_equal(r[s$k + 1L], s$value, tolerance = 1e-12,
                 info = paste("N =", s$N[1L], "lambda =", s$lambda[1L]))
  }
})

test_that("the emission matches the reference", {
  d <- fx("emission.csv")
  for (key in unique(paste(d$N, d$x, d$n))) {
    s <- d[paste(d$N, d$x, d$n) == key, ]
    e <- gc_emission(as.integer(s$N[1L]), as.integer(s$x[1L]),
                     as.integer(s$n[1L]))
    expect_equal(e[s$k + 1L], s$value, tolerance = 1e-12, info = key)
  }
})

test_that("the integer lattice rounds the same way as the reference", {
  d <- fx("nlat.csv")
  expect_equal(gc_nlat(d$z, d$M), as.integer(d$Nlat))
  ## and the floor is real
  expect_equal(gc_nlat(log(1e-9), 16), 2L)
})

test_that("rounding is half-to-even, and the boundary is a real edge", {
  ## constructed so M * exp(z) is exactly a half-integer
  half <- function(v, M) gc_nlat(log(v / M), M)
  expect_equal(half(2.5, 32), 2L)    # to even, not up
  expect_equal(half(3.5, 32), 4L)
  expect_equal(half(4.5, 32), 4L)
  expect_equal(half(5.5, 32), 6L)

  ## The likelihood is piecewise constant in z, so the lattice is a step
  ## function and an ulp either side of a half-integer gives a different
  ## model.  Under a continuous z this has probability zero and cannot
  ## affect inference, but it does mean no fixture may sit on a boundary,
  ## and that bitwise agreement with another implementation is not a
  ## reasonable thing to ask for.
  z <- log(2.5 / 32)
  expect_equal(gc_nlat(z, 32), 2L)
  expect_equal(gc_nlat(z + 1e-12, 32), 3L)
  expect_equal(gc_nlat(z - 1e-12, 32), 2L)
})

## ---- the likelihood and the t* posterior, on the real panel ------------

setting_row <- function(d, nm = "median") {
  s <- fx("settings.csv")
  s <- s[s$setting == nm, ]
  s[match(d$demes, s$deme), ]
}

panel <- function() {
  d <- suppressMessages(gc_example_data(M = 16))
  list(data = d, dd = d$deme_data(),
       settings = fx("settings.csv"))
}

test_that("the per-deme log-likelihood matches the reference", {
  p <- panel()
  ref <- fx("loglik.csv")
  for (nm in unique(ref$setting)) {
    s <- p$settings[p$settings$setting == nm, ]
    s <- s[match(p$data$demes, s$deme), ]
    expect_equal(s$deme, p$data$demes)
    got <- gc_loglik_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda,
                        s$T[1L], p$dd)
    want <- ref$loglik[ref$setting == nm]
    want <- want[match(p$data$demes, ref$deme[ref$setting == nm])]
    expect_equal(got, want, tolerance = 1e-10, info = nm)
  }
})

test_that("the t* posterior matches the reference and is a distribution", {
  p <- panel()
  ref <- fx("tstar.csv")
  for (nm in unique(ref$setting)) {
    s <- p$settings[p$settings$setting == nm, ]
    s <- s[match(p$data$demes, s$deme), ]
    got <- gc_tstar_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda,
                       s$T[1L], p$dd)
    ## D x (T+1): the extra column is p(t* beyond the record), which the
    ## uniform prior fixes at zero
    expect_equal(ncol(got), s$T[1L] + 1L)
    expect_equal(unname(got[, s$T[1L] + 1L]), rep(0, nrow(got)))
    expect_equal(rowSums(got), rep(1, nrow(got)), tolerance = 1e-12)
    got <- got[, seq_len(s$T[1L]), drop = FALSE]
    r <- ref[ref$setting == nm, ]
    for (i in seq_along(p$data$demes)) {
      q <- r[r$deme == p$data$demes[i], ]
      q <- q[order(q$t), ]
      expect_equal(got[i, ], q$prob, tolerance = 1e-10,
                   info = paste(nm, p$data$demes[i]))
    }
  }
})

test_that("t* cannot postdate the first non-zero survey", {
  p <- panel()
  s <- p$settings[p$settings$setting == "median", ]
  s <- s[match(p$data$demes, s$deme), ]
  post <- gc_tstar_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda, s$T[1L],
                      p$dd)[, seq_len(s$T[1L]), drop = FALSE]
  cut <- p$data$cut
  for (i in seq_along(cut)) {
    if (!is.na(cut[i]) && cut[i] < s$T[1L]) {
      expect_equal(sum(post[i, (cut[i] + 1L):s$T[1L]]), 0,
                   info = p$data$demes[i])
      expect_gt(sum(post[i, 1:cut[i]]), 1 - 1e-12)
    } else {
      ## never positive: no constraint, so mass reaches the last step
      expect_gt(post[i, s$T[1L]], 0)
    }
  }
})

test_that("a larger influx raises the likelihood of a late-appearing deme", {
  p <- panel()
  s <- p$settings[p$settings$setting == "median", ]
  s <- s[match(p$data$demes, s$deme), ]
  lo <- gc_loglik_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda * 0.1,
                     s$T[1L], p$dd)
  hi <- gc_loglik_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda * 3,
                     s$T[1L], p$dd)
  late <- which(p$data$demes == "Kasese")     # one non-zero, at the very end
  expect_gt(hi[late], lo[late])
})

test_that("zero influx is impossible, not merely unlikely", {
  p <- panel()
  s <- p$settings[p$settings$setting == "median", ]
  s <- s[match(p$data$demes, s$deme), ]
  ll <- gc_loglik_pt(as.integer(s$Nlat), s$alpha[1L],
                     rep(0, length(s$lambda)), s$T[1L], p$dd)
  ## with no imports an empty deme can still be founded at t*, so the
  ## likelihood is finite -- but every re-entry is shut off
  expect_true(all(is.finite(ll)))
  base <- gc_loglik_pt(as.integer(s$Nlat), s$alpha[1L], s$lambda, s$T[1L],
                       p$dd)
  expect_lt(sum(ll), sum(base))
})


## ---- t* beyond the observation record -----------------------------------

test_that("a geometric t* puts mass beyond the record, a uniform one cannot", {
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting_row(d)
  Nl <- as.integer(s$Nlat); al <- s$alpha[1L]; lam <- s$lambda
  u <- gc_tstar_pt(Nl, al, lam, 112L, d$deme_data())
  expect_true(all(u[, 113L] == 0))

  ## a real hazard: only demes never seen positive can be un-founded, since
  ## a positive count is proof of arrival
  g <- gc_tstar_pt(Nl, al, lam, 112L, d$deme_data(), 2L,
                   -expm1(-0.3 / 16))
  seen <- !is.na(d$cut)
  expect_true(all(g[seen, 113L] == 0))
  expect_true(all(g[!seen, 113L] > 0))
  expect_equal(rowSums(g), rep(1, nrow(g)), tolerance = 1e-12)

  ## a higher hazard leaves less mass un-founded
  g2 <- gc_tstar_pt(Nl, al, lam, 112L, d$deme_data(), 2L,
                    -expm1(-1.0 / 16))
  expect_true(all(g2[!seen, 113L] < g[!seen, 113L]))
})

test_that("the geometric nests the uniform as the hazard vanishes", {
  d <- suppressMessages(gc_example_data(M = 16))
  s <- setting_row(d)
  Nl <- as.integer(s$Nlat); al <- s$alpha[1L]; lam <- s$lambda
  ## both arms must use the SAME founding law, or the comparison is of two
  ## different models rather than of two t* priors
  u <- gc_tstar_pt(Nl, al, lam, 112L, d$deme_data(), 2L)[, 1:112]
  g <- gc_tstar_pt(Nl, al, lam, 112L, d$deme_data(), 2L, -expm1(-1e-6 / 16))
  ## conditional on arrival inside the window the two agree; unconditionally
  ## they cannot, because at a vanishing hazard a deme never seen positive
  ## puts nearly all its mass beyond the record
  gc <- g[, 1:112] / rowSums(g[, 1:112])
  expect_equal(gc, u, tolerance = 1e-6)
  expect_gt(min(g[is.na(d$cut), 113L]), 0.99)
})
