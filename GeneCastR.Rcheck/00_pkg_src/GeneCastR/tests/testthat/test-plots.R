skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")

pj <- local({
  cache <- NULL
  function(n_paths = 0L) {
    key <- paste0("p", n_paths)
    if (is.null(cache[[key]])) {
      d <- suppressMessages(gc_example_data(M = 16))
      fit <- GCFit$new(d, chains = 2L, iter = 120L, warmup = 40L, seed = 5L,
                       quiet = TRUE)
      cache[[key]] <<- list(
        fit = fit,
        pr = fit$project(horizon = 1, n_draws = 6L, n_paths = n_paths))
    }
    cache[[key]]
  }
})

test_that("every plot function returns a ggplot", {
  skip_if_no_ggplot()
  z <- pj(60L)
  for (p in list(gc_plot_frequency(z$pr), gc_plot_present(z$pr),
                 gc_plot_exceed(z$pr), gc_plot_times(z$pr, "cdf"),
                 gc_plot_times(z$pr, "pmf"), gc_plot_joint(z$pr),
                 gc_plot_trace(z$fit)))
    expect_s3_class(p, "ggplot")
  expect_s3_class(plot(z$pr), "ggplot")
  expect_s3_class(plot(z$fit), "ggplot")
})

test_that("plots build without warnings", {
  skip_if_no_ggplot()
  z <- pj(60L)
  expect_silent(ggplot2::ggplot_build(gc_plot_frequency(z$pr)))
  expect_silent(ggplot2::ggplot_build(gc_plot_present(z$pr)))
  expect_silent(ggplot2::ggplot_build(gc_plot_times(z$pr, "cdf")))
})

test_that("deme subsetting selects exactly what was asked for", {
  skip_if_no_ggplot()
  z <- pj(0L)
  p <- gc_plot_frequency(z$pr, demes = c("Agago", "Jinja"))
  expect_setequal(as.character(unique(p$data$deme)), c("Agago", "Jinja"))
  expect_error(gc_plot_frequency(z$pr, demes = "Nowhere"), "no such deme")
  ## the survey layer is subset too, not just the curves
  sv <- p$layers[[length(p$layers) - 1L]]$data
  if (is.data.frame(sv) && "deme" %in% names(sv))
    expect_setequal(as.character(unique(sv$deme)), c("Agago", "Jinja"))
})

test_that("every plot puts a deme in the same panel position", {
  skip_if_no_ggplot()
  z <- pj(60L)
  ord <- lapply(list(gc_plot_frequency(z$pr), gc_plot_present(z$pr),
                     gc_plot_exceed(z$pr), gc_plot_times(z$pr, "cdf"),
                     gc_plot_times(z$pr, "pmf"), gc_plot_joint(z$pr)),
                function(p) levels(p$data$deme))
  expect_length(unique(ord), 1L)
  expect_equal(ord[[1L]], sort(z$pr$demes$deme))
  ## The chains line up with the panels too. Per-deme sizes exist only when
  ## they are fitted per deme, so this part needs the hierarchy: under the
  ## shared default `"N_eff"` is the one size, not sixteen copies of it.
  h <- GCFit$new(z$fit$data, chains = 1L, iter = 60L, warmup = 20L, seed = 5L,
                 control = gc_control(neff = "hierarchical"), quiet = TRUE)
  nm <- levels(gc_plot_trace(h, "N_eff")$data$variable)
  expect_equal(sub("^N_eff\\[(.*)\\]$", "\\1", nm), ord[[1L]])
  expect_equal(levels(gc_plot_hist(h, "N_eff")$data$variable), nm)
})

test_that("the panel-level size is named for the structure that was fitted", {
  skip_if_no_ggplot()
  z <- pj()
  ## shared: one size, and asking for "N_eff" gives that, not sixteen copies
  expect_equal(levels(gc_plot_trace(z$fit, "N_eff")$data$variable), "N_eff")
  expect_true("N_eff" %in% levels(gc_plot_trace(z$fit)$data$variable))
  expect_false("exp_mu" %in% levels(gc_plot_trace(z$fit)$data$variable))
  ## hierarchical: the panel level is a mean of logs, so it stays exp_mu
  h <- GCFit$new(z$fit$data, chains = 1L, iter = 60L, warmup = 20L, seed = 5L,
                 control = gc_control(neff = "hierarchical"), quiet = TRUE)
  expect_true("exp_mu" %in% levels(gc_plot_trace(h)$data$variable))
  expect_false("N_eff" %in% levels(gc_plot_trace(h)$data$variable))
})

test_that("numeric deme names sort numerically, not as text", {
  skip_if_no_ggplot()
  d <- GCData$new(data.frame(deme = rep(c("2", "10", "1"), each = 3L),
                             time = rep(2016:2018, 3L),
                             x = c(0L, 1L, 2L, 0L, 0L, 1L, 0L, 2L, 3L),
                             n = rep(30L, 9L)), M = 4L)
  fit <- GCFit$new(d, chains = 1L, iter = 60L, warmup = 20L, seed = 2L,
                   quiet = TRUE)
  pr <- fit$project(horizon = 0, n_draws = 4L)
  expect_equal(levels(gc_plot_frequency(pr)$data$deme), c("1", "2", "10"))
})

test_that("the sorted orders are still reachable", {
  skip_if_no_ggplot()
  z <- pj(60L)
  dm <- z$pr$demes
  expect_equal(levels(gc_plot_frequency(z$pr, order = "p_present")$data$deme),
               dm$deme[order(-dm$p_present_final)])
  expect_equal(levels(gc_plot_times(z$pr, "cdf", order = "p_alive")$data$deme),
               dm$deme[order(dm$p_alive_T, -dm$n_nonzero)])
})

test_that("T_est plots refuse politely when no paths were sampled", {
  skip_if_no_ggplot()
  z <- pj(0L)
  expect_error(gc_plot_joint(z$pr), "needs sampled paths")
  expect_error(gc_plot_times(z$pr, which = "test"), "needs sampled paths")
  ## but t* alone is always available, and the default degrades to it
  expect_s3_class(gc_plot_times(z$pr, which = "tfirst"), "ggplot")
  p <- gc_plot_times(z$pr)
  expect_setequal(as.character(unique(p$data$what)), "t* = T_first")
})

test_that("exceedance only offers thresholds that were computed", {
  skip_if_no_ggplot()
  z <- pj(0L)
  expect_error(gc_plot_exceed(z$pr, 0.5), "not computed")
  p <- gc_plot_exceed(z$pr, 0.05)
  expect_setequal(as.character(unique(p$data$thr)), "5%")
  bare <- suppressMessages(gc_example_data(M = 16))
  f2 <- GCFit$new(bare, chains = 1L, iter = 60L, warmup = 20L, seed = 6L,
                  quiet = TRUE)
  p2 <- f2$project(horizon = 0, n_draws = 3L, thresholds = numeric(0))
  expect_error(gc_plot_exceed(p2), "no exceedance thresholds")
})

test_that("the pmf scaling puts both times on one axis", {
  skip_if_no_ggplot()
  z <- pj(60L)
  p <- gc_plot_times(z$pr, "pmf")
  ## densities are per unit time.  Under the default geometric prior t* may
  ## postdate the record, so its density integrates to 1 - p_tstar_beyond
  ## rather than to one; the missing mass is not truncation, it is the
  ## probability the deme has not been reached.
  tf <- p$data[p$data$what == "t* = T_first", ]
  by <- as.numeric(tapply(tf$v, tf$deme,
                          function(v) sum(v) / z$pr$settings$M)[z$pr$demes$deme])
  expect_equal(by, 1 - z$pr$demes$p_tstar_beyond, tolerance = 1e-10)
  te <- p$data[p$data$what == "T_est", ]
  by <- as.numeric(tapply(te$v, te$deme,
                          function(v) sum(v) / z$pr$settings$M)[z$pr$demes$deme])
  seen <- !is.na(z$pr$demes$cut_time)
  expect_equal(by[seen], rep(1, sum(seen)), tolerance = 1e-10)
  expect_equal(by[!seen], rep(0, sum(!seen)), tolerance = 1e-10)
})

test_that("the joint plot carries only the supported region", {
  skip_if_no_ggplot()
  z <- pj(60L)
  p <- gc_plot_joint(z$pr)
  ## nothing below the diagonal reaches the plot
  expect_true(all(p$data$te >= p$data$tf - 1e-12))
  expect_true(all(p$data$p > 0))
})

test_that("trace plots carry the diagnostics in the panel labels", {
  skip_if_no_ggplot()
  z <- pj(0L)
  p <- gc_plot_trace(z$fit, vars = c("u", "nu"))
  expect_setequal(as.character(unique(p$data$variable)), c("u", "nu"))
  expect_equal(nlevels(p$data$chain), 2L)
  ## the labeller takes a data frame of the facet variable
  lab <- p$facet$params$labeller(
    data.frame(variable = factor(c("u", "nu"), levels = c("u", "nu"))))
  expect_true(all(grepl("Rhat", lab$variable)))
  expect_true(all(grepl("ESS", lab$variable)))
})
