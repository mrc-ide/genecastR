toy <- function(...) {
  d <- data.frame(
    deme = rep(c("a", "b"), each = 3L),
    time = rep(2016:2018, 2L),
    x    = c(0L, 1L, 3L, 0L, 0L, 0L),
    n    = c(10L, 20L, 20L, 10L, 10L, 10L),
    stringsAsFactors = FALSE
  )
  GCData$new(d, ...)
}

test_that("a clean panel maps onto the lattice as specified", {
  d <- toy(M = 4)
  expect_equal(d$t0, 2015)
  expect_equal(d$M, 4L)
  expect_equal(d$surveys$step, rep(c(4L, 8L, 12L), 2L))
  expect_equal(d$T, 12L)
  expect_equal(d$span, 3)
  expect_equal(d$n_demes, 2L)
  expect_equal(d$n_surveys, 6L)
})

test_that("cut is the first positive step, NA where there is none", {
  d <- toy(M = 4)
  expect_equal(unname(d$cut), c(8L, NA_integer_))
  expect_named(d$cut, c("a", "b"))
})

test_that("a deme that is never positive is kept, not dropped", {
  d <- toy(M = 4)
  expect_true("b" %in% d$demes)
  expect_equal(sum(d$surveys$deme == "b"), 3L)
})

test_that("t0 must be strictly before the first survey", {
  expect_error(toy(M = 4, t0 = 2016), "strictly before")
  expect_error(toy(M = 4, t0 = 2017), "strictly before")
  expect_silent(toy(M = 4, t0 = 2015.75))
})

test_that("non-integral steps are an error, not a silent rounding", {
  expect_error(toy(M = 1, t0 = 2015.5), "not a whole number")
  ## the same t0 is fine once M resolves it
  expect_silent(toy(M = 2, t0 = 2015.5))
})

test_that("off-grid surveys error unless snapping is asked for", {
  d <- data.frame(deme = c("a", "a"), time = c(2016.1, 2017.4),
                  x = c(1L, 2L), n = c(10L, 20L))
  expect_error(GCData$new(d, M = 1, t0 = 2015), "not a whole number")
  snapped <- suppressMessages(GCData$new(d, M = 1, t0 = 2015, snap = TRUE))
  expect_equal(snapped$surveys$step, c(1L, 2L))
  expect_message(GCData$new(d, M = 1, t0 = 2015, snap = TRUE), "snapped")
  expect_silent(GCData$new(d, M = 1, t0 = 2015, snap = TRUE, quiet = TRUE))
})

test_that("snapping can collide, and collisions error before summing", {
  ## 2016.1 and 2016.4 both round to step 1 at M = 1
  d <- data.frame(deme = c("a", "a"), time = c(2016.1, 2016.4),
                  x = c(1L, 2L), n = c(10L, 20L))
  expect_error(
    suppressMessages(GCData$new(d, M = 1, t0 = 2015, snap = TRUE)),
    "more than one survey")
  summed <- suppressMessages(
    GCData$new(d, M = 1, t0 = 2015, snap = TRUE, on_collision = "sum"))
  expect_equal(nrow(summed$surveys), 1L)
  expect_equal(summed$surveys$x, 3)
  expect_equal(summed$surveys$n, 30)
  expect_equal(summed$surveys$time, 2016)
  ## at a resolution that separates them, nothing is summed
  fine <- suppressMessages(GCData$new(d, M = 20, t0 = 2015, snap = TRUE))
  expect_equal(nrow(fine$surveys), 2L)
})

test_that("snapping cannot push a survey off the front of the window", {
  d <- data.frame(deme = "a", time = 2015.2, x = 1L, n = 10L)
  expect_error(suppressMessages(GCData$new(d, M = 1, t0 = 2015, snap = TRUE)),
               "step < 1")
})

test_that("duplicate (deme, time) is always an error", {
  d <- data.frame(deme = c("a", "a"), time = c(2016, 2016),
                  x = c(1L, 2L), n = c(10L, 20L))
  expect_error(GCData$new(d, M = 4), "duplicate")
})

test_that("n == 0 rows are dropped with a message, not an error", {
  d <- data.frame(deme = c("a", "a", "a"), time = 2016:2018,
                  x = c(0L, 1L, 0L), n = c(0L, 20L, 10L))
  expect_message(GCData$new(d, M = 4), "n == 0")
  expect_equal(nrow(suppressMessages(GCData$new(d, M = 4))$surveys), 2L)
  expect_silent(GCData$new(d, M = 4, quiet = TRUE))
})

test_that("a deme whose every row has n == 0 is an error", {
  d <- data.frame(deme = c("a", "a", "b"), time = c(2016, 2017, 2016),
                  x = c(0L, 1L, 0L), n = c(10L, 20L, 0L))
  expect_error(suppressMessages(GCData$new(d, M = 4)), "no usable surveys")
})

test_that("field-level violations are caught", {
  base <- data.frame(deme = "a", time = 2016, x = 1L, n = 10L)
  bad <- base; bad$x <- 11L
  expect_error(GCData$new(bad, M = 4), "0 <= x <= n")
  bad <- base; bad$x <- 1.5
  expect_error(GCData$new(bad, M = 4), "whole numbers")
  bad <- base; bad$n <- NA_integer_
  expect_error(GCData$new(bad, M = 4), "NA")
  bad <- base; bad$time <- Inf
  expect_error(GCData$new(bad, M = 4), "finite|whole number")
})

test_that("missing columns are reported by name", {
  d <- data.frame(district = "a", year = 2016, x = 1L, n = 10L)
  expect_error(GCData$new(d, M = 4), "no column")
  expect_silent(GCData$new(d, M = 4,
                           cols = c(deme = "district", time = "year",
                                    x = "x", n = "n")))
})

test_that("Jeffreys limits bracket the observed frequency", {
  d <- toy(M = 4)
  s <- d$surveys
  expect_true(all(s$lo <= s$f + 1e-12))
  expect_true(all(s$hi >= s$f - 1e-12))
  expect_true(all(s$lo[s$x == 0] == 0))
  expect_true(all(s$lo >= 0 & s$hi <= 1))
})

test_that("deme_data() has one entry per deme, in deme order", {
  d <- toy(M = 4)
  dd <- d$deme_data()
  expect_length(dd, 2L)
  expect_equal(vapply(dd, function(e) e$name, character(1L)), d$demes)
  expect_type(dd[[1L]]$step, "integer")
  expect_equal(dd[[1L]]$cut, 8L)
  expect_true(is.na(dd[[2L]]$cut))
})

## ---- the shipped panel, against the reference Python pipeline -----------

test_that("the Ugandan panel reproduces the reference preparation", {
  d <- suppressMessages(gc_example_data(M = 16))
  expect_equal(d$n_demes, 16L)
  expect_equal(d$n_surveys, 96L)
  expect_equal(d$t0, 2015)
  expect_equal(d$T, 112L)
  expect_equal(d$span, 7)

  expect_equal(d$cut, c(
    Agago = 32L, Amolatar = 32L, Arua = 96L, Hoima = 80L,
    Jinja = NA_integer_, Kaabong = 48L, Kanungu = NA_integer_,
    Kapchorwa = 96L, Kasese = 112L, Katakwi = 48L, Koboko = 80L,
    Kole = 16L, Lamwo = 16L, Mubende = 32L, Rukiga = 112L, Tororo = 96L
  ))

  s <- d$summary()
  expect_equal(s$n_nonzero, c(6L, 1L, 2L, 1L, 0L, 4L, 0L, 2L, 1L, 5L,
                              3L, 5L, 7L, 2L, 1L, 2L))
  expect_equal(sum(s$n_surveys), 96L)
  ## every survey lands on a year boundary at M = 16
  expect_true(all(d$surveys$step %% 16L == 0L))
})

test_that("M is a resolution knob: it changes steps, not the panel", {
  a <- suppressMessages(gc_example_data(M = 16))
  b <- suppressMessages(gc_example_data(M = 32))
  expect_equal(a$n_surveys, b$n_surveys)
  expect_equal(a$demes, b$demes)
  expect_equal(a$span, b$span)
  expect_equal(b$T, 2L * a$T)
  expect_equal(unname(b$cut), unname(2L * a$cut))
  expect_equal(a$surveys$x, b$surveys$x)
})

## ---- dates become decimal years ------------------------------------------

test_that("a Date time column becomes a decimal year, not days since 1970", {
  mk <- function(tt, ...) GCData$new(data.frame(
    deme = rep(c("a", "b"), each = 3L), time = rep(tt, 2L),
    x = c(0L, 1L, 2L, 0L, 0L, 1L), n = rep(30L, 6L)), ...)
  d <- mk(as.Date(c("2016-01-01", "2017-01-01", "2018-01-01")), M = 4L,
          quiet = TRUE)
  expect_equal(sort(unique(d$surveys$time)), c(2016, 2017, 2018))
  ## as.numeric() would have given 16801, 17167, 17532
  expect_false(any(d$surveys$time > 1e4))
  expect_message(mk(as.Date(c("2016-01-01", "2017-01-01", "2018-01-01")),
                    M = 4L), "decimal years")
})

test_that("the fraction is measured against each year's own length", {
  one <- function(s, ...) GCData$new(
    data.frame(deme = "a", time = as.Date(s), x = 0L, n = 30L),
    M = 1L, snap = TRUE, quiet = TRUE)$surveys$time
  ## 2016 is a leap year, so July 1 is day 182 rather than 181
  expect_equal(one("2016-07-01"), 2016 + 182 / 366)
  expect_equal(one("2017-07-01"), 2017 + 181 / 365)
  expect_equal(one("2020-12-31"), 2020 + 365 / 366)
  ## and the elapsed time between two dates is right
  expect_equal(one("2017-01-01") - one("2016-01-01"), 1)
})

test_that("POSIXct is converted too", {
  d <- GCData$new(data.frame(deme = rep("a", 2L),
                             time = as.POSIXct(c("2016-01-01", "2017-01-01"),
                                               tz = "UTC"),
                             x = c(0L, 1L), n = c(30L, 30L)),
                  M = 2L, quiet = TRUE)
  expect_equal(d$surveys$time, c(2016, 2017))
})

test_that("dates that miss a step say so, and snap accepts them", {
  mk <- function(...) GCData$new(data.frame(
    deme = rep(c("a", "b"), each = 3L),
    time = rep(as.Date(c("2016-03-15", "2017-04-02", "2018-05-20")), 2L),
    x = c(0L, 1L, 2L, 0L, 0L, 1L), n = rep(30L, 6L)), ...)
  expect_error(mk(M = 4L, quiet = TRUE), "was a date")
  d <- mk(M = 12L, snap = TRUE, quiet = TRUE)
  expect_true(all(d$surveys$step == round(d$surveys$step)))
})
