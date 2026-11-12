#' Survey panel for the Wright-Fisher introduction model
#'
#' Holds a panel of allele-frequency surveys, validated and mapped onto the
#' lattice the model runs on. Everything downstream -- the likelihood, the
#' smoother, the projections -- is written against this object, so the
#' validation here is the package's contract with its input.
#'
#' @section Time and the lattice:
#' `M` is the number of lattice steps per unit of the `time` column -- not
#' per survey interval, and not per year unless `time` is in years. A survey
#' at `time` sits at step `M * (time - t0)`, which must be a positive whole
#' number, so `M` has to be a multiple of whatever divides the reporting
#' interval: quarterly surveys in year units need `M` divisible by four.
#' Nothing requires the surveys to be equally spaced, or the demes to share
#' a calendar.
#' `t0` must be strictly earlier than the first survey, because the
#' introduction time is supported on steps `1..T` and step 0 is not reachable.
#' The default `t0 = min(time) - 1` gives the introduction a year of support
#' before anything is observed.
#'
#' @section What is rejected, and why:
#' \describe{
#'   \item{`n == 0`}{Dropped with a message rather than an error: survey
#'     extracts routinely carry rows for a district-year that was not
#'     sampled, and they carry no information.}
#'   \item{non-integral steps}{An error. Annual data at `M = 16` is fine;
#'     sub-annual data at small `M` is not, and silently rounding it would
#'     move observations in time.}
#'   \item{two surveys on one step}{Reachable only with `snap = TRUE`, and
#'     then an error by default. Pooling them is often right, but it is a
#'     decision, so it needs `on_collision = "sum"`.}
#'   \item{duplicate `(deme, time)`}{Always an error. This is a defect in the
#'     input, not a resolution question.}
#'   \item{a deme with no usable surveys}{An error, naming it.}
#' }
#' A deme whose surveys are all zero is kept: it is the informative case for
#' the question this package exists to answer.
#'
#' @examples
#' d <- gc_example_data()
#' d
#' head(d$surveys)
#' d$cut
#'
#' @export
GCData <- R6::R6Class(
  "GCData",
  public = list(

    #' @field surveys Validated surveys, one row per observation, with the
    #'   lattice `step`, the observed frequency `f`, and Jeffreys 90% limits
    #'   `lo` and `hi` precomputed for plotting.
    surveys = NULL,

    #' @field demes Deme names, sorted.
    demes = NULL,

    #' @field M Lattice steps per unit time.
    M = NULL,

    #' @field t0 Origin of the model window; step 1 is `t0 + 1/M`.
    t0 = NULL,

    #' @field T Number of lattice steps in the observation window.
    T = NULL,

    #' @field cut Named vector, one entry per deme: the first step carrying a
    #'   non-zero count, or `NA` where the deme was never positive. The
    #'   introduction time is constrained to precede this.
    cut = NULL,

    #' @description Build a panel.
    #' @param surveys A data frame, one row per survey.
    #' @param M Lattice steps per unit time. The published fits use 16.
    #' @param t0 Window origin; defaults to `min(time) - 1`.
    #' @param cols Named character vector mapping the four required fields
    #'   (`deme`, `time`, `x`, `n`) onto columns of `surveys`.
    #' @param snap `FALSE` (default) requires every survey to land exactly on
    #'   a lattice step. `TRUE` rounds each survey to the nearest step and
    #'   reports how far the worst one moved -- the option for dated data,
    #'   where exact alignment is not achievable at any practical `M`.
    #' @param on_collision `"error"` (default) or `"sum"`, for two surveys
    #'   landing on the same lattice step. Reachable only when `snap = TRUE`.
    #' @param quiet Suppress the messages about dropped and moved rows.
    initialize = function(surveys,
                          M = 16,
                          t0 = NULL,
                          cols = c(deme = "deme", time = "time",
                                   x = "x", n = "n"),
                          snap = FALSE,
                          on_collision = c("error", "sum"),
                          quiet = FALSE) {

      on_collision <- match.arg(on_collision)
      stopifnot(isTRUE(snap) || isFALSE(snap))
      stopifnot(is.data.frame(surveys), nrow(surveys) > 0L)

      if (!(length(M) == 1L && is.finite(M) && M >= 1 && M == round(M))) {
        stop("`M` must be a single whole number >= 1.", call. = FALSE)
      }
      M <- as.integer(M)

      need <- c("deme", "time", "x", "n")
      if (!all(need %in% names(cols))) {
        stop("`cols` must name all of: ", paste(need, collapse = ", "),
             call. = FALSE)
      }
      missing_cols <- setdiff(unname(cols[need]), names(surveys))
      if (length(missing_cols)) {
        stop("`surveys` has no column(s): ",
             paste(missing_cols, collapse = ", "),
             ".\nAvailable: ", paste(names(surveys), collapse = ", "),
             call. = FALSE)
      }

      d <- data.frame(
        deme = as.character(surveys[[cols[["deme"]]]]),
        time = as.numeric(surveys[[cols[["time"]]]]),
        x    = as.numeric(surveys[[cols[["x"]]]]),
        n    = as.numeric(surveys[[cols[["n"]]]]),
        stringsAsFactors = FALSE
      )

      ## --- field-level checks -------------------------------------------
      if (anyNA(d)) {
        stop("`surveys` contains NA in deme, time, x or n (rows ",
             private$show_rows(which(!stats::complete.cases(d))), ").",
             call. = FALSE)
      }
      if (any(!is.finite(d$time))) {
        stop("`time` must be finite.", call. = FALSE)
      }
      private$check_integer(d$x, "x")
      private$check_integer(d$n, "n")
      if (any(d$n < 0)) {
        stop("`n` must be non-negative (rows ",
             private$show_rows(which(d$n < 0)), ").", call. = FALSE)
      }
      if (any(d$x < 0 | d$x > d$n)) {
        bad <- which(d$x < 0 | d$x > d$n)
        stop("`x` must satisfy 0 <= x <= n (rows ",
             private$show_rows(bad), ").", call. = FALSE)
      }

      ## --- duplicates, before anything is dropped ------------------------
      key <- paste(d$deme, d$time, sep = "\r")
      if (anyDuplicated(key)) {
        dup <- unique(key[duplicated(key)])
        stop("duplicate (deme, time) rows: ",
             paste(utils::head(gsub("\r", " @ ", dup, fixed = TRUE), 5L),
                   collapse = "; "),
             if (length(dup) > 5L) paste0(" ... and ", length(dup) - 5L,
                                          " more") else "",
             ".\nAggregate them yourself, or remove the duplicates.",
             call. = FALSE)
      }

      ## --- drop uninformative rows ---------------------------------------
      demes_in <- sort(unique(d$deme))
      keep <- d$n > 0
      if (any(!keep)) {
        if (!quiet) {
          message(sum(!keep), " row(s) with n == 0 dropped (no information).")
        }
        d <- d[keep, , drop = FALSE]
      }
      if (nrow(d) == 0L) {
        stop("no surveys left after dropping n == 0.", call. = FALSE)
      }

      ## --- the window ----------------------------------------------------
      if (is.null(t0)) {
        t0 <- min(d$time) - 1
      }
      if (!(length(t0) == 1L && is.finite(t0))) {
        stop("`t0` must be a single finite number.", call. = FALSE)
      }
      if (t0 >= min(d$time)) {
        stop("`t0` (", t0, ") must be strictly before the first survey (",
             min(d$time), "): the introduction time is supported on steps ",
             "1..T, so a survey at step 0 cannot be represented.",
             call. = FALSE)
      }

      step_raw <- M * (d$time - t0)
      off <- abs(step_raw - round(step_raw))
      if (any(off > 1e-8)) {
        bad <- which(off > 1e-8)
        if (!snap) {
          stop("M * (time - t0) is not a whole number for ", length(bad),
               " row(s); e.g. time = ",
               paste(utils::head(unique(d$time[bad]), 4L), collapse = ", "),
               " with M = ", M, " and t0 = ", t0,
               ".\nRaise M, or move t0, so every survey lands on a step -- ",
               "or pass snap = TRUE to round to the nearest step.",
               call. = FALSE)
        }
        if (!quiet) {
          message(length(bad), " survey(s) snapped to the nearest step; ",
                  "the largest move was ",
                  format(max(off) / M, digits = 3), " time units ",
                  "(one step is ", format(1 / M, digits = 3), ").")
        }
      }
      d$step <- as.integer(round(step_raw))
      if (any(d$step < 1L)) {
        stop("survey(s) at step < 1 after rounding; `t0` must leave at least ",
             "one step before the first survey.", call. = FALSE)
      }

      ## --- step collisions ------------------------------------------------
      skey <- paste(d$deme, d$step, sep = "\r")
      if (anyDuplicated(skey)) {
        dup <- unique(skey[duplicated(skey)])
        if (on_collision == "error") {
          stop(length(dup), " lattice step(s) carry more than one survey at ",
               "M = ", M, ", e.g. ",
               paste(utils::head(gsub("\r", " @ step ", dup, fixed = TRUE),
                                 4L), collapse = "; "),
               ".\nRaise M to separate them, or pass on_collision = \"sum\" ",
               "to sum x and n.", call. = FALSE)
        }
        if (!quiet) {
          message(length(dup), " lattice step(s) carried more than one ",
                  "survey; x and n shared.")
        }
        agg <- stats::aggregate(cbind(x, n) ~ deme + step, data = d, FUN = sum)
        agg$time <- t0 + agg$step / M
        d <- agg[, c("deme", "time", "x", "n", "step")]
      }

      ## --- order, derive --------------------------------------------------
      d <- d[order(d$deme, d$step), , drop = FALSE]
      rownames(d) <- NULL
      d$f <- d$x / d$n
      d$lo <- ifelse(d$x == 0, 0,
                     stats::qbeta(0.05, d$x + 0.5, d$n - d$x + 0.5))
      d$hi <- stats::qbeta(0.95, d$x + 0.5, d$n - d$x + 0.5)

      demes <- sort(unique(d$deme))
      empty <- setdiff(demes_in, demes)
      if (length(empty)) {
        stop("deme(s) with no usable surveys (every row had n == 0): ",
             paste(empty, collapse = ", "),
             ".\nDrop them from `surveys`, or supply their counts.",
             call. = FALSE)
      }

      ## T is a whole number of time units, so a horizon in years is exact.
      self$T <- as.integer(M * ceiling(max(d$step) / M))
      self$surveys <- d
      self$demes <- demes
      self$M <- M
      self$t0 <- t0

      cut <- vapply(demes, function(dm) {
        s <- d$step[d$deme == dm & d$x > 0]
        if (length(s)) min(s) else NA_integer_
      }, integer(1L))
      names(cut) <- demes
      self$cut <- cut

      invisible(self)
    },

    #' @description Per-deme view, in the shape the compiled code consumes.
    #' @return A list with one element per deme, each a list of integer
    #'   vectors `step`, `x`, `n`, the scalar `cut` (`NA` where the deme was
    #'   never positive) and the deme's `name`.
    deme_data = function() {
      lapply(self$demes, function(dm) {
        r <- self$surveys[self$surveys$deme == dm, , drop = FALSE]
        list(name = dm,
             step = as.integer(r$step),
             x    = as.integer(r$x),
             n    = as.integer(r$n),
             cut  = unname(self$cut[[dm]]))
      })
    },

    #' @description One row per deme: survey counts and the first positive step.
    summary = function() {
      s <- self$surveys
      data.frame(
        deme       = self$demes,
        n_surveys  = as.integer(table(factor(s$deme, levels = self$demes))),
        n_nonzero  = as.integer(table(factor(s$deme[s$x > 0],
                                             levels = self$demes))),
        cut_step   = unname(self$cut),
        cut_time   = unname(self$t0 + self$cut / self$M),
        stringsAsFactors = FALSE,
        row.names  = NULL
      )
    },

    #' @description Print a compact description.
    #' @param ... Ignored.
    print = function(...) {
      cat("<GCData>\n")
      cat("  demes    : ", length(self$demes), " (",
          paste(utils::head(self$demes, 3L), collapse = ", "),
          if (length(self$demes) > 3L) ", ..." else "", ")\n", sep = "")
      cat("  surveys  : ", nrow(self$surveys), " (",
          sum(self$surveys$x > 0), " non-zero)\n", sep = "")
      cat("  window   : ", format(self$t0), " to ",
          format(self$t0 + self$T / self$M),
          "   M = ", self$M, ", T = ", self$T, " steps\n", sep = "")
      never <- sum(is.na(self$cut))
      if (never) {
        cat("  never positive: ", never, " deme(s) -- ",
            paste(self$demes[is.na(self$cut)], collapse = ", "), "\n", sep = "")
      }
      invisible(self)
    }
  ),

  active = list(
    #' @field n_demes Number of demes.
    n_demes = function(value) {
      if (!missing(value)) stop("read-only", call. = FALSE)
      length(self$demes)
    },
    #' @field n_surveys Number of surveys retained.
    n_surveys = function(value) {
      if (!missing(value)) stop("read-only", call. = FALSE)
      nrow(self$surveys)
    },
    #' @field span Length of the observation window, measured in units of
    #'   the `time` column. `M` divides that same unit, so `span = T / M`
    #'   whatever the column holds -- years, months, or anything else. The
    #'   model never assumes years.
    span = function(value) {
      if (!missing(value)) stop("read-only", call. = FALSE)
      self$T / self$M
    }
  ),

  private = list(
    check_integer = function(v, nm) {
      if (any(abs(v - round(v)) > 1e-8)) {
        bad <- which(abs(v - round(v)) > 1e-8)
        stop("`", nm, "` must be whole numbers (rows ",
             private$show_rows(bad), ").", call. = FALSE)
      }
      invisible(TRUE)
    },
    show_rows = function(i) {
      paste0(paste(utils::head(i, 5L), collapse = ", "),
             if (length(i) > 5L) paste0(", ... and ", length(i) - 5L,
                                        " more") else "")
    }
  )
)


#' The Ugandan C469Y panel
#'
#' Sixteen Ugandan districts surveyed annually from 2016 to 2022 for the
#' PfKelch13 C469Y artemisinin-resistance marker. The raw extract carries
#' rows for district-years that were not sampled; those have `n == 0` and are
#' dropped on construction, leaving 96 surveys.
#'
#' @param M Lattice steps per year.
#' @param ... Passed to [GCData]'s constructor.
#' @return A [GCData] object.
#' @export
gc_example_data <- function(M = 16, ...) {
  f <- system.file("extdata", "uganda_c469y.csv", package = "GeneCastR")
  if (!nzchar(f)) {
    stop("example data not found; is the package installed?", call. = FALSE)
  }
  raw <- utils::read.csv(f, stringsAsFactors = FALSE)
  GCData$new(raw, M = M,
             cols = c(deme = "district", time = "year", x = "x", n = "n"),
             ...)
}


#' The fitted Ugandan panel
#'
#' A production fit of the shipped panel, computed once and stored so the
#' vignette can show diagnostics worth reading: four chains of 4000 with 2000
#' warm-up, and a projection with 200 draws and 400 paths each. Running it
#' takes about four minutes on four cores.
#'
#' See `data-raw/uganda_fit.R` for the call that produced it.
#'
#' @return A list with `fit`, a [GCFit], and `pr`, a [GCProjection].
#' @examples
#' z <- gc_example_fit()
#' z$fit$summary()
#' @export
gc_example_fit <- function() {
  f <- system.file("extdata", "uganda_fit.rds", package = "GeneCastR")
  if (!nzchar(f))
    stop("the stored fit is missing; is the package installed?",
         call. = FALSE)
  readRDS(f)
}
