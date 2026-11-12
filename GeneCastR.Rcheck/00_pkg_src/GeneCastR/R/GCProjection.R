#' Smoothed and projected quantities from a fit
#'
#' Built by [GCFit]`$project()`. Inside the observation window these are
#' *smoothed* -- every survey in the deme, later ones included, informs the
#' state at each step -- and beyond the last survey they are projections from
#' the same kernel, started from the smoothed lattice distribution at the
#' last step.
#'
#' `p_present` can nevertheless step at the join, and the step can be large.
#' A survey at the final step is a hard constraint that the smoother obeys
#' and the projection does not: a deme seen positive there has
#' `p_present = 1` at `T` and less immediately after, because extinction
#' becomes possible again. A deme carrying mass at zero jumps the other way,
#' since a fraction `1 - exp(-lambda)` of that mass -- a quarter to a third
#' at the fitted influx -- arrives in a single step. On the Ugandan panel
#' the jump runs from -0.13 to +0.09. It is a property of the model, not an
#' artefact of stitching two calculations together: the lattice distribution
#' itself is propagated continuously.
#'
#' @section Two kinds of interval:
#' `f_q*` are posterior predictive quantiles of the latent frequency:
#' `theta`, the latent path and the Wright-Fisher noise are all integrated
#' out, so `f_q97.5` is a predictive bound on the frequency itself.
#' `p_present_q*` are quantiles of `P(present | X, theta)` **across draws**:
#' everything but `theta` is integrated out, so they say how well the data
#' pin a probability, not how uncertain the outcome is. Presence is a
#' Bernoulli whose probability is itself uncertain, and
#' `Var(1) = E[p(1-p)] + Var(p)`; the `p_present_q*` columns show only the
#' second term.
#'
#' @section What is unconditional:
#' The `f_*` columns carry the atom at zero: it has mass
#' `1 - p_present`, so any quantile below that level is exactly `0`, and a
#' band resting on zero means absence is the likelier state rather than a
#' frequency near zero.
#'
#' @section The two introduction-time columns:
#' Under `founding = "rate"` the first presence time is no longer `t*`:
#' the deme joins the flow process at `t*` and may sit empty for a while, so
#' `pmf_tfirst` is the posterior for `t*` and `T_first` from sampled paths
#' is a later, separate quantity.
#'
#' `pmf_tfirst` is free and exact. In this model the founding event deposits
#' a copy at `t*` and nothing is present before it, so the first presence
#' time *is* `t*`, which the likelihood already marginalises analytically.
#' `pmf_test` is a path functional -- no marginal of the smoother says when
#' the excursion carrying the first detected copy began -- so it needs
#' `n_paths > 0` and is `NA` otherwise. Walk back from the first survey with
#' a non-zero count until the deme is empty: that is where the lineage which
#' was actually seen began. It is a proper distribution in a deme that was
#' detected and absent in one that was not, there being no observed lineage
#' to date, and it does not condition on survival -- the lineage that was
#' seen need not be the one still present at `T`. It is bounded above by the
#' detection that defines it and below by `t*`.
#' Both are zero on projected rows, not `NA`: `t*` is supported on `1..T`
#' only, and that is the documented limitation of this model rather than
#' missing information.
#'
#' @export
GCProjection <- R6::R6Class(
  "GCProjection",
  public = list(

    #' @field curves One row per deme and step; see the class documentation.
    curves = NULL,
    #' @field demes One row per deme, for what has no time index. Carries
    #'   `p_tstar_beyond`, the posterior probability that the variant had
    #'   not reached THAT deme within the record, given its own surveys --
    #'   zero under the uniform prior on `t*`, and zero for any deme ever
    #'   seen positive. Not to be confused with the fit's
    #'   `p_beyond_baseline`, which is the same probability for a deme with
    #'   no data and is a transform of `theta` alone.
    demes = NULL,
    #' @field surveys The observations, carried through for plotting.
    surveys = NULL,
    #' @field thresholds The exceedance thresholds, so callers select columns
    #'   by lookup rather than by matching column names.
    thresholds = NULL,
    #' @field probs The quantile levels used for `f_q*`.
    probs = NULL,
    #' @field joint `T x T x D` array of `p(T_first = i, T_est = j)`, or
    #'   `NULL` when no paths were sampled. Supported on `j >= i` only, and
    #'   its diagonal is the probability that the lineage which was detected
    #'   is the one that originally arrived -- a quantity neither marginal
    #'   can give, and free, because every sampled path carries both times.
    joint = NULL,
    #' @field settings Horizon, draw count, seed and lattice settings.
    settings = NULL,

    #' @description Assemble the tables. Called by [GCFit]`$project()`.
    #' @param fit The [GCFit].
    #' @param horizon Time units to project beyond the last step.
    #' @param n_draws Posterior draws to average over.
    #' @param thresholds Frequencies at which to report exceedance.
    #' @param probs Quantile levels for the frequency columns.
    #' @param nbins Bins for the frequency density.
    #' @param n_paths Paths sampled per deme per draw for `T_est`. Zero
    #'   skips the path sampler and leaves the `*_test` columns `NA`; that
    #'   is the fast path, and everything else is unaffected.
    initialize = function(fit, horizon = 3, n_draws = 120L,
                          thresholds = c(0.01, 0.05, 0.10),
                          probs = c(.025, .25, .5, .75, .975),
                          nbins = 250L, n_paths = 0L) {
      d <- fit$data
      M <- d$M; T <- d$T
      H <- as.integer(round(horizon * M))
      if (H < 0) stop("`horizon` must be non-negative", call. = FALSE)
      stopifnot(all(thresholds >= 0), all(thresholds < 1),
                all(probs > 0), all(probs < 1))
      thresholds <- sort(unique(as.numeric(thresholds)))
      probs <- sort(unique(as.numeric(probs)))

      m <- posterior::as_draws_matrix(fit$draws)
      idx <- unique(round(seq(1, nrow(m), length.out = n_draws)))
      D <- d$n_demes
      Z <- log(as.matrix(m[idx, seq_len(D), drop = FALSE]))
      dimnames(Z) <- NULL

      fc <- .founding_code(fit$control$founding)
      geom <- identical(fit$control$tstar, "geometric")
      ## the smoother and the path sampler take one hazard for the whole
      ## sweep, so the posterior median stands in; the t* posterior below
      ## uses each draw's own, which is the finer average
      tp <- if (geom) -expm1(-stats::median(as.numeric(m[idx, "theta"])) / M)
            else -1
      r <- gc_project_pt(Z, as.numeric(m[idx, "u"]), as.numeric(m[idx, "nu"]),
                         M = M, T = T, H = H, demes = d$deme_data(),
                         thresholds = thresholds, nbins = nbins,
                         founding = fc, tstar_p = tp)

      TH <- T + H
      fmid <- (seq_len(nbins) - 0.5) / nbins
      step <- rep(seq_len(TH), times = D)
      deme <- rep(d$demes, each = TH)

      ## unconditional frequency quantiles from the mixture over draws, with
      ## the atom at zero carrying mass 1 - p_present
      fq <- matrix(0, D * TH, length(probs))
      for (dd in seq_len(D)) for (t in seq_len(TH)) {
        row <- r$dens[, t, dd]
        p <- r$pres[dd, t]
        rs <- sum(row)
        if (rs > 0 && p > 0) {
          cw <- cumsum(row / rs)
          for (j in seq_along(probs)) {
            q <- probs[j]
            fq[(dd - 1L) * TH + t, j] <-
              if (q <= 1 - p) 0 else
                stats::approx(cw, fmid, xout = (q - (1 - p)) / p,
                              rule = 2, ties = "ordered")$y
          }
        }
      }

      ## P(present | X, theta) quantiles across draws
      pq <- c(.025, .5, .975)
      pqm <- matrix(NA_real_, D * TH, 3L)
      for (dd in seq_len(D)) {
        v <- r$pres_draws[, dd, , drop = FALSE]       # TH x 1 x R
        dim(v) <- c(TH, dim(r$pres_draws)[3L])
        pqm[((dd - 1L) * TH + 1L):(dd * TH), ] <-
          t(apply(v, 1L, stats::quantile, probs = pq, names = FALSE))
      }

      ## t* posterior: analytic, averaged over the same draws.  The pmf
      ## columns then sum to 1 - p_tstar_beyond rather than to one.
      ## as.numeric(): a draws_matrix keeps its dimensions on a single-cell
      ## subset, and a 1 x 1 matrix does not recycle against a deme vector
      uu <- as.numeric(m[idx, "u"]); nn <- as.numeric(m[idx, "nu"])
      ddat <- d$deme_data()
      ## D x (T+1): the last column is p(t* > T), the probability the
      ## variant had not reached the deme within the record at all
      th <- if (geom) as.numeric(m[idx, "theta"]) else rep(NA_real_, nrow(Z))
      ts <- matrix(0, D, T + 1L)
      for (i in seq_len(nrow(Z))) {
        tpi <- if (geom) -expm1(-th[i] / M) else -1
        ts <- ts + gc_tstar_pt(gc_nlat(Z[i, ], M), gc_alpha(uu[i], M),
                               gc_lambda(nn[i], exp(Z[i, ])), T, ddat, fc, tpi)
      }
      ts <- ts / nrow(Z)
      beyond <- ts[, T + 1L]
      ts <- ts[, seq_len(T), drop = FALSE]
      pmf_tf <- numeric(D * TH)
      for (dd in seq_len(D)) pmf_tf[((dd - 1L) * TH + 1L):((dd - 1L) * TH + T)] <-
        ts[dd, ]

      out <- data.frame(
        deme = deme, step = step, time = d$t0 + step / M,
        phase = factor(ifelse(step <= T, "smoothed", "projected"),
                       levels = c("smoothed", "projected")),
        f_mean = as.vector(t(r$fmean)),
        stringsAsFactors = FALSE
      )
      for (j in seq_along(probs))
        out[[private$qname("f_q", probs[j])]] <- fq[, j]
      out$p_present <- as.vector(t(r$pres))
      out$p_present_q2.5 <- pqm[, 1L]
      out$p_present_q50  <- pqm[, 2L]
      out$p_present_q97.5 <- pqm[, 3L]
      for (j in seq_along(thresholds)) {
        v <- r$exceed[j, , ]
        out[[private$qname("p_exceed_", thresholds[j], 100)]] <- as.vector(v)
      }
      out$pmf_tfirst <- pmf_tf
      out$cdf_tfirst <- unlist(lapply(split(pmf_tf, deme)[unique(deme)],
                                      cumsum), use.names = FALSE)
      p_alive <- rep(NA_real_, D)
      if (n_paths > 0L) {
        tm <- gc_times_pt(Z, uu, nn, M = M, T = T, demes = ddat,
                          n_paths = as.integer(n_paths), joint = TRUE,
                          founding = fc, tstar_p = tp)
        self$joint <- tm$joint
        pm <- numeric(D * TH)
        for (dd in seq_len(D))
          pm[((dd - 1L) * TH + 1L):((dd - 1L) * TH + T)] <- tm$pmf_test[dd, ]
        out$pmf_test <- pm
        out$cdf_test <- unlist(lapply(split(pm, deme)[unique(deme)], cumsum),
                               use.names = FALSE)
        p_alive <- tm$p_alive
      } else {
        out$pmf_test <- NA_real_
        out$cdf_test <- NA_real_
      }

      ## per-deme scalars: a different grain, not a schema variant
      ds <- d$summary()
      last <- out[out$step == TH, ]
      last <- last[match(d$demes, last$deme), ]
      self$demes <- data.frame(
        ds,
        N_eff_median = apply(as.matrix(m[idx, seq_len(D), drop = FALSE]), 2L,
                             stats::median),
        p_tstar_beyond = beyond,
        p_alive_T = p_alive,
        p_present_final = last$p_present,
        f_median_final = last[[private$qname("f_q", 0.5)]],
        row.names = NULL, stringsAsFactors = FALSE
      )

      self$curves <- out
      self$surveys <- d$surveys
      self$thresholds <- thresholds
      self$probs <- probs
      self$settings <- list(M = M, T = T, H = H, horizon = horizon,
                            t0 = d$t0, n_draws = length(idx), nbins = nbins,
                            n_paths = n_paths)
      invisible(self)
    },

    #' @description Print a compact description.
    #' @param ... Ignored.
    print = function(...) {
      s <- self$settings
      cat("<GCProjection>\n")
      cat("  ", nrow(self$curves), " rows: ", nrow(self$demes), " demes x ",
          s$T + s$H, " steps (", s$T, " smoothed, ", s$H, " projected)\n",
          sep = "")
      cat("  window   : ", format(s$t0 + 1 / s$M), " to ",
          format(s$t0 + (s$T + s$H) / s$M), "\n", sep = "")
      cat("  draws    : ", s$n_draws,
          if (s$n_paths > 0L) paste0(", ", s$n_paths, " paths each")
          else "  (no path sampling: T_est absent)", "\n", sep = "")
      cat("  exceedance thresholds: ",
          paste(self$thresholds, collapse = ", "), "\n", sep = "")
      if (any(self$demes$p_tstar_beyond > 0))
        cat("  P(t* beyond the record): mean ",
            sprintf("%.3f", mean(self$demes$p_tstar_beyond)), ", max ",
            sprintf("%.3f", max(self$demes$p_tstar_beyond)), "\n", sep = "")
      fin <- self$curves[self$curves$step == s$T + s$H, ]
      cat("\n  at ", format(s$t0 + (s$T + s$H) / s$M),
          ": mean P(present) ", sprintf("%.3f", mean(fin$p_present)),
          ", mean median f ",
          sprintf("%.3f", mean(fin[[private$qname("f_q", 0.5)]])), "\n",
          sep = "")
      invisible(self)
    }
  ),

  private = list(
    ## formatC() pads to a width here, which would put spaces inside column
    ## names; format(trim = TRUE) does not.
    qname = function(prefix, p, scale = 100) {
      paste0(prefix, format(p * scale, trim = TRUE, scientific = FALSE,
                            drop0trailing = TRUE))
    }
  )
)
