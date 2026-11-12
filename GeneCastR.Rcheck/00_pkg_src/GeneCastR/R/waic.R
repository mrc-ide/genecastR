#' Widely applicable information criterion
#'
#' WAIC is computed at the **deme** level, not the survey level. The
#' likelihood marginalises each deme's whole latent path, so it does not
#' factorise over that deme's surveys: there are as many terms as there are
#' demes, and no more. That makes the standard errors crude, and it makes a
#' comparison between two models that are not the best one need its own
#' standard error rather than the difference of two columns.
#'
#' @param fit A [GCFit].
#' @param n_draws Posterior draws to use, thinned evenly.
#' @return A list with `waic`, `elpd`, `lppd`, `p_waic`, the per-deme
#'   `pointwise` contributions and `se`, the standard error of `waic` across
#'   demes.
#' @examples
#' \donttest{
#' z <- gc_example_fit()
#' gc_waic(z$fit)
#' }
#' @export
gc_waic <- function(fit, n_draws = 500L) {
  d <- fit$data
  m <- posterior::as_draws_matrix(fit$draws)
  idx <- unique(round(seq(1, nrow(m), length.out = n_draws)))
  D <- d$n_demes
  Z <- log(as.matrix(m[idx, seq_len(D), drop = FALSE])); dimnames(Z) <- NULL
  ## a per-deme fit carries one u per deme; otherwise the panel shares one.
  ## `%||%` guards fits saved before a control existed, which would
  ## otherwise give a zero-length condition rather than a sensible default.
  sel <- fit$control$selection
  perd <- !is.null(sel) && identical(sel, "hierarchical")
  U <- if (perd) as.matrix(m[idx, paste0("u[", d$demes, "]"), drop = FALSE])
       else matrix(as.numeric(m[idx, "u"]), length(idx), D)
  dimnames(U) <- NULL
  nn <- as.numeric(m[idx, "nu"])
  geom <- identical(fit$control$tstar, "geometric")
  fnd <- fit$control$founding
  th <- if (geom) as.numeric(m[idx, "theta"]) else rep(NA_real_, length(idx))
  fc <- .founding_code(if (is.null(fnd)) "conditional" else fnd)
  ddat <- d$deme_data()
  ## WAIC must score the model that was fitted, so the t* prior each draw
  ## carried has to be rebuilt here too, not left at the default
  nb <- identical(fit$control$tstar, "negbin")
  tset <- if (nb) .tstar_setup(fit$control, d) else NULL
  sh <- if (nb) as.numeric(m[idx, "tstar_shape"]) else rep(NA_real_, length(idx))
  mo <- if (nb && tset$learn_mean)
    as.numeric(m[idx, "tstar_mean_offset"]) * d$M else NULL

  LL <- matrix(NA_real_, length(idx), D)
  for (i in seq_along(idx)) {
    tp <- if (geom) -expm1(-th[i] / d$M) else -1
    lwi <- if (nb) .tstar_logw(tset, d$T, sh[i],
                               if (is.null(mo)) NULL else mo[i]) else NULL
    ## gc_loglik_pt takes one alpha, so a per-deme fit needs one call per
    ## deme; the shared case still does all D in a single call
    LL[i, ] <- if (perd)
      vapply(seq_len(D), function(k)
        gc_loglik_pt(gc_nlat(Z[i, k], d$M), gc_alpha(U[i, k], d$M),
                     gc_lambda(nn[i], exp(Z[i, k])), d$T, ddat[k], fc, tp,
                     if (is.null(lwi)) NULL else lwi[, k, drop = FALSE]),
        numeric(1L))
    else
      gc_loglik_pt(gc_nlat(Z[i, ], d$M), gc_alpha(U[i, 1L], d$M),
                   gc_lambda(nn[i], exp(Z[i, ])), d$T, ddat, fc, tp, lwi)
  }

  lse <- function(v) { mx <- max(v); mx + log(sum(exp(v - mx))) }
  lppd_d <- apply(LL, 2L, function(v) lse(v) - log(length(v)))
  p_d <- apply(LL, 2L, stats::var)
  elpd_d <- lppd_d - p_d

  list(waic = -2 * sum(elpd_d), elpd = sum(elpd_d),
       lppd = sum(lppd_d), p_waic = sum(p_d),
       se = sqrt(D * stats::var(-2 * elpd_d)),
       n_draws = length(idx),
       pointwise = data.frame(deme = d$demes, lppd = lppd_d, p_waic = p_d,
                              elpd = elpd_d, row.names = NULL))
}


#' Compare models by WAIC
#'
#' @param ... Two or more [GCFit] objects, ideally named.
#' @param n_draws Posterior draws to use per fit.
#' @return A data frame ordered best first, with `d_waic` against the best
#'   model and `se_diff`, the standard error of **that pairwise difference**
#'   computed from the per-deme contributions. A difference inside one or
#'   two of its own standard errors is not a separation.
#' @export
gc_waic_compare <- function(..., n_draws = 500L) {
  fits <- list(...)
  nm <- names(fits)
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("model", seq_along(fits))
  w <- lapply(fits, gc_waic, n_draws = n_draws)
  waic <- vapply(w, function(x) x$waic, numeric(1L))
  ord <- order(waic)
  best <- w[[ord[1L]]]
  se_diff <- vapply(w, function(x) {
    dd <- -2 * x$pointwise$elpd - (-2 * best$pointwise$elpd)
    if (all(dd == 0)) 0 else sqrt(length(dd) * stats::var(dd))
  }, numeric(1L))
  out <- data.frame(model = nm, waic = waic,
                    p_waic = vapply(w, function(x) x$p_waic, numeric(1L)),
                    se = vapply(w, function(x) x$se, numeric(1L)),
                    d_waic = waic - waic[ord[1L]], se_diff = se_diff,
                    row.names = NULL)
  out[ord, ]
}
