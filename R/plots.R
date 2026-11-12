## Plot helpers.  Every function returns a ggplot and adds nothing a caller
## cannot remove; none of them print, save or set options.

.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is needed for the plotting functions; install it, or work ",
         "from `pr$curves` directly.", call. = FALSE)
}

.COL <- c(freq = "#c0392b", present = "#2a78d6", ink = "#0b0b0b",
          muted = "#8a8982", tfirst = "#2f7d4f", test = "#7b4fa6",
          shade = "#f6f2ea", hot = "#d62728")

## Contiguous runs where presence exceeds a level, per deme, as rectangles
## spanning the steps they cover.  The half-step padding makes neighbouring
## runs meet rather than leaving hairlines between them.
.above_runs <- function(cv, level) {
  ts <- sort(unique(cv$time))
  half <- if (length(ts) > 1L) stats::median(diff(ts)) / 2 else 0
  out <- do.call(rbind, lapply(split(cv, cv$deme), function(x) {
    x <- x[order(x$time), ]
    a <- !is.na(x$p_present) & x$p_present > level
    if (!any(a)) return(NULL)
    r <- rle(a)
    e <- cumsum(r$lengths); b <- e - r$lengths + 1L
    k <- which(r$values)
    data.frame(deme = as.character(x$deme[1L]),
               xmin = x$time[b[k]] - half, xmax = x$time[e[k]] + half,
               stringsAsFactors = FALSE)
  }))
  if (is.null(out)) return(NULL)
  out$deme <- factor(out$deme, levels = levels(cv$deme))
  out
}

.gc_theme <- function() {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 11.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.4, colour = "grey25"))
}

## Column name for a quantile level, matching GCProjection's own naming.
.qcol <- function(prefix, p, scale = 100) {
  paste0(prefix, format(p * scale, trim = TRUE, scientific = FALSE,
                        drop0trailing = TRUE))
}

## Facet order and labels.  Panels are ordered by something the reader needs
## anyway, so the ordering itself carries information.
## Facet order and labels.  Every plot uses the same order by default, so a
## deme sits in the same panel position whichever figure you are looking at.
## Alphabetical, or numeric where the names are numbers -- "2" before "10",
## which a plain sort gets wrong.
.deme_order <- function(demes) {
  n <- suppressWarnings(as.numeric(demes))
  if (!anyNA(n)) demes[order(n)] else demes[order(demes)]
}

.panels <- function(pr, by = c("name", "p_present", "p_alive"),
                    extra = NULL) {
  by <- match.arg(by)
  dm <- pr$demes
  ord <- switch(by,
    name      = .deme_order(dm$deme),
    p_present = dm$deme[order(-dm$p_present_final)],
    p_alive   = dm$deme[order(dm$p_alive_T, -dm$n_nonzero)])
  lab <- sprintf("%s  (%d/%d non-zero)", dm$deme, dm$n_nonzero, dm$n_surveys)
  if (!is.null(extra)) lab <- paste0(lab, extra)
  list(order = ord, labels = stats::setNames(lab, dm$deme))
}

.subset_demes <- function(pr, demes) {
  if (is.null(demes)) return(pr$demes$deme)
  bad <- setdiff(demes, pr$demes$deme)
  if (length(bad))
    stop("no such deme(s): ", paste(bad, collapse = ", "), call. = FALSE)
  demes
}

.window <- function(pr) {
  s <- pr$settings
  c(s$t0 + 1 / s$M, s$t0 + s$T / s$M, s$t0 + (s$T + s$H) / s$M)
}


## Expand a variable spec: an exact name, or a prefix like "N_eff" standing
## for every indexed variable with that stem.
.expand_vars <- function(fit, vars) {
  have <- posterior::variables(fit$draws)
  out <- unlist(lapply(vars, function(v) {
    if (v %in% have) return(v)
    hit <- grep(paste0("^", v, "\\["), have, value = TRUE)
    if (length(hit)) return(hit)
    stop("no such variable: ", v, call. = FALSE)
  }), use.names = FALSE)
  out <- unique(out)
  ## A structure the fit does not use is reported as NA rather than left out,
  ## so the draws array has one shape whatever `gc_control()` was set to --
  ## `sigma` under `neff = "shared"` is the case that matters, that being the
  ## default.  Those columns are not parameters of this model and there is
  ## nothing to draw for them.
  m <- posterior::as_draws_matrix(
    .exact_draws(fit$draws, out))
  keep <- out[vapply(out, function(v) any(is.finite(m[, v])), TRUE)]
  if (!length(keep)) {
    stop("nothing to plot: ", paste(out, collapse = ", "),
         " are not parameters of this fit.", call. = FALSE)
  }
  keep
}

## Draws as a long data frame, chain kept.
.long_draws <- function(fit, vars) {
  d <- .exact_draws(fit$draws, vars)
  df <- as.data.frame(posterior::as_draws_df(d))
  do.call(rbind, lapply(vars, function(v)
    data.frame(iteration = df$.iteration, chain = factor(df$.chain),
               variable = factor(v, levels = vars), value = df[[v]])))
}


#' Latent frequency with predictive bands
#'
#' The bands are **unconditional**: the atom at zero carries
#' `1 - p_present`, so a band resting on zero means absence is the likelier
#' state rather than a frequency near zero. They are predictive -- the
#' parameters, the latent path and the Wright-Fisher noise are all
#' integrated out.
#'
#' @param pr A [GCProjection].
#' @param demes Demes to show; `NULL` for all.
#' @param order Facet order: `"name"` (the default, alphabetical, or numeric
#'   where the names are numbers) so that a deme keeps the same panel
#'   position across every figure; `"p_present"` or `"p_alive"` to sort by
#'   the quantity instead.
#' @param inner Also draw the inner band, when those quantiles were computed.
#' @param present Overlay `P(present | X)` as a dashed line.
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_frequency <- function(pr, demes = NULL, inner = TRUE,
                              present = TRUE, order = "name", ncol = 4L) {
  .need_ggplot()
  keep <- .subset_demes(pr, demes)
  p <- .panels(pr, by = order)
  cv <- pr$curves[pr$curves$deme %in% keep, ]
  cv$deme <- factor(cv$deme, levels = intersect(p$order, keep))
  sv <- pr$surveys[pr$surveys$deme %in% keep, ]
  sv$deme <- factor(sv$deme, levels = levels(cv$deme))
  w <- .window(pr)
  lo <- .qcol("f_q", min(pr$probs)); hi <- .qcol("f_q", max(pr$probs))
  mid <- if (0.5 %in% pr$probs) .qcol("f_q", 0.5) else NULL

  g <- ggplot2::ggplot(cv, ggplot2::aes(.data$time))
  if (w[3] > w[2])
    g <- g + ggplot2::annotate("rect", xmin = w[2], xmax = w[3],
                               ymin = -Inf, ymax = Inf, fill = .COL[["shade"]])
  g <- g + ggplot2::geom_ribbon(
    ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]),
    fill = .COL[["freq"]], alpha = .13)
  inner_pair <- pr$probs[pr$probs > min(pr$probs) & pr$probs < max(pr$probs) &
                         pr$probs != 0.5]
  if (inner && length(inner_pair) >= 2L) {
    il <- .qcol("f_q", min(inner_pair)); ih <- .qcol("f_q", max(inner_pair))
    g <- g + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data[[il]], ymax = .data[[ih]]),
      fill = .COL[["freq"]], alpha = .20)
  }
  if (!is.null(mid))
    g <- g + ggplot2::geom_line(ggplot2::aes(y = .data[[mid]]),
                                colour = .COL[["freq"]], linewidth = .7)
  if (present)
    g <- g + ggplot2::geom_line(ggplot2::aes(y = .data$p_present),
                                colour = .COL[["present"]], linewidth = .5,
                                linetype = 2)
  g +
    ggplot2::geom_vline(xintercept = w[2], colour = .COL[["ink"]],
                        linewidth = .3, linetype = 3) +
    ggplot2::geom_pointrange(
      data = sv, ggplot2::aes(x = .data$time, y = .data$f,
                              ymin = .data$lo, ymax = .data$hi),
      size = .18, linewidth = .35, colour = .COL[["ink"]],
      inherit.aes = FALSE) +
    ggplot2::facet_wrap(~deme, ncol = ncol,
                        labeller = ggplot2::labeller(deme = p$labels)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    .gc_theme() +
    ggplot2::labs(x = "year", y = "frequency",
      title = "Latent frequency, smoothed and projected")
}


#' Presence probability
#'
#' Two different objects, and the distinction matters. `P(present | X)` is
#' fully integrated. The quantiles of `P(present | X, theta)` across draws
#' hold `theta` out, and say how well the data pin a probability rather than
#' how uncertain the outcome is: presence is a Bernoulli whose probability is
#' itself uncertain, and the across-draw spread is only part of the variance.
#'
#' @param pr A [GCProjection].
#' @param demes Demes to show; `NULL` for all.
#' @param order Facet order; see [gc_plot_frequency].
#' @param bands Draw the across-draw quantiles behind the mean. Off by
#'   default: they say how well the data pin the probability rather than how
#'   uncertain the outcome is, and three overlapping fills are hard to read
#'   at once.
#' @param highlight Shade the stretches where presence exceeds this
#'   probability, or `NULL` for none. The comparison is against the mean, so
#'   a shaded stretch reads as "probably present" rather than "confidently
#'   probably present" -- for the stricter reading, pass `bands = TRUE` and
#'   look at where the lower edge sits.
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_present <- function(pr, demes = NULL, bands = FALSE,
                            highlight = 0.8, order = "name", ncol = 4L) {
  .need_ggplot()
  keep <- .subset_demes(pr, demes)
  p <- .panels(pr, by = order)
  cv <- pr$curves[pr$curves$deme %in% keep, ]
  cv$deme <- factor(cv$deme, levels = intersect(p$order, keep))
  sv <- pr$surveys[pr$surveys$deme %in% keep & pr$surveys$x > 0, ]
  sv$deme <- factor(sv$deme, levels = levels(cv$deme))
  w <- .window(pr)

  g <- ggplot2::ggplot(cv, ggplot2::aes(.data$time))
  if (w[3] > w[2])
    g <- g + ggplot2::annotate("rect", xmin = w[2], xmax = w[3],
                               ymin = -Inf, ymax = Inf, fill = .COL[["shade"]])
  if (!is.null(highlight) && !is.na(highlight)) {
    if (!(length(highlight) == 1L && highlight > 0 && highlight < 1))
      stop("`highlight` must be a single probability in (0, 1), or NULL.",
           call. = FALSE)
    hl <- .above_runs(cv, highlight)
    if (!is.null(hl))
      g <- g + ggplot2::geom_rect(
        data = hl, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                ymin = -Inf, ymax = Inf),
        fill = .COL[["hot"]], alpha = .12, inherit.aes = FALSE)
  }
  if (bands)
    g <- g + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$p_present_q2.5, ymax = .data$p_present_q97.5),
      fill = .COL[["present"]], alpha = .18)
  g +
    ggplot2::geom_line(ggplot2::aes(y = .data$p_present),
                       colour = .COL[["ink"]], linewidth = .8) +
    ggplot2::geom_vline(xintercept = w[2], colour = .COL[["ink"]],
                        linewidth = .3, linetype = 3) +
    ggplot2::geom_rug(data = sv, ggplot2::aes(x = .data$time), sides = "b",
                      colour = .COL[["ink"]], linewidth = .4,
                      inherit.aes = FALSE) +
    ggplot2::facet_wrap(~deme, ncol = ncol,
                        labeller = ggplot2::labeller(deme = p$labels)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    .gc_theme() +
    ggplot2::labs(x = "year", y = "P(present)",
      title = "Presence probability",
      subtitle = if (!is.null(highlight) && !is.na(highlight))
        sprintf("shaded where presence exceeds %.2f", highlight))
}


#' Exceedance probability
#'
#' @param pr A [GCProjection].
#' @param thresholds Thresholds to draw; defaults to all that were computed.
#' @param demes Demes to show; `NULL` for all.
#' @param order Facet order; see [gc_plot_frequency].
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_exceed <- function(pr, thresholds = NULL, demes = NULL,
                           order = "name", ncol = 4L) {
  .need_ggplot()
  if (!length(pr$thresholds))
    stop("no exceedance thresholds were computed; pass `thresholds` to ",
         "project()", call. = FALSE)
  thresholds <- if (is.null(thresholds)) pr$thresholds else thresholds
  bad <- setdiff(thresholds, pr$thresholds)
  if (length(bad))
    stop("threshold(s) not computed: ", paste(bad, collapse = ", "),
         ".  Available: ", paste(pr$thresholds, collapse = ", "),
         call. = FALSE)
  keep <- .subset_demes(pr, demes)
  p <- .panels(pr, by = order)
  cv <- pr$curves[pr$curves$deme %in% keep, ]
  lv <- intersect(p$order, keep)
  w <- .window(pr)

  lab <- paste0(format(thresholds * 100, trim = TRUE, drop0trailing = TRUE),
                "%")
  long <- do.call(rbind, lapply(seq_along(thresholds), function(i)
    data.frame(deme = factor(cv$deme, levels = lv), time = cv$time,
               thr = factor(lab[i], levels = lab),
               p = cv[[.qcol("p_exceed_", thresholds[i])]])))

  g <- ggplot2::ggplot(long, ggplot2::aes(.data$time, .data$p,
                                          colour = .data$thr))
  if (w[3] > w[2])
    g <- g + ggplot2::annotate("rect", xmin = w[2], xmax = w[3],
                               ymin = -Inf, ymax = Inf, fill = .COL[["shade"]])
  g +
    ggplot2::geom_line(linewidth = .7) +
    ggplot2::geom_vline(xintercept = w[2], colour = .COL[["ink"]],
                        linewidth = .3, linetype = 3) +
    ggplot2::scale_colour_viridis_d(option = "plasma", end = .8,
                                    name = "P(f > x)") +
    ggplot2::facet_wrap(~deme, ncol = ncol,
                        labeller = ggplot2::labeller(deme = p$labels)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    .gc_theme() +
    ggplot2::labs(x = "year", y = "probability",
      title = "Exceedance: P(latent frequency above a threshold)")
}


#' Introduction and establishment times
#'
#' `t*` is the founding event. `T_est` is the start of the excursion that
#' produced the first detected count: walk back from the first survey with a
#' non-zero count until the deme is empty. Both integrate to one in a deme
#' that was detected; `T_est` is absent in one that was not, there being no
#' observed lineage to date. `T_est` lies between `t*` and the first
#' detection, and the gap from `t*` is the time the lineage spent
#' undetected -- or the span over which an earlier arrival died out.
#'
#' @param pr A [GCProjection].
#' @param type `"cdf"` or `"pmf"`. Densities are per unit time, so the two
#'   quantities are on one scale whatever the lattice resolution.
#' @param which Which times to draw.
#' @param demes Demes to show; `NULL` for all.
#' @param order Facet order; see [gc_plot_frequency].
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_times <- function(pr, type = c("cdf", "pmf"),
                          which = c("both", "tfirst", "test"),
                          demes = NULL, order = "name", ncol = 4L) {
  .need_ggplot()
  type <- match.arg(type); which <- match.arg(which)
  have_test <- !all(is.na(pr$curves$pmf_test))
  if (which == "test" && !have_test)
    stop("T_est needs sampled paths; call project(n_paths = ...)",
         call. = FALSE)
  keep <- .subset_demes(pr, demes)
  dm <- pr$demes
  ## T_est is the start of the excursion that produced the first detected
  ## count, so the reference a reader needs is that detection, not survival
  ## to T.  A deme that was never detected has no T_est at all.
  extra <- if (have_test && which != "tfirst")
    ifelse(is.na(dm$cut_time), "   never detected",
           sprintf("   first seen %.0f", dm$cut_time)) else NULL
  p <- .panels(pr, by = order, extra = extra)
  cv <- pr$curves[pr$curves$deme %in% keep & pr$curves$phase == "smoothed", ]
  lv <- intersect(p$order, keep)
  M <- pr$settings$M
  sc <- if (type == "pmf") M else 1

  want <- switch(which, both = c("tfirst", "test"), which)
  if (!have_test) want <- setdiff(want, "test")
  nm <- c(tfirst = "t* = T_first", test = "T_est")
  long <- do.call(rbind, lapply(want, function(k)
    data.frame(deme = factor(cv$deme, levels = lv), time = cv$time,
               what = factor(nm[[k]], levels = unname(nm)),
               v = cv[[paste0(type, "_", k)]] * sc)))

  cut <- dm[!is.na(dm$cut_time) & dm$deme %in% keep, c("deme", "cut_time")]
  cut$deme <- factor(cut$deme, levels = lv)

  g <- ggplot2::ggplot(long, ggplot2::aes(.data$time, .data$v,
                                          colour = .data$what))
  g <- g + ggplot2::geom_vline(data = cut,
      ggplot2::aes(xintercept = .data$cut_time), colour = .COL[["muted"]],
      linewidth = .4, linetype = 2, inherit.aes = FALSE)
  g <- g + if (type == "cdf") ggplot2::geom_step(linewidth = .8)
           else ggplot2::geom_line(linewidth = .7)
  g +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c(.COL[["tfirst"]], .COL[["test"]]),
                               unname(nm)),
      name = NULL, drop = TRUE) +
    ggplot2::facet_wrap(~deme, ncol = ncol,
      labeller = ggplot2::labeller(deme = p$labels),
      scales = if (type == "pmf") "free_y" else "fixed") +
    (if (type == "cdf") ggplot2::coord_cartesian(ylim = c(0, 1)) else NULL) +
    .gc_theme() +
    ggplot2::labs(x = "year",
      y = if (type == "cdf") "posterior probability"
          else "posterior density (per unit time)",
      title = "Introduction and establishment times")
}


#' Joint posterior of the two introduction times
#'
#' Every sampled path carries both times, so the pairing costs nothing. Mass
#' lies on and above the diagonal only: the excursion alive at `T` cannot
#' begin before the variant first appears. On the diagonal means the lineage
#' that was detected is the one that originally arrived; above it means the
#' first arrival died out and a later one founded what was eventually seen.
#' Neither marginal can distinguish the two.
#'
#' @param pr A [GCProjection] built with `n_paths > 0`.
#' @param demes Demes to show; `NULL` for all.
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_joint <- function(pr, demes = NULL, ncol = 4L) {
  .need_ggplot()
  if (is.null(pr$joint))
    stop("the joint needs sampled paths; call project(n_paths = ...)",
         call. = FALSE)
  keep <- .subset_demes(pr, demes)
  dm <- pr$demes
  s <- pr$settings
  yr <- s$t0 + seq_len(s$T) / s$M
  idx <- which(dm$deme %in% keep)
  diagonal <- vapply(idx, function(k) sum(diag(pr$joint[, , k])), numeric(1L))
  share <- diagonal / pmax(dm$p_alive_T[idx], 1e-12)

  df <- do.call(rbind, lapply(seq_along(idx), function(i) {
    m <- pr$joint[, , idx[i]]
    nz <- which(m > 0, arr.ind = TRUE)
    if (!nrow(nz)) return(NULL)
    data.frame(deme = dm$deme[idx[i]], tf = yr[nz[, 1L]], te = yr[nz[, 2L]],
               p = m[nz] * s$M * s$M)
  }))
  if (is.null(df)) stop("no sampled paths survived to T", call. = FALSE)
  ord <- .deme_order(dm$deme[idx])
  lab <- stats::setNames(
    sprintf("%s  (%d/%d)   P(alive) %.2f   diagonal %.0f%%",
            dm$deme[idx], dm$n_nonzero[idx], dm$n_surveys[idx],
            dm$p_alive_T[idx], 100 * share), dm$deme[idx])
  df$deme <- factor(df$deme, levels = ord)

  ggplot2::ggplot(df, ggplot2::aes(.data$tf, .data$te, fill = .data$p)) +
    ggplot2::geom_raster() +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "white",
                         linewidth = .35, linetype = 2) +
    ggplot2::facet_wrap(~deme, ncol = ncol,
                        labeller = ggplot2::labeller(deme = lab)) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1,
                                  trans = "sqrt",
                                  name = "density\nper unit time\u00b2") +
    ggplot2::coord_fixed(xlim = range(yr), ylim = range(yr)) +
    .gc_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::labs(x = "t* = T_first  (first appearance)",
      y = "T_est  (founding of the detected lineage)",
      title = "Joint posterior of t* and T_est")
}


#' Trace plots for the shared parameters
#'
#' @param fit A [GCFit].
#' @param vars Variables to show. A bare stem such as `"N_eff"` stands for
#'   every indexed variable with that name. `NULL`, the default, takes the
#'   shared parameters, naming the panel-level size as the fitted structure
#'   does.
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_trace <- function(fit, vars = NULL, ncol = 2L) {
  .need_ggplot()
  if (is.null(vars)) vars <- c("u", "nu", .neff_var(fit$control), "sigma")
  vars <- .expand_vars(fit, vars)
  d <- .exact_draws(fit$draws, vars)
  long <- .long_draws(fit, vars)
  rh <- posterior::summarise_draws(d, "rhat", "ess_bulk")
  lab <- stats::setNames(
    sprintf("%s    Rhat %.3f   ESS %.0f", rh$variable, rh$rhat, rh$ess_bulk),
    rh$variable)

  ggplot2::ggplot(long, ggplot2::aes(.data$iteration, .data$value,
                                     colour = .data$chain)) +
    ggplot2::geom_line(linewidth = .3, alpha = .8) +
    ggplot2::facet_wrap(~variable, ncol = ncol, scales = "free_y",
                        labeller = ggplot2::labeller(variable = lab)) +
    ggplot2::scale_colour_viridis_d(option = "viridis", end = .85,
                                    name = "chain") +
    .gc_theme() +
    ggplot2::labs(x = "iteration (warm-up discarded)", y = NULL,
      title = "Shared parameters")
}



#' Posterior histograms, by chain
#'
#' One panel per variable, with the chains drawn separately rather than
#' shared. Chains that disagree show up here as displaced humps, which a
#' shared histogram hides and a trace plot can make hard to judge.
#'
#' @param fit A [GCFit].
#' @param vars Variables to show. A bare stem such as `"N_eff"` stands for
#'   every indexed variable with that name. `NULL`, the default, takes the
#'   shared parameters, naming the panel-level size as the fitted structure
#'   does.
#' @param bins Histogram bins.
#' @param alpha Opacity of the bars. Kept low on purpose: four overlapping
#'   fills cannot encode chain identity at any opacity, so the outline does
#'   that work and the fill only marks where the mass is.
#' @param linewidth Width of the step outline.
#' @param ncol Facet columns.
#' @return A ggplot.
#' @export
gc_plot_hist <- function(fit, vars = NULL,
                         bins = 35L, alpha = 0.1, linewidth = 0.8,
                         ncol = 2L) {
  .need_ggplot()
  if (is.null(vars))
    vars <- c("u", "A", "nu", .neff_var(fit$control), "sigma")
  vars <- .expand_vars(fit, vars)
  long <- .long_draws(fit, vars)
  rh <- posterior::summarise_draws(
    .exact_draws(fit$draws, vars), "rhat", "ess_bulk")
  lab <- stats::setNames(
    sprintf("%s    Rhat %.3f   ESS %.0f", rh$variable, rh$rhat, rh$ess_bulk),
    rh$variable)

  ## Bars carry almost no alpha and the step outline does the identifying.
  ## Four overlapping fills cannot encode which chain is which at any
  ## opacity, so the fill only marks where the mass is; the outline is what
  ## shows a chain sitting where the others do not.
  ggplot2::ggplot(long, ggplot2::aes(.data$value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(.data$density), fill = .data$chain),
      bins = bins, position = "identity", alpha = alpha, colour = NA) +
    ggplot2::geom_step(
      ggplot2::aes(y = ggplot2::after_stat(.data$density),
                   colour = .data$chain),
      stat = "bin", bins = bins, direction = "mid", linewidth = linewidth) +
    ggplot2::facet_wrap(~variable, ncol = ncol, scales = "free",
                        labeller = ggplot2::labeller(variable = lab)) +
    ggplot2::scale_fill_viridis_d(option = "viridis", end = .85,
                                  name = "chain") +
    ggplot2::scale_colour_viridis_d(option = "viridis", end = .85,
                                    name = "chain") +
    .gc_theme() +
    ggplot2::labs(x = NULL, y = "density", title = "Posterior, by chain")
}


#' Corner plot of the shared parameters
#'
#' The lower triangle of pairwise posteriors, with Pearson correlations. Use
#' it to tell a genuine ridge -- two parameters the data cannot separate --
#' from one parameter simply absorbing signal from another.
#'
#' @param fit A [GCFit].
#' @param vars Variables to pair. Keep the list short; the panel count is
#'   quadratic. `NULL`, the default, takes the shared parameters.
#' @param n Draws to plot, thinned evenly.
#' @return A ggplot.
#' @export
gc_plot_pairs <- function(fit, vars = NULL, n = 2000L) {
  .need_ggplot()
  if (is.null(vars)) vars <- c("u", "nu", .neff_var(fit$control), "sigma")
  vars <- .expand_vars(fit, vars)
  if (length(vars) < 2L) stop("need at least two variables", call. = FALSE)
  m <- posterior::as_draws_matrix(
    .exact_draws(fit$draws, vars))
  idx <- if (nrow(m) > n) unique(round(seq(1, nrow(m), length.out = n)))
         else seq_len(nrow(m))
  ## a draws_matrix keeps its class and its dimensions through `[`, so a
  ## single-column subset is still a matrix; strip it to a plain one
  M <- matrix(as.numeric(m[idx, vars, drop = FALSE]), nrow = length(idx),
              dimnames = list(NULL, vars))

  ## The last variable never appears on x and the first never on y, so
  ## their levels are dropped: otherwise the grid carries an empty row and
  ## an empty column.
  xl <- vars[-length(vars)]
  yl <- rev(vars[-1L])
  pairs <- utils::combn(seq_along(vars), 2L)
  cell <- function(i, j, extra) {
    data.frame(xv = factor(vars[i], levels = xl),
               yv = factor(vars[j], levels = yl), extra)
  }
  df <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1L, k]; j <- pairs[2L, k]     # i before j: x = i, y = j
    cell(i, j, data.frame(x = as.numeric(M[, i]), y = as.numeric(M[, j])))
  }))
  cors <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1L, k]; j <- pairs[2L, k]
    cell(i, j, data.frame(r = stats::cor(M[, i], M[, j])))
  }))

  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
    ## a 2-D histogram rather than points with contours: it reads better at
    ## several thousand draws, and the contour stat warns in the empty
    ## upper-triangle panels that make this a corner plot
    ggplot2::geom_bin2d(bins = 40L) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1,
                                  trans = "sqrt", guide = "none") +
    ggplot2::geom_label(data = cors, ggplot2::aes(label = sprintf("r = %+.2f",
                                                                  .data$r)),
                        x = -Inf, y = Inf, hjust = -0.08, vjust = 1.2,
                        size = 2.6, label.size = 0, alpha = .75,
                        inherit.aes = FALSE) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$yv),
                        cols = ggplot2::vars(.data$xv),
                        scales = "free", switch = "y", drop = FALSE) +
    .gc_theme() +
    ggplot2::theme(legend.position = "none",
                   strip.placement = "outside",
                   panel.spacing = ggplot2::unit(2, "pt")) +
    ggplot2::labs(x = NULL, y = NULL,
      title = "Shared parameters, pairwise")
}

#' @describeIn gc_plot_frequency Default plot for a projection.
#' @param x A [GCProjection].
#' @param ... Passed to [gc_plot_frequency].
#' @export
plot.GCProjection <- function(x, ...) gc_plot_frequency(x, ...)

#' @describeIn gc_plot_trace Default plot for a fit.
#' @param x A [GCFit].
#' @param ... Passed to [gc_plot_trace].
#' @export
plot.GCFit <- function(x, ...) gc_plot_trace(x, ...)
