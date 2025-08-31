
#------------------------------------------------
#' @title Plot posterior prevalence and observed data
#'
#' @description
#' Visualizes posterior estimates of allele prevalence over time, optionally overlaying
#' simulated trajectories and observed data with confidence intervals. Posterior distributions
#' are shown as quantile bands computed from the posterior probability mass at each time point.
#'
#' @param df_post Optional data frame of posterior prevalence distributions as returned by
#' \code{\link{get_posterior_prev}}. Must include columns \code{pop}, \code{t}, \code{p}, and \code{post}.
#' @param df_trajectory Optional data frame of simulated prevalence trajectories (e.g., from \code{sim_WF} or \code{sim_diffusion}),
#' with columns \code{pop}, \code{t}, and \code{p}.
#' @param df_data Optional data frame of observed prevalence data with columns \code{pop}, \code{t},
#' \code{p_est} (estimated prevalence), and confidence interval bounds \code{CI_lower} and \code{CI_upper}.
#' @param CI_type If plotting data, which type of confidence interval to produce (see \code{?get_CIs()}).
#'
#' @return A \code{ggplot} object visualizing the specified components:
#' \itemize{
#'   \item Posterior prevalence quantile bands (if \code{df_post} is provided)
#'   \item Simulated allele frequency trajectories (if \code{df_trajectory} is provided)
#'   \item Observed data with confidence intervals (if \code{df_data} is provided)
#' }
#'
#' @details
#' The plot includes posterior quantile shading based on the cumulative posterior mass, grouped into bands
#' (e.g. 0–20\%, 20–50\%, etc.). Multiple populations are shown using \code{facet_wrap()}.
#' At least one of the three arguments must be supplied.
#'
#' @importFrom ggplot2 ggplot geom_raster geom_line geom_pointrange theme_bw aes facet_wrap scale_fill_brewer scale_y_continuous geom_tile scale_fill_manual
#' @importFrom dplyr group_by arrange mutate filter
#' @export

plot_prev <- function(df_post = NULL, df_trajectory = NULL, df_data = NULL, CI_type = "CP") {

  if (is.null(df_post) & is.null(df_trajectory) & is.null(df_data)) {
    stop("You must specify at least one of the three optional inputs")
  }

  # create initial plotting object
  ret <- ggplot() + theme_bw()

  # plot posterior prevalence
  if (!is.null(df_post)) {

    # funny business to make sure all fill levels are visible in legend
    quant_levels <- c("95%", "80%", "50%", "20%")
    quant_colors <- RColorBrewer::brewer.pal(4, "Blues")

    ret <- ret + geom_tile(aes(x = x, y = y, fill = quant_group), show.legend = TRUE, alpha = 0,
              data = data.frame(x = 1:4, y = 1, quant_group = factor(quant_levels, levels = quant_levels))) +
      scale_fill_manual(values = quant_colors, name = "Credible\nInterval", drop = FALSE,
                        guide = guide_legend(override.aes = list(alpha = 1)))

    # more funny business to expand limits so no part of raster gets cut out
    dx <- df_post |>
      filter(pop == df_post$pop[1]) |>
      filter(week == df_post$week[1]) |>
      pull(p) |>
      diff() |>
      min()
    ret <- ret + scale_y_continuous(limits = c(-dx, 1 + dx)*100, expand = c(0, 0))

    # convert posterior probabilities to quantile groups
    df_quant_group <- df_post |>
      group_by(pop, week) |>
      arrange(-post) |>
      mutate(post_above = cumsum(post),
             quant_group = cut(post_above, breaks = c(0, 0.2, 0.5, 0.8, 0.95)),
             quant_group = 5 - as.numeric(quant_group),
             quant_group = quant_levels[quant_group],
             quant_group = factor(quant_group, levels = quant_levels)) |>
      filter(!is.na(quant_group))

    ret <- ret +
      geom_raster(aes(x = week, y = p*100, fill = quant_group), data = df_quant_group)
  } else {
    ret <- ret + scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
  }

  # overlay trajectories
  if (!is.null(df_trajectory)) {
    ret <- ret +
      geom_line(aes(x = week, y = p*100), col = "firebrick2", data = df_trajectory)
  }

  # overlay raw data CIs
  if (!is.null(df_data)) {
    ret <- ret +
      geom_pointrange(aes(x = week, y = p_est*100, ymin = CI_lower*100, ymax = CI_upper*100),
                      size = 0.2, data = get_CIs(df_data, CI_type = CI_type))
  }

  # final aesthetics
  ret <- ret +
    xlab("Time (weeks)") + ylab("Allele Prevalence (%)") +
    facet_wrap(~pop)

  return(ret)
}

#------------------------------------------------
#' @title Scatter plot of posterior MCMC samples
#'
#' @description
#' Generates a scatter plot of posterior samples for the selection coefficient (\code{s})
#' and diffusion rate (\code{sigma}) from a completed MCMC run. This visualization helps
#' assess the joint posterior distribution and any correlation between parameters.
#'
#' @param mcmc An object returned by \code{run_mcmc()}, containing the MCMC output
#' including both burn-in and sampling phases.
#'
#' @return A \code{ggplot} object showing a scatter plot of posterior samples with \code{s}
#' on the x-axis and \code{sigma} on the y-axis.
#'
#' @importFrom ggplot2 ggplot geom_point theme_bw xlab ylab aes
#' @importFrom dplyr filter
#' @importFrom RColorBrewer brewer.pal
#' @export

plot_mcmc_scatter <- function(mcmc) {

  mcmc$output |>
    filter(phase == "sampling") |>
    ggplot() + theme_bw() +
    geom_point(aes(x = s, y = sigma), col = brewer.pal(n = 3, name = "Reds")[3], alpha = 0.5) +
    xlab("Selection coefficient (s)") + ylab("Diffusion rate (sigma)")
}

#------------------------------------------------
#' @title Density plot of posterior MCMC samples
#'
#' @description
#' Generates a filled 2D density plot of posterior samples for the selection coefficient (\code{s})
#' and diffusion rate (\code{sigma}) from a completed MCMC run. The contours represent different
#' quantile levels of the joint posterior distribution, allowing for visualization of uncertainty
#' and correlation between parameters.
#'
#' @param mcmc An object returned by \code{run_mcmc()}, containing the MCMC output
#' including both burn-in and sampling phases.
#'
#' @return A \code{ggplot} object showing a filled contour plot of the posterior density over
#' the parameter space of \code{s} and \code{sigma}.
#'
#' @importFrom ggplot2 ggplot geom_contour_filled theme_bw xlab ylab aes scale_fill_brewer guides guide_legend
#' @importFrom dplyr filter
#' @importFrom MASS kde2d
#' @export

plot_mcmc_density <- function(mcmc) {

  # extract sampling phase
  df_draws <- mcmc$output |>
    filter(phase == "sampling")

  # compute 2d kde
  kde <- MASS::kde2d(x = df_draws$s, y = df_draws$sigma, n = 64)

  # produce contour plot
  data.frame(x = kde$x, y = rep(kde$y, each = 64), z = as.vector(kde$z)) |>
    mutate(z = 1 - z / max(z)) |>
    ggplot() + theme_bw() +
    geom_contour_filled(aes(x = x, y = y, z = z), breaks = seq(0, 1, 0.1)) +
    scale_fill_manual(values = rev(c('white', brewer.pal(n = 9, name = "Reds"))),
                      labels = sprintf("%s - %s%%", seq(0, 90, 10), seq(10, 100, 10)),
                      name = "Posterior quantile") +
    xlab("Selection coefficient (s)") + ylab("Diffusion rate (sigma)")

}
