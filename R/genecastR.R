#' genecastR: Inference of Selection from Time-Series Allele Frequencies
#'
#' *genecastR* provides tools to model changes in allele prevalence over time
#' from genetic surveillance data and to estimate the strength of selection
#' while accounting for stochastic variation due to genetic drift.
#' The package focuses on infectious-disease applications (e.g., malaria) but is
#' broadly applicable to any system with time-series allele-frequency data.
#'
#' @section Key features:
#' \itemize{
#'   \item \strong{Bayesian inference of selection:} Estimate selection
#'         coefficients \eqn{s} (and diffusion/noise parameters) using a
#'         Wright–Fisher diffusion approximation with selection and drift via MCMC
#'         (see \code{\link{run_mcmc}}).
#'   \item \strong{Data visualization and summaries:} Quick plotting helpers for
#'         observed prevalence with binomial confidence intervals
#'         (see \code{\link{plot_prev}} and \code{\link{get_CIs}}).
#'   \item \strong{Example data:} Included toy dataset for getting started
#'         (see \code{\link{example_data}}).
#' }
#'
#' @section Basic workflow:
#' \preformatted{
#'   library(genecastR)
#'
#'   # Explore example data
#'
#'   # Visualize prevalence over time
#'   plot_prev(df_data = example_data)
#'
#'   # Estimate selection parameters via MCMC
#'   fit <- run_mcmc(df_data = example_data,
#'                   burnin = 1e2,
#'                   samples = 1e3)
#' }
#'
#' @section Getting help:
#' See function reference pages (e.g., \code{\link{run_mcmc}},
#' \code{\link{plot_prev}}, \code{\link{get_CIs}}) and the README for a
#' step-by-step introduction and installation notes.
#'
#' @keywords package
#' @aliases genecastR-package
"_PACKAGE"
