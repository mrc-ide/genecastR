#' Example allele prevalence data
#'
#' A small dataset simulating observed allele frequencies across time for multiple populations.
#' Useful for testing or demonstrating functions such as `estimate_ML()` or `run_mcmc()`.
#'
#' @format A data frame with one row per sampling point and the following columns:
#' \describe{
#'   \item{pop}{Population identifier (integer).}
#'   \item{week}{Sampling week (integer).}
#'   \item{n_samp}{Number of individuals sampled at that time point.}
#'   \item{n_pos}{Number of individuals carrying the allele.}
#' }
#'
#' @examples
#' head(example_data)
#'
"example_data"
