#' @keywords internal
#' @section Pipeline:
#' Three objects, in order. Each is documented on its own page.
#'
#' \describe{
#'   \item{[GCData]}{validates a survey panel onto the model's lattice.
#'     One row per survey: `deme`, `time`, `x`, `n`.}
#'   \item{[GCFit]}{runs blocked Metropolis-within-Gibbs against the
#'     **exact** likelihood. Because the likelihood is evaluated rather than
#'     estimated, this is ordinary Metropolis-within-Gibbs and not a
#'     pseudo-marginal scheme.}
#'   \item{[GCProjection]}{smoothed curves, projections, presence and
#'     exceedance probabilities and the introduction times, as one row per
#'     deme and step. Built by [GCFit]`$project()`.}
#' }
#'
#' Plotting: [gc_plot_frequency], [gc_plot_present], [gc_plot_exceed],
#' [gc_plot_times], [gc_plot_joint], [gc_plot_trace].
#'
#' @section The three introduction times:
#' `t*` is the founding parameter, obtained analytically from a backward row
#' the likelihood already builds, so it carries no Monte-Carlo error.
#' `T_first` is the first step the variant is present, which under the
#' default founding law is identically `t*`. `T_est` is the start of the
#' excursion that produced the first detected count -- walk back from the
#' first non-zero survey until the deme is empty -- so it is bounded by that
#' detection and is undefined where there was none. It is a path functional,
#' needs sampled paths, and is the only expensive output. Their joint
#' distribution has the probability that the detected lineage is the
#' original arrival on its diagonal.
#'
#' @section Caveats worth reading first:
#' `nu` is an arrival rate *conditional on the deme being empty*, not a
#' rate of gene flow. `t*` has no prior support beyond the observation window,
#' so projected presence is biased low in never-detected demes. `nu` and `A`
#' are shared across demes. Nothing here corrects for the ascertainment that
#' selects which demes enter a panel.
"_PACKAGE"

## usethis namespace: start
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel mclapply makePSOCKcluster stopCluster clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom posterior as_draws_array
#' @importFrom rlang .data
#' @importFrom stats qbeta
#' @importFrom utils read.csv
#' @useDynLib GeneCastR, .registration = TRUE
## usethis namespace: end
NULL
