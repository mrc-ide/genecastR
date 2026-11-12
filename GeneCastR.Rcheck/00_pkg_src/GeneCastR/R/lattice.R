#' Lattice size from a log effective size
#'
#' The model is parameterised by `N_eff` per unit time; the lattice holds
#' `N_lat = M * N_eff` gene copies, and it must be a whole number, so the
#' likelihood is piecewise constant in `z`. Rounding is half-to-even, matching
#' the reference implementation; `floor(x + 0.5)` and truncation are each
#' self-consistent but give different lattices, and therefore different
#' likelihoods, at exact half-integers.
#'
#' @param z Log effective size, `log(N_eff)`.
#' @param M Lattice steps per unit time.
#' @return Integer lattice sizes, floored at 2. A lattice of 0 or 1 has no
#'   interior, so there is nothing for selection or drift to act on.
#' @export
gc_nlat <- function(z, M) {
  pmax(2L, as.integer(round(M * exp(z))))
}

#' Per-step selection coefficient from the per-unit-time log growth advantage
#'
#' With `u = log(1 + A)` and `A` the advantage per unit time, a step of
#' `1/M` carries `alpha = exp(u / M) - 1`.
#'
#' @param u Log growth advantage per unit time.
#' @param M Lattice steps per unit time.
#' @return Per-step selection coefficient.
#' @export
gc_alpha <- function(u, M) {
  exp(u / M) - 1
}

#' Per-step Poisson mean from the per-unit-time influx frequency
#'
#' `nu` is an influx in frequency units per unit time, so a deme of size
#' `N_eff` receives `lambda = nu * N_eff` copies per step. Holding this fixed
#' is what keeps the rate of established introductions close to invariant in
#' the resolution `M`: a single copy establishes with probability of order
#' `2 * alpha`, which falls as `1 / M`, while the arrival count per unit time
#' rises as `M`.
#'
#' @param nu Influx in frequency units per unit time.
#' @param neff Effective size, `exp(z)`.
#' @return Per-step Poisson mean, one per deme.
#' @export
gc_lambda <- function(nu, neff) {
  nu * neff
}
