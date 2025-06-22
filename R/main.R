
# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("pop", "n_samp", "p", "n_pos"))
}

#------------------------------------------------
#' @title Simulate from a Wright-Fisher model
#'
#' @description Simulates from Wright-Fisher model with drift and selection.
#'
#' @param s TODO.
#' @param p0 TODO.
#' @param N TODO.
#' @param t_max TODO.
#' @param n_deme TODO.
#'
#' @importFrom stats rbinom
#' @export

sim_WF <- function(s, p0, N, t_max, n_deme) {

  # TODO - checks

  # initialize matrix for storing allele frequencies. Rows are
  # generations, columns are demes
  p_mat <- matrix(0, nrow = t_max, ncol = n_deme)
  p_mat[1,] <- round(p0*N) / N

  # run WF simulation
  for (t in 2:t_max) {
    p <- p_mat[t-1,]

    # apply selection
    p_sel <- p*(1 + s) / (p*(1 + s) + (1 - p))

    # apply drift
    p_mat[t,] <- rbinom(n_deme, N, prob = p_sel) / N
  }

  # convert to long data.frame
  ret <- data.frame(pop = rep(1:n_deme, each = t_max),
                    t = 1:t_max,
                    p = as.vector(p_mat))

  return(ret)
}

#------------------------------------------------
# produces transition matrix for diffusion model given mean and SD
#' @importFrom stats pnorm
#' @noRd
get_trans_mat <- function(s, sigma, dx = 0.01) {

  # generate daily transition matrix
  nx <- 1 / dx + 1
  x <- seq(0, 1, dx)
  trans_mat <- matrix(0, nx, nx)
  trans_mat[1,1] <- trans_mat[nx,nx] <- 1

  for (i in 2:(nx - 1)) {

    # mean and scaled SD of diffusion
    mu <- x[i] + sigma*x[i]*(1 - x[i])
    tau <- sqrt(sigma^2*x[i]*(1 - x[i]))

    # transition prob for intermediate frequencies
    trans_mat[i,] <- pnorm(x + dx/2, mu, tau) - pnorm(x - dx/2, mu, tau)

    # special case at boundaries
    trans_mat[i,1] <- pnorm(x[1] + dx/2, mu, sigma)
    trans_mat[i,nx] <- pnorm(x[nx] - dx/2, mu, sigma, lower.tail = FALSE)
  }

  return(trans_mat)
}

#------------------------------------------------
#' @title Simulate from a diffusion approximation to the Wright-Fisher model
#'
#' @description Simulates from diffusion approximation to Wright-Fisher model.
#'
#' @inheritParams sim_WF
#' @param sigma TODO.
#' @param dx TODO.
#'
#' @export

sim_diffusion <- function(s, p0, sigma, dx = 0.01, t_max, n_deme) {

  # generate transition matrix
  trans_mat <- get_trans_mat(s = s, sigma = sigma, dx = dx)
  nx <- 1 / dx + 1

  # find which frequency bin contains the initial value and initialize
  x_init <- round(p0 / dx) + 1
  p_mat <- matrix(x_init, t_max, n_deme)

  # simulate trajectories by drawing from transition matrix
  for (j in 1:n_deme) {
    for (i in 2:t_max) {
      probs <- trans_mat[p_mat[i-1,j],]
      p_mat[i,j] <- sample.int(nx, 1, prob = probs)
    }
  }

  # convert to long data.frame
  ret <- data.frame(pop = rep(1:n_deme, each = t_max),
                    t = 1:t_max,
                    p = as.vector(p_mat) / nx)

  return(ret)
}

#------------------------------------------------
#' @title Sample from a simulated prevalence trajectory
#'
#' @description TODO
#'
#' @param df_sim TODO.
#' @param df_sample TODO.
#'
#' @importFrom epitools binom.exact
#' @importFrom dplyr inner_join mutate join_by n
#' @export

sample_prev <- function(df_sim, df_sample) {

  # draw from binomial and calculate 95% CIs
  ret <- df_sim |>
    inner_join(df_sample,
               by = join_by(pop, t)) |>
    mutate(n_pos = rbinom(n = n(), size = n_samp, prob = p),
           p_est = n_pos / n_samp,
           CI_lower = epitools::binom.exact(n = n_samp, x = n_pos)$lower,
           CI_upper = epitools::binom.exact(n = n_samp, x = n_pos)$upper)

  return(ret)
}

#------------------------------------------------
#' @title Logistic function
#'
#' @description TODO
#'
#' @inheritParams sim_WF
#' @param t TODO.
#'
#' @export

get_logistic <- function(s, p0, t) {
  1 / (1 + (1 - p0)/p0*exp(-s*t))
}
