
# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("pop", "n_samp", "p", "n_pos", "phase", "s", "p0", "p_logistic", "ll",
                           "CI_upper", "CI_lower", "MOE", "prob", "prob_cum", "w", "x", "p_est",
                           "z_stat", "is_signif", "post_above", "quant_group", "post", "sigma",
                           "week", "y"))
}

#------------------------------------------------
#' @title Simulate from a Wright-Fisher model
#'
#' @description
#' Simulates allele frequency trajectories over time using the Wright-Fisher model with selection and genetic drift.
#' Each trajectory represents a population evolving independently under the same parameters.
#'
#' @param s Selection coefficient. A positive value favours the allele, negative values represent deleterious effects.
#' @param p0 Initial allele frequency (between 0 and 1).
#' @param N Effective population size (number of diploid individuals).
#' @param week_max Number of generations to simulate, where a week is a single generation (valid under a scaled effective population size).
#' @param n_pop Number of independent populations (replicates) to simulate.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{pop}{Population (replicate) ID, from 1 to \code{n_pop}.}
#'   \item{week}{Week number, from 1 to \code{max_weeks}.}
#'   \item{p}{Allele frequency in that population and generation.}
#' }
#'
#' @importFrom stats rbinom
#' @export

sim_WF <- function(s, p0, N, week_max, n_pop) {

  # input checks
  assert_single_numeric(s)
  assert_single_bounded(p0)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos_int(week_max, zero_allowed = FALSE)
  assert_single_pos_int(n_pop, zero_allowed = FALSE)

  # initialize matrix for storing allele frequencies. Rows are
  # generations, columns are populations (demes)
  p_mat <- matrix(0, nrow = week_max, ncol = n_pop)
  p_mat[1,] <- round(p0*N) / N

  # run WF simulation
  for (t in 2:week_max) {
    p <- p_mat[t-1,]

    # apply selection
    p_sel <- p*(1 + s) / (p*(1 + s) + (1 - p))

    # apply drift
    p_mat[t,] <- rbinom(n_pop, N, prob = p_sel) / N
  }

  # convert to long data.frame
  ret <- data.frame(pop = rep(1:n_pop, each = week_max),
                    week = 1:week_max,
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
    mu <- x[i] + s*x[i]*(1 - x[i])
    tau <- sqrt(sigma^2*x[i]*(1 - x[i]))

    # transition prob for intermediate frequencies
    trans_mat[i,] <- pnorm(x + dx/2, mu, tau) - pnorm(x - dx/2, mu, tau)

    # special case at boundaries
    trans_mat[i,1] <- pnorm(dx/2, mu, tau)
    trans_mat[i,nx] <- pnorm(1 - dx/2, mu, tau, lower.tail = FALSE)
  }

  return(trans_mat)
}

#------------------------------------------------
#' @title Simulate from a diffusion approximation to the Wright-Fisher model
#'
#' @description
#' Simulates allele frequency trajectories using a diffusion approximation to the Wright-Fisher model,
#' incorporating selection and drift. Allele frequencies evolve as a discrete approximation to a continuous diffusion process.
#'
#' @inheritParams sim_WF
#' @param sigma Diffusion standard deviation. Controls the strength of genetic drift.
#' @param dx Discretization step size for the allele frequency space (default is 0.01).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{pop}{Population (replicate) ID, from 1 to \code{n_deme}.}
#'   \item{week}{Week number, from 1 to \code{max_weeks}.}
#'   \item{p}{Allele frequency (discretized), between 0 and 1.}
#' }
#'
#' @export

sim_diffusion <- function(s, p0, sigma, dx = 0.01, week_max, n_pop) {

  # input checks
  assert_single_numeric(s)
  assert_single_bounded(p0)
  assert_single_pos(sigma, zero_allowed = FALSE)
  assert_single_bounded(dx, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(week_max, zero_allowed = FALSE)
  assert_single_pos_int(n_pop, zero_allowed = FALSE)

  # dx must generate a valid sequence
  if (round(1/dx) != 1/dx) {
    message(paste(c("dx must split the interval [0,1] into an integer number of sections. For example"),
                  "dx=0.01 is valid (100 sections), but dx=0.03 is not."))
  }

  # generate transition matrix
  trans_mat <- get_trans_mat(s = s, sigma = sigma, dx = dx)
  nx <- 1 / dx + 1

  # find which frequency bin contains the initial value and initialize
  x_init <- round(p0 / dx) + 1
  x_mat <- matrix(x_init, week_max, n_pop)

  # simulate trajectories by drawing from transition matrix
  for (j in 1:n_pop) {
    for (i in 2:week_max) {
      probs <- trans_mat[x_mat[i-1, j],]
      x_mat[i,j] <- sample.int(nx, 1, prob = probs)
    }
  }

  # convert to long data.frame
  ret <- data.frame(pop = rep(1:n_pop, each = week_max),
                    week = 1:week_max,
                    p = as.vector(x_mat - 1) * dx)

  return(ret)
}

#------------------------------------------------
#' @title Sample from a simulated prevalence trajectory
#'
#' @description
#' Samples observed allele prevalence from a simulated trajectory by drawing from a binomial distribution,
#' based on the true simulated allele frequencies and a specified number of sampled individuals.
#'
#' @param df_sim A data frame of simulated allele frequency trajectories, typically from \code{sim_WF()} or \code{sim_diffusion()}.
#' Must include columns \code{pop}, \code{week}, and \code{p} (true frequency).
#' @param df_data A data frame specifying sampling points, with columns \code{pop}, \code{week}, and \code{n_samp}
#' (number of individuals sampled per population and time point).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{pop}{Population ID.}
#'   \item{week}{Week of study.}
#'   \item{n_samp}{Number of individuals sampled at that time point.}
#'   \item{n_pos}{Number of positive (allele-carrying) individuals sampled (binomial draw).}
#' }
#'
#' @importFrom dplyr inner_join mutate join_by n
#' @export

sample_prev <- function(df_sim, df_data) {

  # input checks
  assert_dataframe(df_sim)
  assert_in(c("pop", "week", "p"), names(df_sim))
  assert_dataframe(df_data)
  assert_in(c("pop", "week", "n_samp"), names(df_data))

  # draw from binomial
  df_sim |>
    inner_join(df_data,
               by = join_by(pop, week)) |>
    mutate(n_pos = rbinom(n = n(), size = n_samp, prob = p)) |>
    select(-p)
}

#------------------------------------------------
#' @title Calculate 95\% binomial confidence intervals for allele prevalence
#'
#' @description
#' Given observed allele counts from sampled populations, calculates the estimated allele frequency
#' and 95\% confidence intervals.
#'
#' @param df_data A data frame with columns:
#' \describe{
#'   \item{n_pos}{Number of individuals carrying the allele.}
#'   \item{n_samp}{Total number of individuals sampled.}
#' }
#' @param CI_type Which type of confidence interval to produce. \code{CP} =
#'   Clopper-Pearson interval (default), \code{Wald} = Wald interval. The Clopper-Pearson
#'   interval tends to be wider and more robust, especially for extreme values.
#'
#' @return A data frame identical to the input \code{df_data}, with three additional columns:
#' \describe{
#'   \item{p_est}{Estimated allele frequency (\code{n_pos / n_samp}).}
#'   \item{CI_lower}{Lower bound of the 95\% confidence interval.}
#'   \item{CI_upper}{Upper bound of the 95\% confidence interval.}
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom epitools binom.exact
#' @export

get_CIs <- function(df_data, CI_type = "CP") {

  # input checks
  assert_dataframe(df_data)
  assert_in(c("n_pos", "n_samp"), names(df_data))
  assert_in(CI_type, c("CP", "Wald"))

  if (CI_type == "Wald") {
    ret <- df_data |>
      mutate(p_est = n_pos / n_samp,
             MOE = 1.96 * p_est*(1 - p_est) / sqrt(n_samp),
             CI_lower = p_est - MOE,
             CI_upper = p_est + MOE) |>
      select(-MOE)
  } else if (CI_type == "CP") {
    ret <- df_data |>
      mutate(p_est = n_pos / n_samp,
             CI_lower = epitools::binom.exact(n = n_samp, x = n_pos)$lower,
             CI_upper = epitools::binom.exact(n = n_samp, x = n_pos)$upper)
  }
  return(ret)
}

#------------------------------------------------
#' @title Logistic function
#'
#' @description
#' Computes the deterministic logistic trajectory of allele frequency over time under selection,
#' assuming no genetic drift. This models the expected change in allele frequency due to selection alone.
#'
#' @inheritParams sim_WF
#' @param week A numeric vector of time points (weeks) at which to evaluate the trajectory.
#'
#' @return A numeric vector of allele frequencies at each time point.
#'
#' @export

get_logistic <- function(s, p0, week) {

  # input checks
  assert_vector_numeric(s)
  assert_vector_bounded(p0)
  assert_vector_pos_int(week, zero_allowed = FALSE)

  1 / (1 + (1 - p0)/p0*exp(-s*week))
}

#------------------------------------------------
#' @title Estimate parameters via maximum likelihood
#'
#' @description
#' Estimates the selection coefficient (\code{s}) and initial allele frequency (\code{p0})
#' using maximum likelihood. A grid search is performed over a parameter space, and the
#' log-likelihood is computed based on a deterministic logistic model. Confidence intervals
#' for \code{s} are calculated using a profile likelihood approach. If multiple populations
#' are provided, estimates are computed independently for each one.
#'
#' @param df_data A data frame with observed prevalence data. Must include columns:
#' \code{pop} (population ID), \code{week} (sampling week), \code{n_samp} (sample size), and
#' \code{n_pos} (number of individuals carrying the allele).
#' @param s_range A numeric vector of length 2 specifying the range of selection coefficients
#' to explore (default is \code{c(-0.1, 0.1)}).
#'
#' @return A data frame with one row per population, containing:
#' \describe{
#'   \item{pop}{Population ID.}
#'   \item{p0}{Maximum likelihood estimate of initial allele frequency.}
#'   \item{s}{Maximum likelihood estimate of selection coefficient.}
#'   \item{s_lower}{Lower bound of 95\% confidence interval for \code{s}.}
#'   \item{s_upper}{Upper bound of 95\% confidence interval for \code{s}.}
#' }
#'
#' @importFrom dplyr group_split mutate summarise group_by bind_rows pull filter
#' @importFrom tidyr expand_grid
#' @importFrom stats dbinom qchisq
#' @export

estimate_ML <- function(df_data, s_range = c(-0.1, 0.1)) {

  # input checks
  assert_dataframe(df_data)
  assert_in(c("pop", "week", "n_pos", "n_samp"), names(df_data))
  assert_limit(s_range)

  # create data.frame of parameter combinations to explore
  df_params <- expand_grid(s = seq(s_range[1], s_range[2], l = 101),
                           p0 = seq(0, 1, 0.01))

  # perform calculation separately for each deme
  ret <- df_data |>
    group_split(pop) |>
    lapply(function(df) {

      # calculate log-likelihood of every parameter combination
      df_ll <- df_params |>
        expand_grid(df) |>
        mutate(p_logistic = get_logistic(s = s, p0 = p0, week = week),
               ll = dbinom(x = n_pos, size = n_samp, prob = p_logistic, log = TRUE)) |>
        group_by(s, p0) |>
        summarise(ll = sum(ll), .groups = "drop")

      # calculate global ML estimates
      w <- which.max(df_ll$ll)
      s_ML <- df_params$s[w]
      p0_ML <- df_params$p0[w]

      # calculate 95% CI from profile log-likelihood
      s_CI <- df_ll |>
        group_by(s) |>
        summarise(ll = max(ll)) |>
        filter(ll >= (max(ll) - qchisq(0.95, df = 1) / 2)) |>
        pull(s) |>
        range()

      # return
      data.frame(pop = df$pop[1],
                 p0 = p0_ML,
                 s = s_ML,
                 s_lower = s_CI[1],
                 s_upper = s_CI[2])
    }) |>
    bind_rows()

  ret
}

#------------------------------------------------
# reparameterization of the lognormal distribution in terms of mean and SD
#' @importFrom stats dlnorm
#' @noRd
dlnorm_reparam <- function(x, mean, sd, return_log = FALSE) {
  sigma_sq <- log(sd^2 / mean^2 + 1)
  dlnorm(x, meanlog = log(mean) - sigma_sq / 2, sdlog = sqrt(sigma_sq), log = return_log)
}

#------------------------------------------------
# drjacoby loglikelihood function
#' @noRd
loglike_R <- function(params, data, misc) {

  # calculate log-likelihood over populations
  trans_mat <- get_trans_mat(s = params["s"], sigma = params["sigma"], dx = misc$dx)
  ll <- sapply(data$df, function(df) {
    get_forward_ll(trans_mat = trans_mat,
                   t_samp = df$week,
                   n_samp = df$n_samp,
                   n_pos = df$n_pos,
                   week_min = misc$week_min,
                   x_init = misc$x_init)
  })

  sum(ll)
}

#------------------------------------------------
# drjacoby logprior function
#' @noRd
logprior_R <- function(params, misc) {
  misc$prior_s(params["s"]) + misc$prior_sigma(params["sigma"])
}

#------------------------------------------------
#' @title Estimate parameters via Bayesian MCMC
#'
#' @description
#' Runs MCMC using the \code{drjacoby} package to estimate the selection coefficient (\code{s}) and diffusion parameter (\code{sigma})
#' under a diffusion approximation to the Wright-Fisher model. Allows for multiple populations (each assumed to follow independent trajectories)
#' and returns full posterior samples.
#'
#' @param df_data A data frame of observed prevalence data, with columns:
#' \code{pop} (population ID), \code{t} (generation), \code{n_samp} (sample size), and \code{n_pos} (number of positive observations).
#' @param dx Step size for discretizing allele frequencies (default is \code{0.01}).
#' @param week_min The starting week for likelihood calculation (default is \code{1}).
#' @param burnin Number of burn-in iterations per chain (default is \code{1e2}).
#' @param samples Number of post-burn-in MCMC samples per chain (default is \code{1e3}).
#' @param chains Number of independent MCMC chains to run (default is \code{1}).
#' @param prior_s A function specifying the log-prior for the selection coefficient \code{s}.
#' @param prior_sigma A function specifying the log-prior for the diffusion parameter \code{sigma}.
#' @param s_init Initial value for \code{s} in the MCMC (default is \code{0.01}).
#' @param sigma_init Initial value for \code{sigma} in the MCMC (default is \code{0.01}).
#' @param x_init Initial allele frequency distribution vector (default is \code{seq(0, 1, dx)}).
#' @param pb_markdown Logical; whether to show progress bars suitable for markdown output (default is \code{FALSE}).
#' @param silent Logical; whether to suppress console output during MCMC (default is \code{FALSE}).
#'
#' @return A list containing MCMC results as returned by \code{drjacoby::run_mcmc()}, including posterior samples for all parameters.
#'
#' @importFrom dplyr select group_split
#' @importFrom drjacoby run_mcmc
#' @export

run_mcmc <- function(df_data,
                     dx = 0.01,
                     week_min = 1,
                     burnin = 1e2,
                     samples = 1e3,
                     chains = 1,
                     prior_s = function(s) dlnorm_reparam(s, mean = 0.05, sd = 0.05, return_log = TRUE),
                     prior_sigma = function(sigma) dlnorm_reparam(sigma, mean = 0.1, sd = 0.2, return_log = TRUE),
                     s_init = 0.01,
                     sigma_init = 0.01,
                     x_init = seq(0, 1, dx),
                     pb_markdown = FALSE,
                     silent = FALSE) {

  # input checks
  assert_dataframe(df_data)
  assert_in(c("pop", "week", "n_samp", "n_pos"), names(df_data))
  assert_single_bounded(dx, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_int(week_min)
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_single_numeric(s_init)
  assert_single_pos(sigma_init, zero_allowed = FALSE)
  assert_vector_pos(x_init)
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)

  # dx must generate a valid sequence
  if (round(1/dx) != 1/dx) {
    message(paste(c("dx must split the interval [0,1] into an integer number of sections. For example"),
                  "dx=0.01 is valid (100 sections), but dx=0.03 is not."))
  }

  # split data by pop
  data_list <- df_data |>
    select(pop, week, n_samp, n_pos) |>
    group_split(pop, .keep = FALSE)

  # define parameters data.frame
  df_params <- rbind.data.frame(list(name = "s", min = -Inf, max = Inf, init = s_init),
                                list(name = "sigma", min = 0, max = Inf, init = sigma_init))

  # misc elements
  misc_list <- list(dx = dx,
                    week_min = week_min,
                    prior_s = prior_s,
                    prior_sigma = prior_sigma,
                    x_init = x_init / sum(x_init))

  # run MCMC
  mcmc <- drjacoby::run_mcmc(data = list(df = data_list),
                             df_params = df_params,
                             misc = misc_list,
                             loglike = loglike_R,
                             logprior = logprior_R,
                             burnin = burnin,
                             samples = samples,
                             chains = chains,
                             pb_markdown = pb_markdown,
                             silent = silent)

  return(mcmc)
}

#------------------------------------------------
#' @title Sample from the posterior distribution of a completed MCMC
#'
#' @description
#' Draws samples from the posterior distribution of parameters obtained from a completed MCMC run using the \code{drjacoby} package.
#' Returns a data frame with one row per posterior draw and columns for each sampled parameter.
#'
#' @param mcmc An object returned by \code{drjacoby::run_mcmc()}, containing MCMC output for model parameters.
#' @param n_draws Number of posterior samples to draw from the combined MCMC chains (default is \code{100}).
#'
#' @return A data frame with \code{n_draws} rows and columns:
#' \describe{
#'   \item{sample}{Sample index.}
#'   \item{s}{Sampled value of the selection coefficient.}
#'   \item{sigma}{Sampled value of the diffusion parameter.}
#' }
#'
#' @importFrom drjacoby sample_chains
#' @importFrom dplyr select
#' @export

sample_mcmc <- function(mcmc, n_draws = 1e2) {

  # input checks
  assert_class(mcmc, "drjacoby_output")
  assert_single_pos_int(n_draws, zero_allowed = FALSE)

  ret <- drjacoby::sample_chains(x = mcmc, sample_n = n_draws) |>
    select(sample, s, sigma)
  rownames(ret) <- NULL
  return(ret)
}

#------------------------------------------------
#' @title Get credible interval on s
#'
#' @description
#' Extracts the posterior samples for the selection coefficient \code{s} from a completed MCMC run
#' and returns the 95\% credible interval, marginalized over \code{sigma}. Uses posterior samples
#' from the sampling phase only.
#'
#' @param mcmc An object returned by \code{drjacoby::run_mcmc()}, containing MCMC output for model parameters.
#'
#' @return A data frame with one row and three columns:
#' \describe{
#'   \item{Q2.5}{Lower bound of the 95\% credible interval for \code{s}.}
#'   \item{Q50}{Posterior median of \code{s}.}
#'   \item{Q97.5}{Upper bound of the 95\% credible interval for \code{s}.}
#' }
#'
#' @importFrom dplyr filter pull
#' @importFrom stats quantile
#' @export

estimate_Bayesian <- function(mcmc) {

  # input checks
  assert_class(mcmc, "drjacoby_output")

  s_quant <- mcmc$output |>
    filter(phase == "sampling") |>
    pull(s) |>
    quantile(probs = c(0.025, 0.5, 0.975))
  names(s_quant) <- NULL

  data.frame(Q2.5 = s_quant[1],
             Q50 = s_quant[2],
             Q97.5 = s_quant[3])
}

#------------------------------------------------
# runs forward algorithm in a single population, returns integrated log-likelihood.
#' @noRd
get_forward_ll <- function(trans_mat, t_samp = numeric(), n_samp = numeric(), n_pos = numeric(),
                           week_min, x_init) {

  # quick exit if no data
  if (all(is.na(t_samp))) {
    return(0)
  }

  # set up allele prevalence bins
  nx <- nrow(trans_mat)
  x <- seq(0, 1, l = nx)
  z_forward <- x_init

  # run forward algorithm
  ret <- 0
  for (t in week_min:max(t_samp)) {

    # apply transition matrix
    if (t > week_min) {
      z_forward <- (z_forward %*% trans_mat)[1,]
    }

    # apply emission probability
    if (t %in% t_samp) {
      w <- which(t_samp == t)[1]
      emission <- dbinom(x = n_pos[w], size = n_samp[w], prob = x)
      z_forward <- z_forward * emission
      ret <- ret + log(sum(z_forward))
      if (is.infinite(ret)) {
        return(-1e300)
      }
      z_forward <- z_forward / sum(z_forward)
    }
  }

  return(ret)
}

#------------------------------------------------
# runs forward algorithm in a single population, returns forward matrix.
#' @noRd
get_forward_mat <- function(trans_mat, t_samp = numeric(), n_samp = numeric(), n_pos = numeric(),
                            t_min, t_max, x_init) {

  # set up allele prevalence bins
  nx <- nrow(trans_mat)
  x <- seq(0, 1, l = nx)

  # set up objects to store forward matrix
  nt <- t_max - t_min + 1
  z_forward <- matrix(0, nrow = nt, ncol = nx)
  z_forward[1,] <- x_init

  # run forward algorithm
  for (i in 1:nt) {

    # apply transition matrix
    if (i > 1) {
      z_forward[i,] <- (z_forward[i-1,] %*% trans_mat)[1,]
    }

    # apply emission probability
    t_i <- t_min + i - 1
    if (t_i %in% t_samp) {
      w <- which(t_samp == t_i)[1]
      emission <- dbinom(x = n_pos[w], size = n_samp[w], prob = x)
      z_forward[i,] <- z_forward[i,] * emission
    }

    # normalize
    z_forward[i,] <- z_forward[i,] / sum(z_forward[i,])
  }

  return(z_forward)
}

#------------------------------------------------
# runs backward algorithm in a single population, returns backward matrix.
#' @noRd
get_backward_mat <- function(trans_mat, t_samp = numeric(), n_samp = numeric(), n_pos = numeric(),
                             t_min, t_max) {

  # set up allele prevalence bins
  nx <- nrow(trans_mat)
  x <- seq(0, 1, l = nx)

  # set up objects to store backward matrix
  nt <- t_max - t_min + 1
  z_backward <- matrix(1, nrow = nt, ncol = nx)

  # run backward algorithm
  for (i in (nt - 1):1) {
    z_backward[i,] <- z_backward[i+1,]

    # apply emission probability to next step
    t <- t_min + i
    if (t %in% t_samp) {
      w <- which(t_samp == t)[1]
      emission <- dbinom(x = n_pos[w], size = n_samp[w], prob = x)
      z_backward[i,] <- z_backward[i,] * emission
    }

    # reverse transitions and normalize
    z_backward[i,] <- (trans_mat %*% z_backward[i, ])[,1]
    z_backward[i,] <- z_backward[i,] / sum(z_backward[i,])
  }

  return(z_backward)
}

#------------------------------------------------
#' @title Posterior prevalence from forward-backward algorithm
#'
#' @description
#' Runs the forward-backward algorithm to estimate the posterior distribution of allele prevalence
#' over time, given observed data and model parameters. Supports multiple populations and can marginalize over
#' uncertainty in parameters by averaging over a set of posterior samples.
#'
#' @param df_data A data frame of observed data with columns: \code{pop} (population ID),
#' \code{t} (generation), \code{n_samp} (sample size), and \code{n_pos} (number of positive individuals).
#' @param s A scalar or vector of selection coefficients. If a vector, results are marginalized over the values.
#' @param sigma A scalar or vector of diffusion parameters, corresponding to the values of \code{s}.
#' @param dx Step size for discretizing allele frequencies (default is \code{0.01}).
#' @param week_min The first week to include in the trajectory (default is \code{1}).
#' @param week_max The last week to include in the trajectory (default is \code{52}).
#' @param x_init A vector specifying the initial allele frequency distribution (default is \code{seq(0, 1, dx)}).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{pop}{Population ID.}
#'   \item{p}{Allele frequency bin (from \code{x_init}).}
#'   \item{week}{Time in weeks.}
#'   \item{post}{Posterior probability of allele frequency \code{p} at time \code{t}.}
#' }
#'
#' @importFrom dplyr mutate bind_rows group_split
#' @importFrom tidyr expand_grid
#' @export

get_posterior_prev <- function(df_data,
                               s,
                               sigma,
                               dx = 0.01,
                               week_min = 1,
                               week_max = 52,
                               x_init = seq(0, 1, dx)) {

  # input checks
  assert_dataframe(df_data)
  assert_in(c("pop", "week", "n_samp", "n_pos"), names(df_data))
  assert_vector_numeric(s)
  assert_vector_pos(sigma, zero_allowed = FALSE)
  assert_single_bounded(dx, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_int(week_min)
  assert_single_int(week_max)
  assert_gr(week_max, week_min)

  # dx must generate a valid sequence
  if (round(1/dx) != 1/dx) {
    message(paste(c("dx must split the interval [0,1] into an integer number of sections. For example"),
                  "dx=0.01 is valid (100 sections), but dx=0.03 is not."))
  }

  # loop over pops
  df_data_split <- group_split(df_data, pop)
  df_post <- mapply(function(df, pop) {

    message(sprintf("Population %s:", pop))
    pb <- progress_bar$new(
      format = "  computing [:bar] :percent eta: :eta",
      total = length(s), clear = FALSE, width = 60
    )

    # run forward-backward algorithm for each value of model parameters
    post_mat <- 0
    for (i in seq_along(s)) {
      trans_mat <- get_trans_mat(s = s[i], sigma = sigma[i], dx = dx)
      forward_mat <- get_forward_mat(trans_mat = trans_mat,
                                     t_samp = df$week,
                                     n_samp = df$n_samp,
                                     n_pos = df$n_pos,
                                     t_min = week_min,
                                     t_max = week_max,
                                     x_init = x_init)
      backward_mat <- get_backward_mat(trans_mat = trans_mat,
                                       t_samp = df$week,
                                       n_samp = df$n_samp,
                                       n_pos = df$n_pos,
                                       t_min = week_min,
                                       t_max = week_max)
      tmp <- forward_mat*backward_mat
      post_mat <- post_mat + tmp / rowSums(tmp)

      pb$tick()
    }
    post_mat <- post_mat / rowSums(post_mat)

    # format
    expand_grid(pop = df$pop[1],
                p = seq(0, 1, dx),
                week = week_min:week_max) |>
      mutate(post = as.vector(post_mat))
  }, df_data_split, 1:length(df_data_split), SIMPLIFY = FALSE) |>
    bind_rows()

  return(df_post)
}

#------------------------------------------------
# calculates expected number of positive samples at a given time, integrated
# over the predicted prevalence.
#' @noRd
project_num_pos <- function(pop, week, n_samp, df_post) {

  # rename args so dplyr can filter with same name
  pop_ <- pop
  week_ <- week

  ret <- df_post |>
    filter(pop == pop_ & week == week_) |>
    expand_grid(n_pos = 0:n_samp) |>
    mutate(prob = dbinom(n_pos, n_samp, p)) |>
    group_by(n_pos) |>
    summarise(prob = sum(prob*post),
              .groups = "drop")

  if (nrow(ret) == 0) {
    stop(sprintf("No results corresponding to pop = %s, week = %s", pop, week))
  }

  return(ret)
}

#------------------------------------------------
#' @title Estimate expected margin of error for a future study
#'
#' @description
#' Calculates the expected margin of error (MOE) for a future prevalence study based on the posterior
#' distribution of allele frequencies. Returns the smallest MOE such that there is at least an 80\%
#' chance the observed confidence interval will fall within that margin. Assumes sampling is performed
#' at a single time point.
#'
#' @param pop Integer or character ID specifying the population for which to estimate margin of error.
#' @param week A single time point (week) at which sampling will occur.
#' @param n_samp A single scalar specifying the number of individuals to be sampled.
#' @param df_post A data frame containing the posterior distribution over allele frequencies, as returned
#' by \code{\link{get_posterior_prev}}. Must include columns \code{pop}, \code{week}, \code{p}, and \code{post}.
#' @param CI_type Which type of MOE to use (see \code{?get_CIs()}). If Clopper-Pearson then MOE is taken as the average of the upper and lower CIs, which may be different.
#'
#' @return A single numeric value: the estimated margin of error (half-width of the 95\% confidence interval)
#' that is expected to contain the true prevalence with at least 80\% probability.
#'
#' @importFrom epitools binom.exact
#' @importFrom dplyr mutate summarise arrange pull
#' @export

get_MOE <- function(pop, week, n_samp, df_post, CI_type = "CP") {

  # input checks
  assert_single_pos_int(pop, zero_allowed = FALSE)
  assert_single_int(week)
  assert_single_pos_int(n_samp, zero_allowed = FALSE)
  assert_dataframe(df_post)
  assert_in(c("pop", "p", "week", "post"), names(df_post))
  assert_in(CI_type, c("CP", "Wald"))

  if (CI_type == "CP") {
    df_MOE <- project_num_pos(pop = pop, week = week, n_samp = n_samp, df_post = df_post) |>
      mutate(CI_lower = binom.exact(x = n_pos, n = n_samp)$lower,
             CI_upper = binom.exact(x = n_pos, n = n_samp)$upper,
             MOE = (CI_upper - CI_lower) / 2)
  } else if (CI_type == "Wald") {
    df_MOE <- project_num_pos(pop = pop, week = week, n_samp = n_samp, df_post = df_post) |>
      mutate(p_est = n_pos / n_samp,
             MOE = 1.96*p_est*(1 - p_est) / sqrt(n_samp))
  }

  df_MOE |>
    arrange(MOE) |>
    mutate(prob_cum = cumsum(prob)) |>
    summarise(w = which(prob_cum > 0.8)[1],
              MOE = MOE[w]) |>
    pull(MOE)
}

#------------------------------------------------
#' @title Get minimum sample size based on margin of error
#'
#' @description
#' Calculates the minimum sample size required to ensure at least an 80\% chance that the margin of error (MOE)
#' of a future prevalence estimate falls below a specified threshold. Like \code{\link{get_MOE}}, this function
#' accounts for uncertainty in allele frequency using the posterior distribution. Sampling is assumed to occur
#' at a single time point.
#'
#' @param pop Population ID (integer or character) for which to calculate the required sample size.
#' @param week A single time point (week) at which sampling is planned.
#' @param MOE Desired margin of error (half-width of the 95\% confidence interval).
#' @param n_min Minimum sample size to search (default is \code{10}).
#' @param n_max Maximum sample size to search (default is \code{1000}).
#' @param df_post A data frame of posterior allele frequencies as returned by \code{\link{get_posterior_prev}}.
#' Must include columns \code{pop}, \code{t}, \code{p}, and \code{post}.
#' @param CI_type Which type of MOE to use (see \code{?get_MOE()}).
#'
#' @return A named numeric vector with the minimum required sample size, or \code{NULL} if no sample size within
#' the range \code{n_min} to \code{n_max} satisfies the target margin of error.
#'
#' @importFrom dplyr mutate summarise arrange pull
#' @importFrom progress progress_bar
#' @export

get_sample_size_MOE <- function(pop, week, MOE, n_min = 10, n_max = 1e3, df_post, CI_type = "CP") {

  # input checks
  assert_single_pos_int(pop, zero_allowed = FALSE)
  assert_single_int(week)
  assert_single_bounded(MOE, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(n_min, zero_allowed = FALSE)
  assert_single_pos_int(n_max)
  assert_gr(n_max, n_min)
  assert_dataframe(df_post)
  assert_in(c("pop", "p", "week", "post"), names(df_post))
  assert_in(CI_type, c("CP", "Wald"))

  pb <- progress_bar$new(
    format = "  computing [:bar] :percent eta: :eta",
    total = n_max, clear = FALSE, width = 60
  )

  # calculate MOE for all sample sizes from 1 to n_max
  for (i in n_min:n_max) {
    MOE_i <- get_MOE(pop = pop, week = week, n_samp = i, df_post = df_post, CI_type = CI_type)
    if (MOE_i < MOE) {
      return(c(sample_size = i))
    }
    pb$tick()
  }

  # if reached this point then have failed to find sample size
  message(sprintf(paste("Could not find a sample size that achieves desired MOE within the",
                        "user-defined range (%s:%s)"), n_min, n_max))
  return(NULL)
}

#------------------------------------------------
#' @title Get power under a z-test
#'
#' @description
#' Calculates the statistical power to detect whether allele prevalence exceeds a specified threshold
#' using a one-sample z-test of proportions. The function accounts for uncertainty in the true prevalence using a posterior distribution over allele frequencies. Note that this assumes a one-sided test at \code{alpha = 0.05}, and returns the probability that we correctly conclude that prevalence is \emph{above} this threshold. For example, if there is a very low posterior probability that the prevalence is truly above the threshold then power will be low, not because sample size is insufficient, but because there is a low probability that the null hypothesis is false.
#'
#' @param pop Population ID (integer or character) for which to compute power.
#' @param week A single time point (week) at which sampling will occur.
#' @param n_samp A single scalar specifying the number of individuals to be sampled.
#' @param prev_thresh Prevalence threshold to test against (null hypothesis: prevalence = \code{prev_thresh}).
#' @param df_post A data frame of posterior allele frequencies as returned by \code{\link{get_posterior_prev}}.
#' Must include columns \code{pop}, \code{week}, \code{p}, and \code{post}.
#'
#' @return A data frame with a single value: the estimated power (probability of detecting a significant deviation
#' from \code{prev_thresh} under the posterior predictive distribution).
#'
#' @importFrom dplyr mutate summarise
#' @importFrom stats qnorm
#' @export

get_power_ztest <- function(pop, week, n_samp, prev_thresh = 0.05, df_post) {

  # input checks
  assert_single_pos_int(pop, zero_allowed = FALSE)
  assert_single_int(week)
  assert_single_pos_int(n_samp, zero_allowed = FALSE)
  assert_single_bounded(prev_thresh)
  assert_dataframe(df_post)
  assert_in(c("pop", "p", "week", "post"), names(df_post))

  df_post |>
    filter(p > prev_thresh) |>
    project_num_pos(pop = pop, week = week, n_samp = n_samp) |>
    mutate(p_est = n_pos / n_samp,
           z_stat = abs(p_est - prev_thresh) / sqrt(p_est*(1 - p_est) / n_samp),
           is_signif = abs(z_stat) > qnorm(1 - 0.05),
           is_signif = ifelse(n_pos == 0, FALSE, is_signif)) |>
    summarise(power = sum(prob * is_signif))

}

#------------------------------------------------
#' @title Get minimum sample size based on z-test
#'
#' @description
#' Calculates the minimum sample size required to achieve a desired statistical power for detecting
#' whether prevalence exceeds a specified threshold using a one-sample z-test of proportions.
#' Power is computed based on the posterior predictive distribution of allele frequencies.
#'
#' @param pop Population ID (integer or character) for which to calculate the required sample size.
#' @param week A single time point (week) at which sampling will occur.
#' @param prev_thresh Prevalence threshold to test against (null hypothesis: prevalence = \code{prev_thresh}).
#' @param power Desired statistical power (default is \code{0.8}).
#' @param n_min Minimum sample size to begin the search (default is \code{10}).
#' @param n_max Maximum sample size to search (default is \code{1000}).
#' @param df_post A data frame of posterior allele frequencies as returned by \code{\link{get_posterior_prev}}.
#' Must include columns \code{pop}, \code{week}, \code{p}, and \code{post}.
#'
#' @return A named numeric vector with the minimum required sample size, or \code{NULL} if no sample size
#' within the specified range achieves the desired power.
#'
#' @importFrom dplyr mutate summarise
#' @importFrom progress progress_bar
#' @export

get_sample_size_ztest <- function(pop, week, prev_thresh, power = 0.8, n_min = 10, n_max = 1e3, df_post) {

  # input checks
  assert_single_pos_int(pop, zero_allowed = FALSE)
  assert_single_int(week)
  assert_single_bounded(prev_thresh)
  assert_bounded(power)
  assert_single_pos_int(n_min, zero_allowed = FALSE)
  assert_single_pos_int(n_max)
  assert_gr(n_max, n_min)
  assert_dataframe(df_post)
  assert_in(c("pop", "p", "week", "post"), names(df_post))

  pb <- progress_bar$new(
    format = "  computing [:bar] :percent eta: :eta",
    total = n_max, clear = FALSE, width = 60
  )

  # calculate MOE for all sample sizes from n_min to n_max
  for (i in n_min:n_max) {
    power_i <- get_power_ztest(pop = pop, week = week, n_samp = i, prev_thresh = prev_thresh, df_post = df_post)
    if (power_i >= power) {
      return(c(sample_size = i))
    }
    pb$tick()
  }

  # if reached this point then have failed to find sample size
  message(sprintf(paste("Could not find a sample size that achieves desired power within the",
                        "user-defined range (n_min = %s:n_max = %s)"), n_min, n_max))
  return(NULL)
}

