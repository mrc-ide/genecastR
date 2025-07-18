% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_prev}
\alias{plot_prev}
\title{Plot posterior prevalence and observed data}
\usage{
plot_prev(df_post = NULL, df_trajectory = NULL, df_data = NULL, CI_type = "CP")
}
\arguments{
\item{df_post}{Optional data frame of posterior prevalence distributions as returned by
\code{\link{get_posterior_prev}}. Must include columns \code{pop}, \code{t}, \code{p}, and \code{post}.}

\item{df_trajectory}{Optional data frame of simulated prevalence trajectories (e.g., from \code{sim_WF} or \code{sim_diffusion}),
with columns \code{pop}, \code{t}, and \code{p}.}

\item{df_data}{Optional data frame of observed prevalence data with columns \code{pop}, \code{t},
\code{p_est} (estimated prevalence), and confidence interval bounds \code{CI_lower} and \code{CI_upper}.}

\item{CI_type}{If plotting data, which type of confidence interval to produce (see \code{?get_CIs()}).}
}
\value{
A \code{ggplot} object visualizing the specified components:
\itemize{
  \item Posterior prevalence quantile bands (if \code{df_post} is provided)
  \item Simulated allele frequency trajectories (if \code{df_trajectory} is provided)
  \item Observed data with confidence intervals (if \code{df_data} is provided)
}
}
\description{
Visualizes posterior estimates of allele prevalence over time, optionally overlaying
simulated trajectories and observed data with confidence intervals. Posterior distributions
are shown as quantile bands computed from the posterior probability mass at each time point.
}
\details{
The plot includes posterior quantile shading based on the cumulative posterior mass, grouped into bands
(e.g. 0–20\%, 20–50\%, etc.). Multiple populations are shown using \code{facet_wrap()}.
At least one of the three arguments must be supplied.
}
