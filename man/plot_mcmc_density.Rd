% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_mcmc_density}
\alias{plot_mcmc_density}
\title{Density plot of posterior MCMC samples}
\usage{
plot_mcmc_density(mcmc)
}
\arguments{
\item{mcmc}{An object returned by \code{run_mcmc()}, containing the MCMC output
including both burn-in and sampling phases.}
}
\value{
A \code{ggplot} object showing a filled contour plot of the posterior density over
the parameter space of \code{s} and \code{sigma}.
}
\description{
Generates a filled 2D density plot of posterior samples for the selection coefficient (\code{s})
and diffusion rate (\code{sigma}) from a completed MCMC run. The contours represent different
quantile levels of the joint posterior distribution, allowing for visualization of uncertainty
and correlation between parameters.
}
