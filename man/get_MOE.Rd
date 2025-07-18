% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_MOE}
\alias{get_MOE}
\title{Estimate expected margin of error for a future study}
\usage{
get_MOE(pop, week, n_samp, df_post, CI_type = "CP")
}
\arguments{
\item{pop}{Integer or character ID specifying the population for which to estimate margin of error.}

\item{week}{A single time point (week) at which sampling will occur.}

\item{n_samp}{A single scalar specifying the number of individuals to be sampled.}

\item{df_post}{A data frame containing the posterior distribution over allele frequencies, as returned
by \code{\link{get_posterior_prev}}. Must include columns \code{pop}, \code{week}, \code{p}, and \code{post}.}

\item{CI_type}{Which type of MOE to use (see \code{?get_CIs()}). If Clopper-Pearson then MOE is taken as the average of the upper and lower CIs, which may be different.}
}
\value{
A single numeric value: the estimated margin of error (half-width of the 95\% confidence interval)
that is expected to contain the true prevalence with at least 80\% probability.
}
\description{
Calculates the expected margin of error (MOE) for a future prevalence study based on the posterior
distribution of allele frequencies. Returns the smallest MOE such that there is at least an 80\%
chance the observed confidence interval will fall within that margin. Assumes sampling is performed
at a single time point.
}
