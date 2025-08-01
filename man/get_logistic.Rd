% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{get_logistic}
\alias{get_logistic}
\title{Logistic function}
\usage{
get_logistic(s, p0, week)
}
\arguments{
\item{s}{Selection coefficient. A positive value favours the allele, negative values represent deleterious effects.}

\item{p0}{Initial allele frequency (between 0 and 1).}

\item{week}{A numeric vector of time points (weeks) at which to evaluate the trajectory.}
}
\value{
A numeric vector of allele frequencies at each time point.
}
\description{
Computes the deterministic logistic trajectory of allele frequency over time under selection,
assuming no genetic drift. This models the expected change in allele frequency due to selection alone.
}
