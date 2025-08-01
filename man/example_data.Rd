% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_data}
\alias{example_data}
\title{Example allele prevalence data}
\format{
A data frame with one row per sampling point and the following columns:
\describe{
  \item{pop}{Population identifier (integer).}
  \item{week}{Sampling week (integer).}
  \item{n_samp}{Number of individuals sampled at that time point.}
  \item{n_pos}{Number of individuals carrying the allele.}
}
}
\usage{
example_data
}
\description{
A small dataset simulating observed allele frequencies across time for multiple populations.
Useful for testing or demonstrating functions such as `estimate_ML()` or `run_mcmc()`.
}
\examples{
head(example_data)

}
\keyword{datasets}
