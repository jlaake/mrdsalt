% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu.mr.R
\name{mu.mr}
\alias{mu.mr}
\title{Computes probability of detection by at least one observer for fitted mr model at a range of distances to
integrate pooled detection function}
\usage{
mu.mr(x, dd, models, par, indep)
}
\arguments{
\item{x}{vector of distance values}

\item{dd}{a single pair of detection records for the integration; the pair is replicated length(x) times to compute probability at each x for the integration}

\item{models}{model list of p and delta formulas for detection probabilities}

\item{par}{parameters}

\item{indep}{if TRUE, uses full independence model}
}
\value{
vector of values of probability of detection by at least one observer for fitted bpi model at a set of distance
}
\description{
Computes probability of detection by at least one observer for fitted mr model at a range of distances to
integrate pooled detection function
}
