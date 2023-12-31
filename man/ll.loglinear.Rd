% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ll.loglinear.R
\name{ll.loglinear}
\alias{ll.loglinear}
\title{Compute negative log likelihood log-linear model with possible dependence in detections for mark-recapture distance data.}
\usage{
ll.loglinear(par, cmat, dd, W, pformula, dformula, debug)
}
\arguments{
\item{par}{parameter starting values; if NULL they are computed}

\item{cmat}{capture history matrix with 3 columns (10,01,11) and n rows; contains 1 if that object has the capture history in the column and 0 otherwise}

\item{dd}{dataframe containing double observer data (2 records per observation) with fields observer, detected, distance, ch and any covariate values}

\item{W}{half-width of the transect}

\item{pformula}{formula for detection but excluding the dependence portion}

\item{dformula}{model for dependence due to residual (unmodelled) heterogeneity}

\item{debug}{if TRUE, outputs parameter values and negative log-likelihood as it iterates in optimization}
}
\value{
negative log-likelihood value at those values of par
}
\description{
Compute negative log likelihood log-linear model with possible dependence in detections for mark-recapture distance data.
}
\seealso{
\link{mu.loglinear}, \link{p.loglinear}
}
\author{
Jeff Laake
}
