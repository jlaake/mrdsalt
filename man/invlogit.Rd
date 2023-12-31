% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invlogit.r
\name{invlogit}
\alias{invlogit}
\title{Inverse logit function with and without an offset}
\usage{
invlogit(x,beta,offset=NULL)
}
\arguments{
\item{x}{a design matrix}

\item{beta}{parameters for logit function}

\item{offset}{offset value for logit}
}
\value{
Vector of inverse logit with a value for each row in design matrix x
}
\description{
Computes inverse logit values with and without an offset
}
