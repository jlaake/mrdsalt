% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.general.R
\name{sim.general}
\alias{sim.general}
\title{Simulates some mark-recapture distance data using aspecified detection model.}
\usage{
sim.general(
  x,
  detfct,
  W,
  p0 = c(1, 1),
  scale.formula = ~-1 + observer,
  par,
  ...
)
}
\arguments{
\item{x}{population data with 1 row per observer and must contain distance and observer and object fields at a minimum}

\item{detfct}{a function specifying detection proability; uses distance, scale and optionally covariate cov}

\item{W}{half-width of the transect; observations are out to distance W.}

\item{p0}{probability of detection on the line for the 2 observers; no covariates are used}

\item{scale.formula}{parameter values}

\item{par}{scale parameter values on log scale}

\item{...}{additional parameters passed to detfct}
}
\value{
dataframe containing a pair of records (one for each observer) with the fields to be used for analysis including the capture history. Objects not seen ch="00" are excluded from the dataframe.
}
\description{
Dependence can be induced with a covariate that is not used in the fitted detection model.
}
\examples{
N=1000
W=10
par=c(1,1,1)
p0=c(.9,.9)
x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2),
  distance=rep(runif(N,0,W),each=2),cov=rep(rnorm(N,0,.5),each=2))
dd=sim.general(x,halfnormal,W=10,par=par,p0=p0,
               scale.formula=~-1+observer+cov)
par(mfrow=c(2,1))
hist(dd$distance[dd$observer==1&dd$ch\%in\%c("11","10")],
  main="Observer 1",xlab="Distance")
hist(dd$distance[dd$observer==2&dd$ch\%in\%c("11","01")],
  main="Observer 2",xlab="Distance")
#Fit loglinear model to data
results=ddf(data=dd,meta.data=list(width=W),method="loglinear",
mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
results$criterion
# \%diff
100*(results$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)
#Fit PI bpi model to data
results=ddf(data=dd,meta.data=list(width=W),method="bpi",
 control=list(indep=FALSE,PI=TRUE),
 mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
results$criterion
# \%diff
100*(results$Nhat-N)/N
# unconditional detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)
#Fit full independence bpi model to data
results=ddf(data=dd,meta.data=list(width=W),
  control=list(indep=TRUE,PI=FALSE),method="bpi",
  mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~0))
results$criterion
# \%diff
100*(results$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)
}
\author{
Jeff Laake
}
