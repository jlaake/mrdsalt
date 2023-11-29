#' Integral of probability of detection by at least one observer for fitted bpi model
#'
#' Calls \link{mu.bpi} to compute integrals over detection function values for
#' each observation to integrate average detection probability by at least one observer.
#'
#' @param x observation dataframe
#' @param n number of integrals to compute
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param par parameter values
#' @param width transect half-width
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @return vector of n integral values
#' @author Jeff Laake
#' @export
# compute integral over pooled detection function
integrate_mu.bpi=function(x,n,pformula,dformula,par,width,indep,PI,use.offset,posdep)
{
  # compute integral of detection function over distance x; loop over covariate values
  integrals=vector(mode="numeric",length=n)
  if(all(unique(c(all.vars(pformula),all.vars(dformula)))%in%c("distance","observer")))
  {
    integrals=rep(integrate(mu.bpi,lower=0,upper=width,
              stop.on.error=FALSE,dd=x[1:2,], models=list(pformula=pformula,dformula=dformula),
               par=par,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)$value,n)
  } else
  {
    for (i in 1:n)
      integrals[i]=integrate(mu.bpi,lower=0,upper=width,
                             stop.on.error=FALSE,dd=x[((i-1)*2+1):(i*2),],
                             models=list(pformula=pformula,dformula=dformula),
                             par=par,indep=indep,PI=PI,use.offset=use.offset,posdep=posdep)$value
   }
  return(integrals)
}

