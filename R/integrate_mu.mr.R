#' Integral of probability of detection by at least one observer for fitted mr model
#'
#' Calls \link{mu.mr} to compute integrals over detection function values for
#' each observation to integrate average detection probability by at least one observer.
#'
#' @param x observation dataframe
#' @param n number of integrals to compute
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param par parameter values
#' @param width transect half-width
#' @param indep if TRUE, uses full independence model
#' @return vector of n integral values
#' @author Jeff Laake
#' @export
# compute integral over pooled detection function
integrate_mu.mr=function(x,n,pformula,dformula,par,width,indep,PI,use.offset,posdep)
{
  # compute integral of detection function over distance x; loop over covariate values
  integrals=vector(mode="numeric",length=n)
  if(all(unique(c(all.vars(pformula),all.vars(dformula)))%in%c("distance","observer")))
  {
    integrals=rep(integrate(mu.mr,lower=0,upper=width,
              stop.on.error=FALSE,dd=x[1:2,], models=list(pformula=pformula,dformula=dformula),
               par=par,indep=indep)$value,n)
  } else
  {
    for (i in 1:n)
      integrals[i]=integrate(mu.mr,lower=0,upper=width,
                             stop.on.error=FALSE,dd=x[((i-1)*2+1):(i*2),],
                             models=list(pformula=pformula,dformula=dformula),
                             par=par,indep=indep)$value
   }
  return(integrals)
}

