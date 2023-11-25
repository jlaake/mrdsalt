#' Integral of probability of detection by at least one observer for fitted loglinear model
#'
#' Calls \link{mu.loglinear} to compute integrals over detection function values for
#' each observation to integrate average detection probability by at least one observer.
#'
#' @param x observation dataframe
#' @param n number of integrals to compute
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param par parameters for model
#' @param width transect half-width
#' @return vector of n integral values
#' @author Jeff Laake
#' @export
# compute integral over pooled detection function
integrate_mu.loglinear=function(x,pformula,dformula,par,width)
{
  # get number of observations
  n=nrow(x)/2
  # if formula only uses distance and observer then compute a single effective strip with (integral of detection function mu)
  if(all(unique(c(all.vars(pformula),all.vars(dformula)))%in%c("distance","observer")))
  {
    integral=integrate(mu.loglinear,lower=0,upper=width,par=par,dd=x[1:2,],pformula=pformula,dformula=dformula)$value
    return(rep(integral,n))
  } else
  {
    # if other covariates then compute the integrals of mu for each observation
    for(i in 1:n)
    {
      integrals=integrate(mu.loglinear,lower=0,upper=width,par=par,dd=x[((i-1)*2+1):(i*2),],pformula=pformula,dformula=dformula)$value
      return(integrals)
    }
  }
}

