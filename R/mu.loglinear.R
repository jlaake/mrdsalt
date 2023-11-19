#' Integral of probability of detection by at least one observer for fitted log-linear model
#'
#' Calls \link{p.loglinear} to compute 1-1/K for a range of x (distance) values for
#' each observation to integrate average detection probability by at least one observer.
#'
#' @param x vector of distance values from integrate to compute integral
#' @param par parameter starting values; if NULL they are computed
#' @param dd a single pair of detection records for the integration; the pair is replicated length(x) times to compute probability at each x for the integration
#' @param pformula formula for detection but excluding the dependence portion
#' @param dformula model for dependence due residual (unmodelled) heterogeneity
#' @return vector of values at different x (distance) values for integrate
#' @author Jeff Laake
#' @seealso \link{p.loglinear}
#' @export
# compute integral over pooled detection function
mu.loglinear=function(x,par,dd,pformula,dformula)
{
  # replicate the pair of observation length(x) times
  dd=dd[rep(1:2,length(x)),]
  # create distance field using input x values
  dd$distance=rep(x,each=2)
  # create object numbers
  dd$object=rep(1:length(x),each=2)
  # compute K which is the sum of the 3 multinomial logit values
  K=p.loglinear(par,dd,pformula=pformula,dformula=dformula)$K
  # return vector of probabilities that at least one observer detected the object for each distance
  return(1-1/K)
}
