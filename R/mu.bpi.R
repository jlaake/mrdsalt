#' Integral of probability of detection by at least one observer for fitted bpi model
#'
#' Calls \link{integratelogistic} to compute integrals over detection function values for
#' each observation to integrate average detection probability by at least one observer.
#'
#' @param x observation dataframe
#' @param n number of integrals to compute
#' @param p.formula formula for detection probabilities
#' @param delta.formula formula for delta dependence function
#' @param beta parameters for p.formula
#' @param gamma parameters for delta.formula
#' @param width transect half-width
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @return vector of n integral values
#' @author Jeff Laake
#' @export
# compute integral over pooled detection function
mu.bpi=function(x,n,p.formula,delta.formula,beta,gamma,width,indep,PI,use.offset)
{
  # compute integral of detection function over distance x; loop over covariate values
  integrals=vector(mode="numeric",length=n)
  if(all(unique(c(all.vars(p.formula),all.vars(delta.formula)))%in%c("distance","observer")))
  {
    integrals=rep(integratelogistic(x[x$observer==1,][1,],x[x$observer==2,][1,],
                                    list(p.formula=p.formula,delta.formula=delta.formula),
                                    beta,gamma,width,indep,PI,use.offset),n)
  } else
  {
    for (i in 1:n)
      integrals[i]=integratelogistic(x[x$observer==1,][i,],x[x$observer==2,][i,],
                                     list(p.formula=p.formula,delta.formula=delta.formula),
                                     beta,gamma,width,indep,PI,use.offset)
  }
  return(integrals)
}

