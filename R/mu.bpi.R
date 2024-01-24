#' Computes probability of detection by at least one observer for fitted bpi model at a range of distances to
#' integrate pooled detection function
#'
#' @param x vector of distance values
#' @param dd a single pair of detection records for the integration; the pair is replicated length(x) times to compute probability at each x for the integration
#' @param models model list of p and delta formulas for detection probabilities
#' @param par parameters for p
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @param posdep if TRUE enforces positive dependence in delta
#' @export
#' @return vector of values of probability of detection by at least one observer for fitted bpi model at a set of distance
mu.bpi=function (x, dd, models,par,indep,PI,use.offset,posdep)
{
  # replicate the pair of observation length(x) times
  dd=dd[rep(1:2,length(x)),]
  # create distance field using input x values
  dd$distance=rep(x,each=2)
  # create object numbers
  dd$object=rep(1:length(x),each=2)
  plist=p.bpi(par,dd,models$pformula,models$dformula,indep,PI,use.offset,posdep)
  p1=plist$p1
  p2=plist$p2
  delta.values=plist$delta.values
  xx=p1+p2-delta.values*p1*p2
  xx[xx<0]=0
  xx[xx>1]=1
  return(xx)
}

