#' Computes probability of detection by at least one observer for fitted mr model at a range of distances to
#' integrate pooled detection function
#'
#' @param x vector of distance values
#' @param dd a single pair of detection records for the integration; the pair is replicated length(x) times to compute probability at each x for the integration
#' @param models model list of p and delta formulas for detection probabilities
#' @param par parameters
#' @param indep if TRUE, uses full independence model
#' @export
#' @return vector of values of probability of detection by at least one observer for fitted bpi model at a set of distance
mu.mr=function (x, dd, models,par,indep)
{
  # replicate the pair of observation length(x) times
  dd=dd[rep(1:2,length(x)),]
  # create distance field using input x values
  dd$distance=rep(x,each=2)
  # create object numbers
  dd$object=rep(1:length(x),each=2)
  plist=p.mr(par,dd,models$pformula,models$dformula,indep)
  p1=plist$p1
  p20=plist$p20
  p21=plist$p21
  xx=p1+p20*(1-p1)
  return(xx)
}

