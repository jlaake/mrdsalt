#' Compute negative log likelihood log-linear model with possible dependence in detections for mark-recapture distance data.
#' @param par parameter starting values; if NULL they are computed
#' @param cmat capture history matrix with 3 columns (10,01,11) and n rows; contains 1 if that object has the capture history in the column and 0 otherwise
#' @param dd dataframe containing double observer data (2 records per observation) with fields observer, detected, distance, ch and any covariate values
#' @param W half-width of the transect
#' @param pformula formula for detection but excluding the dependence portion
#' @param dformula model for dependence due to residual (unmodelled) heterogeneity
#' @param debug if TRUE, outputs parameter values and negative log-likelihood as it iterates in optimization
#' @return negative log-likelihood value at those values of par
#' @seealso \link{mu.loglinear}, \link{p.loglinear}
#' @author Jeff Laake
#' @export
#'
ll.loglinear=function(par,cmat,dd,W,pformula,dformula,debug)
{
  # get number of observations
  n=nrow(dd)/2
  # compute multinomial probability values at current set of par values
  prob=p.loglinear(par=par,x=dd,pformula=pformula,dformula=dformula)$prob
  # compute integrals
  esw=integrate_mu.loglinear(dd,pformula,dformula,par,width=W)
  #compute negative log-likelihood
  nll=-sum(log(rowSums(prob*cmat)))+sum(log(esw))
  # if debug, print out results for this iteration and then return nll value to optimx
  if(debug) cat("Par = ",par,"\n","-ll = ",nll)
  return(nll)
}
