#' Log likelihood calculation for mark-recapture type Huggins model
#'
#' @param par model parameter values
#' @param x observation dataframe
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param width transect half-width
#' @param debug if TRUE, outputs parameter and log likelihood values
#' @param indep if TRUE, uses full independence model
#' @export
#' @return negative log-likelihood value at values of par
ll.mr=function(par,x,pformula,dformula,width,debug,indep)
{
  # create design matrix for p1,p2 and delta
  if(debug)cat("\npar=",par,"\n")
  # compute probabilities and delta values
  p.list=p.mr(par,x,pformula,dformula,indep=indep)
  p1=p.list$p1
  p20=p.list$p20
  p21=p.list$p21
  #  11
  p11=p1*p21
  #  10
  p10=p1*(1-p21)
  #  01
  p01=(1-p1)*p20
  # compute integrals
  esw=integrate_mu.mr(x,n=nrow(x)/2,pformula,dformula,par,width,indep)
  # compute negative log-likelihood
  lnl=sum((x$detected[x$observer==1]*(1-x$detected[x$observer==2]))*log(p10))+
    sum((x$detected[x$observer==2]*(1-x$detected[x$observer==1]))*log(p01))+
    sum((x$detected[x$observer==1]*x$detected[x$observer==2])*log(p11))    -
    sum(log(esw))
  if(debug)cat("\n-lnl=",-lnl)
  return(-lnl)
}
