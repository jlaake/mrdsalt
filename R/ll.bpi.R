#' Log likelihood calculation for Buckland et al 2010 model
#'
#' For a set of observations, computes observer specific probabilities, integrals and delta values See Buckland et al, (2010) for
#' definitions
#'
#' @param par model parameter values
#' @param x observation dataframe
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param width transect half-width
#' @param debug if TRUE, outputs parameter and log likelihood values
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @param posdep if TRUE enforces positive dependence in delta
#' @export
#' @return negative log-likelihood value at values of par
ll.bpi=function(par,x,pformula,dformula,width,debug,indep,PI,
                 use.offset,posdep)
{
  # create design matrix for p1,p2 and delta
  if(debug)cat("\npar=",par,"\n")
  # compute probabilities and delta values
  p.list=p.bpi(par,x,pformula,dformula,indep=indep,PI=PI,
               use.offset=use.offset,posdep=posdep)
  p1=p.list$p1
  p2=p.list$p2
  delta.values=p.list$delta.values
  if(debug)print(summary(data.frame(p1=p1)))
  if(debug)print(summary(data.frame(p2=p2)))
  if(debug&!indep)print(summary(data.frame(delta=delta.values)))
  #  11
  p11=delta.values*p1*p2
  #  10
  p10=p1-p11
  #  01
  p01=p2-p11
  # compute integrals
  esw=integrate_mu.bpi(x,n=nrow(x)/2,pformula,dformula,par,width,indep,PI,use.offset,posdep)
  # compute negative log-likelihood
  lnl=sum((x$detected[x$observer==1]*(1-x$detected[x$observer==2]))*log(p10))+
    sum((x$detected[x$observer==2]*(1-x$detected[x$observer==1]))*log(p01))+
    sum((x$detected[x$observer==1]*x$detected[x$observer==2])*log(p11))    -
    sum(log(esw))
  if(debug)cat("\n-lnl=",-lnl)
  return(-lnl)
}
