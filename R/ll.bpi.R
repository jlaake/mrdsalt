#' Log likelihood calculation for Buckland et al 2010 model
#'
#' For a set of observations, computes observer specific probabilities, integrals and delta values See Buckland et al, (2010) for
#' definitions
#'
#' @usage ll.bpi(par,x,p.formula,delta.formula,width,debug=FALSE,indep=FALSE,PI=FALSE,use.offset=FALSE,posdep=FALSE)
#' @param par model parameter values
#' @param x observation dataframe
#' @param p.formula formula for detection probabilities
#' @param delta.formula formula for delta dependence function
#' @param width transect half-width
#' @param debug if TRUE, outputs parameter and log likelihood values
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @param posdep if TRUE enforces positive dependence in delta
#' @export
#' @return negative log-likelihood value at values of par
ll.bpi=function(par,x,p.formula,delta.formula,width,debug=FALSE,indep=FALSE,PI=FALSE,
                 use.offset=FALSE,posdep=FALSE)
{
  # create design matrix for p1,p2 and delta
  if(debug)cat("\npar=",par,"\n")
  # compute probabilities and delta values
  p.list=p.bpi(par,x,p.formula,delta.formula,width,indep=indep,PI=PI,
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
  lnl=sum((x$detected[x$observer==1]*(1-x$detected[x$observer==2]))*log(p10))+
    sum((x$detected[x$observer==2]*(1-x$detected[x$observer==1]))*log(p01))+
    sum((x$detected[x$observer==1]*x$detected[x$observer==2])*log(p11))    -
    sum(log(p.list$mudot))
  if(debug)cat("\n-lnl=",-lnl)
  return(-lnl)
}
