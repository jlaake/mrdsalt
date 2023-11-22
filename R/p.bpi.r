#' Probability calculations for Buckland et al 2010 model
#'
#' For a set of observations, computes observer specific probabilities, integrals and delta values See Buckland et al, (2010) for
#' definitions
#'
#' @param par model parameter values
#' @param x observation dataframe
#' @param p.formula formula for detection probabilities
#' @param delta.formula formula for delta dependence function
#' @param width transect half-width
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @param posdep if TRUE enforces positive dependence in delta
#' @export
#' @return list of vectors for p1,p2, integrals and delta values
p.bpi=function(par,x,p.formula,delta.formula,width,indep=FALSE,PI=FALSE,
               use.offset=FALSE,posdep=FALSE)
{
  # create design matrix for p and delta
  xmat=model.matrix(p.formula,x)
  #xmat2=model.matrix(p.formula,x[x$observer==2,])

  dmat=model.matrix(delta.formula,x[x$observer==1,])
  parnames=c(paste("p:",colnames(xmat)),paste("delta:",colnames(dmat)))
  if(PI & all(dmat[,1]==1)) stop("\nError: No intercept allowed with PI model\n")
  # extract parameter vectors: beta for detection and gamma for delta
  beta=par[1:dim(xmat)[2]]
  if(posdep)
    gamma=exp(par[(1+dim(xmat)[2]):length(par)])
  else
    gamma=par[(1+dim(xmat)[2]):length(par)]
  # compute probabilities and delta values
  px=invlogit(xmat,beta)
  p1=px[seq(1,length(px),2)]
  p2=px[seq(2,length(px),2)]
  if(indep)
    delta.values=1
  else
  {
    delta.values=delta(dmat,p1,p2,gamma,PI,use.offset)
  }
  # compute integral of detection function over distance x; loop over covariate values
  n=length(p1)
  integrals=mu.bpi(x,n,p.formula,delta.formula,beta,gamma,width,indep,PI,use.offset)
  return(list(p1=p1,p2=p2,mudot=integrals,delta.values=delta.values,parnames=parnames))
}
