#' Probability calculations for Buckland et al 2010 model
#'
#' For a set of observations, computes observer specific probabilities, integrals and delta values See Buckland et al, (2010) for
#' definitions
#'
#' @usage p.bpi(par,x,p.formula,delta.formula,width,indep=FALSE,PI=FALSE,use.offset=FALSE,posdep=FALSE)
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
  # create design matrix for p1,p2 and delta
  xmat1=model.matrix(p.formula,x[x$observer==1,])
  xmat2=model.matrix(p.formula,x[x$observer==2,])
  dmat=model.matrix(delta.formula,x[x$observer==1,])
  parnames=c(paste("p:",colnames(xmat1)),paste("delta:",colnames(dmat)))
  if(PI & all(dmat[,1]==1)) stop("\nError: No intercept allowed with PI model\n")
  # extract parameter vectors: beta for detection and gamma for delta
  beta=par[1:dim(xmat1)[2]]
  if(posdep)
    gamma=exp(par[(1+dim(xmat1)[2]):length(par)])
  else
    gamma=par[(1+dim(xmat1)[2]):length(par)]
  # compute probabilities and delta values
  p1=invlogit(xmat1,beta)
  p2=invlogit(xmat2,beta)
  if(indep)
    delta.values=1
  else
  {
    delta.values=delta(dmat,p1,p2,gamma,PI,use.offset)
  }
  # compute integral of detection function over distance x; loop over covariate values
  n=dim(xmat1)[1]
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
  return(list(p1=p1,p2=p2,mudot=integrals,delta.values=delta.values,parnames=parnames))
}
