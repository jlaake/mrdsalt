#' Generate a set of double observer detection histories from population data and model
#'
#' The return value is the dataframe of the observations seen by at least
#' one observer and the field detected is added to indicate whether the observer
#' detected the observation.  Those missed by both observers are not included.
#' @usage gen.history(x,par,p.formula,delta.formula,indep=FALSE,PI=FALSE,debug=FALSE,use.offset=FALSE)
#' @param x population data - 1 row per observer and must contain and observer and object fields at a minimum
#' @param par model parameter values
#' @param p.formula formula for true detection probabilities
#' @param delta.formula formula for delta dependence function
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param debug plots p1, conditional p1 and delta values
#' @param use.offset if TRUE use offset value for logit
#' @export
#' @return simulation observation dataframe;The return value is the dataframe of the observations
#' seen by at least one observer and the field detected is added to indicate whether the observer
#'  detected the observation.  Those missed by both observers are not included.
gen.history=function(x,par,p.formula,delta.formula,indep=FALSE,PI=FALSE,debug=FALSE,use.offset=FALSE)
{
  #
  #  Setup p's and deltas
  #
  xmat1=model.matrix(p.formula,x[x$observer==1,])
  xmat2=model.matrix(p.formula,x[x$observer==2,])
  N=dim(xmat1)[1]
  beta=par[1:dim(xmat1)[2]]
  p1=invlogit(xmat1,beta)
  p2=invlogit(xmat2,beta)
  if(indep)
    delta.values=1
  else
  {
    dmat=model.matrix(delta.formula,x[x$observer==1,])
    gamma=par[(1+dim(xmat1)[2]):length(par)]
    delta.values=delta(dmat,p1,p2,gamma,PI,use.offset=use.offset)
  }
  if(debug)
  {
    par(mfrow=c(2,2))
    plot(x$distance[x$observer==1],p1,ylim=c(0,1),xlab="Distance",ylab="Single obs p(x)")
    plot(x$distance[x$observer==2],p1*delta.values,ylim=c(0,1),xlab="Distance",ylab="Single obs conditional p(x)")
    plot(x$distance[x$observer==2],delta.values,xlab="Distance",ylab="Delta",ylim=c(0,max(delta.values)))
  }
  # generate observations for platform1
  x1=rbinom(N,1,p1)
  # compute conditional probabilities for platform2 and observations
  p2c=delta.values*p2*x1+(1-x1)*(p2-delta.values*p1*p2)/(1-p1)
  x2=rbinom(N,1,p2c)
  # return the set of x which were seen by at least one platform
  seen=x1+x2>0
  data1=x[x$observer==1,][seen,]
  data2=x[x$observer==2,][seen,]
  data1$detected=x1[seen]
  data2$detected=x2[seen]
  data=rbind(data1,data2)
  if(debug)
  {
    ncov=100
    xx=hist(data$distance[data$observer==1],breaks=0:10/10,plot=FALSE)
    barplot(xx$count/(N/10),names.arg=1:10/10,xlab="Distance",main="",ylab="Combined detection probability",axis.lty=1,ylim=c(0,1),xlim=c(0,1),width=.1,space=0)
    if(length(beta)<=2)
      xx=cbind(rep(1,101),(0:100)/100)
    else
      xx=cbind(rep(rep(1,101),ncov),rep((0:100)/100,ncov),rep(0:(ncov-1)/(ncov-1),each=101))
    p=invlogit(xx,beta=beta)
    if(length(beta)<=2)
    {
      if(delta.formula==~1)
        xx=matrix(1,nrow=101,ncol=1)
      else
        xx=cbind(rep(1,101),(0:100)/100)
      if(indep)
        delta.values=1
      else
        delta.values=delta(xx,p,p,gamma,PI=FALSE,use.offset=use.offset)
      lines((0:100)/100,2*p-delta.values*p^2,ylim=c(0,1))
    } else
    {
      if(delta.formula==~1)
        xx=matrix(1,nrow=101*ncov,ncol=1)
      else
        xx=cbind(rep(rep(1,101),ncov),rep((0:100)/100,ncov))
      if(indep)
        delta.values=1
      else
        delta.values=delta(xx,p,p,gamma,PI=FALSE,use.offset=use.offset)
      ave.p=by(2*p-delta.values*p^2,factor(rep(0:100,times=ncov)),mean)
      lines((0:100)/100,ave.p,ylim=c(0,1))
    }
  }
  return(data[order(data$object,data$observer),])
}
