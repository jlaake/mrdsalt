#' Generate a set of double observer detection histories from population data and bpi model
#'
#' The return value is the dataframe of the observations seen by at least
#' one observer and the field detected is added to indicate whether the observer
#' detected the observation.  Those missed by both observers are not included.
#' @param x population data - 1 row per observer and must contain distance and observer and object fields at a minimum
#' @param par model parameter values
#' @param pformula formula for true detection probabilities
#' @param dformula formula for delta dependence function
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param debug plots p1, conditional p1 and delta values
#' @param use.offset if TRUE use offset value for logit
#' @export
#' @return simulation observation dataframe;The return value is the dataframe of the observations seen by at least one observer and the field detected is added to indicate whether the observer detected the observation.  Those missed by both observers are not included.
#' @examples
#' N=1000
#' W=10
#' x=sim.bpi(x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2),
#'  distance=rep(runif(N/2,0,W),each=2)),
#'  par=c( 1.517424, 1.663175, -0.28, -0.5, 0.1311178),
#'  pformula=~observer*distance,dformula=~-1+distance,PI=TRUE)
#'
sim.bpi=function(x,par,pformula,dformula,indep=FALSE,PI=TRUE,
                     debug=FALSE,use.offset=FALSE)
{
  #
  #  Setup p's and deltas
  #
  xmat=model.matrix(pformula,x)
  N=nrow(x)/2
  beta=par[1:ncol(xmat)]
  px=plogis(xmat%*%beta)
  p1=px[seq(1,nrow(x),2)]
  p2=px[seq(2,nrow(x),2)]
  if(indep)
    delta.values=1
  else
  {
    dmat=model.matrix(dformula,x[x$observer==1,])
    gamma=par[(1+length(beta)):length(par)]
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
  data$ch=paste(data1$detected,data2$detected,sep="")
  if(debug)
  {
    ncov=100
    xx=hist(data$distance[data$observer==1],breaks=0:10/10,plot=FALSE)
    barplot(xx$count/(N/10),names.arg=1:10/10,xlab="Distance",main="",ylab="Combined detection probability",axis.lty=1,ylim=c(0,1),xlim=c(0,1),width=.1,space=0)
    if(length(beta)<=2)
      xx=cbind(rep(1,101),(0:100)/100)
    else
      xx=cbind(rep(rep(1,101),ncov),rep((0:100)/100,ncov),rep(0:(ncov-1)/(ncov-1),each=101))
    p=plogis(xx%*%beta)
    if(length(beta)<=2)
    {
      if(dformula==~1)
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
      if(dformula==~1)
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
