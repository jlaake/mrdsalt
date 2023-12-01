#' Simulates some mark-recapture distance data using a a specified detection model.
#'
#' Dependence can be induced with a covariate that is not used in the fitted detection model.
#'
#' @param x population data - 1 row per observer and must contain distance and observer and object fields at a minimum
#' @param detfct a function specifying detection proability; uses distance, scale and optionally ...
#' @param W half-width of the transect; observations are out to distance W.
#' @param p0 probability of detection on the line for the 2 observers
#' @param scale.formula parameter values
#' @param par scale parameter values on log scale
#' @param ... additional parameters passed to detfct
#' @return dataframe containing a pair of records (one for each observer) with the fields to be used for analysis including the capture history. Objects not seen ch="00" are excluded from the dataframe.
#' @author Jeff Laake
#' @examples
#' N=1000
#' W=10
#' par=c(1,1,1)
#' p0=c(.9,.9)
#' x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2), distance=rep(runif(N,0,W),each=2),cov=rep(rnorm(N,0,.5),each=2))
#' dd=sim.general(x,halfnormal,W=10,par=par,p0=p0,scale.formula=~-1+observer+cov)
#' par(mfrow=c(2,1))
#' hist(dd$distance[dd$observer==1&dd$ch%in%c("11","10")],main="Observer 1",xlab="Distance")
#' hist(dd$distance[dd$observer==2&dd$ch%in%c("11","01")],main="Observer 2",xlab="Distance")
#' #Fit loglinear model to data
#' results=ddf(data=dd,meta.data=list(width=W),method="loglinear",mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
#' results$criterion
#' # %diff
#' 100*(results$Nhat-N)/N
#' # unconditonal detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#' #Fit PI bpi model to data
#' results=ddf(data=dd,meta.data=list(width=W),method="bpi",control=list(indep=FALSE,PI=TRUE),mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
#' results$criterion
#' # %diff
#' 100*(results$Nhat-N)/N
#' # unconditional detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#' #Fit full independence bpi model to data
#' results=ddf(data=dd,meta.data=list(width=W),control=list(indep=TRUE,PI=FALSE),method="bpi",mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~0))
#' results$criterion
#' # %diff
#' 100*(results$Nhat-N)/N
#' # unconditonal detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#' @export
sim.general=function(x,detfct,W,p0=c(1,1),scale.formula=~-1+observer,par,...)
{
  if(any(!c("distance","observer")%in%colnames(x))) stop("x must contain distance and observer fields at a minimum")
  if(!is.factor(x$observer))
  {
    cat("\nforcing observer to be a factor variable")
    x$observer=factor(x$observer)
  }
  if(!any(levels(x$observer) %in% c("1", "2"))) stop("\n observer must have values 1 and 2 only")
  if(length(p0)!=2)stop("a value of p0 must be given for each of the 2 observers")
  N=nrow(x)/2
  # compute scale values for each value in x; can be observer dependent
  scale.mat=model.matrix(scale.formula,x)
  if(ncol(scale.mat)!=length(par)) stop("length of parameters doesn't match formula")
  scale=scale.mat%*%par
  # compute detection probabilities for each row in x
  px=rep(p0,times=N)*detfct(x$distance,scale,...)
  # generate uniform random variates between 0 and 1 to determine detection
  seen=as.numeric(runif(2*N)<px)
  # add capture history to x
  x$detected=seen
  x$ch=rep(paste(seen[seq(1,nrow(x),2)],seen[seq(2,nrow(x),2)],sep=""),each=2)
  x=x[x$ch!="00",]
  x$object=rep(1:(nrow(x)/2),each=2)
  return(x)
}
#' Half-normal detection function
#'
#' @param distance vector of distances
#' @param scale vector of scale values on log scale for sigma
#' @return vector of detection probabilities
#' @author Jeff Laake
#' @export
halfnormal=function(distance,scale)
  return(exp(-(distance/exp(scale))^2/2))

