#' Simulates some mark-recapture distance data using a log-linear model with possible dependence in detections.
#' data.
#' @param N abundance of population which is sampled by dual observers
#' @param W half-width of the transect; observations are out to distance W.
#' @param par parameter values
#' @param cov a function that creates a single covariate value named cov to be used in the data simulation
#' @param pformula formula for detection but excluding the dependence portion
#' @param dformula model for dependence due residual (unmodelled) heterogeneity
#' @return dataframe containing a pair of records (one for each observer) with the fields to be used for analysis including the capture history.
#' Objects not seen ch="00" are excluded from the dataframe.
#' @author Jeff Laake
#' @export
sim_data=function(N,W,par,cov=NULL,pformula=~-1+observer+observer:distance,dformula=~-1+distance)
{
  # generate N uniform random variates between 0 and W
  distance=runif(N,0,W)
  # if no cov function, create dataframe with object number and distance
  # if there is a cov function add the cov field as well to be used in simulation
  if(is.null(cov))
    x=data.frame(object=1:N,distance=distance)
  else
    x=data.frame(object=1:N,distance=distance,cov=cov(N))
  # create pairs of records (one for each observer)
  dd=make.design.data(x)
  # generate the 4 multinomial cell probabilities for each of the N objects
  prob=p.ll(par,dd,pformula=pformula,dformula=dformula)$prob
  prob=cbind(1-rowSums(prob),prob)
  all_ch=c("00","10","01","11")
  # generate N multinomial random variables 1 to 4 and with that value lookup ch value
  x$ch=all_ch[colSums(apply(prob,1,rmultinom,n=1,size=1)*(1:4))]
  # toss out those not seen by either observer
  x=x[x$ch!="00",]
  # renumber object #s and return simulated data
  x$object=1:nrow(x)
  return(x)
}

