#' Compute average detection probabilities for log-linear model with possible dependence in detections for mark-recapture distance data.
#' @param object fitted loglinear model
#' @param ...	additional arguments affecting the predictions produced
#' @return vector of predicted average detection probabilities (integrated over x) for each observation.
#' @seealso \link{mu.loglinear}, \link{p.loglinear}
#' @author Jeff Laake
#' @export
#
predict.loglinear=function(object,...)
{
  data=object$data
  # get number of observations
  n=nrow(data)/2
  pformula=object$mrmodel$pformula
  dformula=object$mrmodel$dformula
  par=object$par
  W=object$meta.data$width
  # compute multinomial probability values at current set of par values
  prob=p.loglinear(par=par,x=data,pformula=pformula,dformula=dformula)$prob
  # if formula only uses distance and observer then compute a single effective strip with (integral of detection function mu)
  if(all(unique(c(all.vars(pformula),all.vars(dformula)))%in%c("distance","observer")))
  {
    avep=integrate(mu.loglinear,lower=0,upper=W,par=par,dd=data[1:2,],pformula=pformula,dformula=dformula)$value/W
    return(rep(avep,n))
  } else
  {
    # if other covariates then compute the integrals of mu
    for(i in 1:n)
    {
      avep=integrate(mu.loglinear,lower=0,upper=W,par=par,dd=data[((i-1)*2+1):(i*2),],pformula=pformula,dformula=dformula)$value/W
      return(avep)
    }
  }
}
