#' Compute average detection probabilities for Buckland et al model with possible dependence in detections for mark-recapture distance data.
#'
#' @param object fitted bpi model
#' @param ...	additional arguments affecting the predictions produced
#' @return vector of predicted average detection probabilities (integrated over x) for each observation.
#' @seealso \link{integrate_mu.bpi}
#' @author Jeff Laake
#' @export
#
predict.bpi=function(object,...)
{
  W=object$meta.data$width
  # compute integrals
  esw=integrate_mu.bpi(x=object$data,n=nrow(object$data)/2,pformula=object$mrmodel$pformula,
                       dformula=object$mrmodel$dformula,
                       object$par,width=object$meta.data$width,
                       indep=object$control$indep,PI=object$control$PI,
                       use.offset=object$control$use.offset,posdep=object$control$posdep)
  #return integral divided by width
  return(esw/W)
}
