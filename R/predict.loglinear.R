#' Compute average detection probabilities for log-linear model with possible dependence in detections for mark-recapture distance data.
#'
#' @param object fitted loglinear model
#' @param ...	additional arguments affecting the predictions produced
#' @return vector of predicted average detection probabilities (integrated over x) for each observation.
#' @seealso \link{integrate_mu.loglinear}
#' @author Jeff Laake
#' @export
#
predict.loglinear=function(object,...)
{
  W=object$meta.data$width
  # compute integrals
  esw=integrate_mu.loglinear(x=object$data,pformula=object$mrmodel$pformula,
                  dformula=object$mrmodel$dformula,par=object$par,
                  width=W)
  #return integral divided by width
  return(esw/W)
}
