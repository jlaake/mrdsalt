#' Compute average detection probabilities for mark-recapture with possible dependence in detections for mark-recapture distance data.
#'
#' @param object fitted mr model
#' @param ...	additional arguments affecting the predictions produced
#' @return vector of predicted average detection probabilities (integrated over x) for each observation.
#' @seealso \link{integrate_mu.mr}
#' @author Jeff Laake
#' @export
#
predict.mr=function(object,...)
{
  W=object$meta.data$width
  # compute integrals
  esw=integrate_mu.mr(x=object$data,n=nrow(object$data)/2,pformula=object$mrmodel$pformula,
                       dformula=object$mrmodel$dformula,
                       object$par,width=object$meta.data$width,
                       indep=object$control$indep)
  #return integral divided by width
  return(esw/W)
}
