#' Summary of bpi model fit
#'
#' Creates a summary object of a ddf loglinear model fitted to double observer data in
#' an independent observer configuration
#'
#' @usage \method{summary}{bpi}(object,se=TRUE,...)
#' @param object a ddf model object
#' @param se if TRUE, computes std errors
#' @param ... additional non-specified argument for S3 generic function
#' @export
#' @return A list with each of the summarized objects that depends on the model
summary.bpi=function (object, se = TRUE, ...)
{
  model <- object
  n <- nrow(model$xmat)/2
  ans <- list(par=model$par, vc=solve(model$hessian), Nhat = model$Nhat, AIC = model$criterion,
              average.p = n/model$Nhat, convergence=model$mod$convcode)

  # if (se) {
  #   se.obj <- calc.se.Np(model, avgp, n, ans$average.p)
  #   ans$average.p.se <- se.obj$average.p.se
  #   ans$Nhat.se <- se.obj$Nhat.se
  #   ans$cv <- ans$Nhat.se/model$Nhat
  # }
  class(ans) <- "summary.bpi"
  return(ans)
}
