#' Distance Detection Function Fitting
#'
#' Generic function for fitting detection functions for distance sampling with
#' single and double observer configurations. Independent observer, trial and
#' dependent observer (removal) configurations are included. This is a generic
#' function which does little other than to validate the calling arguments and
#' methods and then calls the appropriate \code{method} specific function to do
#' the analysis. This is a temp version to replace what is in mrds. Not all functionality
#' is included here.
#'
#' @param dsmodel distance sampling model specification
#' @param mrmodel mark-recapture model specification
#' @param data dataframe containing data to be analyzed
#' @param method analysis method
#' @param meta.data list containing settings controlling data structure
#' @param control list containing settings controlling model fitting
#' @param call not implemented for top level ddf function, this is set by ddf as it is passed to the other ddf generics.
#' @return model object of class=(method, "ddf")
#' @export
#' @import optimx
#' @author Jeff Laake
#' # load data
ddf=function (dsmodel = call(), mrmodel = call(), data, method = "ds",
          meta.data = list(), control = list(), call = NULL)
{
  save.options <- options()
  options(contrasts = c("contr.treatment", "contr.poly"))
  if (!is.null(data$observer) & !is.null(data$object))
    data <- data[order(data$object, data$observer), ]
  method <- match.arg(method, c("ds", "io", "io.fi", "trial",
                                "trial.fi", "rem", "rem.fi","loglinear","bpi"))
  if (method %in% c("ds", "io", "trial", "rem")) {
    if (missing(dsmodel)) {
      stop("For method=", method, ", dsmodel must be specified")
    }
  }
  if (method != "ds") {
    if (missing(mrmodel)) {
      stop("For method=", method, ", mrmodel must be specified")
    }
  }
  ptm=proc.time()
  result <- switch(method,
                   ds = mrds:::ddf.ds(dsmodel = dsmodel, data = data,
                                       meta.data = meta.data, control = control, call = match.call()),
                   io = mrds:::ddf.io(dsmodel = dsmodel, mrmodel = mrmodel, data = data,
                               meta.data = meta.data, control = control, call = match.call()),
                   io.fi = mrds:::ddf.io.fi(mrmodel = mrmodel, data = data, meta.data = meta.data,
                                     control = control, call = match.call(), method = method),
                   trial = mrds:::ddf.trial(dsmodel = dsmodel, mrmodel = mrmodel,
                                     data = data, meta.data = meta.data, control = control,
                                     call = match.call()),
                   trial.fi = mrds:::ddf.trial.fi(mrmodel = mrmodel, data = data, meta.data = meta.data, control = control,
                                     call = match.call(), method = "trial.fi"),
                   rem = mrds:::ddf.rem(dsmodel = dsmodel, mrmodel = mrmodel, data = data, meta.data = meta.data, control = control,
                                     call = match.call()),
                   rem.fi = mrds:::ddf.rem.fi(mrmodel = mrmodel, data = data, meta.data = meta.data, control = control,
                                      call = match.call(), method = method),
                   loglinear=ddf.loglinear(mrmodel = mrmodel, data = data, meta.data = meta.data, control = control,
                                           call = match.call()),
                   bpi=ddf.bpi(mrmodel = mrmodel, data = data, meta.data = meta.data, control = control,
                                           call = match.call()))
  options(save.options)
  result$runtime=(proc.time()-ptm)["elapsed"]
  return(result)
}
