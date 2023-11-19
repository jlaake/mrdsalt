#' Plot some detection function values. This is very simple at present.
#' @param mod fitted result from fit_mrdsalt
#' @return None
#' @author Jeff Laake
#' @export plot_det plot_uncond plot_cond plot.loglinear
#' @import mrds
#' @importFrom graphics barplot hist lines par  points title
#' @importFrom stats as.formula glm integrate model.matrix predict rbinom rmultinom rnorm runif
plot_det=function(mod)
{
  if(is.null(mod$mrmodel$pformula))stop("full model result with formulas not used in call")
  # show plots of detection function for observer 1,2 duplicates, pooled and conditional
  prob=p.loglinear(par=mod$par,mod$data,pformula=mod$mrmodel$pformula,dformula=mod$mrmodel$dformula)$prob
  delta=prob[,3]/((prob[,1]+prob[,3])*(prob[,2]+prob[,3]))

  par(mfrow=c(2,2))
  plot(mod$data$distance[seq(1,nrow(mod$data),2)],prob[,1]+prob[,3],ylim=c(0,1),xlab="Distance",ylab="Detection probability",main="Observer 1")
  points(mod$data$distance[seq(1,nrow(mod$data),2)],delta*(prob[,1]+prob[,3]),ylim=c(0,1),xlab="Distance",ylab="Detection probability",main="Observer 1")
  plot(mod$data$distance[seq(1,nrow(mod$data),2)],prob[,2]+prob[,3],ylim=c(0,1),xlab="Distance",ylab="Detection probability",main="Observer 2",pch="c")
  points(mod$data$distance[seq(1,nrow(mod$data),2)],delta*(prob[,2]+prob[,3]),ylim=c(0,1),xlab="Distance",ylab="Detection probability",main="Observer 2",pch="c")
  plot(mod$data$distance[seq(1,nrow(mod$data),2)],prob[,2]+prob[,3],ylim=c(0,1),xlab="Distance",ylab="Detection probability",main="Duplicate")

  pooled=prob[,1]+prob[,2]+prob[,3]
  plot(mod$data$distance[seq(1,nrow(mod$data),2)],pooled,xlab="Distance",ylim=c(0,1),ylab="Detection probability",pch="p",main="Pooled and Conditional")
  points(mod$data$distance[seq(1,nrow(mod$data),2)],pooled*delta,pch="c")
  return()
}

#' Plot conditional detection function from distance sampling model
#'
#' Plot proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data. Internal function called by \code{plot} methods.
#'
#' @aliases plot_cond
#' @param obs observer code
#' @param xmat processed data
#' @param gxvalues detection function values for each observation
#' @param model fitted model from \code{ddf}
#' @param nc number of equal-width bins for histogram
#' @param breaks user define breakpoints
#' @param finebr fine break values over which line is averaged
#' @param showpoints logical variable; if \code{TRUE} plots predicted value
#'   for each observation
#' @param showlines logical variable; if \code{TRUE} plots average predicted
#'   value line
#' @param maintitle main title line for each plot
#' @param ylim range of y axis (default \code{c(0,1)})
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
plot_cond=function (obs, xmat, gxvalues, model, nc, breaks, finebr, showpoints,
          showlines, maintitle, ylim, angle = -45, density = 20, col = "black",
          jitter = NULL, xlab = "Distance", ylab = "Detection probability",
          subtitle = TRUE, ...)
{
  selection <- xmat$detected[xmat$observer != obs] == 1
  selmat <- (xmat[xmat$observer == obs, ])[selection, ]
  inside_breaks <- selmat$distance >= min(breaks) & selmat$distance <=
    max(breaks)
  selmat <- selmat[inside_breaks, ]
  gxvalues <- gxvalues[inside_breaks]
  shist <- hist(xmat$distance[xmat$observer != obs & xmat$detected ==
                                1], breaks = breaks, plot = FALSE)
  mhist <- hist(xmat$distance[xmat$timesdetected == 2 & xmat$observer ==
                                obs], breaks = breaks, plot = FALSE)
  if (length(mhist$counts) < length(shist$counts)) {
    prop <- c(mhist$counts/shist$counts[1:length(mhist$counts)],
              rep(0, (length(shist$counts) - length(mhist$counts))))
  }
  else {
    prop <- mhist$counts/shist$counts
  }
  mhist$density <- prop
  mhist$equidist <- FALSE
  mhist$intensities <- mhist$density
  mrds:::histline(mhist$density, breaks = breaks, lineonly = FALSE,
           xlab = xlab, ylab = ylab, ylim = ylim, angle = angle,
           density = density, col = col, det.plot = TRUE, ...)
  if (showlines) {
    if(model$method=="loglinear")
      line <- average.line.cond.ll(finebr, obs, model)
    else
      line <- mrds:::average.line.cond(finebr, obs, model)
    linevalues <- line$values
    xgrid <- line$xgrid
    lines(xgrid, linevalues, ...)
  }
  if (showpoints) {
    ifelse(is.null(jitter), jitter.p <- 1, jitter.p <- rnorm(length(gxvalues),
                                                             1, jitter))
    points(selmat$distance, gxvalues * jitter.p, ...)
  }
  if (maintitle != "") {
    if (subtitle) {
      title(paste(maintitle, "\nConditional detection probability\nObserver=",
                  obs, " | Observer = ", 3 - obs), ...)
    }
    else {
      title(maintitle)
    }
  }
  else {
    if (subtitle) {
      title(paste("Conditional detection probability\nObserver=",
                  obs, " | Observer = ", 3 - obs), ...)
    }
  }
}

#' Plot unconditional detection function from distance sampling model
#'
#' Plots unconditional detection function for observer=obs observations
#' overlays histogram, average detection function and values for individual
#' observations data. Internal function called by \code{plot} methods.
#'
#' @aliases plot_uncond
#' @param model fitted model from \code{ddf}
#' @param obs value of observer for plot
#' @param xmat processed data
#' @param gxvalues detection function values for each observation
#' @param nc number of equal-width bins for histogram
#' @param finebr fine break values over which line is averaged
#' @param breaks user define breakpoints
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE plots average predicted value line
#' @param maintitle main title line for each plot
#' @param ylim range of y axis; defaults to (0,1)
#' @param return.lines if TRUE, returns values for line
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return if \code{return.lines==TRUE} returns dataframe \code{average.line}
#'  otherwise just plots
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
plot_uncond=function (model, obs, xmat, gxvalues, nc, finebr, breaks, showpoints,
          showlines, maintitle, ylim, return.lines = FALSE, angle = -45,
          density = 20, col = "black", jitter = NULL, xlab = "Distance",
          ylab = "Detection probability", subtitle = TRUE, ...)
{
  if (obs <= 2) {
    n <- length(xmat$distance[xmat$observer == obs & xmat$detected ==
                                1])
    selmat <- xmat[xmat$observer == obs, ]
    det.detected <- xmat$detected[xmat$observer == obs] ==
      1
  }
  else if (obs == 3) {
    n <- length(xmat$distance[xmat$observer == 1])
    selmat <- xmat[xmat$observer == 1, ]
    det.detected <- selmat$observer == 1
  }
  else if(obs==4){
    n <- length(xmat$distance[xmat$observer == 1 & xmat$timesdetected ==
                                2])
    selmat <- xmat[xmat$observer == 1, ]
    det.detected <- selmat$timesdetected == 2
  } else if(obs==5){
    n <- length(xmat$distance[xmat$observer == 1 & xmat$detected ==
                                1 & xmat$timesdetected==1])
    selmat <- xmat[xmat$observer == 1,]
    det.detected <- selmat$detected==1 & selmat$timesdetected ==1
  } else{
    n <- length(xmat$distance[xmat$observer == 2 & xmat$detected ==
                                              1 & xmat$timesdetected==1])
    selmat <- xmat[xmat$observer == 2, ]
    det.detected <- selmat$detected==1 & selmat$timesdetected==1
  }
  inside_breaks <- selmat$distance >= min(breaks) & selmat$distance <=
    max(breaks)
  selmat <- selmat[inside_breaks, ]
  gxvalues <- gxvalues[inside_breaks]
  hist.obj <- hist(selmat$distance[det.detected], breaks = breaks,
                   plot = FALSE)
  if (!model$meta.data$point) {
    expected.counts <- (breaks[2:(nc + 1)] - breaks[1:nc]) *
      (model$Nhat/breaks[nc + 1])
  }
  else {
    expected.counts <- -apply(matrix(c(breaks[2:(nc + 1)]^2,
                                       breaks[1:nc]^2), ncol = 2, nrow = nc), 1, diff) *
      (model$Nhat/breaks[nc + 1]^2)
  }
  hist.obj$density <- hist.obj$counts/(expected.counts)
  hist.obj$intensities <- hist.obj$density
  freq <- hist.obj$density
  hist.obj$equidist <- FALSE
  ylim <- c(0, max(ylim, hist.obj$density))
  mrds:::histline(hist.obj$density, breaks = breaks, lineonly = FALSE,
           xlab = xlab, ylab = ylab, ylim = ylim, angle = angle,
           density = density, col = col, det.plot = TRUE, ...)
  if (showlines)
  {
    if(model$method=="loglinear")
      line <- average.line.ll(finebr, obs, model)
    else
      line <- mrds:::average.line(finebr, obs, model)
    linevalues <- line$values
    xgrid <- line$xgrid
    lines(xgrid, linevalues, ...)
  }
  if (showpoints) {
    ifelse(is.null(jitter), jitter.p <- 1, jitter.p <- rnorm(length(gxvalues),
                                                             1, jitter))
    points(selmat$distance, gxvalues * jitter.p, ...)
  }
  if (!subtitle) {
    if (maintitle != "")
      maintitle <- title(maintitle)
  }
  else {
    if (maintitle != "")
      maintitle <- paste(maintitle, "\n", sep = "")
    if (obs <= 2) {
      title(paste(maintitle, "Observer = ", obs, " detections"),
            ...)
    }
    else {
      if (obs == 3) {
        title(paste(maintitle, "Pooled detections"),
              ...)
      }
      else {
        if(obs==4){
           title(paste(maintitle, "Duplicate detections"),
              ...)
        } else
          if (obs==5){
            title(paste(maintitle, "First observer only"),
                 ...)
           } else{
               title(paste(maintitle, "Second observer only"),
               ...)
           }
        }
     }
  }
  if (return.lines)
    invisible(line)
  else invisible(NULL)
}

plot.io=function (x, which = 1:6, breaks = NULL, nc = NULL, maintitle = "",
                  showlines = TRUE, showpoints = TRUE, ylim = c(0, 1), angle = NULL,
                  density = NULL, col = "lightgrey", jitter = NULL, divisions = 25,
                  pages = 0, xlab = "Distance", ylab = "Detection probability",
                  subtitle = TRUE, ...)
{
  model <- x
  which <- sort(which)
  xmat.p0 <- model$mr$mr$data
  xmat.p0$offsetvalue <- 0
  xmat.p0$distance <- 0
  ddfobj <- model$ds$ds$aux$ddfobj
  if (ddfobj$type == "gamma") {
    xmat.p0$distance <- rep(mrds:::apex.gamma(ddfobj), 2)
  }
  p0 <- predict(model$mr, newdata = xmat.p0, integrate = FALSE)$fitted
  xmat <- model$mr$mr$data
  cond.det <-predict(model$mr, newdata = xmat, integrate = FALSE)
  width <- model$meta.data$width
  left <- model$meta.data$left
  detfct.pooled.values <- mrds:::detfct(xmat$distance[xmat$observer ==
                                                 1], ddfobj, width = width - left)
  delta <- cond.det$fitted/(p0 * detfct.pooled.values)
  p1 <- cond.det$p1
  p2 <- cond.det$p2
  gxlist <- list(p1/delta, p2/delta, (p1 + p2 - p1 * p2)/delta,
                 p1 * p2/delta, p1*(1-p2)/delta,p2*(1-p1)/delta)
  if (is.null(nc)) {
    nc <- round(sqrt(min(length(xmat$distance[xmat$observer ==
                                                1 & xmat$detected == 1]), length(xmat$distance[xmat$observer ==
                                                                                                 2 & xmat$detected == 1]), length(xmat$distance[xmat$observer ==
                                                                                                                                                  1 & xmat$timesdetected == 2]))), 0)
  }
  if (model$meta.data$binned) {
    breaks <- model$meta.data$breaks
    nc <- length(breaks) - 1
  }
  else {
    if (is.null(breaks)) {
      breaks <- left + ((width - left)/nc) * (0:nc)
    }
    else {
      nc <- length(breaks) - 1
    }
  }
  #oask <- plot_layout(which, pages)
  #on.exit(devAskNewPage(oask))
  for (wh in which[which < 7]) {
    plot_uncond(model, wh, xmat, gxvalues = gxlist[[wh]],
                nc, finebr = (width/divisions) * (0:divisions), breaks,
                showpoints, showlines, maintitle, ylim, angle = angle,
                density = density, col = col, jitter = jitter, xlab = xlab,
                ylab = ylab, subtitle = subtitle, ...)
  }
  data <- model$mr$mr$data
  data$offsetvalue <- 0
  if (is.element(7, which)) {
    gxvalues <- p1[xmat$detected[xmat$observer == 2] == 1]
    plot_cond(1, data, gxvalues, model, nc, breaks, finebr = (width/divisions) *
                (0:divisions), showpoints, showlines, maintitle,
              ylim, angle = angle, density = density, col = col,
              jitter = jitter, xlab = xlab, ylab = ylab, subtitle = subtitle,
              ...)
  }
  if (is.element(8, which)) {
    gxvalues <- p2[xmat$detected[xmat$observer == 1] == 1]
    plot_cond(2, data, gxvalues, model, nc, breaks, finebr = (width/divisions) *
                (0:divisions), showpoints, showlines, maintitle,
              ylim, angle = angle, density = density, col = col,
              jitter = jitter, xlab = xlab, ylab = ylab, subtitle = subtitle,
              ...)
  }
  invisible(NULL)
}

average.line<-function (finebr, obs, model)
{
  xgrid <- NULL
  linevalues <- NULL
  if (model$method == "io") {
    newdat <- model$mr$mr$data
  }
  else {
    if (model$method == "trial" | model$method == "trial.fi") {
      newdat <- mrds::process.data(model$data, model$meta.data)$xmat
      newdat <- newdat[newdat$observer == obs & newdat$detected ==
                         1, ]
    }
    else {
      if (model$method == "rem.fi") {
        newdat <- model$data
      }
      else {
        newdat <- model$mr$data
      }
    }
  }
  newdat$offsetvalue <- rep(0, dim(newdat)[1])
  if (model$method == "io" | model$method == "trial" | model$method ==
      "rem") {
    prob.det <- predict(model)$fitted
    newdat$distance <- 0
    ddfobj <- model$ds$ds$aux$ddfobj
    if (ddfobj$type == "gamma") {
      if (model$method == "io") {
        newdat$distance <- rep(mrds:::apex.gamma(ddfobj), 2)
      }
      else {
        newdat$distance <- as.vector(mrds:::apex.gamma(ddfobj))
      }
    }
    if (model$method == "trial") {
      g0 <- predict(model$mr$mr, newdat, type = "response")
    }
    else {
      g0 <- predict(model$mr, newdat, integrate = FALSE)$fitted
    }
  }
  else {
    prob.det <- predict(model, newdat, integrate = TRUE)$fitted
    newdat$distance <- 0
    g0 <- predict(model, newdat, integrate = FALSE)$fitted
  }
  for (i in 1:(length(finebr) - 1)) {
    x <- (finebr[i] + finebr[i + 1])/2
    xgrid <- c(xgrid, x)
    newdat$distance <- rep(x, dim(newdat)[1])
    if (model$method != "io" & model$method != "rem") {
      cond.det <- predict(model, newdata = newdat, integrate = FALSE)
    }
    else {
      cond.det <- predict(model$mr, newdata = newdat, integrate = FALSE)
    }
    if (model$method == "io" | model$method == "io.fi" |
        model$method == "rem" | model$method == "rem.fi") {
      p1 <- cond.det$p1
      p2 <- cond.det$p2
    }
    else {
      p1 <- cond.det$fitted
    }
    par <- model$ds$par
    if (model$method == "io" | model$method == "trial" |
        model$method == "rem") {
      detfct.pooled.values <- mrds:::detfct(newdat$distance[newdat$observer ==
                                                       1], ddfobj, width = model$meta.data$width - model$meta.data$left)
      deltax <- detfct.pooled.values/(cond.det$fitted/g0)
    }
    else {
      detfct.pooled.values <- cond.det$fitted/g0
      deltax <- rep(1, length(detfct.pooled.values))
    }
    if (obs == 1) {
      linevalues <- c(linevalues, sum(p1 * deltax/prob.det)/sum(1/prob.det))
    }
    else if (obs == 2) {
      linevalues <- c(linevalues, sum(p2 * deltax/prob.det)/sum(1/prob.det))
    }
    else if (obs == 3) {
      linevalues <- c(linevalues, sum(g0 * detfct.pooled.values/prob.det)/sum(1/prob.det))
    }
    else if(obs==4) {
      linevalues <- c(linevalues, sum(p1 * p2 * deltax/prob.det)/sum(1/prob.det))
    }
    else if(obs==5) {
      linevalues <- c(linevalues, sum(p1 * (1-p2) * deltax/prob.det)/sum(1/prob.det))
    } else
      linevalues <- c(linevalues, sum((1-p1)* p2 * deltax/prob.det)/sum(1/prob.det))
  }
  return(list(xgrid = xgrid, values = linevalues))
}

average.line.ll<-function (finebr, obs, model)
{
  xgrid <- NULL
  linevalues <- NULL
  newdat <- model$data
  prob.det=predict.loglinear(model)
  for (i in 1:(length(finebr) - 1)) {
    x <- (finebr[i] + finebr[i + 1])/2
    xgrid <- c(xgrid, x)
    newdat$distance <- rep(x, dim(newdat)[1])
    prob=p.loglinear(par=model$par,newdat,pformula=model$mrmodel$pformula,dformula=model$mrmodel$dformula)$prob
    p1 = prob[,3]/(prob[,2]+prob[,3]) # conditional probability 1|2
    p2 = prob[,3]/(prob[,1]+prob[,3]) # conditional probability 2|1
    delta=prob[,3]/((prob[,1]+prob[,3])*(prob[,2]+prob[,3]))
    width <- model$meta.data$width
    left <- model$meta.data$left
    if (obs == 1) {
      linevalues <- c(linevalues, sum(p1/delta/prob.det)/sum(1/prob.det))
    }
    else if (obs == 2) {
      linevalues <- c(linevalues, sum(p2/delta/prob.det)/sum(1/prob.det))
    }
    else if (obs == 3) {
      linevalues <- c(linevalues, sum((p1/delta + p2/delta - p1 * p2/delta^2)/prob.det)/sum(1/prob.det))
    }
    else if(obs==4) {
      linevalues <- c(linevalues, sum((p1 * p2)/delta^2/prob.det)/sum(1/prob.det))
    }
    else if(obs==5) {
      linevalues <- c(linevalues, sum(p1/delta*(1-p2)/prob.det)/sum(1/prob.det))
    } else
      linevalues <- c(linevalues, sum(p2/delta*(1-p1)/prob.det)/sum(1/prob.det))
  }
  return(list(xgrid = xgrid, values = linevalues))
}

#' Average conditional detection function line for plotting
#'
#' For models with covariates the detection probability for each observation
#' can vary.  This function computes an average value for a set of distances to
#' plot an average line to graphically represent the fitted model in plots that
#' compare histograms and the scatter of individual estimated detection
#' probabilities.
#'
#' @param finebr set of fine breaks in distance over which detection function
#'   values are averaged and plotted
#' @param obs value of observer for averaging (1-2 individual observers)
#' @param model ddf model object
#' @return list with 2 elements:
#'   \tabular{ll}{\code{xgrid} \tab vector of gridded distance values \cr
#'   \code{values} \tab vector of average detection function values at
#'    the \code{xgrid} values\cr}
#' @note Internal function called from plot functions for \code{ddf} objects
#' @author Jeff Laake
#' @keywords utility
average.line.cond.ll <- function(finebr,obs,model){
  xgrid <- NULL
  linevalues <- NULL

  # Depending on the type of model, setup data to be used for prediction
  # selecting the specified observer and only those detected (detected=1).
  if(model$method%in%c("io","trial")){
    newdat <- model$mr$mr$data
  }else{
    # if(model$method=="trial" | model$method=="trial.fi"){
    #   newdat=process.data(model$data,model$meta.data)$xmat
    #   newdat=newdat[newdat$observer!=obs & newdat$detected==1,]
    # }else
    if(model$method=="rem.fi"){
      newdat <- model$data
    }else{
      newdat <- model$mr$data
    }
  }

  newdat$offsetvalue <- rep(0,dim(newdat)[1])

  # For each element in the finebr grid, compute a value for x and the
  # appropriate p(x) averaged over all the covariate values
  for (i in 1:(length(finebr)-1)){
    # Compute x as the midpoint of the breaks and store as distance in dataframe
    x <- (finebr[i]+finebr[i+1])/2
    xgrid <- c(xgrid,x)
    newdat$distance <- rep(x,dim(newdat)[1])

    # Based on model compute p(x) from conditional detection function
    if(model$method=="io" | model$method=="trial"|model$method=="rem"){
      cond.det <- predict(model$mr,newdata=newdat,integrate=FALSE)
    }else{
      cond.det <- predict(model,newdata=newdat,integrate=FALSE)
    }

    if(model$method=="io" | model$method=="io.fi"|
       model$method=="rem" | model$method=="rem.fi"){
      p1 <- cond.det$p1
      p2 <- cond.det$p2
    }else{
      p1 <- cond.det$fitted
    }

    # Depending on observer compute average values for all observed
    # covariate values at specified distance (x)
    if(obs==1){
      linevalues <- c(linevalues,mean(p1))
    }else{
      linevalues <- c(linevalues,mean(p2))
    }
  }
  return(list(xgrid=xgrid,
              values=linevalues))
}

