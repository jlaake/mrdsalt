#' Plot fit of detection functions and histograms of data from distance
#' sampling independent observer (\code{mr}) model
#'
#' Plots the fitted detection functions for a distance sampling model and
#' histograms of the distances (for unconditional detection functions) or
#' proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.
#'
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}.  The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#'
#' It is not intended for the user to call \code{plot.bpi} but its arguments
#' are documented here. Instead the generic \code{plot} command should be used
#' and it will call the appropriate function based on the class of the
#' \code{ddf} object.
#'
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced.
#'  \tabular{ll}{1 \tab Plot primary unconditional detection function \cr
#'               2 \tab Plot secondary unconditional detection function \cr
#'               3 \tab Plot pooled unconditional detection function \cr
#'               4 \tab Plot duplicate unconditional detection function \cr
#'               5 \tab Plot unconditional for seen by primary only\cr
#'               6 \tab Plot unconditional for seen by secondary only \cr
#'               7 \tab Plot primary conditional detection function\cr
#'               8 \tab Plot secondary conditional detection function \cr}
#'  Note that the order of which is ignored and plots are produced in the above
#'  order.
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param maintitle main title line for each plot
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param ylim range of vertical axis; defaults to (0,1)
#' @param angle shading angle for histogram bars.
#' @param density shading density for histogram bars.
#' @param col colour for histogram bars.
#' @param jitter scaling option for plotting points. Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
#' @param divisions number of divisions for averaging line values; default = 25
#' @param pages the number of pages over which to spread the plots. For
#'  example, if \code{pages=1} then all plots will be displayed on one page.
#'  Default is 0, which prompts the user for the next plot to be displayed.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return Just plots
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
plot.mr=function (x, which = 1:6, breaks = NULL, nc = NULL, maintitle = "",
                         showlines = TRUE, showpoints = TRUE, ylim = c(0, 1), angle = NULL,
                         density = NULL, col = "lightgrey", jitter = NULL, divisions = 25,
                         pages = 0, xlab = "Distance", ylab = "Detection probability",
                         subtitle = TRUE, ...)
{
  model <- x
  width <- model$meta.data$width
  left <- model$meta.data$left
  which <- sort(which)
  xmat <- model$data
  p.list=p.mr(par=model$par,xmat,model$mrmodel$pformula,model$mrmodel$dformula,indep=model$control$indep)
  p1=p.list$p1
  p20=p.list$p20
  p21=p.list$p21
  #  11
  p11=p1*p21
  #  10
  p10=p1*(1-p21)
  #  01
  p01=(1-p1)*p20
  delta.values=p21/p20
  gxlist <- list(p1, p20+p21, (p10+ p01 + p11),
                 p11,p10,p01)
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
  #  oask <- plot_layout(which, pages)
  #  on.exit(devAskNewPage(oask))
  for (wh in which[which < 7]) {
    plot_uncond(model, obs=wh, xmat=xmat, gxvalues = gxlist[[wh]],
                nc, finebr = (width/divisions) * (0:divisions), breaks,
                showpoints, showlines, maintitle, ylim, angle = angle,
                density = density, col = col, jitter = jitter, xlab = xlab,
                ylab = ylab, subtitle = subtitle, ...)
  }
  #data <- model$mr$mr$data
  #data$offsetvalue <- 0
  if (is.element(7, which)) {
    gxvalues <- (p11/(p01+p11))[xmat$detected[xmat$observer == 2] == 1]
    plot_cond(1, xmat, gxvalues, model, nc, breaks, finebr = (width/divisions) *
                (0:divisions), showpoints, showlines, maintitle,
              ylim, angle = angle, density = density, col = col,
              jitter = jitter, xlab = xlab, ylab = ylab, subtitle = subtitle,
              ...)
  }
  if (is.element(8, which)) {
    gxvalues <- (p21)[xmat$detected[xmat$observer == 1] == 1]
    plot_cond(2, xmat, gxvalues, model, nc, breaks, finebr = (width/divisions) *
                (0:divisions), showpoints, showlines, maintitle,
              ylim, angle = angle, density = density, col = col,
              jitter = jitter, xlab = xlab, ylab = ylab, subtitle = subtitle,
              ...)
  }
  invisible(NULL)
}
