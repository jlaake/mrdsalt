#' Golf tee data examples
#' @name tees
#' @author Jeff Laake
#' @examples
#' # Analysis of golf tee data from mrds package useing 4 alternative models
#' N=250
#'data(book.tee.data,package="mrds")
#'x=book.tee.data$book.tee.dataframe
#'results=ddf(data=x,mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance),
#'            meta.data=list(width=4),method="loglinear")
#'results$criterion
#'# %diff
#'100*(results$Nhat-N)/N
#'# unconditional detection fcts
#'par(mfrow=c(3,2))
#'plot(results,1:6,showpoints=FALSE)
#'# conditional detection functions
#'par(mfrow=c(1,2))
#'plot(results,7:8,showlines=FALSE)
#'
#'results=ddf(data=x,mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
#'            meta.data=list(width=4),method="bpi")
#'results$criterion
#'results$Nhat
#'# %diff
#'100*(results$Nhat-N)/N
#'# unconditional detection fcts
#'par(mfrow=c(3,2))
#'plot(results,1:6,showpoints=FALSE)
#'# conditional detection functions
#'par(mfrow=c(1,2))
#'plot(results,7:8,showlines=FALSE)
#'
#'results=ddf(data=x,mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
#'            meta.data=list(width=4),method="mr")
#'results$criterion
#'results$Nhat
#'# %diff
#'100*(results$Nhat-N)/N
#'# unconditional detection fcts
#'par(mfrow=c(3,2))
#'plot(results,1:6,showpoints=FALSE)
#'# conditional detection functions
#'par(mfrow=c(1,2))
#'plot(results,7:8,showlines=FALSE)
#'
NULL
