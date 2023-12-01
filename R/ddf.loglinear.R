#' Fits a log-linear model with possible dependence in detection for mark-recapture distance
#' data.
#' @param data dataframe containing double observer data (2 records per observation) with fields observer, detected, distance and possibly ch (capture history=10,01 or 11 for pair of observations
#' @param mrmodel list containing pformula (for detection) and dformula (for dependence)
#' @param meta.data list containing width and other fields not currently used
#' @param control list containing method (for optimization), debug (TRUE or FALSE), initial (optional starting values) and other fields not currently used
#' @param call call syntax to ddf
#' @return list containing: fitted parameter values (par), AIC, optimx result with avep and NHat added, pformula and dformula.
#' @seealso ll.loglinear
#' @author Jeff Laake
#' @examples
#' # example code
#' par=c( 1.517424, 1.663175, -0.28, -0.5, 0.1311178)
#' W=10
#' N=1000
#' set.seed(1083821)
#' x=sim.loglinear(N,W,par=par)
#' #create data structure with 2 records for each observation
#' dd=make.design.data(x)
#' #fit model
#' results=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
#' data=dd,method="loglinear",meta.data=list(width=W),
#' control=list(method="nlminb",debug=FALSE))
#' results$criterion
#' # %diff
#' 100*(results$Nhat-N)/N
#' # unconditonal detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#'
#' set.seed(1083821)
#' x=sim.bpi(x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2),
#'  distance=rep(runif(N/2,0,W),each=2)),
#'  par=c( 1.517424, 1.663175, -0.28, -0.5, 0.1311178),
#'  pformula=~observer*distance,dformula=~-1+distance,PI=TRUE)
#'  #' #fit model
#' results=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
#' data=x,method="loglinear",meta.data=list(width=W),
#' control=list(method="nlminb",debug=FALSE))
#' results$criterion
#' # %diff
#' 100*(results$Nhat-N)/N
#' # unconditonal detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#'
# fit mrds data with log linear model
ddf.loglinear=function(mrmodel,data,meta.data,control,call)
{
  save.options <- options()
  options(contrasts = c("contr.treatment", "contr.poly"))
  meta.data <- mrds:::assign.default.values(meta.data, left = 0, width = NA,
                                     binned = FALSE, int.range = NA, point = FALSE)
  control <- mrds:::assign.default.values(control, showit = 0, estimate = TRUE,
                                   refit = TRUE, nrefits = 25, initial = NA, lowerbounds = NA,
                                   upperbounds = NA, mono.points = 20,par=NULL,method="nlminb",debug=FALSE)
  if(!is.na(control$initial))
    par=control$initial
  else
    par=NULL
  data.list <- mrds:::process.data(data, meta.data)
  meta.data <- data.list$meta.data
  result <- list(call = call, data = data.list$xmat, mrmodel = mrmodel,
                 meta.data = meta.data, control = control,
                 method = "loglinear")
  pformula=mrmodel$pformula
  dformula=mrmodel$dformula
  # check to make sure all variables in formula are in the data
  vars=unique(c(all.vars(pformula),all.vars(dformula)))
  if(any(!vars%in%colnames(data)))
    stop("variables missing from dataframe: ",vars[!vars%in%colnames(data)])
  # no parameter values given use glm to create initial values of pformula and use 0 for each column in dformula
  if(is.null(par))
    par=c(glm(as.formula(paste("detected~",as.character(pformula)[2],sep="")),family="binomial",data=data)$coefficients,
          rep(0.1,ncol(model.matrix(dformula,data[seq(1,nrow(data),2),]))))
  names(par)=NULL
  # construct 0/1 matrix with column 1 being 10 capture history, 2 being 01 and 3 being 11; it is 1 if that
  # object has the particular capture history (ch); only a single 1 in each row
  #
  # if ch is not in data create it
  if(is.null(data$ch))data$ch=paste(data$detected[seq(1,nrow(data),2)],data$detected[seq(2,nrow(data),2)],sep="")
  # if ch is 00 remove record
  if(any(data$ch=="00"))cat("00 records")
  data=data[data$ch!="00",]
  # only need one record from each data pair
  x=data[seq(1,nrow(data),2),]
  n=nrow(x)
  cmat=matrix(0,ncol=3,nrow=n)
  cmat[,1]=as.numeric(x$ch=="10")
  cmat[,2]=as.numeric(x$ch=="01")
  cmat[,3]=as.numeric(x$ch=="11")
  # call detprobs just to to check if length par is ok; only really needed if user provided par vector
  xx=p.loglinear(par,data,pformula,dformula)
  # fit model using optimx
  mod=optimx(par=par,fn=ll.loglinear,cmat=cmat,dd=data,W=meta.data$width,method=control$method,
             pformula=pformula,dformula=dformula,debug=control$debug,control=control$control)
  if(mod$convcode!=0)warning("\nmodel did not converge")
   # extract parameters from model output
  npar=attributes(mod)$npar
  par=as.vector(unlist(mod[paste("p",1:npar,sep="")]))
  if(any(is.na(par)))stop("invalid parameter values")
  names(par)=xx$parnames
  # return result with parameter values, AIC, optimx object, formula and data that was used to fit model
  result$hessian=attributes(mod)$details[[3]]
  colnames(result$hessian)=names(par)
  rownames(result$hessian)=names(par)
  result$par=par
  result$lnl=-mod$value
  result$criterion=2*mod$value+2*length(par)
  class(result) <- c("loglinear", "ddf")
  result$fitted=predict(result)
  result$Nhat=sum(1/result$fitted)
  result$avep=n/result$Nhat
  result$mod=mod
  return(result)
}

