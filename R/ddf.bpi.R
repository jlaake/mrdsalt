#' Fits a Buckland et al (2010) model with possible dependence in detection for mark-recapture distance
#' data.
#' @param data dataframe containing double observer data (2 records per observation) with fields observer, detected, distance and possibly ch (capture history=10,01 or 11 for pair of observations
#' @param mrmodel list containing pformula (for detection) and dformula (for dependence)
#' @param meta.data list containing width and other fields not currently used
#' @param control list containing method (for optimization), debug (TRUE or FALSE), initial (optional starting values), indepm(TRUE if full independence), PI point independence if TRUE, use.offset and posdep.
#' @return list containing: fitted parameter values (par), AIC, optimx result with avep and NHat added, pformula and dformula.
#' @seealso ll.bpi
#' @author Jeff Laake
#' @export
#' @examples
#' # example code
#' par=c( 1.517424, 1.663175, -0.28, -0.5, 0.1311178)
#' W=10
#' N=1000
#' set.seed(1083821)
#' x=sim_data(N,W,par=par)
#' #create data structure with 2 records for each observation
#' dd=make.design.data(x)
#' #fit model
#' results=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),data=dd,method="bpi",meta.data=list(width=W),control=list(method="nlminb",debug=FALSE))
#' #N hat
#' Nhat=nrow(x)/results$mod$avep
#' Nhat
#' # %diff
#' 100*(N-Nhat)/N
#' # unconditonal detection fcts
#' par(mfrow=c(3,2))
#' plot(results,1:6,showpoints=FALSE)
#' # conditional detection functions
#' par(mfrow=c(1,2))
#' plot(results,7:8,showlines=FALSE)
#'

ddf.bpi=function(mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance),data,meta.data=meta.data,control=control,call="")
{
  save.options <- options()
  options(contrasts = c("contr.treatment", "contr.poly"))
  meta.data <- mrds:::assign.default.values(meta.data, left = 0, width = NA,
                                     binned = FALSE, int.range = NA, point = FALSE)
  control <- mrds:::assign.default.values(control, showit = 0, estimate = TRUE,
                                   refit = TRUE, nrefits = 25, initial = NA, lowerbounds = NA,
                                   upperbounds = NA, mono.points = 20,par=NULL,method="nlminb",debug=FALSE,
                                   indep=FALSE,PI=FALSE,use.offset=FALSE,posdep=FALSE)
  par=control$par
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
  # call p.bpi just to to check if length par is ok and to get parnames
  xx=p.bpi(par,data,p.formula=pformula,delta.formula=dformula,width=meta.data$width,indep=control$indep,PI=control$PI,use.offset=control$use.offset,
           posdep=control$posdep)
  # fit model using optimx
  mod=optimx(par=par,fn=ll.bpi,x=data,width=meta.data$width,method=control$method,
             p.formula=pformula,delta.formula=dformula,debug=control$debug,indep=control$indep,PI=control$PI,use.offset=control$use.offset,
             control=list(maxit=5000),posdep=control$posdep)
  # extract parameters from model output
  npar=attributes(mod)$npar
  par=as.vector(unlist(mod[paste("p",1:npar,sep="")]))
  names(par)=xx$parnames
  # compute average p and Nhat
  n=nrow(x)
  #without individual covariates
  if(all(unique(c(all.vars(pformula),all.vars(dformula)))%in%c("distance","observer")))
  {
    fitted=integrate(mu,lower=0,upper=meta.data$width,dd=data[1:2,],par=par,pformula=pformula,dformula=dformula)$value/meta.data$width
    mod$avep=fitted
    fitted=rep(fitted,n)
    mod$Nhat=n/mod$avep
  }
  else
  {
    # with individual covariates
    mod$Nhat=0
    fitted=NULL
    for(i in 1:n)
      fitted=c(fitted,integrate(mu,lower=0,upper=meta.data$width,par=par,dd=data[((i-1)*2+1):(i*2),],pformula=pformula,dformula=dformula)$value/meta.data$width)
    mod$Nhat=sum(1/fitted)
    mod$avep=n/mod$Nhat
  }
  # return result with parameter values, AIC, optimx object, formula and data that was used to fit model
  result$hessian=attributes(mod)$details[[3]]
  colnames(result$hessian)=names(par)
  rownames(result$hessian)=names(par)
  result$par=par
  result$lnl=-mod$value
  result$criterion=2*mod$value+2*length(par)
  result$fitted=fitted
  result$Nhat=mod$Nhat
  result$mod=mod
  class(result) <- c("loglinear", "ddf")
  return(result)
}


