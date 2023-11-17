#' Compute detection probabilities for each observable outcome (10,01,11) and K for log-linear model with possible dependence in detections for mark-recapture distance data
#' data.
#' @param par parameter values
#' @param x dataframe containing double observer data (2 records per observation) with fields observer, detected, distance, ch and any covariate values
#' @param pformula formula for detection but excluding the dependence portion
#' @param dformula model for dependence due to residual (unmodelled) heterogeneity
#' @return list with prob: probability matrix with 3 columns for 10,01,11 at each of n rows of observations and K the denominator
#' @author Jeff Laake
#' @export
p.ll=function(par,x,pformula,dformula)
{
  # compute design matrix for observer probability part of model
  mm=model.matrix(pformula,x)
  # compute design matrix for dependence part of model but for only 1 value for each pair
  if(dformula==as.formula("~0"))
  {
    dm=matrix(1,nrow=nrow(x[seq(1,nrow(x),2),]),ncol=1)
    colnames(dm)="0"
  } else
    dm=model.matrix(dformula,x[seq(1,nrow(x),2),])
  # produce error if length of par does not match formulas
  if(dformula==as.formula("~0"))
  {
    if(length(par)!=ncol(mm))stop("Mismatch between length of par and formulas")
  } else
    if(length(par)!=ncol(mm)+ncol(dm))stop("Mismatch between length of par and formulas")
  # extract the theta parameter values
  theta=par[1:ncol(mm)]
  # the remainder are for the dependence part of the model
  if(dformula==as.formula("~0"))
    delta=0
  else
    delta=par[(length(theta)+1):length(par)]
  # compute numerator of 10,01,11 probabilities
  prob=cbind(exp(mm[seq(1,nrow(x),2),]%*%theta),exp(mm[seq(2,nrow(x),2),]%*%theta),
             exp(sapply(split(mm%*%theta,x$object), sum)+dm%*%delta))
  # compute denominator return probabilities and denominator
  K=1+rowSums(prob)
  return(list(prob=prob/K,K=K,parnames=c(paste("p:",colnames(mm)),paste("delta:",colnames(dm)))))
}


