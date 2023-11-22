#' Computes values of logistic at a range of x for integration
#'
#' For a given z it computes logistic function with those covariate values for a set of distances
#'
#' @usage logisticbyx(distance, x1, x2, models, beta, gamma,indep=FALSE,PI=FALSE,use.offset=FALSE)
#' @param distance vector of distance values
#' @param x1 covariate data for observer1
#' @param x2 covariate data for observer2
#' @param models model list of p and delta formulas for detection probabilities
#' @param beta parameters for p
#' @param gamma parameters for delta dependence function
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @export
#' @return vector of values of logistic with specific covariate values over a set of distances
logisticbyx=function (distance, x1, x2, models, beta, gamma,indep=FALSE,PI=FALSE,use.offset=FALSE)
{
  xlist1 <- as.list(x1)
  xlist1$distance <-distance
  xlist2 <- as.list(x2)
  xlist2$distance <-distance
  x1 <- expand.grid(xlist1)
  x2 <- expand.grid(xlist2)
  xmat1=model.matrix(models$p.formula,x1)
  xmat2=model.matrix(models$p.formula,x2)
  p1=invlogit(xmat1,beta)
  p2=invlogit(xmat2,beta)
  if(indep)
    delta.values=1
  else
  {
    dmat=model.matrix(models$delta.formula,x1)
    delta.values=delta(dmat,p1,p2,gamma,PI,use.offset)
  }
  xx=p1+p2-delta.values*p1*p2
  xx[xx<0]=0
  xx[xx>1]=1
  return(xx)
}

#' Computes integral of logistic function from 0 to width
#'
#' Uses logisticbyx to compute integral from 0 to width for a specific set of covariate values
#'
#' @usage integratelogistic(x1, x2, models, beta, gamma, width,indep,PI,use.offset)
#' @param x1 covariate data for observer1
#' @param x2 covariate data for observer2
#' @param models model list of p and delta formulas for detection probabilities
#' @param beta parameters for p
#' @param gamma parameters for delta dependence function
#' @param width transect half-width
#' @param indep if TRUE, uses full independence model
#' @param PI if TRUE, uses point independence model
#' @param use.offset if TRUE use offset value for logit
#' @export
#' @return integral value
integratelogistic=function (x1, x2, models, beta, gamma, width,indep,PI,use.offset)
{
  integrate(logisticbyx,lower=0,upper=width,subdivisions=10, abs.tol=0.01,rel.tol=0.01,
            stop.on.error=FALSE,x1=x1,x2=x2, models=models, beta=beta, gamma=gamma,
            indep=indep,PI=PI,use.offset=use.offset)$value
}
