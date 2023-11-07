#' Delta computation for Buckland et al 2010 model
#'
#' Computes delta values for all values in design matrix, See Buckland et al, (2010) for
#' definition
#'
#' @usage delta(x,p1,p2,gamma,PI=FALSE,use.offset=FALSE)
#' @param x a design matrix for delta
#' @param p1 detection probabilities for observer 1
#' @param p2 detection probabilities for observer 2
#' @param PI TRUE if this is a point independence model
#' @param use.offset if TRUE use offset value for logit
#' @export
#' @return Vector of delta values for each row in design matrix x
delta=function(x,p1,p2,gamma,PI=FALSE,use.offset=FALSE)
{
  lower=apply(cbind(p1,p2),1,function(x)return(max(0,(x[1]+x[2]-1)/(x[1]*x[2]))))
  upper=apply(cbind(p1,p2),1,function(x)return(min(1/x[1],1/x[2])))
  if(PI | use.offset)
  {
    delta.values=rep(1,length(p1))
    offset=rep(0,length(p1))
    ix=(upper<1.0000001 & upper>.9999999) | ( lower<1.0000001 & lower>.9999999)
    offset[!ix]=log((1-lower[!ix])/(upper[!ix]-1))
    xx=(upper-lower)*invlogit(x,gamma,offset)+lower
    delta.values[!ix]=xx[!ix]
    return(delta.values)
  }
  else
  {
    offset=NULL
    return((upper-lower)*invlogit(x,gamma,offset)+lower)
  }
}
