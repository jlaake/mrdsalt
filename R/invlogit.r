#' Inverse logit function with and without an offset
#'
#' Computes inverse logit values with and without an offset
#'
#' @usage invlogit(x,beta,offset=NULL)
#' @param x a design matrix
#' @param beta parameters for logit function
#' @param offset offset value for logit
#' @export
#' @return Vector of inverse logit with a value for each row in design matrix x
invlogit=function(x,beta,offset=NULL)
{
  if(is.null(offset))
    return(1/(1+exp(-x%*%beta)))
  else
    return(1/(1+exp(-x%*%beta-offset)))
}
