#' Probability calculations for mark-recapture type Huggins model
#'
#' For a set of observations, computes observer specific probabilities, integrals and delta values See Buckland et al, (2010) for
#' definitions
#'
#' @param par model parameter values
#' @param x observation dataframe
#' @param pformula formula for detection probabilities
#' @param dformula formula for delta dependence function
#' @param indep if TRUE, uses full independence model
#' @export
#' @return list of vectors for p1,p2, integrals and delta values
p.mr=function(par,x,pformula,dformula,indep)
{
  # create design matrix for p and delta
  xmat=model.matrix(pformula,x)
  dmat=model.matrix(dformula,x[x$observer==1,])
  parnames=c(paste("p:",colnames(xmat)),paste("delta:",colnames(dmat)))
  if(colnames(dmat)[1]=="(Intercept)") stop("\nError: No intercept allowed with dependence model\n")
  # extract parameter vectors: beta for detection and gamma for delta
  beta=par[1:dim(xmat)[2]]
  # compute probabilities and delta values
  px=plogis(xmat%*%beta)
  p1=px[seq(1,length(px),2)] # probability first observer makes a detection
  p20=px[seq(2,length(px),2)] # probability second observer detects when first observer did not
  if(indep)
  {
    p21=p20  # they are the same under independence
  } else
  {
    gamma=par[(1+dim(xmat)[2]):length(par)]
    p21=plogis(log(p20/(1-p20))+dmat%*%gamma)  # probability second observer detects when first observer did
                                               # this is c in MARK lingo - the recapture probability
  }
  return(list(p1=p1,p20=p20,p21=p21,parnames=parnames))
}

