#' Create a dataframe with 2 records (one fore each observer) for each observation from a dataframe with a single record per observation contianing a capture history (ch), distance and
#' any covariate values.
#'
#' @param x a dataframe with a single record per observation
#' @param persons a dataframe containing a single record per observation and 2 columns indicating which person was in observer (position) 1 and 2 which provides the opportnity to estimate of personnel effects on detection
#' @return a dataframe with 2 records for each observation (format of mrds package)
#' @author Jeff Laake
#' @export
make.design.data=function(x,persons=NULL)
{
  dd=rbind(cbind(observer=rep(1,nrow(x)),x),cbind(observer=rep(2,nrow(x)),x))
  dd=dd[order(dd$object,dd$observer),]
  dd$observer=factor(dd$observer)
  if(!is.null(x$ch))dd$detected=as.numeric(as.vector(sapply(strsplit(x$ch,""),function(x) as.vector(x))))
  if(!is.null(persons))
  {
    if(is.dataframe(persons) & ncol(persons)==2&nrow(persons)==nrow(x))
    {
      dd$person=as.vector(t(person))
      dd$person=factor(dd$person)
    } else
      stop("Incorrect structure for persons. Must be dataframe with 2 columns and same number of rows in data")
  }
  rownames(dd)=1:nrow(dd)
  # this little tweak handles whether called with a dataframe with detected and another for simulation (no detected field) from sim_data
  if("detected"%in%names(dd))
      dd=dd[,c("object","observer","detected",names(dd)[!names(dd)%in%c("object","observer","detected")])]
  else
      dd=dd[,c("object","observer",names(dd)[!names(dd)%in%c("object","observer")])]
  return(dd)
}

