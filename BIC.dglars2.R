BIC.dglars2 <- function(object,complexity=c("df","gdf"), ...){
  if (!inherits(object,"dglars2"))
    stop("This function only works for objects of class 'dglars2'")
  complexity <- match.arg(complexity)
  dev <- deviance(object)
  if (complexity=="df")  compl <- object$df
  else compl <- gdf(object)
  if (object$family=="binomial" | !is.null(object$control$dispersion)) gof <- dev + log(dim(object$X)[1]) * compl
  else gof <- dev + log(dim(object$X)[1]) * (compl+1)
  return(gof)
}
