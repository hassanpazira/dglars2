summary.dglars2 <- function(object,k=c("BIC","AIC"),complexity=c("df","gdf"),digits = max(3, getOption("digits") - 3),...){
  if (!inherits(object,"dglars2"))
    stop("This function only works for objects of class 'dglars2'")
  if(is.numeric(k)){
    if(k<=0) stop("k must be greater than zero")
    knm <- "GoF"
  }
  else{
    knm <- match.arg(k)
    k <- ifelse(knm == "BIC",log(dim(object$X)[1]),2)
  }
  complexity <- match.arg(complexity)
  tbl <- make_summary_table(object,k,complexity)
  names(tbl$table)[5] <- complexity
  names(tbl$table)[6] <- knm
  action <- object$action
  id <- which(action!="")
  n.tbl <- dim(tbl$table)[1]
  n.space <- length(id)
  id.tbl <- vector(mode="numeric",length=n.tbl+n.space)
  id.space <- id + seq(1,n.space)
  id.tbl[-id.space] <- seq(1:n.tbl)
  id.tbl[id.space] <- id
  tbl.format <- format(tbl$table[id.tbl,], digits = digits)
  tbl.format[id.space-1,1] <- ""
  tbl.format[id.space,-1] <- ""
  b.gof.names <- names(tbl$b.gof)[-1]
  if(!length(b.gof.names)) b.gof.names <- "1"
  best.model.formula <- paste("y",paste(b.gof.names,collapse=" + "),sep=" ~ ")
  cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  print.data.frame(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
  cat("\n=====================================================\n")
  cat("\nBest model identified by",knm,"criterion ( k =",k,"and complexity =",complexity,"):\n\n", best.model.formula,"\n")
  cat("\nCoefficients:\n\n")
  print.default(format(tbl$b.gof, digits = digits), print.gap = 2, quote = FALSE,...)
  cat("\n",knm,": ",format(min(tbl$table[6]), digits = digits))
  cat("\n\n=====================================\n\nMethod    :  \"",object$control$method,"\"")
  if (object$conv==0) cat("\n\nAlgorithm :  \"",object$control$algorithm,"\"  ( converged! )\n\n")
  else cat("\n\nAlgorithm :  \"",object$control$algorithm,"\"  ( conv =",object$conv,")\n\n")
  invisible(tbl)
}
