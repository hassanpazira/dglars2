print.dglars2 <- function (x,digits = max(3, getOption("digits") - 3), ...){
  if (!inherits(x,"dglars2"))
    stop("This function only works for objects of class 'dglars2'")
  action <- x$action
  g <- x$g
  #dispersion <- disper.hat(object)
  res.dev <- x$res.dev
  dev.ratio <- 1 - res.dev/res.dev[1]
  df <- x$df
  #tbl <- data.frame(action,g,dispersion,res.dev,dev.ratio,df)
  tbl <- data.frame(action,g,res.dev,dev.ratio,df)
  names(tbl) <- c("Sequence","g","Res.Dev","%Dev","df")
  id <- which(action!="")
  n.tbl <- dim(tbl)[1]
  n.space <- length(id)
  id.tbl <- vector(mode="numeric",length=n.tbl+n.space)
  id.space <- id + seq(1,n.space)
  id.tbl[-id.space] <- seq(1:n.tbl)
  id.tbl[id.space] <- id
  tbl.format <- format(tbl[id.tbl,], digits = digits)
  tbl.format[id.space-1,1] <- ""
  tbl.format[id.space,-1] <- ""
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  print.data.frame(tbl.format,print.gap = 2, quote = FALSE,row.names=FALSE, ...)
  cat("\n\n=====================================\n\nMethod    :  \"",x$control$method,"\"")
  if (x$conv==0) cat("\n\nAlgorithm :  \"",x$control$algorithm,"\"  ( converged! )\n\n")
  else cat("\n\nAlgorithm :  \"",x$control$algorithm,"\"  ( conv =",x$conv,")\n\n")
  invisible(tbl)
}
