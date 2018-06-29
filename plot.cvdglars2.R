plot.cvdglars2 <- function(x,...){
  if (!inherits(x,"cvdglars2"))
    stop("This function only works for objects of class 'cvdglars2'")
  dev_m <- x$dev_m
  dev_v <- x$dev_v
  nfold <- x$control$nfold
  k <- qt(0.975,nfold-1)
  dev_up <- dev_m + k * sqrt(dev_v/nfold)
  dev_low <- dev_m - k * sqrt(dev_v/nfold)
  ng <- x$control$ng
  g <- seq(x$g_max,x$g0,length=ng)
  g_hat <- g[which.min(dev_m)]
  df <- x$df  # sum(abs(x$beta)>0)
  if (any(is.na(dev_up)) | any(is.na(dev_low)) | any(dev_up==Inf) | any(dev_low==Inf)) {
    len_err<-c(sum(is.na(dev_up),na.rm=T),sum(is.na(dev_low),na.rm=T),sum(dev_up==Inf,na.rm=T),sum(dev_low==Inf,na.rm=T))
    dev_up <- dev_up[-c((length(dev_up)-max(len_err)+1):length(dev_up))]
    dev_low <- dev_low[-c((length(dev_low)-max(len_err)+1):length(dev_low))]
    gg <- g[-c((length(g)-max(len_err)+1):length(g))]
  } else {
    gg <- g
  }
  plot(g,dev_m,xlab=expression(gamma),ylab="Deviance",ylim=c(min(dev_low),max(dev_up)),pch=20,type="n",main="Cross-Validation Deviance",...)
  segments(x0=gg,y0=dev_low,y1=dev_up,col=8,lty=2)
  points(gg,dev_low,pch="-",col="deepskyblue")
  points(g,dev_m,pch=20,col="blue")
  points(gg,dev_up,pch="-",col="deepskyblue")
  abline(v=g_hat,col=2,lty=2,lwd=2)
  axis(3,at=g_hat,labels=paste("df = ",df,sep=""),padj=1)
}
