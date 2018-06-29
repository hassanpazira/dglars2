plot.dglars2 <- function(x,k=c("AIC & BIC","BIC","AIC"),complexity=c("df","gdf"),g.gof=NULL,...){
  if (!inherits(x,"dglars2"))
    stop("This function only works for objects of class 'dglars2'")
  if(is.numeric(k)){
    if(k <= 0) stop("k must be greater than zero")
    knm <- "GoF"
    k1 <- k2 <- k
  }
  else{
    knm <- match.arg(k)
    k1 <- 2
    k2 <- log(dim(x$X)[1])
  }
  complexity <- match.arg(complexity)
  n <- dim(x$X)[1]
  beta <- x$beta
  rao <- x$rao
  dev <- deviance(x)
  g <- x$g
  g.action <- g[x$action!=""]
  if(is.null(g.gof)){
    if(complexity == "df")  {
      df <- x$df
    } else {
      df <- gdf(x)
    }
    gof1 <- dev + k1 * df
    gof2 <- dev + k2 * df
    g.gof2 <- g[which.min(gof2)]
    g.gof1 <- g[which.min(gof1)]
    if (knm=="AIC & BIC") {
      plot(g,gof1,ylim=c(min(gof1),max(gof2)),xlab=expression(gamma),ylab=knm,type="n",main="Model Selection Criteria")  # ylab="AIC (Blue) & BIC (Red)"
      abline(v=g.action,lty=3,col=8)
      if (g.gof1==g.gof2) {
        axis(3,g.gof1,"minAIC=minBIC",padj=1)
        abline(v=g.gof1,lty=2,col='purple',lwd=2)
      } else {
        axis(3,g.gof1,"AIC",padj=1)  # expression(paste("AIC",circ))
        axis(3,g.gof2,"BIC",padj=1)  # expression(paste("BIC",bullet))
        abline(v=g.gof2,lty=2,col=2,lwd=2)
        abline(v=g.gof1,lty=2,col='blue',lwd=2)
      }
      points(g,gof1,type="o",pch=21,lty=2,...)
      points(g,gof2,type="o",pch=20,lty=1,...)
      #legend(mean(g),max(gof2),"bot", "(x,y)", pch = 1)
    }
    if (knm=="AIC") {
      plot(g,gof1,xlab=expression(gamma),ylab=knm,type="n",main="Model Selection Criterion")
      axis(3,g.gof1,knm,padj=1)
      abline(v=g.action,lty=2,col=8)
      abline(v=g.gof1,lty=2,col=2,lwd=2)
      points(g,gof1,type="o",pch=20,lty=1,...)
    }
    if (knm=="BIC" | knm=="GoF") {
      plot(g,gof2,xlab=expression(gamma),ylab=knm,type="n",main="Model Selection Criterion")
      if (knm=="GoF") {
        axis(3,g.gof2,paste("GoF ( k =",k,")"),padj=1)
      } else {
        axis(3,g.gof2,knm,padj=1)
      }
      abline(v=g.action,lty=2,col=8)
      abline(v=g.gof2,lty=2,col=2,lwd=2)
      points(g,gof2,type="o",pch=20,lty=1,...)
    }
    op <- par(ask=dev.interactive())
  } else knm <- "g.gof"
  matplot(g,t(beta[-1,]),col=1,type="n",xlab=expression(gamma),ylab="Regression Coefficients",main="Coefficients Path")
  abline(v=g.action,lty=3,col=8)
  if (is.null(g.gof)) {
    if (knm=="GoF") {
      abline(v=g.gof2,lty=2,col=2,lwd=2)
      axis(3,g.gof2,knm,padj=1)
    } else {
      if (g.gof1==g.gof2) {
        axis(3,g.gof1,"minAIC=minBIC",padj=1)
        abline(v=g.gof1,lty=2,col='purple',lwd=2)
      } else {
        axis(3,g.gof1,"AIC",padj=1)
        axis(3,g.gof2,"BIC",padj=1)
        abline(v=g.gof2,lty=2,col=2,lwd=2)
        abline(v=g.gof1,lty=2,col='blue',lwd=2)
      }
    }
  } else {
    abline(v=g.gof,lty=2,col=2,lwd=2)
    axis(3,g.gof,paste(knm,"=",g.gof),padj=1)
  }
  matpoints(g,t(beta[-1,]),col=1,type="l",lty=1,...)
  if(!is.null(g.gof)) op <- par(ask=dev.interactive())
  if(x$control$algorithm %in% c("pc","ipc")){
    matplot(g,t(rao),col=1,type="n",xlab=expression(gamma),ylab="| Rao Score Statistics |",main="Rao Score Path")
    abline(v=g.action,lty=3,col=8)
    if (is.null(g.gof)) {
      if (knm=="GoF") {
        abline(v=g.gof2,lty=2,col=2,lwd=2)
        axis(3,g.gof2,knm,padj=1)
      } else {
        if (g.gof1==g.gof2) {
          axis(3,g.gof1,"minAIC=minBIC",padj=1)
          abline(v=g.gof1,lty=2,col='purple',lwd=2)
        } else {
          axis(3,g.gof1,"AIC",padj=1)
          axis(3,g.gof2,"BIC",padj=1)
          abline(v=g.gof2,lty=2,col=2,lwd=2)
          abline(v=g.gof1,lty=2,col='blue',lwd=2)
        }
      }
    } else {
      abline(v=g.gof,lty=2,col=2)
      axis(3,g.gof,paste(knm,"=",g.gof),padj=1)
    }
    matpoints(g,t(rao),col=1,type="l",lty=1,...)
  }
  par(op)
}
