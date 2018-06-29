##
### Pearson's Estimator of Dispersion parameter
##
disper.hat <- function(object) {
  family <- object$family
  link <- object$link
  np <- object$np
  dispersion <- object$control$dispersion
  b_mat <- object$beta
  rownames(b_mat) <- NULL
  c1X <- cbind(1,object$X)
  n <- dim(c1X)[1]
  y <-  object$y
  lac <- apply(abs(b_mat)>0,2,sum)
  eta <- tcrossprod(c1X,t(b_mat))
  dispers <- c()
  for (i in 1:np) {
    mu <- linkinv(drop(eta[,i]),family,link)
    varf_mu <- varf(mu,family)
    if (family=="binomial") {
      disper <- 1
    } else {
      if (is.null(dispersion)) {
        disper <- sum((y-mu)^2/varf_mu)/(n-ifelse((i==1 & b_mat[1,1]==0),0,lac[i]))  # Pearson's X2 statistic
      } else disper <- dispersion
    }
    dispers <- c(dispers,disper)
  }
  return(dispers)
}
