deviance <- function(object) {
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
  devi <- c()
  for (i in 1:np) {
    mu <- linkinv(drop(eta[,i]),family,link)
    varf_mu <- varf(mu,family)
    if (family=="binomial") {
      dev <- -2*sum(dbinom(round(y),1,mu,log=TRUE))  # This deviance is: -2*loglikelihood, also we can use:  dev <- -2*(sum(log(mu[y>0.5]))+sum(log(1-mu[y<0.5]))). They are similar because dispersion=1 (!!!).
    } else {
      if (is.null(dispersion)) {
        disper <- sum((y-mu)^2/varf_mu)/(n-ifelse((i==1 & b_mat[1,1]==0),0,lac[i]))  # Pearson's X2 statistic
      } else disper <- dispersion
      if (family=="Gamma") {
        dev <- -2*sum(dgamma(y,1/disper,scale=mu*disper,log=TRUE))  # 'dev' means '-2*logliklihood'.
      }
      if (family=="inverse.gaussian") {
        dev <- n*log(disper*2*pi)+3*sum(log(y))-2*sum(1/mu)/disper+sum(y/(mu^2))/disper+sum(1/y)/disper   # 'dev' means '-2*logliklihood'
      }
    }
    devi <- c(devi,dev)
  }
  return(devi)
}
