gdf <- function(object){
  if (!inherits(object,"dglars2"))
    stop("This function only works for objects of class 'dglars2'")
  family <- object$family
  link <- object$link
  n <- dim(object$X)[1]
  y <- object$y
  X <- cbind(rep(1,n),object$X) # with intercept
  beta <- object$beta
  np <- object$np
  gdf_n <- vector(mode="numeric",length=np)
  if (family=="Gamma") {
    out.glm <- switch(link,
                      inverse=glm(y~X[,-1],family=Gamma("inverse")),
                      log=glm(y~X[,-1],family=Gamma("log")),
                      identity=glm(y~X[,-1],family=Gamma("identity"))
    )
  }
  if (family=="inverse.gaussian") {
    out.glm <- switch(link,
                      "1/mu^2"=glm(y~X[,-1],family=inverse.gaussian("1/mu^2")),
                      inverse=glm(y~X[,-1],family=inverse.gaussian("inverse")),
                      log=glm(y~X[,-1],family=inverse.gaussian("log")),
                      identity=glm(y~X[,-1],family=inverse.gaussian("identity"))
    )
    
  }
  if (family=="binomial") {
    out.glm <- switch(link,
                      logit=glm(y~X[,-1],family=binomial("logit")),
                      probit=glm(y~X[,-1],family=binomial("probit")),
                      cloglog=glm(y~X[,-1],family=binomial("cloglog")),
                      log=glm(y~X[,-1],family=binomial("log"))
    )
  }
  mu_for_MLE <- out.glm$fit
  V_for_MLE <- varf(mu_for_MLE,family)
  if(!out.glm$converged){
    warning("complexity was set equal to 'df'")
    gdf_n <- colSums(abs(beta) > 0)  # Or beta!=0
  } else {
    for(i in 1:np){
      b <- beta[,i,drop=TRUE]
      ac <- which(abs(b)>0) #
      eta <- drop(tcrossprod(X[,ac],t(b[ac])))
      mu <- linkinv(eta,family,link)
      V <- varf(mu,family)
      mu.eta_mu <- mu.eta.mu(mu,family,link)
      r <- y-mu
      w.for.Imn <- V_for_MLE*((mu.eta_mu/V)^2)
      Imn <- crossprod(sqrt(w.for.Imn)*X[,ac])       #  Or  t(X[,ac])%*%diag(w.for.Imn)%*%X[,ac]
      w.for.dmnl <- ((d2teta_dmu2(mu,family,link)*((mu.eta_mu)^2)*V+d2mu_deta2.mu(mu,family,link))*r-(mu.eta_mu)^2)/V
      dmnl <- t(X[,ac])%*%diag(w.for.dmnl)%*%X[,ac]  #  != crossprod(sqrt(w.for.dmnl)*X[,ac])  because 'w.for.dmnl' can be negative!
      gdf_n[i] <- drop(sum(solve(-dmnl)*Imn))        #  Or  drop(sum(diag(crossprod(t(inv.dmnl),Imn))))  Or drop(crossprod(as.vector(inv.dmnl),as.vector(Imn)))	
    }
  }
  return(gdf_n)
}
