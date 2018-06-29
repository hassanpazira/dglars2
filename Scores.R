Scores <- function(X,y,X2,b,ac,dispersion,n,family,link) {
  #if (length(ac)==1) {eta <- rep(b[1],n);mu <- rep(mean(y),n)}
  eta <- drop(tcrossprod(X[,ac],t(b[ac])))  # Or drop(X[,ac]%*%as.matrix(b[ac]))   Or   drop(crossprod(t(X[,ac]),b[ac]))
  mu <- linkinv(eta,family,link)
  if (family=="binomial") {
    if (!all(is.finite(mu)) | !all(mu>0 & mu<1)) {
      #cat(paste("A",paste(ac[-1],collapse=" , "),sep=" = "));cat("\n")
      #cat(paste("mu",paste(mu,collapse=" , "),sep=" = "))
      stop("Invalid 'mu'! Valid mu: all(is.finite(mu)) & all(mu>0 & mu<1) !")
    }
  }
  if (family=="Gamma") {
    if (!all(is.finite(mu)) | !all(mu>0)) {
      #cat(paste("A",paste(ac[-1],collapse=" , "),sep=" = "));cat("\n")
      #cat(paste("mu",paste(mu,collapse=" , "),sep=" = "))
      stop("Invalid 'mu'! Valid mu: all(is.finite(mu)) & all(mu>0) !")
    }
  }
  mu.eta_mu <- mu.eta.mu(mu,family,link)            # Or we can use 'mu.eta.eta' function which is based on 'eta'
  varf_mu <- varf(mu,family)
  d2mu_deta2_mu <- d2mu_deta2.mu(mu,family,link)    # Or we can use 'd2mu_deta2.eta' function which is based on 'eta'
  d2teta_dmu2_mu <- d2teta_dmu2(mu,family,link)
  if (family=="binomial") {
    dispersion <- 1
    #dev <- -2*sum(dbinom(round(y),1,mu,log=TRUE))  # This deviance is: -2*loglikelihood, also we can use:  dev <- -2*(sum(log(mu[y>0.5]))+sum(log(1-mu[y<0.5]))). They are similar because dispersion=1 (!!!).
    res.dev <- -2*sum(log(mu[y>0.5]))-2*sum(log(1-mu[y<0.5]))  # res.dev and the scaled deviance (which is res.dev/dispersion) are the same because dispersion=1.
  } else {
    if (is.null(dispersion)) dispersion <- 1  # can be omitted !!!
    if (family=="Gamma") {
      #dev <- -2*sum(dgamma(y,1/dispersion,scale=mu*dispersion,log=TRUE))  # 'dev' means '-2*logliklihood'.
      res.dev <- 2*sum(log(mu/y)+(y-mu)/mu)             #  Residual Deviance means '-2*phi*(l(mu_hat)-l(y))', which is phi*scaled deviance.
    }
    if (family=="inverse.gaussian") {
      #dev <- n*log(dispersion*2*pi)+3*sum(log(y))-2*sum(1/mu)/dispersion+sum(y/(mu^2))/dispersion+sum(1/y)/dispersion  # 'dev' means '-2*logliklihood'
      res.dev <- sum(((y-mu)^2)/((mu^2)*y))             #  Residual Deviance means '-2*phi*(l(mu_hat)-l(y))', which is phi*scaled deviance.
    }
  }
  w.for.dl <- mu.eta_mu/(varf_mu*dispersion) 
  dl <- as.numeric(apply(X*w.for.dl*(y-mu),2,sum))         # dl <- drop(t(X) %*% diag(w.for.dl) %*% (y-mu))    Or    dl <- Scorefun(X,y,family,b)
  #dla <- as.numeric(apply(as.matrix(X[,ac])*w.for.dl*(y-mu),2,sum))
  w.for.Im <- ((mu.eta_mu)^2)/(varf_mu*dispersion)
  Imm <- as.numeric(apply(w.for.Im*X2,2,sum))              # Imm <- diag(t(X) %*% diag(w.for.Im) %*% X)   Or  diag(crossprod(sqrt(w.for.Im)*X))   Or   Im=XTWX=t(glm()$R)%*%glm()$R
  #Imma <- as.numeric(apply(w.for.Im*as.matrix(X2[,ac]),2,sum))
  w.for.dmnl <- ((d2teta_dmu2_mu*((mu.eta_mu)^2)*varf_mu+d2mu_deta2_mu)*(y-mu)-(mu.eta_mu)^2)/(varf_mu*dispersion)
  w.for.dmIn <- (d2teta_dmu2_mu*((mu.eta_mu)^3)+2*(mu.eta_mu)*(d2mu_deta2_mu)/varf_mu)/dispersion
  rao <- dl/sqrt(Imm)
  rao.and.rest <- list(eta=eta,mu=mu,mu.eta_mu=mu.eta_mu,varf_mu=varf_mu,d2mu_deta2_mu=d2mu_deta2_mu,
                       d2teta_dmu2_mu=d2teta_dmu2_mu,dispersion=dispersion,rao=rao,dl=dl,res.dev=res.dev,
                       w.for.dl=w.for.dl,Imm=Imm,w.for.Im=w.for.Im,w.for.dmnl=w.for.dmnl,w.for.dmIn=w.for.dmIn)
}
