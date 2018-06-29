cvdglars2.fit <- function(X,y,family=c("Gamma","inverse.gaussian","binomial"),link=c("inverse","1/mu^2","log","identity","logit","probit","cloglog"),control=list()){
  m.call <- match.call()
  family <- match.arg(family)
  link <- match.arg(link)
  if (!family %in% c("Gamma","inverse.gaussian","binomial")) {ifelse (family %in% c("poisson"), warning("For Poisson distributions use 'dglars' package!"), warning("This distribution is not 'Gamma' or 'Inverse Gaussian'!"))}
  if(is.data.frame(X)) X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if(is.null(colnames(X))) colnames(X) <- paste("x",1:p,sep="")
  min_np <- min(n-1,p)
  if (NCOL(y) != 1) stop("for the 'binomial' family, y must be a vector, not a matrix")
  if (length(y)!=n) stop("length of 'y' != sample size 'n'")
  if (family=="binomial") {
    if (is.character(y)) y <- factor(y)
    if (is.factor(y)){
      if (nlevels(y)!=2) stop("only factors with two levels are allowed")
      yy <- as.numeric(y) - 1
    } else {
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      if (any(abs(y - round(y)) > 0.001)) warning("non-integer #successes in a binomial glm!") #
      yy <- y
    }
  } else {
    if (is.character(y)) stop("It's not allowed to use a character or a factor as a response vector for this family (0 < y < Inf)!")
    if (any(y<=0 | is.na(y))) stop("Non-Positive values (y <= 0) and NA are not allowed for this family!")
    yy <- y
  }
  if(!is.list(control)) stop("control is not a list")
  setting <- list(algorithm="ipc",method="dgLARS",nfold=10,foldid=NULL,ng=100,dispersion=NULL,g_hat=0,trace=FALSE,
                  max_var=min_np,np=NULL,g0=ifelse(p<(n-1),1.0e-04,1.0e-03),dg_max=0,nNR=50,NReps=1.0e-06,
                  ncrct=50,cf=0.5,eps=1.0e-05)
  nmsSetting <- names(setting)
  setting[(nms <- names(control))] <- control  
  if(length(noNms <- nms[!nms %in% nmsSetting]))
    warning("unknown names in control: ",paste(noNms,collapse =", "))
  if (family=="Gamma" & (link=="1/mu^2" | link=="logit" | link=="probit" | link=="cloglog")) {
    stop("This link function, ",link,", not available for 'Gamma' family; available links are: 'inverse', 'log' and 'identity'!")
  }
  if (family=="inverse.gaussian" & (link=="logit" | link=="probit" | link=="cloglog")) {
    stop("This link function, ",link,", not available for 'inverse.gaussian' family; available links are: '1/mu^2', 'inverse', 'log' and 'identity'!")
  }
  if (family=="binomial" & (link=="1/mu^2" | link=="identity" | link=="inverse")) {
    stop("This link function, ",link,", not available for 'binomial' family; available links are: 'logit', 'probit', 'cloglog' and 'log'!")
  }
  if(!setting$algorithm %in% c("pc","ipc"))
    stop("'algorithm' should be one of \"pc\" and \"ipc\"")
  if(!setting$method %in% c("dgLARS","dgLASSO"))
    stop("'method' should be \"dgLARS\" or \"dgLASSO\"")
  if(setting$nfold<1 | setting$nfold>n)
    stop("'nfold' should be an integer between 1 and n")
  if(is.null(setting$foldid)) {
    setting$foldid<-sample(n)
  } else{
    if(min(setting$foldid)<1 | max(setting$foldid)>n)
      stop("Read the documentation for 'foldid' more details.")
  }
  if(setting$ng<1)
    stop("'ng' should be a non-negative integer")
  if(setting$max_var<1 | setting$max_var>min_np) 
    stop("'max_var' should be an integer between 1 and min(n-1,p)")
  if(is.null(setting$np))
    setting$np <- ifelse(setting$algorithm %in% c("pc","ipc"),min_np*50,100L)
  if(setting$np<=0)
    stop("'np' should be a non-negative integer. Read the documentation for more details.")
  if(setting$g0<0)
    stop("'g0' should be a non-negative value. Read the documentation for more details.")
  if(setting$dg_max<0)
    stop("'dg_max' should be a non-negative value. Read the documentation for more details.")
  if(!is.null(setting$dispersion)) {
    if(setting$dispersion<=0)
      stop("'dispersion' should be a non-negative value. Read the documentation for more details.")
  }# else {setting$dispersion <- 1}  # Without Dispersion
  if(setting$eps<=0)
    stop("'eps' should be a non-negative value. Read the documentation for more details.")
  if(setting$ncrct<=0)
    stop("'ncrct' should be a non-negative value")
  if(setting$NReps<=0)
    stop("'NReps' should be a non-negative value")
  if(setting$nNR<=0)
    stop("'nNR' should be a non-negative integer")
  if(setting$cf<0 | setting$cf>1)
    stop("'cf' should be a value in the interval (0,1)")
  if(!is.logical(setting$trace))
    stop("'trace' should be one of \"TRUE\" and \"FALSE\"")
  fit <- cvdglars2.G.IG.B(X,yy,family,link,n,p,setting)
  fit <- make_cvdglars2(fit,setting)
  fit$call <- fit$dglars2class$call <- m.call
  fit$family <- fit$dglars2class$family <- family
  fit$link <- fit$dglars2class$link <- link
  fit$y <- y
  if(fit$conv!=0) warning("dgLARS with the ",fit$control$algorithm," algorithm does not converge, conv = ",fit$conv) # or warning
  fit
}
