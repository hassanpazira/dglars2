cvdglars2 <- function(formula,family=c("Gamma","inverse.gaussian","binomial"),link=c("1/mu^2","inverse","log","identity","logit","probit","cloglog"),data,subset,contrast=NULL,control=list()){
  if (missing(data)) data <- environment(formula)
  m.call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(attr(mt,"intercept") == 0) stop("Models without intercept are not allowed in this version of the package")
  y <- model.response(mf, "any")
  X <- if (!is.empty.model(mt)) model.matrix(mt,mf,contrasts)
  else stop("Model matrix is empty")
  X <- X[,-1,drop=FALSE]
  if (!is.null(colnames(mf[,2]))) colnames(X) <- colnames(mf[,2])
  fit <- cvdglars2.fit(X=X,y=y,family=family,link=link,control=control)
  fit$call <- m.call 
  fit$predictors_cv <- names(fit$beta[abs(fit$beta) > 0])
  fit$formula_cv <- update(formula, as.formula(paste(" ~ ", paste(fit$var_cv,collapse = " + "))))
  fit
}
