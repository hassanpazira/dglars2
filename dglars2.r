dglars2 <- function(formula,family=c("Gamma","inverse.gaussian","binomial"),link=c("log","1/mu^2","inverse","identity","logit","probit","cloglog"),data,subset,contrast=NULL,control=list()){
  if (missing(data))  data <- environment(formula)
  m.call <- match.call()
  m.frame <- match.call(expand.dots = FALSE)
  # If set to FALSE lets all '...' arguments be collected as a single argument with the tag...'.
  m <- match(c("formula", "data", "subset"), names(m.frame), nomatch=0L) # 0L == as.integer(0)
  m.frame <- m.frame[c(1L, m)]
  m.frame$drop.unused.levels <- TRUE   
  # To add an argument one uses tagged list assignment.
  m.frame[[1L]] <- quote(model.frame)
  # To change the name of the function called, assign to the first element of the list and make sure that
  # the value is a name, either using the as.name("model.frame") construction here or quote(model.frame).
  m.frame <- eval(m.frame, envir=parent.frame())
  # If you want to re-run code captured with match.call(), you also need to capture the environment
  # in which it was evaluated, usually the parent.frame().
  m.terms <- attr(m.frame, "terms")
  if(attr(m.terms,"intercept") == 0) stop("Models without intercept are not allowed in this version of the package!")
  y <- model.extract(m.frame, "response")
  # model.extract(m.frame, "response") == model.response(m.frame, type="any")
  X <- if (!is.empty.model(m.terms)) model.matrix(m.terms,m.frame,contrasts)
  else stop("Design (or Model) matrix is empty")
  X <- X[,-1,drop=FALSE]            # Without Intercept
  if (!is.null(colnames(m.frame[,2]))) {
    colnames(X)<- colnames(m.frame[,2])
  }
  fit <- dglars2.fit(X=X,y=y,family=family,link=link,control=control)
  fit$call <- m.call
  fit
}
