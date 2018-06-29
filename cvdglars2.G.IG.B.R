cvdglars2.G.IG.B <- function(X,y,family,link,n,p,setting) {
  nfold <- setting$nfold
  foldid <- setting$foldid
  dispersion <- setting$dispersion  # can be omitted !
  ng <- setting$ng
  g0 <- setting$g0
  g <- ((ng-1):0)/(ng-1)
  lfold <- floor(n/nfold)
  b <- rep(0,p+1)
  dev_m <- dev_v <- double(ng)
  dev_cv <- c()
  for (i in 1:nfold) {
    cat(paste(" CV Fold: ",i,"\n"))
    setting$g_hat <- 2
    call_dgl <- dglars2.G.IG.B(X[foldid[(lfold+1):n],,drop=FALSE],y[foldid[(lfold+1):n]],family,link,n-lfold,p,setting) ## Based on 'training' datase '
    conv <- call_dgl$conv
    if (conv!=0) {
      warning("dgLARS in nfold= ",i," does not converge, conv = ",conv) ##
      break
    }
    b_cv <- call_dgl$b
    np_cv <- call_dgl$np
    g_seq <- call_dgl$g_seq
    #dispersion <- call_dgl$dispersion[np_cv]  # It's not requared anymore! # 'dispersion' can be omitted!
    call_pre_b <- predict_b(X[foldid[1:lfold],,drop=FALSE],y[foldid[1:lfold]],lfold,p,ng,g,g_seq,b_cv,dispersion,family,link) # Based on 'validation' dataset 
    dev_cv <- cbind(dev_cv,call_pre_b$dev_ng)
    foldid <- c(tail(foldid,n=-lfold),head(foldid,n=lfold))   # == cshift(foldid,lfold) in Fortran
  }
  if (conv==0) {
    dev_m <- apply(dev_cv,1,mean)
    dev_v <- apply(dev_cv^2,1,mean)-dev_m^2
    dev_v <- dev_v*nfold/(nfold-1)
    g_id <- which.min(dev_m)
    setting$g_hat <- g[g_id]
    call_2_dgl <- dglars2.G.IG.B(X,y,family,link,n,p,setting) 
    dglars2class <- make_dglars2(call_2_dgl,setting)
    dglars2class$family <- family
    dglars2class$link <- link
    conv <- call_2_dgl$conv
    if (conv==0) {
      b <- call_2_dgl$b[,call_2_dgl$np]
      g0 <- call_2_dgl$g0    #  g_hat <- g0
      g[1] <- call_2_dgl$g_seq[1]
      dispersion <- call_2_dgl$dispersion[call_2_dgl$np]
      if (is.null(setting$dispersion)) dispersion <- disper.hat(dglars2class)[dglars2class$np]
    }
  } else {
    stop("dgLARS in nfold= ",i," does not converge, conv = ",conv) ##
  }
  outputcv <- list(X=X, y=y, b=b, g=g, g0=g0, last_dispersion=dispersion, dev_m=dev_m, dev_v=dev_v, dglars2class=dglars2class, conv=conv)
  return(outputcv)
}
