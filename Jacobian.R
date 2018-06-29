Jacobian <- function(X,X2,ac,n,Score_s) {
  nac <- length(ac)
  Jacob <- matrix(0, nrow=nac, ncol=nac)
  w_for_dmnl <- Score_s$w.for.dmnl
  w_for_dmIn <- Score_s$w.for.dmIn
  Im_for_here <- Score_s$Imm[ac]
  rao_for_here <- Score_s$rao[ac]
  rao_over_Im <- rao_for_here/Im_for_here
  Jacob[1,1] <- apply(matrix(w_for_dmnl*X[,ac[1]]*X[,ac[1]],nrow=n),2,sum)
  for (h in 2:nac) {
    Jacob[1,h] <- apply(matrix(w_for_dmnl*X[,ac[1]]*X[,ac[h]],nrow=n),2,sum)
    Jacob[h,1] <- Jacob[1,h]/sqrt(Im_for_here[h])-0.5*rao_over_Im[h]*apply(matrix(w_for_dmIn*X[,ac[1]]*X2[,ac[h]],nrow=n),2,sum)
  }
  if (nac > 2) {
    for (k in 2:(nac-1)) {
      Jacob[k,k] <- apply(matrix(w_for_dmnl*X[,ac[k]]*X[,ac[k]],nrow=n),2,sum)/sqrt(Im_for_here[k])-0.5*rao_over_Im[k]*apply(matrix(w_for_dmIn*X[,ac[k]]*X2[,ac[k]],nrow=n),2,sum)
      for (h in (k+1):nac) {
        xkxh <- apply(matrix(w_for_dmnl*X[,ac[k]]*X[,ac[h]],nrow=n),2,sum)
        Jacob[k,h] <- xkxh/sqrt(Im_for_here[k])-0.5*rao_over_Im[k]*apply(matrix(w_for_dmIn*X[,ac[h]]*X2[,ac[k]],nrow=n),2,sum)
        Jacob[h,k] <- xkxh/sqrt(Im_for_here[h])-0.5*rao_over_Im[h]*apply(matrix(w_for_dmIn*X[,ac[k]]*X2[,ac[h]],nrow=n),2,sum)
      }
    }
  }
  Jacob[nac,nac] <- apply(matrix(w_for_dmnl*X[,ac[nac]]*X[,ac[nac]],nrow=n),2,sum)/sqrt(Im_for_here[nac])-0.5*rao_over_Im[nac]*apply(matrix(w_for_dmIn*X[,ac[nac]]*X2[,ac[nac]],nrow=n),2,sum)
  return(Jacob)
}
