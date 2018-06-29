dglars2.G.IG.B <- function(X,y,family,link,n,p,setting) {
  if (all(as.numeric(X[,1])==rep(1,n))) warning("Look at 'X'; It should not be with 'Intercept' !")
  X <- cbind(rep(1,n),X)
  colnames(X)[1] <- "Int."
  p <- p+1
  X2 <- X*X
  dispersion <- setting$dispersion
  step.size.method <- setting$algorithm
  method <- setting$method
  trace <- setting$trace
  max_var <- setting$max_var
  dg_max <- setting$dg_max
  eps <- setting$eps
  g0 <- setting$g0
  g_hat <- setting$g_hat
  ncrct <- setting$ncrct
  nNR <- setting$nNR
  NReps <- setting$NReps
  cf <- setting$cf
  np <- setting$np
  dzero <- 1.0e-8
  zero <- 1.0e-5
  b <- rao <- vector(mode="numeric",length=p)
  nstp <- 0
  ac <- c(1)
  cac <- c(1:(p))[-ac]
  if(trace) cat(paste("Step: ",nstp,"\n"))
  if(trace) cat(paste("Length( A = {Int.} ): ",length(ac),"\n"))
  mu <- mean(y)                                      # E(y)=mu
  eta <- linkf(mu,family,link)                       # We can also use: b0<-as.numeric(glm(y~1,family,data=data.frame(wafer))$coef), because eta=family$linkfun(mu)=b0, because of (b1,...,bp)=(0,...,0).
  b0 <- eta
  b[ac] <- b0
  df <- c(c(),length(ac))
  if (b0==0) {
    df[1] <- 0
    #b[1] <- dzero
  }
  Matrix_B.Hat <- cbind(c(),b)
  nstp <- 1
  Score_s <- Scores(X,y,X2,b,ac,dispersion,n,family,link)
  rao <- Score_s$rao
  r_at_g <- cbind(c(),abs(rao[-1]))
  #dev <- Score_s$dev
  res.dev <- Score_s$res.dev
  disper <- Score_s$dispersion
  b_pre <- b
  ac_pre <- ac
  ac <- c(ac,as.numeric(which.max(abs(rao[-1])))+1)              # The active set with the intercept
  ac.order <- c(ac[-1])
  cac <- c(1:(p))[-ac]
  gmax <- max(abs(rao[-1]))                         # The first value of the tuning parameter
  g <- gmax                                         # Or  g <- abs(rao[ac.order[1]])
  g_seq <- g
  if(trace) cat(paste("Step: ",nstp,"\n"))
  if(trace) cat(paste("g: ",g,"\n"))
  if(trace) cat(paste("Length(A): ",length(ac),"\n"))
  conv <- 0
  n.rep <- 1
  dg_zero_break <- final <- The.end <- gleg0 <- length_eqcac <- FALSE  #  dg_zero
  if (g <= g0 | g_hat == 1) {
    gleg0 <- TRUE
    if (g_hat == 1) g0 <- g
  }
  if (g_hat != 2 & !gleg0) g0 <- g0 + g_hat*(g-g0)
  repeat{
    if (gleg0) break
    n.rep <- n.rep+1
    if (n.rep >= np) {
      conv <- 3
      # warning("Maximum number of iterations is been reached") #
      break
    }
    if(length(ac[-1]) >= max_var) final <- TRUE   # length(ac[-1]) == nav  and  max_var==nv
    Ja <- Jacobian(X,X2,ac,n,Score_s)
    if(abs(det(Ja))==0 | any(det(Ja)=="NaN")){
      conv <- 2
      # warning("No inverse for Jacobian in 'b' : 1") #
      break
    } else {
      dba <- try(solve(Ja,as.matrix(c(0,sign(as.numeric(Score_s$rao[ac[-1]]))))),TRUE)
      if(any(is.nan(dba)) | !is.null(attr(dba,"class"))) {
        conv <- 2
        break
      }
    }
    if (final) {
      if (dg_max > 0 ) {
        dg <- min(dg_max,g-g0)
      } else {
        dg <- g-g0
      }
      ai <- 0
    } else {
      dgopt <- dg.opt(X,y,X2,g,g0,ac,n,p,dg_max,dba,Score_s)
      dg <- dgopt$dg_opt
      ai <- dgopt$ai
    }
    if (method=="dgLASSO") {
      dg_out <- b[ac[-1]]/dba[-1]
      for (ii in 1:length(dg_out)) {
        if (dg_out[ii] > 0 & dg_out[ii] <= dg) {
          dg <- dg_out[ii]
          ai <- -ii
        }
      }
    }
    pre_g <- g
    g <- g-dg
    pre_b <- b
    for (n.rep.b in 1:ncrct) {
      b <- pre_b
      b[ac] <- b[ac]-dg*dba
      sign_rao_g <- c(0,sign(rao[ac[-1]])*g)
      Newt.list <- Newton.N.L(X,y,X2,b,ac,g,nNR,NReps,eps,dispersion,n,sign_rao_g,family,link) # dispersion should not be disper[length(disper)] ?
      if (Newt.list$conv!=0) {
        dg_old <- dg
        dg <- dg * cf
        if (dg <= dzero) {
          dg_zero_break <- TRUE
          break
        }
        g <- g+(dg_old-dg)
      } else {
        b[ac] <- Newt.list$bac
        Score_s <- Scores(X,y,X2,b,ac,dispersion,n,family,link)
        rao <- Score_s$rao
        # if (final & method=="dgLARS") break ###
        if (!final) {
          d_rao.c_g <- abs(rao[cac])-g
          len_wcac <- any(d_rao.c_g > eps)
          wcac <- cac[which(d_rao.c_g > eps)]
        } else {
          len_wcac <- FALSE
          wcac <- c()
        }
        if (method=="dgLASSO") {
          wac_b_r <- (b[ac]*rao[ac])[-1]
          wac <- any(wac_b_r<0)
        } else wac <- FALSE
        # if ((!length(wcac) & !wac)) {break}      ####  !length(wcac) = length(wcac)==0
        if (len_wcac | wac) {
          rao.base.a <- (rao)
          Score_s <- Scores(X,y,X2,b_pre,ac_pre,dispersion,n,family,link)
          rao <- Score_s$rao
          if (step.size.method=="pc" | (!length(wcac))) {       #  !length(wcac) = length(wcac)==0
            dg_old <- dg
            dg <- dg * cf
            g <- g+(dg_old-dg)
          } else {
            rao.base.b <- (rao)  # c(0,abs(r_at_g[,dim(r_at_g)[2]])) # c(0,abs(r_at_g[,nstp]))
            rao.a <- rao.base.a[wcac]
            rao.b <- rao.base.b[wcac]
            g.a <- g
            g.b <- pre_g
            go.back <- as.numeric(((g.a*rao.b)-(g.b*rao.a))/((rao.b-rao.a)+(g.a-g.b)*(sign(rao.a))))
            go.back_max <- max(go.back)
            if (g.b-go.back_max < 0) go.back_max<-g.b
            dg <- as.numeric(g.b-go.back_max)
            g <- go.back_max
          }
          if (dg <= dzero) {
            # dg_zero <- TRUE
            g <- pre_g
            break
          }
        } else {break}
      }
    }
    if (dg_zero_break) {
      conv <- 2
      #warning("dg <= dzero ; dg_zero_break") #
      break
    }
    if (n.rep.b >= ncrct) {
      conv <- 3
      #warning("Maximum number of iterations, in 'ncrct'")
      break
    }
    if (abs(abs(rao[ac.order[1]])-g0) <= eps) The.end <- TRUE
    if (len_wcac | wac) { # dg_zero
      # dg_zero <- FALSE
      if (len_wcac) {  # !length(wcac)  # !?
        cac_fit <- which(d_rao.c_g >= eps)
        for (aci in cac_fit) {
          if (d_rao.c_g[aci] >= eps) { # just for sure  # it can be omit
            new_ac <- cac[aci]
            #new_ac <- cac[which(d_rao.c_g==d_rao.c_g[which(d_rao.c_g >= eps)][which.min(d_rao.c_g[which(d_rao.c_g >= eps)])])] #
            ac <- sort(c(ac,as.numeric(new_ac)))                # The updated Active set
            ac.order<-c(ac.order,as.numeric(new_ac))
            if(length(ac[-1]) >= max_var) {
              final <- TRUE
              break
            }
          }
        }
        cac <- c(1:p)[-ac]
      }
      if (wac) {
        ai_w <- which(wac_b_r<0)
        for (wac_i in ai_w) {
          ac.order[which(ac.order==ac[wac_i+1])] <- ac.order[length(ac.order)]
          ac.order <- ac.order[-length(ac.order)]
          ac <- sort(c(1,ac.order))
          cac <- c(1:p)[-ac]
          final <- FALSE
        }
      }
    } else {
      nstp <- nstp+1
      b[cac] <- 0
      b_pre <- b
      ac_pre <- ac
      Score_s <- Scores(X,y,X2,b,ac,dispersion,n,family,link)
      rao <- Score_s$rao
      g <- abs(rao[ac.order[1]])
      g_seq <- c(g_seq,g)
      r_at_g <- cbind(r_at_g,abs(rao[-1]))
      Matrix_B.Hat <- cbind(Matrix_B.Hat,b)
      #dev <- c(dev,Score_s$dev)
      res.dev <- c(res.dev,Score_s$res.dev)
      df <- c(df,length(ac))
      disper <- c(disper,Score_s$dispersion)
      if (any(abs(b[ac[-1]])<=zero) & method=="dgLASSO" & ai<0) {
        ai <- abs(ai)
        rao[ac[ai+1]] <- g
        old_ac <- ac
        ac.order[which(ac.order==ac[ai+1])] <- ac.order[length(ac.order)]
        ac.order <- ac.order[-length(ac.order)]
        ac <- sort(c(1,ac.order))
        cac <- c(1:p)[-ac]
        if (abs(b[old_ac[ai+1]])==0) df[length(df)] <- length(ac)
        final <- FALSE
      }
      if (!final & ai>0) {
        eqcac <- cac[which(abs(d_rao.c_g) <= eps)]
        if (!length(eqcac)) {
          length_eqcac<-FALSE
        } else {
          length_eqcac<-TRUE
        }
        if(length_eqcac) {
          if (length(eqcac)==1) {
            new_eqcac <- eqcac
          } else {
            #warning("We have more than one covariate which is in (-eps,eps) !")
            new_eqcac <- cac[which(d_rao.c_g==d_rao.c_g[which(abs(d_rao.c_g) <= eps)][which.max((d_rao.c_g[which(abs(d_rao.c_g) <= eps)]))])] # NOT cac[...][which.max(abs(d_rao.c_g[which(abs(d_rao.c_g) <= eps)]))])]
          }
          ac <- sort(c(ac,as.numeric(new_eqcac)))                # The updated Active set
          ac.order<-c(ac.order,as.numeric(new_eqcac))
          cac <- c(1:(p))[-ac]                                   # The complementary of the Active set
          if(trace) cat(paste("Step: ",nstp,"\n"))
          if(trace) cat(paste("g: ",g,"\n"))
          if(trace) cat(paste("Length(A): ",length(ac),"\n"))
        }
      }
    }
    b[cac] <- 0
    if (The.end) break
  }
  X <- X[,-1,drop=FALSE]
  colnames(Matrix_B.Hat) <- NULL
  #rownames(Matrix_B.Hat) <- colnames(X)
  np <- length(g_seq)  # or np <- nstp   because np==nstp
  order <- ac.order-1
  output <- list(X=X, y=y, b=Matrix_B.Hat, rao=r_at_g, df=df, np=np, g_seq=g_seq, g0=g0, res.dev=res.dev, dispersion=disper, A_1.nav=order, conv=conv)
  return(output)
}
