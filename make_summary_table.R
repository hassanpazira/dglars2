make_summary_table <- function(object,k,complexity){
  n <- dim(object$X)[1]
  np <- object$np
  action <- object$action
  b <- object$beta 
  dispersion <- disper.hat(object)
  dev <- deviance(object)
  res.dev <- object$res.dev
  g <- object$g
  if(complexity=="df")  compl <- object$df
  else compl <- gdf(object)
  if (object$family=="binomial" | !is.null(object$control$dispersion)) gof <- dev + k * compl
  else gof <- dev + k * (compl+1)
  rank.gof <- rank(gof)
  best <- rank.gof==1
  b.gof <- b[,best]
  b.gof <- b.gof[abs(b.gof)>0]
  mark <- rep("   ",np)
  mark[best] <- "<--"
  rank.gof <- paste(rank.gof,mark)
  tbl <- data.frame(Sequence=action,g=g,Disp=dispersion,Res.Dev=res.dev,Complexity=compl,gof=gof,Rank=rank.gof)
  names(tbl)[names(tbl)=="Rank"]<-"Rank    "
  list(table=tbl,b.gof=b.gof)
}
