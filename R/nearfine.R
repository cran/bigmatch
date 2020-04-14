nearfine<-function(z,fine,dist,dat,ncontrol=1,penalty=1000,max.cost=penalty/10,nearexPenalty=max.cost,subX=NULL){
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(length(z)==length(fine))
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  stopifnot(length(z)==(dim(dat)[1]))

  if (!requireNamespace("optmatch",quietly=TRUE)){
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
  }

  nobs<-length(z)
  #Must have treated first
  if(!(min(z[1:(nobs-1)]-z[2:nobs])>=0)){
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
    fine<-fine[o]
    if (!is.null(subX)) subX<-subX[o]
  }

  timeinnet<-proc.time()
  net<-netfine(z,fine,dist,ncontrol,penalty,max.cost,nearexPenalty,subX)
  timeinnet<-proc.time()-timeinnet
  timeinrelax<-proc.time()
  output<-rcbalance::callrelax(net)
  timeinrelax<-proc.time()-timeinrelax

  timeinmatch<-proc.time()
  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    list(feasible=output$feasible,timeinrelax=timeinrelax,d=NULL)
  }else{
    x<-output$x[1:net$tcarcs]
    treated<-net$startn[1:net$tcarcs]
    control<-net$endn[1:net$tcarcs]
    match.df<-data.frame('treat'=treated,'x'=x,'control'=control)
    treated=treated[which(x==1)]
    control=control[which(x==1)]
    match.df=data.frame('treat'=treated,'control'=control)
    match.df$treat<-as.factor(as.character(match.df$treat))
    matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                   function(treat.edges) treat.edges$control,.drop_o=FALSE))
    #matched.or.not<-plyr::daply(match.df,plyr::.(match.df$treat),
    #                              function(treat.edges) c(as.numeric(as.character(treat.edges$treat[1])),sum(treat.edges$x)),.drop_o=FALSE)
    #if(any(matched.or.not[,2]==0)){
    #  match.df<-match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,2]==0),1]),]
    #}
    #match.df$treat<-as.factor(as.character(match.df$treat))
    #matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
    #                                 function(treat.edges) treat.edges$control[treat.edges$x==1],.drop_o=FALSE))
    n<-length(z)
    ntreat<-sum(z)
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    matchid<-matrix(c(id1[as.numeric(row.names(matches))],id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
    matchid<-as.vector(t(matchid))
    dat1<-dat[matchid,]
    mset<-rep(1:length(matches),each=ncontrol+1)
    #dat1<-cbind(mset,dat1)
    dat1$mset<-mset
    timeinmatch<-proc.time()-timeinmatch
    list(feasible=output$feasible,timeinrelax=timeinrelax,timeinnet=timeinnet,timeinmatch=timeinmatch,d=dat1,number=net$tcarcs)
  }
}

