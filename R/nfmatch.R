nfmatch<-function(z,p,fine,X,dat,caliper,constant=NULL,ncontrol=1,rank=T,exact=NULL,penalty=1000,max.cost=penalty/10,nearexact=NULL,nearexPenalty=max.cost,Xextra=NULL,weight=NULL,sub=F,subX=NULL){

  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  if (is.factor(fine)){
    levels(fine)=1:nlevels(fine)
    fine<-as.integer(fine)
  }
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(length(unique(fine))>=2)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  stopifnot(length(z)==length(p))
  stopifnot(length(z)==length(fine))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  if (is.data.frame(X)) X<-as.matrix(X)
  if (is.vector(X)) X<-matrix(X,length(X),1)
  if (is.data.frame(nearexact)) nearexact<-as.matrix(nearexact)
  if (is.factor(nearexact)){
    levels(nearexact)=1:nlevels(nearexact)
    nearexact=as.numeric(nearexact)
  }
  if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
  if (is.data.frame(Xextra)) Xextra<-as.matrix(Xextra)
  if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(length(z)==(dim(dat)[1]))
  stopifnot(caliper>=0)
  if (!is.null(constant)) stopifnot(constant>=0)
  if (!sub) stopifnot(is.null(subX))
  if (sub) stopifnot(!is.null(subX))
  if (!is.null(subX)){
    if (is.factor(subX)){
      levels(subX)=1:nlevels(subX)
      subX<-as.integer(subX)
    }
    stopifnot(is.vector(subX))
  }

  if (!is.null(exact)){
    if (is.factor(exact)){
      levels(exact)=1:nlevels(exact)
      exact<-as.integer(exact)
    }
  }

  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  or <- rank(1-p,ties.method = 'min')
  z<-z[o]
  p<-p[o]
  or <- or[o]
  fine<-fine[o]
  X<-X[o,]
  dat<-dat[o,]
  if (!is.null(nearexact)) nearexact=nearexact[o,]
  if (!is.null(Xextra)) Xextra=Xextra[o,]
  if (!is.null(subX)) subX=subX[o]

  #do match
  timeind=proc.time()
  if (rank){
    dist<-smahal(z,or,X,caliper,constant,exact,nearexact,Xextra,weight)
  }else{
    dist<-smahal(z,p,X,caliper,constant,exact,nearexact,Xextra,weight)
  }
  timeind=proc.time()-timeind
  m<-nearfine(z,fine,dist,dat,X,ncontrol,penalty,max.cost,nearexPenalty,sub,subX)
  if(m[[1]]==0) {
    if (!is.null(exact)) warning("The match you requested is infeasible.  Reconsider caliper, ncontrol and exact.")
    else warning("The match you requested is infeasible, perhaps because the caliper is too small.")
  }else{
    list(data=m$d,timeinrelax=m$timeinrelax,edgenum=m$number,timeind=timeind,timeinnet=m$timeinnet,timeinmatch=m$timeinmatch)
  }
}
