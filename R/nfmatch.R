nfmatch<-function(z,p,fine=rep(1,length(z)),X,dat,caliper,constant=NULL,ncontrol=1,rank=TRUE,exact=NULL,penalty=1000,max.cost=penalty/10,nearexact=NULL,nearexPenalty=max.cost,Xextra=NULL,weight=NULL,subX=NULL,ties.all=TRUE,seed=1){

  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  if (is.factor(fine)){
    levels(fine)<-1:nlevels(fine)
    fine<-as.integer(fine)
  }
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  stopifnot(length(z)==length(p))
  stopifnot(length(z)==length(fine))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  if (is.null(subX)) stopifnot(ncontr>=(ncontrol*ntreat))

  if (is.vector(X)) X<-matrix(X,nrow=length(z))
  if (is.data.frame(nearexact)) nearexact<-as.matrix(nearexact)
  if (is.factor(nearexact)){
    levels(nearexact)<-1:nlevels(nearexact)
    nearexact<-as.numeric(nearexact)
  }
  if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)

  if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(length(z)==(dim(dat)[1]))
  stopifnot(caliper>=0)
  if (!is.null(constant)) stopifnot(constant>=1)
  if (!is.null(subX)){
    if (is.factor(subX)){
      levels(subX)<-1:nlevels(subX)
      subX<-as.integer(subX)
    }
    stopifnot(is.vector(subX))
  }

  if (!is.null(exact)){
    if (is.factor(exact)){
      levels(exact)<-1:nlevels(exact)
      exact<-as.integer(exact)
    }
  }

  if (!ties.all){
    set.seed(seed)
    ra<-sample(1:nobs,nobs)
    z<-z[ra]
    p<-p[ra]
    exact<-exact[ra]
    fine<-fine[ra]
    X<-X[ra,]
    dat<-dat[ra,]
    if (!is.null(nearexact)) nearexact<-nearexact[ra,]
    if (!is.null(Xextra)) Xextra<-Xextra[ra,]
    if (!is.null(subX)) subX<-subX[ra]
    if (is.vector(X)) X<-matrix(X,length(X),1)
    if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
    if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)
  }

  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  or<-rank(p,ties.method='min')
  z<-z[o]
  p<-p[o]
  or<-or[o]
  fine<-fine[o]
  X<-X[o,]
  dat<-dat[o,]
  if (!is.null(nearexact)) nearexact<-nearexact[o,]
  if (!is.null(Xextra)) Xextra<-Xextra[o,]
  if (!is.null(subX)) subX<-subX[o]
  if (is.vector(X)) X<-matrix(X,length(X),1)
  if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
  if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)

  #Must have treated first
  if(!(min(z[1:(nobs-1)]-z[2:nobs])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    or<-or[o]
    X<-X[o,]
    dat<-dat[o,]
    fine<-fine[o]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(nearexact)) nearexact<-nearexact[o,]
    if (!is.null(Xextra)) Xextra<-Xextra[o,]
    if (!is.null(subX)) subX<-subX[o]
    if (is.vector(X)) X<-matrix(X,length(X),1)
    if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
    if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)
  }

  #do match
  timeind<-proc.time()
  if (rank){
    dist<-smahal(z,or,X,caliper,constant,ncontrol,exact,nearexact,nearexPenalty,Xextra,weight,subX,ties.all)
  }else{
    dist<-smahal(z,p,X,caliper,constant,ncontrol,exact,nearexact,nearexPenalty,Xextra,weight,subX,ties.all)
  }
  timeind<-proc.time()-timeind
  m<-nearfine(z,fine,dist,dat,ncontrol,penalty,max.cost,nearexPenalty,subX)
  if(m[[1]]==0) {
    if (!is.null(exact)) warning("The match you requested is infeasible.  Reconsider caliper, ncontrol and exact.")
    else warning("The match you requested is infeasible, perhaps because the caliper is too small.")
  }else{
    list(data=m$d,timeinrelax=m$timeinrelax,edgenum=m$number,timeind=timeind,timeinnet=m$timeinnet,timeinmatch=m$timeinmatch)
  }
}
