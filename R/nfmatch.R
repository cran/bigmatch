nfmatch<-function(z,p,fine,X,caliper,dat,constant=NULL,exact=NULL,nearexact=T,rank=T,ncontrol=1,penalty=1000,max.cost=penalty/10,Xextra=NULL,weight=1,sub=F,subX=NULL){

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
    if (is.vector(exact)) exact=as.matrix(exact,ncol=1)
    EE=findexact(z,exact,ncontrol)
    Ex=EE$NewExact
    missingm=NULL
    if (is.null(Ex)){
      print("Exact matching for any subset of exact variables is infeasible for every caliper.")
      missingm=EE$miss
    }else if (dim(EE$miss)[2]>0){
      print("An exact match for all exact variables is infeasible for every caliper.")
      print("We will select as many important variables as we can.")
      missingm=EE$miss
    }
  }else{
    Ex=NULL
  }

  #sort input
  if (is.null(Ex)){
    o<-order(1-p)
  }else{
    o<-order(Ex,1-p)
    Ex<-Ex[o]
    if (!is.null(missingm)){
      missingm=missingm[o,]
    }
  }

  or <- rank(1-p,ties.method = 'min')
  z<-z[o]
  p<-p[o]
  or <- or[o]
  fine<-fine[o]
  X<-X[o,]
  dat<-dat[o,]
  if (!is.null(Xextra)) Xextra=Xextra[o,]
  if (!is.null(subX)) subX=subX[o]

  #do match
  timeind=proc.time()
  if (rank){
    dist<-smahal(z,or,X,caliper,constant,Ex,nearexact,Xextra, weight)
  }else{
    dist<-smahal(z,p,X,caliper,constant,Ex,nearexact,Xextra, weight)
  }
  timeind=proc.time()-timeind
  m<-nearfine(z,fine,dist,dat,X,ncontrol,penalty,max.cost,sub,subX)
  if(m[[1]]==0) {
    if (!is.null(exact)) warning("The match you requested is infeasible.  Reconsider caliper, ncontrol and exact.")
    else warning("The match you requested is infeasible, perhaps because the caliper is too small.")
  }else{
    list(data=m$d,timeinrelax=m$timeinrelax,edgenum=m$number,timeind=timeind,timeinnet=m$timeinnet,timeinmatch=m$timeinmatch)
  }
}
