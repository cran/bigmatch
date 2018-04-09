optcal<-function(z,p,exact=NULL,ncontrol=1,tol=NULL,rank=T){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(length(z)==length(p))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(sum(z)*ncontrol<=sum(1-z))

  if (is.null(tol)){
    if (rank) tol=1
    else tol=0.01
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

  or<-rank(1-p,ties.method = 'min')
  #sort input
  if (is.null(Ex)){
    o<-order(1-p)
  }else{
    o<-order(Ex,1-p)
    Ex<-Ex[o]
  }
  z<-z[o]
  p<-p[o]
  or<-or[o]

  # Define some variables and functions
  if (rank){
    use1<-rep(or[z==1],each=ncontrol)
    use0<-or[z==0]
  }else{
    use1<-rep(p[z==1],each=ncontrol)
    use0<-p[z==0]
  }
  if (!is.null(Ex)){
    ex1<-rep(Ex[z==1],each=ncontrol)
    ex0<-Ex[z==0]
  }
  nt<-sum(z)*ncontrol
  nc<-sum(1-z)
  cid<-1:nc

  leftrightexact<-function(caliper){
    left<-rep(NA,nt)
    right<-left
    for (i in 1:nt){
      close<-cid[(ex1[i]==ex0)&(abs(use1[i]-use0)<=caliper)]
      if (length(close)==0) return(NULL)
      left[i]<-min(close)
      right[i]<-max(close)
    }
    list(left=left,right=right)
  }

  leftright<-function(caliper){
    left<-rep(NA,nt)
    right<-left
    for (i in 1:nt){
      close<-cid[(abs(use1[i]-use0)<=caliper)]
      if (length(close)==0) return(NULL)
      left[i]<-min(close)
      right[i]<-max(close)
    }
    list(left=left,right=right)
  }

  if (rank){
    highc<-max(or)-1
  }else{
    highc<-max(p)-min(p)
  }
  lowc<-0

  while ((highc-lowc)>tol){
    midc<-(highc+lowc)/2
    if (!is.null(exact)) lr<-leftrightexact(midc)
    else lr<-leftright(midc)
    if (!is.null(lr)){
      res<-glover(lr$left,lr$right)}
    if (is.null(lr)||is.null(res)){
      lowc<-midc
    }
    else highc<-midc
  }
  list(caliper=highc,interval=c(lowc,highc),interval.length=highc-lowc)
}
