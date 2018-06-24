optcal<-function(z,p,exact=NULL,ncontrol=1,tol=NULL,rank=T){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(length(z)==length(p))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(sum(z)*ncontrol<=sum(1-z))

  if (is.null(tol)){
    if (rank) tol<-1
    else tol<-0.01
  }

  if (!is.null(exact)){
    exact<-as.factor(exact)
    nexactlevels<-nlevels(exact)
    levels(exact)<-1:nexactlevels
    tb<-table(z,exact)
    if (!all(tb[2,]<=tb[1,])){
      stop("An exact match for exact is infeasible for every caliper.")
    }
    ratio<-tb[1,]/tb[2,]
    order_ratio<-order(ratio)
    exactlevels<-(1:nexactlevels)[order_ratio]
  }

  or<-rank(1-p,ties.method='min')
  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
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
  if (!is.null(exact)){
    ex1<-rep(exact[z==1],each=ncontrol)
    ex0<-exact[z==0]
  }

  optcalone<-function(u1,u0,ub,lb,tol){
    nt<-length(u1)
    nc<-length(u0)
    cid<-1:nc

    leftright<-function(caliper){
      left<-rep(NA,nt)
      right<-left
      for (i in 1:nt){
        close<-cid[(abs(u1[i]-u0)<=caliper)]
        if (length(close)==0) return(NULL)
        left[i]<-min(close)
        right[i]<-max(close)
      }
      list(left=left,right=right)
    }

    highc<-ub
    lowc<-lb

    while ((highc-lowc)>tol){
      midc<-(highc+lowc)/2
      lr<-leftright(midc)
      if (!is.null(lr)){
        res<-glover(lr$left,lr$right)}
      if (is.null(lr)||is.null(res)){
        lowc<-midc
      }
      else highc<-midc
    }
    list(caliper=highc,interval=c(lowc,highc),interval.length=highc-lowc)
  }

  if (is.null(exact)){
    use<-c(use1,use0)
    ub<-max(use)-min(use)
    lb<-0
    optcalone(use1,use0,ub,lb,tol)
  }else{
    cals<-0
    for (ee in 1:nexactlevels){
      ei<-exactlevels[ee]
      use1e<-use1[ex1==ei]
      use0e<-use0[ex0==ei]
      use<-c(use1e,use0e)
      ub<-max(use)-min(use)
      lb<-max(cals-tol,0)
      if (ub>=lb){
        ree<-optcalone(use1e,use0e,ub,lb,tol)
        if (cals<ree$caliper){
          result<-ree
          cals<-ree$caliper
        }
      }
    }
    result
  }
}
