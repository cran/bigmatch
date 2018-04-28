optconstant<-function(z,p,caliper=NULL,exact=NULL,ncontrol=1,tol=1,rank=T){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(length(z)==length(p))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(sum(z)*ncontrol<=sum(1-z))

  if (!is.null(exact)){
    tb<-table(z,exact)
    if (!all(tb[2,]<=tb[1,])){
      stop("An exact match for exact is infeasible for every caliper.")
    }
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
  nt<-sum(z)*ncontrol
  nc<-sum(1-z)
  cid<-1:nc

  leftrightexactc<-function(constant){
    if (!is.integer(constant)) constant=floor(constant)
    left<-rep(NA,nt)
    right<-left
    for (i in 1:nt){
      dn<-abs(use1[i]-use0)
      if (is.null(caliper)){
        close<-(ex1[i]==ex0)
      }else{
        close<-(ex1[i]==ex0)&(dn<=caliper)
      }
      if (sum(close)==0) return(NULL)
      if (sum(close)>constant){
        closei<-which(close)[order(dn[close])[1:constant]]
        if (max(closei)-min(closei)+1>constant){
          closea=closei[which(dn[closei]!=max(dn[closei]))]
          closeb=which(close)[which(dn[close]==max(dn[closei]))]
          closeb=closeb[order(pmin(abs(closeb-min(closea)),abs(closeb-max(closea))))]
          closeb=closeb[1:(constant-length(closea))]
          close=c(closea,closeb)
        }else{
          close=closei
        }
      }
      left[i]<-min(cid[close])
      right[i]<-max(cid[close])
    }
    list(left=left,right=right)
  }

  leftrightc<-function(constant){
    if (!is.integer(constant)) constant=floor(constant)
    left<-rep(NA,nt)
    right<-left
    for (i in 1:nt){
      dn<-abs(use1[i]-use0)
      if (is.null(caliper)){
        close=rep(T,nc)
      }else{
        close<-(dn<=caliper)
      }
      if (sum(close)==0) return(NULL)
      if (sum(close)>constant){
        closei<-which(close)[order(dn[close])[1:constant]]
        if (max(closei)-min(closei)+1>constant){
          closea=closei[which(dn[closei]!=max(dn[closei]))]
          closeb=which(close)[which(dn[close]==max(dn[closei]))]
          closeb=closeb[order(pmin(abs(closeb-min(closea)),abs(closeb-max(closea))))]
          closeb=closeb[1:(constant-length(closea))]
          close=c(closea,closeb)
        }else{
          close=closei
        }
      }
      left[i]<-min(cid[close])
      right[i]<-max(cid[close])
    }
    list(left=left,right=right)
  }

  highc<-nc
  lowc<-0

  while (((highc-lowc)>tol) && (highc>=1)){
    midc<-(highc+lowc)/2
    if (midc<1) highc=1
    else{
      if (!is.null(exact)) lr<-leftrightexactc(midc)
      else lr<-leftrightc(midc)
      if (!is.null(lr)){
        res<-glover(lr$left,lr$right)
      }
      if (is.null(lr)||is.null(res)){
        lowc<-midc
      }else highc<-midc
    }
  }

  list(constant=highc,interval=c(lowc,highc),interval.length=highc-lowc)
}
