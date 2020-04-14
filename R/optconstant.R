optconstant<-function(z,p,caliper=NULL,exact=NULL,ncontrol=1,tol=1,rank=TRUE,subX=NULL,ties.all=TRUE,seed=1){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(length(z)==length(p))
  stopifnot(all((z==1)|(z==0)))
  if (is.null(subX)) stopifnot(sum(z)*ncontrol<=sum(1-z))

  if (!is.null(subX)){
    subX<-as.factor(subX)
    nsubXlevels<-nlevels(subX)
    levels(subX)<-1:nsubXlevels
  }

  TF<-F
  if (!is.null(exact)){
    exact<-as.factor(exact)
    nexactlevels<-nlevels(exact)
    levels(exact)<-1:nexactlevels
    tb<-table(z,exact)

    if (!is.null(subX)) TF<-(nlevels(exact)!=nlevels(subX)) ||any(exact!=subX)
    if (!all(tb[2 ,]<=tb[1,])){
      if (is.null(subX) || TF){
        stop("An exact match for exact is infeasible for every caliper.")
      }
    }
    ratio<-tb[1,]/tb[2,]
    order_ratio<-order(ratio)
    exactlevels<-(1:nexactlevels)[order_ratio]
  }

  if (!ties.all){
    set.seed(seed)
    ra<-sample(1:length(z),length(z))
    z<-z[ra]
    p<-p[ra]
    exact<-exact[ra]
  }
  or<-rank(p,ties.method='min')
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

  TF<-TF||is.null(exact)||is.null(subX)
  optconstantone<-function(u1,u0,caliper=NULL,ub,lb,tol){
    nt<-length(u1)
    nc<-length(u0)
    cid<-1:nc

    leftrightc<-function(constant){
      if (!is.integer(constant)) constant<-floor(constant)
      left<-rep(NA,nt)
      right<-left
      for (i in 1:nt){
        dn<-abs(u1[i]-u0)
        if (is.null(caliper)){
          close<-rep(T,nc)
        }else{
          close<-(dn<=caliper)
        }
        nexti<-ncontrol
        #if ((sum(close)==0) & TF) return(NULL)
        if (sum(close)==0) return(NULL)
        if (sum(close)>constant){
          closei<-which(close)[rank(dn[close],ties.method='min')<=constant]
          if (ties.all || length(closei)==constant) close<-closei
          else if (constant<=ncontrol){
            closei<-which(close)[order(dn[close])[1:ncontrol]]
            if ((cid[closei][1] %in% left) && (cid[closei][constant] %in% right)){
              nexti<-nexti+1
              if (nexti+constant-1>nc) return(NULL)
              else{
                nextclose<-all(close[order(dn[close])[1:constant]]==close[order(dn[close])[nexti:(nexti+constant-1)]])
                closei<-which(close)[order(dn[close])[nexti:(nexti+constant-1)]]
                while(nextclose && (cid[closei][1] %in% left) && (cid[closei][constant] %in% right) && (nexti+constant-1<=nc)){
                  nexti<-nexti+1
                  if (nexti+constant-1>nc) return(NULL)
                  closei<-which(close)[order(dn[close])[nexti:(nexti+constant-1)]]
                }
              }
            }
            close<-closei
          }else{
            if(any(dn[closei]!=max(dn[closei]))){
              closea<-closei[which(dn[closei]!=max(dn[closei]))]
              closeb<-which(close)[which(dn[close]==max(dn[closei]))]
              closeb<-closeb[order(pmin(abs(closeb-min(closea)),abs(closeb-max(closea))))]
              closeb<-closeb[1:(constant-length(closea))]
              close<-c(closea,closeb)
            }else{
              close<-which(close)[1:constant]
            }
          }
        }
        if (sum(close)>0){
          left[i]<-min(cid[close])
          right[i]<-max(cid[close])
        }
      }

      list(left=left,right=right)
    }

    highc<-floor(ub)
    lowc<-lb

    lr<-leftrightc(lowc)
    res<-NULL
    if (!is.null(lr)&& any(!is.na(lr$left))){
      res<-glover(lr$left[!is.na(lr$left)],lr$right[!is.na(lr$right)])
    }
    if ((!is.null(lr)) && any(!is.na(lr$left)) && (res>=min(nt,nc)/nt)){
      return(list(constant=floor(lowc),interval=c(floor(lowc),floor(lowc)),interval.length=0))
    }

    lr<-leftrightc(highc)
    if ((!is.null(lr)) && any(!is.na(lr$left))){
      res<-glover(lr$left[!is.na(lr$left)],lr$right[!is.na(lr$right)])
    }
    if (is.null(lr)||((!is.null(res))&&(ceiling(res*sum(!is.na(lr$left)))<min(nt,nc)))){
      stop('The caliper itself is not feasible.')
    }

    while (((highc-lowc)>tol) && (highc>=1)){
      midc<-(highc+lowc)/2
      if (midc<1) highc<-1
      else{
        lr<-leftrightc(midc)
        if ((!is.null(lr)) && any(!is.na(lr$left))){
          res<-glover(lr$left[!is.na(lr$left)],lr$right[!is.na(lr$right)])
        }
        if (is.null(lr)||((!is.null(res))&&(round(res*sum(!is.na(lr$left)))<min(nt,nc)))){
          lowc<-midc
        }else highc<-floor(midc)
      }
    }

    list(constant=highc,interval=c(lowc,highc),interval.length=highc-lowc)
  }

  if (is.null(exact)){
    use<-c(use1,use0)
    ub<-length(use0)
    lb<-ncontrol
    optconstantone(use1,use0,caliper,ub,lb,tol)
  }else{
    cons<-0
    for (ee in 1:nexactlevels){
      ei<-exactlevels[ee]
      use1e<-use1[ex1==ei]
      use0e<-use0[ex0==ei]
      use<-c(use1e,use0e)
      if ((length(use1e)>0) & (length(use0e)>0)){
        ub<-length(use0e)
        lb<-max(ncontrol,cons-tol)
        if (ub>=lb){
          ree<-optconstantone(use1e,use0e,caliper,ub,lb,tol)
          if (cons<ree$constant){
            result<-ree
            cons<-ree$constant
          }
        }
      }
    }
    result
  }
}
