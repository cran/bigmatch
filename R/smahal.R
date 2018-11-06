smahal<-function(z,p,X,caliper,constant=NULL,ncontrol=1,exact=NULL,nearexact=NULL,Xextra=NULL,weight=NULL,subX=NULL,ties.all=T){

  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x)) x<-matrix(x,nrow=length(z))

    if(is.data.frame(x) || is.character(x)){
      if(!is.data.frame(x)) x <- as.data.frame(x)
      X.chars <- which(plyr::laply(x, function(y) 'character' %in% class(y)))
      if(length(X.chars) > 0){
        for(i in X.chars){
          x[,i] <- factor(x[,i])

        }
      }
      #if some variables are factors convert to dummies
      X.factors <-  which(plyr::laply(x, function(y) 'factor' %in% class(y)))

      #handle missing data
      for(i in which(plyr::laply(x, function(y) any(is.na(y))))){
        if(i %in% X.factors){
          #for factors, make NA a new factor level
          x[,i] <- addNA(x[,i])
        }else{
          #for numeric/logical, impute means and add a new indicator for missingness
          x[[paste(colnames(x)[i],'NA', sep = '')]] <- is.na(x[,i])
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }
      for(i in rev(X.factors)){
        dummyXi <- model.matrix(as.formula(
          paste('~',colnames(x)[i], '-1')),data=x)
        x <- cbind(x[,-i], dummyXi)
      }

    }else{
      #handle missing data
      for(i in c(1:ncol(x))){
        if(any(is.na(x[,i]))){
          x <- cbind(x,is.na(X[,i]))
          colnames(x)[ncol(x)] <- paste(colnames(X)[i],'NA', sep = '')
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }

    }

    #get rid of columns that do not vary
    varying <- apply(x,2, function(y) length(unique(y)) > 1)
    x <- x[,which(varying),drop = FALSE]

    as.matrix(x)
  }

  X<-Xmatrix(X)
  if (is.data.frame(nearexact)) nearexact<-data.matrix(nearexact)
  if (is.vector(nearexact)) nearexact<-as.matrix(nearexact,ncol=1)
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(caliper>=0)
  if (!is.null(constant)){
    if (!is.integer(constant)) constant<-floor(constant)
  }

  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
    tb<-table(z,exact)
    if (!all(tb[2,]<=tb[1,])){
      if (is.null(subX) || any(exact!=subX)){
        stop("An exact match for exact is infeasible for every caliper.")
      }
    }
  }

  if (!is.null(Xextra)) Xextra<-Xmatrix(Xextra)
  if (!is.null(Xextra) & is.null(weight)) weight<-(dim(matrix(X,nrow=length(z)))[2])/(dim(matrix(Xextra,nrow=length(z)))[2])

  n<-dim(X)[1]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    X<-X[o,]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(nearexact)) nearexact<-nearexact[o,]
    if (!is.null(Xextra)) Xextra<-Xextra[o,]
  }

  if (is.vector(X)) X<-matrix(X,ncol=1)
  if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
  if (is.vector(Xextra)) Xextra<-matrix(Xextra,length(Xextra),1)

  ids<-1:n
  k<-dim(X)[2]
  m<-sum(z)

  for (j in 1:k) X[, j]<-rank(X[, j])
  cv<-stats::cov(X)
  vuntied<-stats::var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  LL<-try(chol(cv), silent=T)
  if(is(LL,"try-error")) {
    LL<-chol(cv+10^{-10}*diag(nrow(cv)))
  }

  if (!is.null(Xextra)){
    if (is.vector(Xextra)) Xextra<-matrix(Xextra,ncol=1)
    for (j in 1:(dim(Xextra)[2])) Xextra[, j]<-rank(Xextra[, j])
    cv_ex<-stats::cov(Xextra)
    rat_ex<-sqrt(vuntied/diag(cv_ex))
    cv_ex<-diag(rat_ex)%*%cv_ex%*%diag(rat_ex)
    LL_ex<-try(chol(cv_ex), silent=T)
    if(is(LL_ex,"try-error")) {
      LL_ex<-chol(cv_ex+10^{-10}*diag(nrow(cv_ex)))
    }
  }

  if (is.vector(X)) X<-matrix(X,ncol=1)
  X0<-X[z==0,]
  X1<-X[z==1,]
  if (is.vector(X0)) X0<-matrix(X0,ncol=1)
  if (is.vector(X1)) X1<-matrix(X1,ncol=1)
  cids <-ids[z==0]
  p0<-p[z==0]
  p1<-p[z==1]
  if (!is.null(exact)){
    exact0<-exact[z==0]
    exact1<-exact[z==1]
  }
  if (!is.null(nearexact)){
    if (is.vector(nearexact)) nearexact<-matrix(nearexact,ncol=1)
    nearex0<-nearexact[z==0,]
    nearex1<-nearexact[z==1,]
    if (is.vector(nearex0)) nearex0<-matrix(nearex0,ncol=1)
    if (is.vector(nearex1)) nearex1<-matrix(nearex1,ncol=1)
  }
  if (!is.null(Xextra)){
    if (is.vector(Xextra)) Xextra<-matrix(Xextra,ncol=1)
    Xextra0<-Xextra[z==0,]
    Xextra1<-Xextra[z==1,]
    if (is.vector(Xextra0)) Xextra0<-matrix(Xextra0,ncol=1)
    if (is.vector(Xextra1)) Xextra1<-matrix(Xextra1,ncol=1)
  }

  distance<-c()
  start_node<-c()
  end_node<-c()
  nearex<-c()


  for (i in 1:m){
    #use caliper
    d<-abs(p1[i]-p0)
    who<-d<=caliper
    #use exact
    if (!is.null(exact)) who<-who&(exact1[i]==exact0)
    if (!any(who) & is.null(subX)) return(NULL)
    num<-sum(who)
    #use constant
    nexti<-ncontrol
    if (!is.null(constant)){
      if(num>constant){
        #        who<-which(who)[order(d[who])[1:constant]]
        #        num<- constant
        whoi<-which(who)[rank(d[who],ties.method='min')<=constant]
        if (ties.all || length(whoi)==constant){
          who<-whoi
          num<-length(who)
        }else if (constant==ncontrol){
          whoi<-which(who)[order(d[who])[1:ncontrol]]
          if ((cids[whoi][1] %in% left) && (cids[whoi][constant] %in% right)){
            nexti<-nexti+1
            if (nexti+constant-1>n-m) return(NULL)
            else{
              nextclose<-all(who[order(d[who])[1:constant]]==who[order(d[who])[nexti:(nexti+constant-1)]])
              whoi<-which(who)[order(d[who])[nexti:(nexti+constant-1)]]
              while(nextclose && (cids[whoi][1] %in% left) && (cids[whoi][constant] %in% right) && (nexti+constant-1<=n-m)){
                nexti<-nexti+1
                if (nexti+constant-1>n-m) return(NULL)
                whoi<-which(who)[order(d[who])[nexti:(nexti+constant-1)]]
              }
            }
          }
          who<-whoi
          num<- constant
          left[i]<-min(cids[who])
          right[i]<-max(cids[who])
        }else if(any(d[whoi]!=max(d[whoi]))){
          closea<-whoi[which(d[whoi]!=max(d[whoi]))]
          closeb<-which(who)[which(d[who]==max(d[whoi]))]
          closeb<-closeb[order(pmin(abs(closeb-min(closea)),abs(closeb-max(closea))))]
          closeb<-closeb[1:(constant-length(closea))]
          who<-c(closea,closeb)
          num<- constant
        }else{
          who<-which(who)[1:constant]
          num<- constant
        }
      }
    }
    cc<-X0[who,]
    if (num==1) cc<-matrix(cc,nrow=1)
    else if (is.vector(cc)) cc<-matrix(cc,nrow=num)
    tt<-t(as.matrix(X1[i,]))
    if (!is.null(Xextra)){
      cc_ex<-Xextra0[who,]
      if (num==1) cc_ex<-matrix(cc_ex,nrow=1)
      else if (is.vector(cc_ex)) cc_ex<-matrix(cc_ex,nrow=num)
      tt_ex<-t(as.matrix(Xextra1[i,]))
      distancei<-mvnfast::maha(cc,tt,LL,isChol=TRUE)+(mvnfast::maha(cc_ex,tt_ex,LL_ex,isChol=TRUE))*weight
    }else{
      distancei<-mvnfast::maha(cc,tt,LL,isChol=TRUE)
    }

    distance<-c(distance, distancei)
    start_node<-c(start_node, rep(i, num))
    end_node<-c(end_node,cids[who])

    #use nearexact
    if (!is.null(nearexact)){
      if (dim(nearexact)[2]==1) nearex<-c(nearex, nearex1[i,]!=nearex0[who,])
      else nearex<-rbind(nearex, nearex1[i,]!=nearex0[who,])
    }
  }

  out<-list(d=distance,start=start_node,end=end_node,nearex=nearex)
  out
}
