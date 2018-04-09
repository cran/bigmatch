smahal<-function (z, p, X, caliper, constant=NULL, exact=NULL, nearexact=T, Xextra=NULL, weight=1){
  if (is.data.frame(X)) X<-as.matrix(X)
  if (is.vector(X)) X<-matrix(X,length(X),1)
  mode(X)<-'numeric'
  if (is.data.frame(Xextra)) Xextra<-as.matrix(Xextra)
  if (is.vector(Xextra)) Xextra<-matrix(Xextra,ncol=1)
  stopifnot(is.matrix(X))
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(caliper>0)
  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
    tb<-table(z,exact)
    if (!all(tb[2,]<=tb[1,])){
      warning("An exact match for exact is infeasible for every caliper.")
      if (!nearexact) return(NULL)
    }
  }
  if (!is.null(constant)){
    if (!is.integer(constant)) constant=floor(constant)
  }
  if (!is.null(Xextra)){
    weight=(dim(X)[2])/(dim(Xextra)[2])
  }

  n <- dim(X)[1]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    X<-X[o,]
    if (!is.null(exact)) exact<-exact[o]
    if (!is.null(Xextra)) Xextra<-Xextra[o,]
  }

  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)

  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- stats::cov(X)
  vuntied <- stats::var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  LL <- chol(cv)
  if (!is.null(Xextra)){
    for (j in 1:(dim(Xextra)[2])) Xextra[, j] <- rank(Xextra[, j])
    cv_ex <- stats::cov(Xextra)
    rat_ex <- sqrt(vuntied/diag(cv_ex))
    cv_ex <- diag(rat_ex) %*% cv_ex %*% diag(rat_ex)
    LL_ex <- chol(cv_ex)
  }

  X0 <- X[z == 0, ]
  X1 <- X[z == 1, ]
  p0 <- p[z == 0]
  p1 <- p[z == 1]
  if (!is.null(exact)){
    exact0 <- exact[z == 0]
    exact1 <- exact[z == 1]
  }
  if (!is.null(Xextra)){
    if (is.vector(Xextra)) Xextra<-matrix(Xextra,ncol=1)
    Xextra0 <- Xextra[z == 0,]
    Xextra1 <- Xextra[z == 1,]
    if (is.vector(Xextra0)) Xextra0<-matrix(Xextra0,ncol=1)
    if (is.vector(Xextra1)) Xextra1<-matrix(Xextra1,ncol=1)
  }

  distance <- c()
  start_node <- c()
  end_node <- c()
  nearex <- c()
  for (i in 1:m){
    #use caliper
    d<-abs(p1[i]-p0)
    who<-d<=caliper
    #use exact
    if ((!is.null(exact)) & (!nearexact)) who <- who & (exact1[i]==exact0)
    if (!any(who)) return(NULL)
    num <- sum(who)
    #use constant
    if (!is.null(constant)){
      if(num>constant){
        who<-which(who)[order(d[who])[1:constant]]
        num<- constant
      }
    }
    cc=X0[who,]
    if (num==1) cc<-matrix(cc,nrow=1)
    tt=t(as.matrix(X1[i,]))
    if (!is.null(Xextra)){
      cc_ex=Xextra0[who,]
      if (num==1) cc_ex<-matrix(cc_ex,nrow=1)
      else if (is.vector(cc_ex)) cc_ex=matrix(cc_ex,nrow=num)
      tt_ex=t(as.matrix(Xextra1[i,]))
      distancei=mvnfast::maha(cc,tt,LL,isChol = TRUE)+(mvnfast::maha(cc_ex,tt_ex,LL_ex,isChol = TRUE))*weight
    }else{
      distancei=mvnfast::maha(cc,tt,LL,isChol = TRUE)
    }

    #use nearexact
    if ((!is.null(exact)) & (nearexact)){
      nearex=c(nearex, exact1[i]!=exact0[who])
    }

    distance <- c(distance, distancei)
    start_node <- c(start_node, rep(i, num))
    end_node <- c(end_node,as.numeric(row.names(X0)[who]))
  }
  out <- list(d=distance,start=start_node,end=end_node,nearex=nearex)
  out
}
