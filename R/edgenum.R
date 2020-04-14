edgenum<-function(z,p,caliper,constant=NULL,exact=NULL,ties.all=TRUE){

  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(caliper>=0)
  if (!is.null(constant)){
    if (!is.integer(constant)) constant<-floor(constant)
  }

  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
  }

  n<-length(z)

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    p<-p[o]
    if (!is.null(exact)) exact<-exact[o]
  }

  ids<-1:n
  m<-sum(z)

  cids <-ids[z==0]
  p0<-p[z==0]
  p1<-p[z==1]
  if (!is.null(exact)){
    exact0<-exact[z==0]
    exact1<-exact[z==1]
  }


  result=0

  for (i in 1:m){
    #use caliper
    d<-abs(p1[i]-p0)
    who<-d<=caliper
    #use exact
    if (!is.null(exact)) who<-who&(exact1[i]==exact0)
    num<-sum(who)
    #use constant
    if (!is.null(constant)){
      if(num>constant){
        #        who<-which(who)[order(d[who])[1:constant]]
        #        num<- constant
        whoi<-which(who)[rank(d[who],ties.method='min')<=constant]
        if (ties.all || length(whoi)==constant){
          who<-whoi
          num<-length(who)
        }else{
          num<- constant
        }
      }
    }
    result=result+num
  }
  result
}
