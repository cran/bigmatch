glover<-function(left,right,mm=F){
  #This is Glover's algorithm as given by Katriel (2008) Matchings in node-weighted convex bipartitie graphs.
  #INFORMS Journal on Computing 20(2):205-211, Figure 1

  nx<-length(left)
  ny<-max(right)
  stopifnot(length(left)==length(right))
  stopifnot(all(left<=right))

  mx<-rep(NA,nx)
  my<-rep(NA,nx)
  qu <- liqueueR::PriorityQueue$new()
  for (i in 1:ny) {
    w<-which(left==i)
    if (any(w)) {
      for (j in w) qu$push(j,priority=ny-right[j])
    }
    if (qu$size()>0){
      xj<-qu$pop()
      if (right[xj]>=i){
        mx[xj]<-xj
        my[xj]<-i
      }else if(!mm){
        return(NULL)
      }
    }
  }
  if (mm) return(data.frame(mx,my))
  else if (qu$size()==0) return(data.frame(mx,my))
  else return(NULL)
}
