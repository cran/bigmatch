glover<-function(left,right){

  nx<-length(left)
  ny<-max(right)
  stopifnot(length(left)==nx)
  stopifnot(length(right)==nx)
  stopifnot(all(left<=right))

  oright<-order(right)
  lefto<-left[oright]
  righto<-right[oright]
  oleft<-order(lefto)

  mx<-rep(NA,nx)
  my<-rep(NA,nx)
  qu<-liqueueR::PriorityQueue$new()
  nb<-1
  ne<-1
  for (i in 1:ny) {
    while ((nb<=nx) & (lefto[oleft[nb]]==i)){
      qu$push(oleft[nb],priority=-oleft[nb])
      nb<-nb+1
    }
    if (qu$size()>0){
      xj<-qu$pop()
      mx[xj]<-oright[xj]
      my[xj]<-i

    }
    while ((ne<=nx)&(righto[ne]==i)){
      if ((qu$size()>0) & (ne>xj)) qu$pop()
      ne<-ne+1
    }
  }
  return(1-sum(is.na(mx))/nx)
}
