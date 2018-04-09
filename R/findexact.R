findexact<-function(z,E,ncontrol=1){
  if (is.factor(E)){
    levels(E)=1:nlevels(E)
    E<-as.integer(E)
  }
  if (is.vector(E)) E=as.matrix(E,ncol=1)
  k<-dim(E)[2]
  r<-rep(F,k)
  for (i in 1:k){
    r[i]=T
    Esub=E[,r]
    dsub=as.data.frame(cbind(Esub,z))
    colnames(dsub)=c(colnames(E)[r],'z')
    t=plyr::count(dsub,c(colnames(E)[r],'z'))
    nlevel=nrow(t)
    control=t[seq(1,nlevel,2),(sum(r)+2)]
    treated=t[seq(2,nlevel,2),(sum(r)+2)]
    if (!((length(control)==length(treated))&&all(control>=ncontrol*treated))){
      r[i]=F
    }
  }
  Esub=as.matrix(E[,r])
  Emiss=as.matrix(E[,!r])
  Esub=as.data.frame(Esub)
  if (sum(r)==0) list(miss=E,variables=NULL,NewExact=NULL)
  else{
    variables=transform(Esub, ClusterID = as.numeric(interaction(Esub, drop=TRUE)))
    list(miss=Emiss,variables=variables[,1:(dim(variables)[2]-1)], NewExact=variables$ClusterID)
  }
}
