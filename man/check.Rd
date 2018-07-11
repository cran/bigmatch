\name{check}
\alias{check}
\title{Check SMDs of the matched data set.}
\usage{
check(fdata, mdata, fz, mz)
}
\arguments{
\item{fdata}{ A full data frame with length(fz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{mdata}{ A matched data frame with length(mz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{fz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the full data set.}
\item{mz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the matched data set.}
}
\value{ A matrix with a row for each variable and five columns being the mean of treated group, mean of matched control group, mean of full control group, SMD of matched control group and SMD of full control group.
}
\description{
The function is used to create a table of mean and SMDs to check the balance before and after matching.
}
\examples{
# To run this example, you must load the optmatch package.

# Caliper of .2 on the propensity score, near fine balance of
# education, a robust Mahalanobis distrance for X.
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr,bmi)
m<-nfmatch(z=z,p=propens,fine=education,X=X,caliper=.2,dat=nh0506,rank=FALSE)
matcheddata<-m$data
Xm<-subset(matcheddata, select=c('female','age','black','hispanic','education','povertyr','bmi'))
check(X,Xm,z,matcheddata$z)
detach(nh0506)
}
