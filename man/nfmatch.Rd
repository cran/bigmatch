\name{nfmatch}
\alias{nfmatch}
\title{Minimum-distance near-fine matching.}
\usage{
nfmatch(z, p, fine=rep(1,length(z)), X, dat, caliper, constant=NULL, ncontrol=1,
rank=TRUE, exact=NULL, penalty=1000, max.cost=penalty/10, nearexact=NULL,
nearexPenalty=max.cost, Xextra=NULL, weight=NULL, subX=NULL, ties.all=TRUE, seed=1)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{p}{ A vector of with length(z)=length(p) giving the numeric values of the variable to be matched with a caliper.  Typically, p is the propensity score.  If p takes a few levels, exact matching for p is attempted if caliper=0. If no caliper is intended to apply, set p=rep(1,length(z)) and caliper=1.}
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{X}{ A matrix of covariates used to create a robust Mahalanobis distance.  X must have length(z) rows.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair.  If caliper is too small, the match may be infeasible. If no caliper is intended to apply, set p=rep(1,length(z)) and caliper=1.}
\item{constant}{ If there are more than constant controls for a treated differ on p within caliper, we select the constant closest controls.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{rank}{ An indicator of whether we want a caliper on rank of p or p.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{max.cost}{ The maximum cost for the each pair of treated and control while rounding the cost.}
\item{nearexact}{ If not NULL, then a vector of length length(z) or matrix with length(z) rows giving variables that need to be exactly matched. If it is not possible to exactly match all variables, we will exactly match as many variables as we can.}
\item{nearexPenalty}{ The penalty for a mismatch on nearexact. If it is a number, then use the same penalty for all nearexact variables. Otherwise, it should be a vector of length the same as number of nearexact variables, indicating the penalty for mismatch on each nearexact variable. The larger nearexPenalty is, the more priorty the variable get in near-exact match.}
\item{Xextra}{ If not NULL, another robust Mahalanobis distance based on Xextra  is calculated. The distance between treated-control pair is a weighted sum of the two distances.}
\item{weight}{ The weight for Mahalanobis distance of Xextra.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
\item{ties.all}{ If ties.all is True, include all ties while choosing nearest neighbors. In this case, some treated may have more than constant controls. Otherwise, randomly select one or several controls to make sure there are not more than constant controls for each treated.}
\item{seed}{ When ties.all is False, seed for randomly select one or several controls when there are ties.}
}

\value{
If the match is infeasible, a warning is issued.  Otherwise, a list of results is returned.

A match may be infeasible if the caliper is too small, or ncontrol is too large, or if exact matching for exact is impossible.
\item{data}{The matched sample.  Selected rows of dat.}
\item{timeinrelax}{Time in RELAX IV spent computing the minimum cost flow.}
\item{edgenum}{ Number of edges between the treated subjects and controls in
the reduced network.}
\item{timeind}{Time in calculating robust Mahalanobis distance between connected pairs.}
\item{timeinnet}{Time in constructing the network.}
\item{timeinmatch}{Time in constructing the matched dataset.}
}

\description{
The program finds an optimal near-fine match with a given caliper on p or rank of p.
}


\details{
The match minimizes the total distance between treated subjects and their matched controls
subject to a near-fine balance constraint imposed as a penalty on imbalances.

For discussion of networks for fine-balance, see Rosenbaum (1989, Section 3) and Rosenbaum (2010).
For near-fine balannce balance, see Yang et al. (2012).

You MUST install and load the optmatch package to use nearfine.
}
\examples{
# To run this example, you must load the optmatch package.

# Caliper of .3 on the propensity score, near fine balance of
# education, a robust Mahalanobis distrance for X.
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr,bmi)
m<-nfmatch(z=z,p=propens,fine=education,X=X,caliper=.3,dat=nh0506,rank=FALSE)
matcheddata=m$data
table(matcheddata$z,matcheddata$education)
head(matcheddata)
detach(nh0506)
\donttest{
  #finds the optimal caliper for the propensity score while exact matching on female
  #near fine balance for education and hispanic jointly.
  data(nh0506)
  attach(nh0506)
  X<-cbind(female,age,black,hispanic,education,povertyr,bmi)

  oc<-optcal(z,propens,exact=female,tol=0.1,rank=FALSE)
  oc
  oco<-optconstant(z,propens,oc$caliper,exact=female,rank=FALSE)
  oco
  m2<-nfmatch(z,propens,factor(hispanic):factor(education),X,nh0506,oc$caliper,oco$constant,
              exact=female,rank=FALSE)

  matcheddata2=m2$data
  table(matcheddata2$z,matcheddata2$female)
  table(matcheddata2$z,matcheddata2$education)
  table(matcheddata2$z,matcheddata2$education,matcheddata2$hispanic)

  #finds the optimal caliper for the propensity score while exact matching on female
  #nearexact on quantiles of povertyr and bmi
  #near fine balance for education and hispanic jointly.
  pq=cut(povertyr,c(-0.1,1,2,3,4,5))
  bq=cut(bmi,(0:7)*20)
  #first assume povertyr and bmi are of the same importance
  m3<-nfmatch(z,propens,factor(hispanic):factor(education),X,nh0506,oc$caliper,oco$constant,
              exact=female,nearexact=cbind(pq,bq),rank=FALSE)
  matcheddata3=m3$data
  head(matcheddata3)

  #then assume povertyr is more important than bmi
  m4<-nfmatch(z,propens,factor(hispanic):factor(education),X,nh0506,oc$caliper,oco$constant,
              exact=female,nearexact=cbind(pq,bq),nearexPenalty=c(100,50),rank=FALSE)
  matcheddata4=m4$data
  head(matcheddata4)
  detach(nh0506)
}
}
\references{
Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.  Fortran and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the optmatch package.

Rosenbaum, P.R. (1989) Optimal matching in observational studies. Journal of the American Statistical Association, 84, 1024-1032.

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

Yang, D., Small, D. S., Silber, J. H., and Rosenbaum, P. R. (2012) Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes.  Biometrics, 68, 628-636.
}
