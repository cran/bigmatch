\name{nfmatch}
\alias{nfmatch}
\title{Minimum-distance near-fine matching.}
\usage{
nfmatch(z, p, fine, X, caliper, dat, constant=NULL, exact=NULL, nearexact=T,
rank=T, ncontrol=1, penalty=1000, max.cost=penalty/10, Xextra=NULL, weight=1,
sub=F, subX=NULL)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{p}{ A vector of with length(z)=length(p) giving the numeric values of the variable to be matched with a caliper.  Typically, p is the propensity score.  If p takes a few levels, exact matching for p is attempted if caliper=0.  }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{X}{ A matrix of covariates used to create a robust Mahalanobis distance.  X must have length(z) rows.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair.  If caliper is too small, the match may be infeasible.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the
matched portion of dat is returned with additional columns that define the match.}
\item{constant}{ If there are more than constant controls for a treated differ on p within caliper, we select the constant closest controls.}
\item{exact}{ If not NULL, then a vector or matrix with length(z) rows giving variables (ordered by the user according to their importance) that need to be exactly matched. If it is not possible to exactly match all variables, we will select as many important variables as we can.}
\item{nearexact}{ An indicator for whether nearexact matching is allowed.}
\item{rank}{ An indicator of whether we want a caliper on rank of p or p.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{max.cost}{ The maximum cost for the each pair of treated and control while rounding the cost.}
\item{Xextra}{ If not NULL, another robust Mahalanobis distance based on Xextra  is calculated. The distance between treated-control pair is a weighted sum of the two distances.}
\item{weight}{ The weight for Mahalanobis distance of Xextra.}
\item{sub}{ An indicator of whether optimal subset matching is required.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on, i.e. only minimum of number of treated and controls will be matched.}
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
The program finds an optimal near-fine match.
}


\details{

The match minimizes the total distance between treated subjects and their matched controls
subject to a near-fine balance constraint imposed as a penalty on imbalances.

For discussion of networks for fine-balance, see Rosenbaum (1989, Section 3) and Rosenbaum (2010).
For near-fine balannce balance, see Yang et al. (2012).

You MUST install and load the rcbalance package to use nearfine.
}
\examples{

# To run this example, you must load the rcbalance package.

# Caliper of .3 on the propensity score, near fine balance of
# education, a robust Mahalanobis distrance for X.
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr,bmi)
m<-nfmatch(z=z,p=propens,fine=education,X=X,caliper=.3,dat=nh0506,rank=FALSE)
matcheddata=m$data
table(matcheddata$z,matcheddata$education)
max(abs(matcheddata$propens[matcheddata$z==1]-matcheddata$propens[matcheddata$z==0]))
detach(nh0506)
\donttest{
#Builds an exact match variable, then finds the
#optimal caliper for the propensity score.  Near fine balance for
#education and female jointly.
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr,bmi)

oc<-optcal(z,propens,exact=X,tol=0.1,rank=F)
oc
oco<-optconstant(z,propens,oc$caliper,exact=X,rank=F)
oco
m2<-nfmatch(z,propens,factor(female):factor(education),X,oc$caliper,
nh0506,oco$constant,exact=X,rank=F)
detach(nh0506)

matcheddata2=m2$data
table(matcheddata2$z,matcheddata2$cexact)
table(matcheddata2$z,matcheddata2$female)
table(matcheddata2$z,matcheddata2$education)
table(matcheddata2$z,matcheddata2$education,matcheddata2$female)
}
}
\references{

Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.  Fortran and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the optmatch package.

Rosenbaum, P.R. (1989) Optimal matching in observational studies. Journal of the American Statistical Association, 84, 1024-1032.

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

Yang, D., Small, D. S., Silber, J. H., and Rosenbaum, P. R. (2012) Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes.  Biometrics, 68, 628-636.

}
