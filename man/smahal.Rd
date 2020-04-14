\name{smahal}
\alias{smahal}
\title{Creates a robust Mahalanobis distance for matching.}
\usage{
smahal(z, p, X, caliper, constant=NULL, ncontrol=1, exact=NULL,
nearexact=NULL, nearexPenalty=100, Xextra=NULL, weight=NULL, subX=NULL, ties.all=TRUE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{p}{ A vector of length(z)=length(p) giving the variable used to define the caliper.  Typically, p is the propensity score or its rank.}
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair.}
\item{constant}{ If the number of pairs within a caliper is greater than constant, we will select the constant closest ones.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) or matrix with length(z) rows giving variables that need to be exactly matched. If it is not possible to exactly match all variables, we will exactly match as many variables as we can.}
\item{nearexPenalty}{ Penalty for mismatch on nearexact if nearexact is not NULL.}
\item{Xextra}{ If not NULL, another robust Mahalanobis distance based on Xextra is calculated. The distance between treated-control pair is a weighted sum of the two distances.}
\item{weight}{ The weight for Mahalanobis distance of Xextra.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
\item{ties.all}{ If ties.all is True, include all ties while choosing nearest neighbors. In this case, some treated may have more than constant controls. Otherwise, randomly select one or several controls to make sure there are not more than constant controls for each treated.}
}

\value{
\item{d}{ A distance list for each pair within the caliper distance and constant constraint.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
\item{nearex}{ A vector or matrix with the same dimesion as nearexact. Its entry is TRUE for a connected pair whose nearexact variable values are different and is FALSE otherwise. }
}
\description{
Computes a robust Mahalanobis distance list for use in matching.

This function and its use are discussed in Rosenbaum (2010).  The robust
Mahalanobis distance in described in Chapter 8 of Rosenbaum (2010).

}
\details{
The usual Mahalanobis distance works well for multivariate Normal covariates, but can exhibit odd behavior with typical covariates. Long tails or an outlier in a covariate can yield a large estimated
variance, so the usual Mahalanobis distance pays little attention to large differences in this covariate. Rare binary covariates have a small variance, so a mismatch on a rare binary covariate is
viewed by the usual Mahalanobis distance as extremely important. If you were matching for binary covariates indicating US state of residence, the usual Mahalanobis distance would regard a
mismatch for Wyoming as much worse than a mismatch for California.

The robust Mahalanobis distance uses ranks of covariates rather than the covariates themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a variable more important. Binary covariates are, of course, heavily tied.
}
\examples{
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr)
dist<-smahal(z,propens,X,0.2)
dist$d[1:10]

dist2<-smahal(z,propens,X,0.2,exact=female,Xextra=hispanic,nearexact=bmi)
dist2$nearex[1:10]

detach(nh0506)
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
