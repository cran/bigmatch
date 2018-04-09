\name{optconstant}
\alias{optconstant}
\title{Finds the optimal constant.}
\usage{
optconstant(z, p, caliper=NULL,exact=NULL,ncontrol=1,tol=1,rank=T)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{p}{ A vector with the same length as z.  The best constant for p is found. Often p is the propensity score.}
\item{caliper}{ If two individuals differ on p by more than caliper, there is no edge between this pair.  If caliper is too small, the match may be infeasible.}
\item{exact}{ If exact is NULL, then there is no exact matching, and the constant refers to p alone.  Otherwise, exact is a vector or matrix of important covariates that need to be exactly matched. We will select as many important covariates as we can and construct a new vector Ex of the same length as z for exact matching, such that two individuals, i and j, can be matched only if they have the same value of Ex, Ex[i]=Ex[j].  In this case, the constant is best among calipers that permit exact matching for exact.  Typically, Ex has a moderate number of possible values, far fewer than length(z).}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{tol}{ The tolerance.  The optimal constant is determined with an error of at most tol.}
\item{rank}{ An indicator of whether we want a constant on rank of p or p.}
}
\value{
\item{constant}{The optimal constant, with an error of at most tol.  This constant is a little too large, at most tol too large, but because its error is on the high side, a match with this constant does exist.}
\item{interval}{An interval that contains the best constant.  The upper bound of the interval was returned as constant above.}
\item{interval.length}{The length of interval.  By definition, length.interval<=tol.}
}
\description{
Finds the smallest constant on variable p or rank of p such that a treated-control matching with that constant exists.  If exact is not NULL, then finds the smallest constant on p or rank of p such that a treated-control matching with that constant exists while also matching exact for the variable exact.

}
\details{
The method uses binary search to find the optimal constant.  At each step in the search, it applies Glover's algorithm to determine whether a proposed constant is feasible.

Often, we need a small and feasible constant, but we do not need to determine the optimal constant very precisely.  Making tol larger will reduce the number of steps in the binary search.
}
\examples{
\donttest{
data(nh0506)
attach(nh0506)

#optimal constant using the propensity score alone
optconstant(z,propens,rank=FALSE)

#optimal constant for the propensity score while requiring
#an exact match for female
optconstant(z,propens,exact=female,rank=FALSE)

#optimal constant for the propensity score given a caliper
oc=optcal(z,propens,tol=0.1,rank=FALSE)
optconstant(z,propens,caliper=oc$caliper,rank=FALSE)
detach(nh0506)
}
}
\references{

Glover, F. (1967). Maximum matching in convex bipartite graphs.  Naval Research Logistics Quarterly,
14, 313-316.

Katriel, I. (2008). Matchings in node-weighted convex bipartite graphs.  INFORMS Journal on
Computing, 20, 2, 205-211.

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}