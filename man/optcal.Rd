\name{optcal}
\alias{optcal}
\title{Finds the optimal caliper width.}
\usage{
optcal(z, p, exact=NULL, ncontrol=1, tol=NULL, rank=TRUE, subX=NULL)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{p}{ A vector with the same length as z. The best caliper for p is found. Often p is the propensity score.}
\item{exact}{ If exact is NULL, then there is no exact matching, and the caliper refers to p alone.  Otherwise, exact is a vector of the same length as z for exact matching, such that two individuals, i and j, can be matched only if they have the same value of exact, exact[i]=exact[j].  In this case, the caliper is best among calipers that permit exact matching for exact.  Typically, exact has a moderate number of possible values, far fewer than length(z).}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{tol}{ The tolerance.  The optimal caliper is determined with an error of at most tol.  tol=0.01 might be used for the propensity score, as it takes values between 0 and 1, whereas tol=1/2 for p=age would mean that the caliper for age errs by at most half a year.}
\item{rank}{ An indicator of whether we want a caliper on rank of p or p.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
}
\value{
\item{caliper}{The optimal caliper, with an error of at most tol.  This caliper is a little too large, at most tol too large, but because its error is on the high side, a match with this caliper does exist.}
\item{interval}{An interval that contains the best caliper. The upper bound of the interval was returned as caliper above.}
\item{interval.length}{The length of interval. By definition, length.interval<=tol.}
}
\description{
Finds the smallest caliper on variable p or rank of p such that a treated-control matching with that caliper exists.  If exact is not NULL, then finds the smallest caliper on p or rank of p such that a treated-control matching with that caliper exists while also matching exactly for the variable exact.
}
\details{
The method uses binary search to find the optimal caliper.  At each step in the search, it applies Glover's algorithm to determine whether a proposed caliper is feasible.

Often, we need a small and feasible caliper, but we do not need to determine the optimal caliper very precisely.  Making tol larger will reduce the number of steps in the binary search.
}
\examples{
data(nh0506)
attach(nh0506)

#optimal caliper using the propensity score alone
optcal(z,propens,tol=0.1,rank=FALSE)

#optimal caliper using the rank of propensity score
#and match each treated subject with two controls
\donttest{
optcal(z,propens,ncontrol=2,rank=TRUE)
}

#optimal caliper for the propensity score while requiring
#an exact match for female
optcal(z,propens,exact=female,tol=0.1,rank=FALSE)

detach(nh0506)

}
\references{
Glover, F. (1967). Maximum Matching In Convex Bipartite Graphs.  Naval Research Logistics Quarterly, 14, pp 313-316.

Lipski, W., Jr, and Preparata, F. P. (1981). Efficient Algorithms For Finding Maximum Matchings In Convex Bipartite Graphs And Related Problems.  Journal Acta Informatica, 15, 4, pp 329-346.
}
