\name{glover}
\alias{glover}
\title{Maximum matching in a convex bipartite graph.}
\usage{
glover(left, right)
}
\arguments{
\item{left}{ Treated person i may be matched to controls j with left[i] <= j <= right[i]. There are length(left) treated individuals and length(left)=length(right). Must have left[i]<=right[i] for every i.  Define maxc = max(right).  Then there are maxc potential control, labeled 1, 2,..., maxc.  The values in left and right are these labels for controls.}
\item{right}{ See left.  }
}
\value{ The maximum matching ratio. A perfect matching has every treated individual i matched to a different control j with left[i] <= j <= right[i].  A perfect matching may not exist. A number smaller than 1 is returned if no perfect matching exists.  Otherwise, 1 is returned.}
\examples{
# A perfect matching exists, and glover produces one.
left<-c(2,1,1,4,5)
right<-c(4,3,1,5,5)
glover(left,right)

# No perfect matching exists, and glover returns maximum matching ratio.
# Treated 4 and treated 5 can only be matched to control 5.

left<-c(2,1,1,5,5)
right<-c(4,3,1,5,5)
glover(left,right)
}
\description{
Uses Glover's (1967) algorithm to find a maximum matching in a doubly convex bipartite graph. The implementation uses a priority queue, not used by Glover, as in Lipski and Preparata (1981).
Of limited interest to most users; function glover() would typically be called by other functions.
}

\details{
The match produced by glover is rarely useful in observational studies, because it finds a match, not the closest match, not the minimum distance match.

The glover algorithm may be used to find the smallest feasible caliper on a certain criterion.  Each caliper on the criterion creates a different convex bipartite graph. Use glover to check whether a particular caliper is feasible. Use bisection search to find the smallest caliper that permits a match. There are many variations on this theme.

The glover algorithm is much faster than optimal matching.  Iterative use of glover's algorithm is often faster than a single minimum distance match.
}

\references{
Glover, F. (1967). Maximum Matching In Convex Bipartite Graphs.  Naval Research Logistics Quarterly, 14, pp 313-316.

Lipski, W., Jr, and Preparata, F. P. (1981). Efficient Algorithms For Finding Maximum Matchings In Convex Bipartite Graphs And Related Problems.  Journal Acta Informatica, 15, 4, pp 329-346.
}
