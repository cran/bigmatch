\name{glover}
\alias{glover}
\title{Perfect matching and maximum matching in a convex bipartite graph.}
\usage{
glover(left, right,mm=F)
}
\arguments{
\item{left}{ Treated person i may be matched to controls j with left[i] <= j <= right[i]. There are length(left) treated individuals and length(left)=length(right). Must have left[i]<=right[i] for every i.  Define maxc = max(right).  Then there are maxc potential control, labeled 1, 2,..., maxc.  The values in left and right are these labels for controls.}
\item{right}{ See left.  }
\item{mm}{ An indicator for whether maximum matching is required. With mm=F, the function tries to find an X-perfect matching. Otherwise, it returns a maximum matching.  }
}
\value{ A perfect matching has every treated individual i matched to a different control j with left[i] <= j <= right[i].  A perfect matching may not exist.  If asked to produce a perfect matching but no perfect matching exists, NULL is returned.  Otherwise, a perfect matching or maximum matching is returned in the form of a data.frame with two columns, mx and my, where mx[i] is the ith treated individual and my[i] is the matched control.}
\examples{
# A perfect matching exists, and glover produces one.
left<-c(2,1,1,4,5)
right<-c(4,3,1,5,5)
glover(left,right)

# No perfect matching exists.
# Treated 4 and treated 5 can only be matched to control 5.
# glover returns NULL for perfect matching.
# glover returns a maximum matching of size 4 with mm=T.

left<-c(2,1,1,5,5)
right<-c(4,3,1,5,5)
glover(left,right)
glover(left,right,mm=TRUE)
}
\description{
Uses Glover's (1967) algorithm to find a perfect matching in a convex bipartite graph.
The implementation uses a priority queue, not used by Glover, as in Figure 1 of Katriel (2008).
Of limited interest to most users; function glover() would typically be called by other functions.
}

\details{

The match produced by glover is rarely useful in observational studies, because it finds a match, not the closest match, not the minimum distance match.

The glover algorithm may be used to find the smallest feasible caliper on a certain criterion.  Each caliper on the criterion creates a different convex bipartite graph.  Use glover to check whether a particular caliper is feasible.  Use binary search to find the smallest caliper that permits a match.  There are many variations on this theme.

The glover algorithm is much faster than optimal matching.  Iterative use of glover's algorithm is often faster than a single minimum distance match.
}

\references{

Glover, F. (1967). Maximum matching in convex bipartite graphs.  Naval Research Logistics Quarterly, 14, 313-316.

Katriel, I. (2008). Matchings in node-weighted convex bipartite graphs.  INFORMS Journal on Computing, 20, 2, 205-211.



}
