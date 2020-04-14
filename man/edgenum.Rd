\name{edgenum}
\alias{edgenum}
\title{Computes the number of edges in the reduced bipartite.}
\usage{
edgenum(z, p, caliper, constant=NULL, exact=NULL,
ties.all=TRUE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{p}{ A vector of length(z)=length(p) giving the variable used to define the caliper.  Typically, p is the propensity score or its rank.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair.}
\item{constant}{ If the number of pairs within a caliper is greater than constant, we will select the constant closest ones.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{ties.all}{ If ties.all is True, include all ties while choosing nearest neighbors. In this case, some treated may have more than constant controls. Otherwise, randomly select one or several controls to make sure there are not more than constant controls for each treated.}
}

\value{
Number of edges in the reduced bipartite graph with the constraints on caliper and number of nearest neighbors (constant).
}
\description{
Computes the number of edges in the reduced bipartite graph after applying the caliper and number of nearest neighbors (constant).

This function can provide users some idea of the required computation time. Smaller caliper and constant removes more edges, hence accelarates computation, but risks infeasibility.

}
\details{
A given choice of caliper and number of nearest neighbors (constant) removes candidate pairs, so there exists a corresponding reduced bipartite graph.

Smaller caliper and constant removes more edges from the original dense graph, hence the computation is faster. However, this risks infeasibility. A smallest caliper that permits a feasible match and its corresponding smallest number of nearest neighbors can be computed by functions optcal() and optconstant().
}
\examples{
data(nh0506)
attach(nh0506)
edgenum(z,propens,0.2)

edgenum(z,propens,0.2,10,exact=female)

detach(nh0506)
}

