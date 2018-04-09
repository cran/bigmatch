\name{nearfine}
\alias{nearfine}
\title{Minimum-distance near-fine matching.}
\usage{
nearfine(z, fine, dist, dat, X, ncontrol = 1, penalty=1000,
max.cost=penalty/10, sub=F, subX)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for
a control. Must have treated subjects (z=1) before controls (z=0).}
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels
that are to be nearly-finely balanced.}
\item{dist}{ A distance list with the starting node (treated subjec), ending
node (control), the distance between them and whether nearexact is needed for each pair.}
\item{dat}{ A data frame with length(z) rows.  One output of the program is
a data frame with some of the rows of dat for the matched sample, together
with additional columns describing the match.}
\item{X}{ A matrix of covariates used to create a robust Mahalanobis distance.
X must have length(z) rows.}
\item{ncontrol}{ A positive integer giving the number of controls to be
matched to each treated subject.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{max.cost}{ The maximum cost for the each pair of treated and control
while rounding the cost. 2 times it is the cost for nearexact matching}
\item{sub}{ An indicator of whether optimal subset matching is required.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on, i.e. only minimum of number of treated and controls will be matched.}
}

\value{
\item{feasible}{feasible=1 if the match is feasible or feasible=0 if the match
is infeasible.}
\item{timeinrelax}{Time in RELAX IV spent computing the minimum cost flow.}
\item{timeinnet}{Time in constructing the network.}
\item{timeinmatch}{Time in constructing the matched dataset.}
\item{d}{The matched sample.  Selected rows of dat.}
\item{number}{ Number of edges between the treated subjects and controls in
the reduced network.}
}

\description{
The program finds an optimal near-fine match.
}


\details{

The match minimizes the total distance between treated subjects and their
matched controls subject to a near-fine balance constraint imposed as a
penalty on imbalances.

The distance list only includes pairs closed based on the caliper, i.e. some
edges are removed from the network.  Because of this, the match may be
infeasible.  This is reported in feasible.

For discussion of networks for fine-balance, see Rosenbaum (1989, Section 3)
and Rosenbaum (2010).

For near-fine balannce balance, see Yang et al. (2012).

You MUST install and load the rcbalance package to use nearfinelate.
}
\examples{
\dontrun{
# To run this example, you must load the optmatch package.
data(nh0506)
attach(nh0506)
X<-cbind(female,age,black,hispanic,education,povertyr)
dist<-smahal(z,propens,X,0.3)
fine<-education
detach(nh0506)

m<-nearfine(nh0506$z,fine,dist,nh0506,X)
}
}
\references{

Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost
network flow problems. Annals of Operations Research, 13, 125-190.  Fortran
and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the
optmatch package.

Rosenbaum, P.R. (1989) Optimal matching in observational studies. Journal of
the American Statistical Association, 84, 1024-1032.

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

Yang, D., Small, D. S., Silber, J. H., and Rosenbaum, P. R. (2012) Optimal
matching with minimal deviation from fine balance in a study of obesity and
surgical outcomes.  Biometrics, 68, 628-636.

}
