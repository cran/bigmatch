\name{netfine}
\alias{netfine}
\title{Optimal near-fine match from a distance matrix.}
\usage{
netfine(z, fine, dist, ncontrol=1, penalty=1000, max.cost=penalty/10,
nearexPenalty=max.cost, subX=NULL)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for
a control. }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels
that are to be nearly-finely balanced.}
\item{dist}{ A distance list with the starting node (treated subjec), ending
node (control), the distance between them and whether nearexact is needed for each pair.}
\item{ncontrol}{ A positive integer giving the number of controls to be
matched to each treated subject.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{max.cost}{ The maximum cost for the each pair of treated and control
while rounding the cost. }
\item{nearexPenalty}{ The penalty for a mismatch on nearexact. If it is a number, then use the same penalty for all nearexact variables. Otherwise, it should be a vector of length the same as number of nearexact variables, indicating the penalty for mismatch on each nearexact variable. The larger nearexPenalty is, the more priorty the variable get in near-exact match.}
\item{subX}{ If a subset matching is required, the variable that the subset matching is based on. That is, for each level of subX, extra treated will be discarded in order to have the number of matched treated subjects being the minimum size of treated and control groups. If exact matching on a variable x is desired and discarding extra treated is fine if there are more treated than controls for a certain level k, set exact=x, subX=x.}
}
\value{ A network for optimal near-fine matching.}
\description{
The function creates the network for optimal near-fine matching to be passed via callrelax to the Fortran code for Bertsekas and Tseng's (1988) Relax IV.

Of limited interest to most users; function netfine() would typically be called by some other functions.
}
\details{
The network contains a bipartite graph for treated and control subjects plus
additional nodes for fine balance categories, plus additional nodes accept
needed deviations from fine balance yielding near-fine balance.

For discussion of fine-balance, see Rosenbaum (1989, Section 3) and Rosenbaum (2010).
For near-fine balance balance, see Yang et al. (2012).
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
