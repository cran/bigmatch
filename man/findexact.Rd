\name{findexact}
\alias{findexact}
\title{Build an exact match variable given a matrix of covarites ordered by the user according to their importance.}
\usage{
findexact(z,E,ncontrol)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{E}{ A matrix or a vector with length(z) rows giving the covariates that need to be exactly matched.}
\item{ncontrol}{ The number of controls to be matched to each treated individual.}
}
\value{
\item{miss}{The covariates that cannot be exactly matched}
\item{variables}{The covariates that can be exactly matched}
\item{NewExact}{The constructed exact match variable.}
}
\description{
The function is used to build an exact match variable given a matrix of covarites ordered by the user according to their importance. It will select as many important covariates as possible.
}
\examples{
data(nh0506)
attach(nh0506)
# The following example uses all of the variables in the propensity score,propens
ex<-findexact(z,cbind(female,age,education,black,hispanic,povertyr))
head(ex$miss)
head(ex$variables)
table(z,ex$NewExact)
detach(nh0506)
}
