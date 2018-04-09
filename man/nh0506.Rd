\name{nh0506}
\alias{nh0506}
\docType{data}
\title{
Smoking and homocysteine levels in NHANES 2005-2006.
}
\description{
Bazzano et al. (2003) noted higher homocysteine levels in smokers than in nonsmokers.  See also Pimentel et al. (2016) for a related analysis.
}
\usage{data("nh0506")}
\format{
  A data frame with 2475 observations on the following 32 variables.
  \describe{
    \item{\code{Row}}{a numeric vector}
    \item{\code{SEQN}}{NHANES id number}
    \item{\code{female}}{1 if female, 0 if male}
    \item{\code{age}}{age in years, >=20}
    \item{\code{black}}{1 if black, 0 otherwise}
    \item{\code{hispanic}}{1 if hispanic, 0 otherwise}
    \item{\code{education}}{Education}
    \item{\code{povertyr}}{Ratio of family income to the poverty level, capped at 5x}
    \item{\code{creactiveprotein}}{creactive protein}
    \item{\code{homocysteine}}{homocysteine}
    \item{\code{cotinine}}{cotinine in blood}
    \item{\code{cadmium}}{cadmium in blood}
    \item{\code{lead}}{lead in blood}
    \item{\code{bmi}}{Body mass index}
    \item{\code{cigs100life}}{1 if smoked more than 100 cigarettes in lifetime, 0 otherwise}
    \item{\code{smokenow}}{1 if smokes now, 0 otherwise}
    \item{\code{cigsdays30}}{Days smoked in last 30 days: 0 if never smoker, 30 if daily smoker}
    \item{\code{cigsperday30}}{Cigarettes smoked per day in last 30 days}
    \item{\code{tobacco5days}}{1 =  used tobacco in the last 30 days, 0 otherwise}
    \item{\code{dailysmoker}}{1 = daily smoker, 0 = never smoker}
    \item{\code{neversmoker}}{1 = never smoker, 0 = daily smoker}
    \item{\code{z}}{1 if daily smoker, 0 if never smoker}
    \item{\code{propens}}{Estimated propensity score.  The score was formed
    by logit regression of z on female, age, education, black, hispanic, povertyr, and bmi.}
    \item{\code{pstrat}}{Propensity score strata: \code{(0,0.0733]} \code{(0.0733,0.131]} \code{(0.131,0.204]} \code{(0.204,0.33]} \code{(0.33,1]}}
    \item{\code{age3}}{Age in 3 categories}
    \item{\code{ed3}}{Education in 3 categories}
    \item{\code{bmi3}}{BMI in 3 categories}
    \item{\code{pov2}}{Income above 2 times poverty, TRUE or FALSE}
    \item{\code{stf}}{A factor defining strata using female, age3, ed3, bmi3 pov2. }
    \item{\code{st}}{A numeric version of stf}
    \item{\code{stfp}}{A factor defining strata using stf and pstrat}
    \item{\code{stp}}{A numeric version of stp}
  }
}
\details{
Data from NHANES 2005-2006 concerning homocysteine levels in daily smokers (z=1) and never smokers (z=0), aged 20 and older.  Daily smokers smoked every day for the last 30 days, smoking an average of at least 10 cigarettes per day.  Never smokers smoked fewer than 100 cigarettes in their lives, do not smoke now, and had no tobacco use in the previous 5 days.
}
\source{
NHANES, the US National Health and Nutrition Examination Survey, 2005-2006.
}
\references{
Bazzano, L. A., He, J., Muntner, P., Vupputuri, S. and Whelton, P. K. (2003) Relationship between cigarette smoking and novel risk factors for cardiovascular disease in the United States.  Annals of Internal Medicine, 138, 891-897.

Pimentel, S. D., Small, D. S. and Rosenbaum, P. R. (2016) Constructed second control groups and attenuation of unmeasured biases.  Journal of the American Statistical Association, 111, 1157-1167.
}
\examples{
data(nh0506)
}
\keyword{datasets}
