\name{outlierTest}
\alias{outlierTest}
\alias{outlierTest.robmixglm}
\title{
Test for the presence of outliers.}
\description{
Uses the parametric bootstrap to test for the presence of outliers.}
\usage{
outlierTest(object, R=999, showProgress = TRUE)
}
\arguments{
  \item{object}{
A robmixglm object with a mixture (robust) random effects distribution.
}
  \item{R}{number of bootstrap replications}
  \item{showProgress}{Show progress bar. Should be set to FALSE if output is sent to a file.}
}
\value{
An outlierTest object.
}
\details{
Performs a parametric bootstrap to compare models with and without outliers.
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
hospcosts.robustmix <- robmixglm(costs~adm+age+dest+ins+loglos+sex,family="gamma",data=hospcosts)
summary(hospcosts.robustmix)
summary(outlierTest(hospcosts.robustmix))
}}
\keyword{method}
