\name{summary.robmixglm}
\alias{summary.robmixglm}
\title{summaryficients for robmixglm object}
\description{
Returns summary for a robmixglm object.
}
\usage{
\method{summary}{robmixglm}(object, \ldots)
}
\arguments{
  \item{object}{robmixglm object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
summary
}
\examples{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres,data=forbes)
summary(forbes.robustmix)
}
\author{Ken Beath}
\keyword{methods}