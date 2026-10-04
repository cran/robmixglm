\name{print.outlierTestMix}
\alias{print.outlierTestMix}
\alias{summary.outlierTestMix}
\alias{print.summary.outlierTestMix}
\title{Print an outlierTestMix object}
\description{
Print an outlierTestMix object.
}
\usage{
\method{print}{outlierTestMix}(x,  \ldots)
}
\arguments{
  \item{x}{outlierTestMix object}
  \item{\ldots}{further arguments (not currently used)}
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
summary(forbes.robustmix)
print(outlierTestMix(forbes.robustmix,  cores  =  1))
}
}
\author{Ken Beath}
\keyword{methods}