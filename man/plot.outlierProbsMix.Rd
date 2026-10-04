\name{plot.outlierProbsMix}
\alias{plot.outlierProbsMix}
\title{
Plot outlier probabilities.
}
\description{
Plots the outlier probability for each observation,  from an outlierProbs object.
}
\usage{
\method{plot}{outlierProbsMix}(x,  ...)
}
\arguments{
  \item{x}{
outlierProbsMix object to be plotted
}
 \item{\ldots}{
additional parameters to plot}
}
\value{
Plot
}
\author{
Ken Beath <ken.beath@mq.edu.au>
}
\examples{
\donttest{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(bp~pres, data = forbes, cores = 1)
plot(outlierProbsMix(forbes.robustmix))
}}
\keyword{method}
