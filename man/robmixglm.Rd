\name{robmixglm}
\alias{robmixglm}
\title{Fits a Robust Generalized Linear Model and Variants}
\description{
Fits robust generalized linear models and variants described in Beath (2018).
}
\usage{
robmixglm(formula,family=c("gaussian","binomial","poisson",
"gamma","truncpoisson"),data,offset=NULL,
quadpoints=21,notrials=20,EMTol=1.0e-4,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula with columns containing outcomes on the left and covariates on the right.}
  \item{family}{Distribution of response}
  \item{data}{Data frame from which variables are obtained}
  \item{offset}{Offset to be incorporated in the linear predictor.}
  \item{quadpoints}{Number of quadrature points used in the Gauss-Hermite integration.}
  \item{notrials}{Number of random starting values to be used for EM}
  \item{EMTol}{Relative change in likelihood for completion of EM algorithm before switching to quasi-Newton}
  \item{verbose}{Print out diagnostic information? This includes the likelihood and parameter estimates for each EM run.}
}
\value{
  robmixglm object. This contains
\item{call}{Call to function}
\item{family}{Family of model to be fitted}
\item{model}{List of model components. These are X, the response, Y, the predictors, and offset, the offset}
\item{quadpoints}{Number of quadrature points.}
\item{notrials}{Number of trials}
\item{EMTol}{EM tolerance used}
\item{verbose}{Was verbose output requested?}
}
\references{Beath, K. J. A mixture-based approach to robust analysis of generalised linear models, Journal of Applied Statistics, 45(12), 2256-2268 (2018) DOI: 10.1080/02664763.2017.1414164}
\author{Ken Beath}
\examples{
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(100*log10(pres)~bp,data=forbes)
}
\keyword{multivariate}
