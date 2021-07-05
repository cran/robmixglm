\name{robmixglm}
\alias{robmixglm}
\title{Fits a Robust Generalized Linear Model and Variants}
\description{
Fits robust generalized linear models and variants described in Beath (2018).
}
\usage{
robmixglm(formula, family = c("gaussian", "binomial", "poisson", 
"gamma", "truncpoisson", "nbinom"), data, offset = NULL, 
quadpoints = 21, notrials = 50, EMTol = 1.0e-4, cores  =  max(detectCores() \%/\% 2,  1), 
verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{Model formula}
  \item{family}{Distribution of response}
  \item{data}{Data frame from which variables are obtained}
  \item{offset}{Offset to be incorporated in the linear predictor.}
  \item{quadpoints}{Number of quadrature points used in the Gauss-Hermite integration.}
  \item{notrials}{Number of random starting values to be used for EM}
  \item{EMTol}{Relative change in likelihood for completion of EM algorithm before switching to quasi-Newton}
  \item{cores}{Number of cores to be used for parallel evaluation of starting values}
  \item{verbose}{Print out diagnostic information? This includes the likelihood and parameter estimates for each EM run.}
}
\value{
  robmixglm object. This contains
\item{fit}{Final model fit from quasi-Newton}
\item{prop}{Posterior probability of observation in each class}
\item{logLik}{final log likelihood}
\item{np}{Number of parameters}
\item{nobs}{Number of observations}
\item{coef.names}{Coefficient names}
\item{call}{Call to function}
\item{family}{Family of model to be fitted}
  \item{model}{model}
  \item{terms}{terms}
  \item{xlevels}{Levels for factors.}
  \item{quadpoints}{Number of quadrature points used in the Gauss-Hermite integration.}
  \item{notrials}{Number of random starting values to be used for EM}
  \item{EMTol}{Relative change in likelihood for completion of EM algorithm before switching to quasi-Newton}
\item{verbose}{Was verbose output requested?}
}
\details{
Fits robust generalized models assuming that data is a mixture of standard observations and outlier abservations, which belong to an overdispersed model (Beath, 2018). For binomial, Poisson, truncated Poisson and gamma, the overdispersed component achieved through including a random effect as part of the linear predictor, as described by Aitkin (1996).  For gaussian and negative binomial data the outlier component is also a gaussian and negative binomial model, respectively but with a higher dispersion. For gaussian this corresponds to a higher value of \eqn{\sigma^2} but for negative binomial this is a lower value of \eqn{\theta}.

The method used is a generalised EM. Random starting values are determined by randomly allocating observations to either the standard or outlier class for the first iteration of the EM. The EM is then run to completion for all sets of starting values. The best set of starting values is then used to obtain the final results using a quasi-Newton method. Where the overdispersed data is obtained using a random effect, the likelihood is obtained by integrating out the random effect using Gauss-Hermite quadrature.
}
\references{Beath, K. J. A mixture-based approach to robust analysis of generalised linear models,  Journal of Applied Statistics, 45(12), 2256-2268 (2018) DOI: 10.1080/02664763.2017.1414164

Aitkin, M. (1996). A general maximum likelihood analysis of overdispersion in generalized linear models. Statistics and Computing, 6, 251262. DOI: 10.1007/BF00140869}
\author{Ken Beath}
\examples{
\donttest{
if (requireNamespace("MASS", quietly = TRUE)) {
library(MASS)
data(forbes)
forbes.robustmix <- robmixglm(100*log10(pres)~bp, data = forbes, cores = 1)
}
}}
\keyword{multivariate}
