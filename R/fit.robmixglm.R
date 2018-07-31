fit.robmixglm <- function(x,y,family,offset,gh=norm.gauss.hermite(21),notrials,EMTol,verbose) {

 out <- switch(family,
          gaussian = gaussian.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL),
          binomial = binomial.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL),
          poisson = poisson.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL),
         gamma = gamma.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL),
         truncpoisson = truncpoisson.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL),
         negbinom = negbinom.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,verbose,starting.values=NULL)
 )
  out
}
