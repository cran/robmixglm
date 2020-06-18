gaussian.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  gaussian.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    sd <- sqrt(mle$deviance/mle$df.residual)
    data$Y <- rnorm(length(data$Y),mean=lp,sd=sd)
    return(data)
  }
  
  gaussian.fun <- function(thedata) {
    gaussian.mle2 <- glm.fit(thedata$X, thedata$Y, 
                             offset = thedata$offset, family = gaussian())
    thesd <- sqrt(gaussian.mle2$deviance/gaussian.mle2$df.residual)
    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(gaussian.mle2)[1:length(coef(gaussian.mle2))],log(0.2/(1-0.2)),1,thesd)
    if (fitno==1) mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(gaussian.mle2$coefficients,ncol=1) + thedata$offset
    gaussian.loglik <- sum(dnorm(thedata$Y,mean=lp,sd=thesd,log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gaussian.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gaussian.mle <- glm.fit(object$X, object$Y, 
                          offset = object$offset, family = gaussian())
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=gaussian.mle)
  fitno <- 1
  gaussian.boot <- boot(thedata, gaussian.fun, R, sim="parametric",
                        ran.gen = gaussian.rg, mle=gaussian.mle,parallel=parallel,
                        ncpus = cores)
  return(list(pos=sum(gaussian.boot$t[,1]>gaussian.boot$t0[1]),R=R,nullstat=gaussian.boot$t[,1]))
}

nbinom.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  nbinom.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    theta <- mle$theta
    data$Y <- rnbinom(length(data$Y),mu=exp(lp),size=theta)
    return(data)
  }
  
  nbinom.fun <- function(data) {
    thedata <- data.frame(data$X,data$Y)
    
    if (dim(object$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$X)[2]],"+",sep='',collapse="")
    else theformula <- "1+"
    
    if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$offset)
    
    theformula <- substr(theformula,1,nchar(theformula)-1)
    
    if(!is.null(object$offset)) donegbin <- paste("nbinom.mle <- glm.nb(object.Y~",theformula,"+offset(offset),",
                                                        "data=thedata)",sep='',collapse='')
    else donegbin <- paste("nbinom.mle <- glm.nb(object.Y~",theformula,",",
                           "data=thedata)",sep='',collapse='')
    
    eval(parse(text=donegbin))
    thetheta <-  nbinom.mle$theta

    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(nbinom.mle)[1:length(coef(nbinom.mle))],log(0.2/(1-0.2)),1,thetheta)
    if (fitno==1) mixfit <- suppressWarnings(nbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(nbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=starting.values, cores=1))
    lp <- data$X %*% matrix(nbinom.mle2$coefficients,ncol=1) + data$offset
    nbinom.loglik <- sum(dnbinom(thedata$Y,mu=lp,size=thetheta,log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-nbinom.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  #stop("Not implemented yet")
  # fit glm to data
  thedata <- data.frame(object$X,object$Y)
  
  if (dim(object$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$X)[2]],"+",sep='',collapse="")
  else theformula <- "1+"
  
  if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$offset)
  
  theformula <- substr(theformula,1,nchar(theformula)-1)
  
  if(!is.null(object$offset)) donegbin <- paste("nbinom.mle <- glm.nb(object.Y~",theformula,"+offset(offset),",
                                                     "data=thedata)",sep='',collapse='')
  else donegbin <- paste("nbinom.mle <- glm.nb(object.Y~",theformula,",",
                        "data=thedata)",sep='',collapse='')
  
  eval(parse(text=donegbin))
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=negbin.mle)
  fitno <- 1
  nbinom.boot <- boot(thedata, nbinom.fun, R, sim="parametric",
                        ran.gen = nbinom.rg, mle=nbinom.mle,parallel=parallel,
                      ncpus = cores)
  return(list(pos=sum(nbinom.boot$t[,1]>nbinom.boot$t0[1]),R=R,nullstat=nbinom.boot$t[,1]))
}


gamma.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  gamma.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    phi <- sum(mle$residuals^2)/mle$df.residual
     data$Y <- rgamma(length(data$Y),shape=1.0/phi,rate=1.0/(phi*exp(lp)))
    return(data)
  }
  
  gamma.fun <- function(thedata) {
    gamma.mle2 <- glm.fit(thedata$X, thedata$Y, 
                          offset = thedata$offset, family = Gamma(link="log"))
    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(gamma.mle2)[1:length(coef(gamma.mle2))],log(0.5/(1-0.5)),1,sum(gamma.mle2$residuals^2)/gamma.mle2$df.residual)
    if (fitno==1)  mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=NULL, cores=1))
     else mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(gamma.mle2$coefficients,ncol=1) + thedata$offset
    phi <- sum(gamma.mle2$residuals^2)/gamma.mle2$df.residual
    gamma.loglik <- sum(dgamma(thedata$Y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gamma.loglik)
     fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gamma.mle <- glm.fit(object$X, object$Y, 
                       offset = object$offset, family = Gamma(link="log"))
  
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=gamma.mle)
  fitno <- 1
  gamma.boot <- boot(thedata, gamma.fun, R, sim="parametric",
                        ran.gen = gamma.rg, mle=gamma.mle, parallel=parallel,
                     ncpus = cores)
  return(list(pos=sum(gamma.boot$t[,1]>gamma.boot$t0[1]),R=R,nullstat=gamma.boot$t[,1]))
}

truncpoisson.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  truncpoisson.rg <- function(data, mle) {
    
    lp <- data$X %*% matrix(mle@coefficients,ncol=1) + data$offset
    data$Y <- actuar::rztpois(length(data$Y),exp(lp))
    return(data)
  }
  
  truncpoisson.fun <- function(data) {
 
    thedata <- data.frame(data$X,data$Y)
    
    if (dim(data$X)[2]>1) theformula <- paste(names(thedata)[2:(dim(data$X)[2])],"+",sep='',collapse="")
    else theformula <- "1+"
    
    if(!is.null(data$offset)) thedata <- data.frame(thedata,offset=data$offset)
    
    theformula <- substr(theformula,1,nchar(theformula)-1)
    
    if(!is.null(object$offset)) dotrunc <- paste("truncpoisson.mle2 <- vglm(data.Y~",theformula,",",
                                                       "family=pospoisson, offset=offset,",
                                                       "data=thedata)",sep='',collapse='')
    else dotrunc <- paste("truncpoisson.mle2 <- vglm(data.Y~",theformula,",",
                          "family=pospoisson,",
                          "data=thedata)",sep='',collapse='')

    eval(parse(text=dotrunc))
      
    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(truncpoisson.mle2)[1:length(coef(truncpoisson.mle2))],log(0.5/(1-0.5)),1)
    if (fitno==1) mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=NULL, cores=1))
    else mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=starting.values, cores=1))
    lp <- data$X %*% matrix(truncpoisson.mle2@coefficients,ncol=1) + data$offset
    truncpoisson.loglik <- sum(dztpois(data$Y,exp(lp),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-truncpoisson.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  thedata <- data.frame(object$X,object$Y)

  if (dim(object$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$X)[2]],"+",sep='',collapse="")
  else theformula <- "1+"
  
  if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$offset)
  
  theformula <- substr(theformula,1,nchar(theformula)-1)
  
  if(!is.null(object$offset)) dotrunc <- paste("truncpoisson.mle2 <- vglm(object.Y~",theformula,",",
                                        "family=pospoisson, offset=offset,",
                                        "data=thedata)",sep='',collapse='')
  else dotrunc <- paste("truncpoisson.mle2 <- vglm(object.Y~",theformula,",",
                        "family=pospoisson,",
                        "data=thedata)",sep='',collapse='')
  
  eval(parse(text=dotrunc))
  fitno <- 1
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=truncpoisson.mle2)
  
  truncpoisson.boot <- boot(thedata, truncpoisson.fun, R, sim="parametric",
                            ran.gen = truncpoisson.rg, mle=truncpoisson.mle2, parallel=parallel,
                            ncpus = cores)
  return(list(pos=sum(truncpoisson.boot$t[,1]>truncpoisson.boot$t0[1]),R=R,nullstat=truncpoisson.boot$t[,1]))
}



poisson.outlierTest.robmixglm <- function(object, R, parallel, cores) {

    poisson.rg <- function(data, mle) {
      
      # calculate linear predictor
      lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
      data$Y <- rpois(length(data$Y),exp(lp))
      return(data)
    }

    poisson.fun <- function(thedata) {
      poisson.mle2 <- glm.fit(thedata$X, thedata$Y, 
                              offset = thedata$offset, family = poisson())
      # starting values assume 50% outliers and tau^2 of 1
      starting.values <- c(coef(poisson.mle2)[1:length(coef(poisson.mle2))],log(0.5/(1-0.5)),1)
      if (fitno==1) mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=TRUE,
                                        verbose=object$verbose, starting.values=NULL, cores=1))
      else mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                        verbose=object$verbose, starting.values=starting.values, cores=1))
      lp <- thedata$X %*% matrix(poisson.mle2$coefficients,ncol=1) + thedata$offset
      poisson.loglik <- sum(dpois(thedata$Y,exp(lp),log=TRUE))
      sim.chisq <- 2*(mixfit$logLik-poisson.loglik)
    	fitno <<- fitno+1
      return(sim.chisq)
    }
    
    # fit glm to data
  poisson.mle <- glm.fit(object$X, object$Y, 
                         offset = object$offset, family = poisson())
  
 
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=poisson.mle)
  fitno <- 1
  poisson.boot <- boot(thedata, poisson.fun, R, sim="parametric",
                       ran.gen = poisson.rg, mle=poisson.mle, parallel=parallel,
                       ncpus = cores)
  return(list(pos=sum(poisson.boot$t[,1]>poisson.boot$t0[1]),R=R,nullstat=poisson.boot$t[,1]))
}

binomial.outlierTest.robmixglm <- function(object, R, parallel, cores) {
  
  binomial.rg <- function(data, mle) {
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    n <- data$Y[,1]+data$Y[,2]
    data$Y[,1] <- rbinom(dim(data$Y)[1],data$Y[,1]+data$Y[,2],1.0/(1.0+exp(-lp)))
    data$Y[,2] <- n-data$Y[,1]
    return(data)
  }
  
  binomial.fun <- function(thedata) {
    binomial.mle2 <- glm.fit(thedata$X, thedata$Y,
                             offset = thedata$offset, family = binomial())
    starting.values <- c(coef(binomial.mle2),log(0.5/(1-0.5)),1)
    if (fitno==1) mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL, cores=1))
    else mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values, cores=1))
    lp <- thedata$X %*% matrix(binomial.mle2$coefficients,ncol=1) + thedata$offset
    binomial.loglik <- sum(dbinom(thedata$Y[,1],thedata$Y[,1]+thedata$Y[,2],1.0/(1.0+exp(-lp)),log=TRUE))
     sim.chisq <- 2*(mixfit$logLik-binomial.loglik)
     fitno <<- fitno+1
     return(sim.chisq)
  }
  
  # fit glm to data
  binomial.mle <- glm.fit(object$X, object$Y,
                          offset = object$offset, family = binomial())
  
  
  thedata <- list(X=object$X,Y=object$Y,offset=object$offset,mle=binomial.mle)
  fitno <- 1
  binomial.boot <- boot(thedata, binomial.fun, R, sim="parametric",
                       ran.gen = binomial.rg, mle=binomial.mle, parallel=parallel,
                       ncpus = cores)
  return(list(pos=sum(binomial.boot$t[,1]>binomial.boot$t0[1]),R=R,nullstat=binomial.boot$t[,1]))
}

outlierTest <-
  ## Short form for generic function
  function(object, R = 999, cores = max(detectCores() - 1, 1))
    UseMethod("outlierTest")

outlierTest.robmixglm <- function(object, R = 999, cores = max(detectCores() - 1, 1)) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  if (cores>1) {
    if(.Platform$OS.type=="unix") parallel <- "multicore"
    else parallel <- "snow"
  } else parallel <- "no"
  out <- switch (
    object$family,
    gaussian = gaussian.outlierTest.robmixglm(object, R, parallel, cores),
    binomial = binomial.outlierTest.robmixglm(object, R, parallel, cores),
    poisson = poisson.outlierTest.robmixglm(object, R, parallel, cores),
    gamma = gamma.outlierTest.robmixglm(object, R, parallel, cores),
    truncpoisson = truncpoisson.outlierTest.robmixglm(object, R, parallel, cores),
    nbinom = nbinom.outlierTest.robmixglm(object, R, parallel, cores)
  )
  class(out) <- "outlierTest"
  return(out)
}
