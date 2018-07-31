gaussian.outlierTest.robmixglm <- function(object, R, showProgress) {
  
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
    if (showProgress) {
      pinterval <- R %/% 40
      if (pinterval < 1)  pinterval <- 1
      if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
    }
    if (fitno==1) mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL))
     else mixfit <- suppressWarnings(gaussian.fit.robmixglm(thedata$X, thedata$Y, 
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values))
    lp <- thedata$X %*% matrix(gaussian.mle2$coefficients,ncol=1) + thedata$offset
    gaussian.loglik <- sum(dnorm(thedata$Y,mean=lp,sd=thesd,log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gaussian.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gaussian.mle <- glm.fit(object$model$X, object$model$Y, 
                          offset = object$model$offset, family = gaussian())
  
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=gaussian.mle)
  fitno <- 1
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  gaussian.boot <- boot(thedata, gaussian.fun, R, sim="parametric",
                        ran.gen = gaussian.rg, mle=gaussian.mle)
  if (showProgress)  close(pb)
  return(list(pos=sum(gaussian.boot$t[,1]>gaussian.boot$t0[1]),R=R,nullstat=gaussian.boot$t[,1]))
}

negbinom.outlierTest.robmixglm <- function(object, R, showProgress) {
  
  negbinom.rg <- function(data, mle) {
    
    # calculate linear predictor
    lp <- data$X %*% matrix(mle$coefficients,ncol=1) + data$offset
    theta <- mle$theta
    data$Y <- rnbinom(length(data$Y),mu=exp(lp),size=theta)
    return(data)
  }
  
  negbinom.fun <- function(data) {
    thedata <- data.frame(data$X,data$Y)
    
    if (dim(object$model$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$model$X)[2]],"+",sep='',collapse="")
    else theformula <- "1+"
    
    if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$model$offset)
    
    theformula <- substr(theformula,1,nchar(theformula)-1)
    
    if(!is.null(object$model$offset)) donegbin <- paste("negbinom.mle <- glm.nb(object.model.Y~",theformula,"+offset(offset),",
                                                        "data=thedata)",sep='',collapse='')
    else donegbin <- paste("negbinom.mle <- glm.nb(object.model.Y~",theformula,",",
                           "data=thedata)",sep='',collapse='')
    
    eval(parse(text=donegbin))
    thetheta <-  negbinom.mle$theta

    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(negbinom.mle)[1:length(coef(negbinom.mle))],log(0.2/(1-0.2)),1,thetheta)
    if (showProgress) {
      pinterval <- R %/% 40
      if (pinterval < 1)  pinterval <- 1
      if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
    }
    if (fitno==1) mixfit <- suppressWarnings(negbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=NULL))
     else mixfit <- suppressWarnings(negbinom.fit.robmixglm(data$X, data$Y, 
                                       offset = data$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=object$verbose, starting.values=starting.values))
    lp <- data$X %*% matrix(negbinom.mle2$coefficients,ncol=1) + data$offset
    negbinom.loglik <- sum(dnbinom(thedata$Y,mu=lp,size=thetheta,log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-negbinom.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  #stop("Not implemented yet")
  # fit glm to data
  thedata <- data.frame(object$model$X,object$model$Y)
  
  if (dim(object$model$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$model$X)[2]],"+",sep='',collapse="")
  else theformula <- "1+"
  
  if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$model$offset)
  
  theformula <- substr(theformula,1,nchar(theformula)-1)
  
  if(!is.null(object$model$offset)) donegbin <- paste("negbinom.mle <- glm.nb(object.model.Y~",theformula,"+offset(offset),",
                                                     "data=thedata)",sep='',collapse='')
  else donegbin <- paste("negbinom.mle <- glm.nb(object.model.Y~",theformula,",",
                        "data=thedata)",sep='',collapse='')
  
  eval(parse(text=donegbin))
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=negbin.mle)
  fitno <- 1
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  negbinom.boot <- boot(thedata, negbinom.fun, R, sim="parametric",
                        ran.gen = negbinom.rg, mle=negbinom.mle)
  if (showProgress) close(pb)
  return(list(pos=sum(negbinom.boot$t[,1]>negbinom.boot$t0[1]),R=R,nullstat=negbinom.boot$t[,1]))
}


gamma.outlierTest.robmixglm <- function(object, R, showProgress) {
  
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
    if (showProgress) {
      pinterval <- R %/% 40
      if (pinterval < 1)  pinterval <- 1
      if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
    }
    if (fitno==1)  mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=NULL))
     else mixfit <- suppressWarnings(gamma.fit.robmixglm(thedata$X, thedata$Y, 
                                    offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                    notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                    verbose=FALSE, starting.values=starting.values))
    lp <- thedata$X %*% matrix(gamma.mle2$coefficients,ncol=1) + thedata$offset
    phi <- sum(gamma.mle2$residuals^2)/gamma.mle2$df.residual
    gamma.loglik <- sum(dgamma(thedata$Y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-gamma.loglik)
     fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  gamma.mle <- glm.fit(object$model$X, object$model$Y, 
                       offset = object$model$offset, family = Gamma(link="log"))
  
  
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=gamma.mle)
  fitno <- 1
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  gamma.boot <- boot(thedata, gamma.fun, R, sim="parametric",
                        ran.gen = gamma.rg, mle=gamma.mle)
  if (showProgress) close(pb)
  return(list(pos=sum(gamma.boot$t[,1]>gamma.boot$t0[1]),R=R,nullstat=gamma.boot$t[,1]))
}

truncpoisson.outlierTest.robmixglm <- function(object, R, showProgress) {
  
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
    
    if(!is.null(object$model$offset)) dotrunc <- paste("truncpoisson.mle2 <- vglm(data.Y~",theformula,",",
                                                       "family=pospoisson, offset=offset,",
                                                       "data=thedata)",sep='',collapse='')
    else dotrunc <- paste("truncpoisson.mle2 <- vglm(data.Y~",theformula,",",
                          "family=pospoisson,",
                          "data=thedata)",sep='',collapse='')

    eval(parse(text=dotrunc))
      
    # starting values assume 50% outliers and tau^2 of 1
    starting.values <- c(coef(truncpoisson.mle2)[1:length(coef(truncpoisson.mle2))],log(0.5/(1-0.5)),1)
    if (showProgress) {
      pinterval <- R %% 100
      if (pinterval < 1)  pinterval <- 1
      if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
    }
    if (fitno==1) mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=NULL))
    else mixfit <- suppressWarnings(truncpoisson.fit.robmixglm(data$X, data$Y, 
                                           offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                           notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                           verbose=object$verbose, starting.values=starting.values))
    lp <- data$X %*% matrix(truncpoisson.mle2@coefficients,ncol=1) + data$offset
    truncpoisson.loglik <- sum(dztpois(data$Y,exp(lp),log=TRUE))
    sim.chisq <- 2*(mixfit$logLik-truncpoisson.loglik)
    fitno <<- fitno+1
    return(sim.chisq)
  }
  
  # fit glm to data
  thedata <- data.frame(object$model$X,object$model$Y)

  if (dim(object$model$X)[2]>2) theformula <- paste(names(thedata)[2:dim(object$model$X)[2]],"+",sep='',collapse="")
  else theformula <- "1+"
  
  if(!is.null(offset)) thedata <- data.frame(thedata,offset=object$model$offset)
  
  theformula <- substr(theformula,1,nchar(theformula)-1)
  
  if(!is.null(object$model$offset)) dotrunc <- paste("truncpoisson.mle2 <- vglm(object.model.Y~",theformula,",",
                                        "family=pospoisson, offset=offset,",
                                        "data=thedata)",sep='',collapse='')
  else dotrunc <- paste("truncpoisson.mle2 <- vglm(object.model.Y~",theformula,",",
                        "family=pospoisson,",
                        "data=thedata)",sep='',collapse='')
  
  eval(parse(text=dotrunc))
  fitno <- 1
  
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=truncpoisson.mle2)
  
  truncpoisson.boot <- boot(thedata, truncpoisson.fun, R, sim="parametric",
                            ran.gen = truncpoisson.rg, mle=truncpoisson.mle2)
  if (showProgress) close(pb)
  return(list(pos=sum(truncpoisson.boot$t[,1]>truncpoisson.boot$t0[1]),R=R,nullstat=truncpoisson.boot$t[,1]))
}



poisson.outlierTest.robmixglm <- function(object, R, showProgress) {

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
      if (showProgress) {
        pinterval <- R %/% 40
        if (pinterval < 1)  pinterval <- 1
        if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
      }
      if (fitno==1) mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=TRUE,
                                        verbose=object$verbose, starting.values=NULL))
      else mixfit <- suppressWarnings(poisson.fit.robmixglm(thedata$X, thedata$Y, 
                                        offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                        notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                        verbose=object$verbose, starting.values=starting.values))
      lp <- thedata$X %*% matrix(poisson.mle2$coefficients,ncol=1) + thedata$offset
      poisson.loglik <- sum(dpois(thedata$Y,exp(lp),log=TRUE))
      sim.chisq <- 2*(mixfit$logLik-poisson.loglik)
    	fitno <<- fitno+1
      return(sim.chisq)
    }
    
    # fit glm to data
  poisson.mle <- glm.fit(object$model$X, object$model$Y, 
                         offset = object$model$offset, family = poisson())
  
 
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=poisson.mle)
  fitno <- 1
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  poisson.boot <- boot(thedata, poisson.fun, R, sim="parametric",
                       ran.gen = poisson.rg, mle=poisson.mle)
  if (showProgress) close(pb)
  return(list(pos=sum(poisson.boot$t[,1]>poisson.boot$t0[1]),R=R,nullstat=poisson.boot$t[,1]))
}

binomial.outlierTest.robmixglm <- function(object, R, showProgress) {
  
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
     if (showProgress) {
       pinterval <- R %/% 40
       if (pinterval < 1)  pinterval <- 1
       if (((fitno %% pinterval)==0) | (fitno==R)) setTxtProgressBar(pb, fitno)
     }
    starting.values <- c(coef(binomial.mle2),log(0.5/(1-0.5)),1)
    if (fitno==1) mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=NULL))
    else mixfit <- suppressWarnings(binomial.fit.robmixglm(thedata$X, thedata$Y,
                                       offset = thedata$offset,gh = norm.gauss.hermite(object$quadpoints),
                                       notrials=object$notrials, EMTol = object$EMTol, calcHessian=FALSE,
                                       verbose=FALSE, starting.values=starting.values))
    lp <- thedata$X %*% matrix(binomial.mle2$coefficients,ncol=1) + thedata$offset
    binomial.loglik <- sum(dbinom(thedata$Y[,1],thedata$Y[,1]+thedata$Y[,2],1.0/(1.0+exp(-lp)),log=TRUE))
     sim.chisq <- 2*(mixfit$logLik-binomial.loglik)
     fitno <<- fitno+1
     return(sim.chisq)
  }
  
  # fit glm to data
  binomial.mle <- glm.fit(object$model$X, object$model$Y,
                          offset = object$model$offset, family = binomial())
  
  
  thedata <- list(X=object$model$X,Y=object$model$Y,offset=object$model$offset,mle=binomial.mle)
  fitno <- 1
  if (showProgress) pb <- txtProgressBar(width=40,min=0,max=R,style=3)
  binomial.boot <- boot(thedata, binomial.fun, R, sim="parametric",
                       ran.gen = binomial.rg, mle=binomial.mle)
  if (showProgress) close(pb)
  return(list(pos=sum(binomial.boot$t[,1]>binomial.boot$t0[1]),R=R,nullstat=binomial.boot$t[,1]))
}

outlierTest <-
  ## Short form for generic function
  function(object, R = 999, showProgress = TRUE)
    UseMethod("outlierTest")

outlierTest.robmixglm <- function(object, R = 999, showProgress = TRUE) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  out <- switch (
    object$family,
    gaussian = gaussian.outlierTest.robmixglm(object, R, showProgress),
    binomial = binomial.outlierTest.robmixglm(object, R, showProgress),
    poisson = poisson.outlierTest.robmixglm(object, R, showProgress),
    gamma = gamma.outlierTest.robmixglm(object, R, showProgress),
    truncpoisson = truncpoisson.outlierTest.robmixglm(object, R, showProgress),
    negbinom = negbinom.outlierTest.robmixglm(object, R, showProgress)
  )
  class(out) <- "outlierTest"
  return(out)
}
