binomial.fit.robmixglm <- function(x,y,offset,gh,notrials,EMTol, calcHessian=TRUE, cores, verbose,   starting.values=NULL) {
  
#  browser()
  
  if (!is.null(starting.values)) {
    notrials <- 1
    cores <- 1
  }
  
  logit <- function(x){
    return(1.0/(1.0+exp(-x)))
  }
  
  llrandbinom <- function(y, lp, tau2, gh){

    dointegrate <- function(oney,onelp,tau2) {

      thefun <- function(nu){
        return(dbinom(oney[1],oney[1]+oney[2],logit(onelp+nu*sqrt(tau2))))
      }
      theint <- sum(thefun(gh[,1])*gh[,2])
      if (is.nan(theint)) theint <- NA
      return(log(theint))
    }

    calconell <- function(x) {
      oney <- x[1:2]
      onelp <- x[3]
      if (tau2<0.0) return(NA)
      if (tau2==0) thell <- dbinom(oney[1],oney[1]+oney[2],logit(onelp),log=TRUE) else
        thell <- dointegrate(oney=oney,onelp=onelp,tau2=tau2)
      return(thell)
    }
    if (is.na(tau2)) stop("tau2 is NA")
    ll <- apply(cbind(y,lp),1,calconell)
    return(ll)
  }
  

  thestatistic <- function(data, ...) {
    
    rlreg <- function(xcoef,lpoutlier,tau2,x,y,offset,prop) {
      poutlier <- 1.0/(1+exp(-lpoutlier))

      lp <- as.vector(x %*% xcoef)+offset

      lp1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)
      lp2 <- llrandbinomcpp(y, lp, tau2, gh)

      if (!missing(prop)) {
        ll <- prop*cbind(lp1,lp2)
        negll <- -sum(apply(ll,1,sum))
      } else {
        l <- exp(cbind(lp1+log(1-poutlier),lp2+log(poutlier)))
        negll <- -sum(log(apply(l,1,sum)))
      }
      if (is.nan(negll)) negll <- 1e108
      if (!is.finite(negll)) negll <- 1e108
      return(negll)
    }

    optimrlreg <- function(p,lpoutlier,x,y,offset,prop) {
      return(rlreg(matrix(p[1:(length(p)-1)],ncol=1),lpoutlier,p[length(p)],x,y,offset,prop))
    }

    
    y <- data[,1:2]
    x <- data[,4:dim(data)[2]]
    if (dim(data)[2] >=4) offset <- data[,3]
    else offset <- NULL
    
    
   tryCatch({
     if (is.null(starting.values)) {
       noutliers <- max(1,round(dim(data)[1]*0.2))
       outliers <- sample(c(rep(1,noutliers),rep(0,dim(data)[1]-noutliers)))
       
       robust.binomial.prefit <- glm(y~-1+x,
                                    family=binomial, offset=offset,subset=(outliers!=1))

      prefit.coef <- coef(robust.binomial.prefit)
      # assume 20% outliers as a starting point
      currlpoutlier <- log(0.2/(1-0.2))
      currxcoef <- matrix(prefit.coef[1:(length(prefit.coef))],ncol=1)
      currxcoef <- ifelse(is.na(currxcoef),0,currxcoef)
      currtau2 <- min(rgamma(1,1.5,1),2)
    } else {
      currxcoef <- matrix(starting.values[1:(length(starting.values)-2)],ncol=1)
      currlpoutlier <- starting.values[length(starting.values)-1]
      currtau2 <- starting.values[length(starting.values)]
    }

    currll <- -1.0e100
    nem <- 0

    repeat {
      nem <- nem+1
      # expectation step
      currpoutlier <- 1.0/(1+exp(-currlpoutlier))

      lp <- as.vector(x %*% currxcoef)+offset

      ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-currpoutlier)
      ll2 <- llrandbinomcpp(y, lp, currtau2, gh)+log(currpoutlier)

      ll <- cbind(ll1,ll2)
      prop <- t(apply(ll,1,function(x) {
        x <- x-max(x)
        x <- ifelse(x==-Inf,-1e100,x)
        return(exp(x)/sum(exp(x)))
      }))
      # calculate outlier proportion
      poutlier <- sum(prop[,2])/dim(prop)[1]
      currlpoutlier <- log(poutlier/(1-poutlier))

      if (is.na(poutlier)) stop("problem with poutlier")

      startvals <- c(currxcoef,currtau2)
      names(startvals) <- c(dimnames(x)[[2]],"tau2")

      results.nlm <- suppressWarnings(nlminb(startvals,optimrlreg,
                                             lower=c(rep(-Inf,length(currxcoef)),0),
                                             #control=list(trace=1,iter.max=5),
                                             control=list(iter.max=5),
                                             lpoutlier=currlpoutlier,
                                             prop=prop,
                                             y=y,x=x,offset=offset))
      currxcoef <- matrix(as.numeric(results.nlm$par)[1:(length(results.nlm$par)-1)],ncol=1)
      currtau2 <- as.numeric(results.nlm$par)[length(results.nlm$par)]

      lastll <- currll
      currll <- -rlreg(currxcoef,currlpoutlier,currtau2,x,y,offset)
      if (verbose) print(sprintf("Likelihood at end of EM step %.4f", currll))
      if (abs((lastll-currll)/currll)<EMTol) break()
      if (nem >100) break()
    }
  return(c(currll,currxcoef,currlpoutlier,currtau2))
   },
  error=function(e) return(rep(NA,5))
   )
  }
  
    ll.robustbinomial <- function(p){

      xcoef <- p[1:(length(p)-2)]
      lpoutlier <- p[length(p)-1]
      tau2 <- p[length(p)]

      poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))

      lp <- as.vector(x %*% xcoef)+offset

      ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-poutlier)
      ll2 <- llrandbinomcpp(y, lp, tau2, gh)+log(poutlier)

      l <- exp(cbind(ll1,ll2))
      negll <- -sum(log(apply(l,1,sum)))
      if (is.nan(negll)) negll <- NA
      if (!is.finite(negll)) negll <- NA      #print(negll)
      return(negll)
    }

    if (is.null(starting.values)) {
      bootdata <- cbind(y, offset, x)

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)


    if (cores>1) {
      if(.Platform$OS.type=="unix")  parallel <- "multicore"
      else parallel <- "snow" }
    else parallel <- "no"
    #browser()
    if (notrials == 1) {
      theboot <- matrix(thestatistic(bootdata), nrow=1)
    } else 
      theboot <- boot(bootdata, thestatistic, R=notrials, sim = "parametric",
                      ran.gen = function(d, p) d,
                      parallel = parallel,
                      ncpus = cores)$t
 
    res <- vector("list", length = notrials)
    for (i in 1:notrials) {
      res[[i]] <- list(logLik=theboot[i,1],start.val=theboot[i,2:dim(theboot)[2]])
    }

    maxll <- -Inf
    nfails <- 0

    for (i in 1:notrials) {
      if (verbose) cat(c(res[[i]]$logLik,res[[i]]$start.val),"\n")
      if (is.na(res[[i]]$logLik)) nfails <- nfails+1
      else {
        if (res[[i]]$logLik>maxll) {
          maxll <- res[[i]]$logLik
          start.val <- res[[i]]$start.val
        }
      }
      if (nfails > 0) warning(sprintf("Failed to obtain starting values for %i starting sets", nfails))
    }
    } else {
      start.val <- starting.values
    }
    
    if(is.null(start.val)) stop("Cannot find valid starting values")

    thenames <- c(dimnames(x)[[2]],"lpoutlier","tau2")

    names(start.val) <- thenames

    parnames(ll.robustbinomial) <- names(start.val)

    lower.val <- c(rep(-Inf,length(start.val)-2),-Inf,0)

    names(lower.val) <- names(start.val)

    if(verbose) thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000,trace=5)
    else thecontrol <- list(eval.max=1000,iter.max=1000,eval.max=2000)

    ll.robustbinomial.fit <- mle2(ll.robustbinomial,start=start.val,vecpar=TRUE,
                                  optimizer="user",optimfun=myoptim,
                                  data=list(y=y,x=x,offset=offset),
                                  skip.hessian=!calcHessian,trace=verbose,
                                  lower=lower.val,
                                  control=thecontrol)

  if (calcHessian) {
    thecoef <- coef(ll.robustbinomial.fit)
    ncoef <- length(thecoef)
    if (thecoef[ncoef]<0.0001) thecoef[ncoef] <- 0.0001
    ll.robustbinomial.fit@details$hessian <- optimHess(thecoef,ll.robustbinomial,control=list(ndeps=c(rep(1.0e-5,length(thecoef)))))
    ll.robustbinomial.fit@vcov <- ginv(ll.robustbinomial.fit@details$hessian)
  }
 
  xcoef <- matrix(coef(ll.robustbinomial.fit)[1:(length(coef(ll.robustbinomial.fit))-2)],ncol=1)
  lpoutlier <- coef(ll.robustbinomial.fit)[length(coef(ll.robustbinomial.fit))-1]
  poutlier <- 1.0/(1+exp(-lpoutlier))
  tau2 <- coef(ll.robustbinomial.fit)[length(coef(ll.robustbinomial.fit))]


  lp <- as.vector(x %*% xcoef)+offset

  ll1 <- dbinom(y[,1], y[,1]+y[,2], logit(lp), log = TRUE)+log(1-poutlier)
  ll2 <- llrandbinomcpp(y, lp, tau2, gh)+log(poutlier)

  ll <- cbind(ll1,ll2)
  prop <- t(apply(ll,1,function(x) {
    x <- x-max(x)
    x <- ifelse(x==-Inf,-1e100,x)
    return(exp(x)/sum(exp(x)))
  }))

  coef.names <- c(dimnames(x)[[2]],"Outlier p.","Tau-sq")
  return(list(fit=ll.robustbinomial.fit,prop=prop,logLik=-ll.robustbinomial.fit@min,np=length(coef.names),nobs=dim(x)[1],coef.names=coef.names))
}
  
  