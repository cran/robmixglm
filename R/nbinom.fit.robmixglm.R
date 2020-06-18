nbinom.fit.robmixglm <- function(x,y,offset,gh,notrials,EMTol,  calcHessian=TRUE, cores, verbose,  starting.values=NULL) {
  
  if (!is.null(starting.values)) notrials <- 1

  fitonemlreg <- function(y,outliers,x=NULL,offset=NULL,fixed) {
    
    rlreg <- function(xcoef,lpoutlier,tau2,theta,x,y,offset,prop) {
      poutlier <- 1.0/(1+exp(-lpoutlier))
      
      lp <- as.vector(x %*% xcoef)
      
      if(!is.null(offset)) lp <- lp+offset

      lp1 <-  dnbinom(y, mu=exp(lp), size=theta, log=TRUE)
      lp2 <- llrandnbinomcpp(y, lp, tau2, theta, gh)
      
      if (!missing(prop)) {
        ll <- prop*cbind(lp1,lp2)
        negll <- -sum(apply(ll,1,sum))
      } else {
        l <- exp(cbind(lp1+log(1-poutlier),lp2+log(poutlier)))
        negll <- -sum(log(apply(l,1,sum)))
      }
      if (is.nan(negll)) negll <- NA
      if (!is.finite(negll)) negll <- NA
      return(negll)
    }
    
    optimrlreg <- function(p,lpoutlier,x,y,offset,prop,fixed) {
      p[names(fixed)] <- fixed
      return(rlreg(matrix(p[1:(length(p)-2)],ncol=1),lpoutlier,p[length(p)-1],p[length(p)],x,y,offset,prop))
    }
    
    thedata <- data.frame(y,x)
    
    lastdata=length(names(thedata))
    if(!is.null(offset)) thedata <- data.frame(thedata,offset=offset)
    
    # fit first with glm
    if (is.null(starting.values)) {
      
      if (lastdata>2) theformula <- paste(names(thedata)[3:lastdata],"+",sep='',collapse="")
      else theformula <- "1+"
      
      theformula <- substr(theformula,1,nchar(theformula)-1)
      
      
      doprefit <- paste("suppressWarnings(robust.nbinom.prefit <- glm.nb(y~",theformula,"+offset(offset),",
                        "subset=(outliers!=1),data=thedata))",sep='',collapse='')
      
      eval(parse(text=doprefit))
      
      # get the starting values
      prefit.coef <- coef(robust.nbinom.prefit)
      # assume 20% outliers as a starting point
      currlpoutlier <- log(0.2/(1-0.2))
      currxcoef <- matrix(prefit.coef[1:(length(prefit.coef))],ncol=1)
      currxcoef <- ifelse(is.na(currxcoef),0,currxcoef)
      currtau2 <- min(rgamma(1,2,1),5)
      currtheta <- robust.nbinom.prefit$theta
    } else
    {
      currxcoef <- matrix(starting.values[1:(length(starting.values)-3)],ncol=1)
      currlpoutlier <- starting.values[length(starting.values)-2]
      currtau2 <- starting.values[length(starting.values)-1]
      currtheta <- starting.values[length(starting.values)]
    }
    
    currll <- -1.0e100
    nem <- 0
    
    repeat {
      nem <- nem+1
      # expectation step
      currpoutlier <- 1.0/(1+exp(-currlpoutlier))
      
      lp <- as.vector(x %*% currxcoef)
      
      if (!is.null(offset)) lp <- lp+offset
      
      ll1 <- dnbinom(y, mu=exp(lp), size=currtheta, log=TRUE)+log(1-currpoutlier) 
      ll2 <- llrandnbinomcpp(y, lp, currtau2, currtheta, gh)+log(currpoutlier)
      
      ll <- cbind(ll1,ll2)
      prop <- t(apply(ll,1,function(x) {
        x <- x-max(x)
        x <- ifelse(x==-Inf,-1e100,x)
        return(exp(x)/sum(exp(x)))
      }))
      # calculate outlier proportion
      poutlier <- sum(prop[,2])/dim(prop)[1]
      currlpoutlier <- log(poutlier/(1-poutlier))
      
      if (is.na(poutlier)) stop()
      
      startvals <- c(currxcoef,currtau2,currtheta)
      names(startvals) <- c(dimnames(x)[[2]],"tau2","theta")
      
      #browser()
      
      results.nlm <- suppressWarnings(nlminb(startvals,optimrlreg,
                            lower=c(rep(-Inf,length(currxcoef)),0,0),
                            #control=list(trace=1,iter.max=10),
                            control=list(iter.max=5),
                            lpoutlier=currlpoutlier,
                            prop=prop,
                            y=y,x=x,offset=offset,
                            fixed=fixed))
      
      #browser()
      
      currxcoef <- matrix(as.numeric(results.nlm$par)[1:(length(results.nlm$par)-2)],ncol=1)
      currtau2 <- as.numeric(results.nlm$par)[length(results.nlm$par)-1]
      currtheta <- as.numeric(results.nlm$par)[length(results.nlm$par)]
      
      lastll <- currll
      currll <- -rlreg(currxcoef,currlpoutlier,currtau2,currtheta,x,y,offset)
      if (abs((lastll-currll)/currll)<EMTol) break()
      if (nem >100) break()
    }
    return(list(ll=currll,start.val=c(currxcoef,currlpoutlier,currtau2,currtheta)))    
  }
  
  ll.robustnbinom <- function(p){
    
    xcoef <- p[1:(length(p)-3)]
    lpoutlier <- p[length(p)-2]
    tau2 <- p[length(p)-1]
    theta <- p[length(p)]
    
    poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
    
    lp <- as.vector(x %*% xcoef)
    
    if (!is.null(offset)) lp <- lp+offset
    
    ll1 <- dnbinom(y, mu=exp(lp), size=theta, log=TRUE)+log(1-poutlier)
    ll2 <- llrandnbinomcpp(y, lp, tau2, theta, gh)+log(poutlier)
    
    l <- exp(cbind(ll1,ll2))
    negll <- -sum(log(apply(l,1,sum)))
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    return(negll)
  }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  if (is.null(starting.values)) {
    if (cores > 1) {
      cl = parallel::makeCluster(cores, setup_strategy = "sequential")
      doParallel::registerDoParallel(cl)
      res = foreach(i = 1:notrials, 
                    .options.RNG=seed[1]) %dorng% {
                      noutliers <- max(1,round(dim(x)[1]*0.2))
                      outliers <- sample(c(rep(1,noutliers),rep(0,dim(x)[1]-noutliers)),dim(x)[1])
                      fitonemlreg(y,outliers,x,offset,fixed=NULL)}
      parallel::stopCluster(cl)
      maxll <- -Inf
      for (i in 1:notrials) {
        if (verbose) print(c(res[[i]]$ll,res[[i]]$start.val))
        if (res[[i]]$ll>maxll) {
          maxll <- res[[i]]$ll
          start.val <- res[[i]]$start.val
        }
      }
    } else {
      maxll <- -Inf
      for (i in 1:notrials) {
        try({
          noutliers <- max(1,round(dim(x)[1]*0.2))
          outliers <- sample(c(rep(1,noutliers),rep(0,dim(x)[1]-noutliers)),dim(x)[1])
          thefit <- fitonemlreg(y,outliers,x,offset,fixed=NULL)
          if (verbose) print(c(thefit$ll,thefit$start.val))
          if (thefit$ll>maxll) {
            maxll <- thefit$ll
            start.val <- thefit$start.val
          }
        }
        )
      }
    }
  } else {
    start.val <- starting.values
  }
  
  thenames <- c(dimnames(x)[[2]],"lpoutlier","tau2", "theta")
  
  names(start.val) <- thenames
  
  parnames(ll.robustnbinom) <- names(start.val)
  
  lower.val <- c(rep(-Inf,length(start.val)-3),-Inf,0,0)
  
  names(lower.val) <- names(start.val)
  
  
  robustnbinom.fit <- mle2(ll.robustnbinom,start=start.val,vecpar=TRUE,
                          optimizer="user",optimfun=myoptim,
                          data=list(y=y,x=x,offset=offset),
                          skip.hessian=TRUE,trace=verbose,
                          lower=lower.val,
                          control=if (verbose) list(eval.max=1000,iter.max=1000,trace=5) else list(eval.max=1000,iter.max=1000))
  
  if (calcHessian) {
    thecoef <- coef(robustnbinom.fit)
    ncoef <- length(thecoef)
    if (thecoef[ncoef-1]<0.0001) thecoef[ncoef-1] <- 0.0001
    robustnbinom.fit@details$hessian <- hessian(ll.robustnbinom,thecoef)
    robustnbinom.fit@vcov <- ginv(robustnbinom.fit@details$hessian)
  }
  
  xcoef <- matrix(coef(robustnbinom.fit)[1:(length(coef(robustnbinom.fit))-3)],ncol=1)
  lpoutlier <- coef(robustnbinom.fit)[length(coef(robustnbinom.fit))-2]
  poutlier <- 1.0/(1+exp(-lpoutlier))
  tau2 <- coef(robustnbinom.fit)[length(coef(robustnbinom.fit))-1]
  theta <- coef(robustnbinom.fit)[length(coef(robustnbinom.fit))]
  
  
  lp <- as.vector(x %*% xcoef)
  
  if (!is.null(offset)) lp <- lp+offset
  
  ll1 <- dnbinom(y, mu=exp(lp), size=theta, log=TRUE)+log(1-poutlier)  
  ll2 <- llrandnbinomcpp(y, lp, tau2, theta, gh)+log(poutlier)
  
  ll <- cbind(ll1,ll2)
  prop <- t(apply(ll,1,function(x) {
    x <- x-max(x)
    x <- ifelse(x==-Inf,-1e100,x)
    return(exp(x)/sum(exp(x)))
  }))
  
  coef.names <- c(dimnames(x)[[2]],"Outlier p.","Tau-sq","theta")
  return(list(fit=robustnbinom.fit,prop=prop,logLik=-robustnbinom.fit@min,np=length(coef.names),nobs=dim(x)[1],coef.names=coef.names))  
}

