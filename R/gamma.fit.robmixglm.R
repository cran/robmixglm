gamma.fit.robmixglm <- function(x,y,offset,gh,notrials,EMTol,  calcHessian=TRUE, verbose, starting.values=NULL) {

    if (!is.null(starting.values)) notrials <- 1

  fitonemlreg <- function(y,outliers,x=NULL,offset=NULL,fixed) {
    
    rlreg <- function(xcoef,lpoutlier,tau2,phi,x,y,offset,prop) {
      poutlier <- 1.0/(1+exp(-lpoutlier))
      
      lp <- as.vector(x %*% xcoef)
      
      if(!is.null(offset)) lp <- lp+offset
      
      lp1 <-  suppressWarnings(dgamma(y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE))
      lp2 <- llrandgammacpp(y, lp, tau2, phi, gh)
       
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
    
    
    doprefit <- paste("robust.gamma.prefit <- glm(y~",theformula,",",
           "family=Gamma(link='log'), offset=offset,subset=(outliers!=1),data=thedata)",sep='',collapse='')
    
    eval(parse(text=doprefit))
    
    # get the starting values
    prefit.coef <- coef(robust.gamma.prefit)
    # assume 50% outliers as a starting point
    currlpoutlier <- log(0.5/(1-0.5))
    currxcoef <- matrix(prefit.coef[1:(length(prefit.coef))],ncol=1)
    currxcoef <- ifelse(is.na(currxcoef),0,currxcoef)
    currtau2 <- min(rgamma(1,2,1),3)
    currphi <- summary(robust.gamma.prefit)$dispersion
    } else
    {
      currxcoef <- matrix(starting.values[1:(length(starting.values)-3)],ncol=1)
      currlpoutlier <- starting.values[length(starting.values)-2]
      currtau2 <- starting.values[length(starting.values)-1]
      currphi <- starting.values[length(starting.values)]
    }
        
    currll <- -1.0e100
    nem <- 0
   
    repeat {
      nem <- nem+1
      # expectation step
      currpoutlier <- 1.0/(1+exp(-currlpoutlier))
      
      lp <- as.vector(x %*% currxcoef)
      
      if (!is.null(offset)) lp <- lp+offset
      
      ll1 <- dgamma(y, shape=1.0/currphi,rate=1.0/(currphi*exp(lp)), log = TRUE)+log(1-currpoutlier) 
      ll2 <- llrandgammacpp(y, lp, currtau2, currphi, gh)+log(currpoutlier)
      
      l <- exp(cbind(ll1,ll2))
      prop <- l/apply(l,1,sum)
      # calculate outlier proportion
      poutlier <- sum(prop[,2])/dim(prop)[1]
      currlpoutlier <- log(poutlier/(1-poutlier))
      
      
      if (is.na(poutlier)) stop("Problem with poutlier")
      
      startvals <- c(currxcoef,currtau2,currphi)
      names(startvals) <- c(dimnames(x)[[2]],"tau2","phi")
      
       results.nlm <- suppressWarnings(nlminb(startvals,optimrlreg,
                            lower=c(rep(-Inf,length(currxcoef)),0,0),
                            #control=list(trace=1,iter.max=10),
                            control=list(iter.max=5),
                            lpoutlier=currlpoutlier,
                            prop=prop,
                            y=y,x=x,offset=offset,
                            fixed=fixed))
      currxcoef <- matrix(as.numeric(results.nlm$par)[1:(length(results.nlm$par)-2)],ncol=1)
      currtau2 <- as.numeric(results.nlm$par)[length(results.nlm$par)-1]
      currphi <- as.numeric(results.nlm$par)[length(results.nlm$par)]
      
      lastll <- currll
      currll <- -rlreg(currxcoef,currlpoutlier,currtau2,currphi,x,y,offset)
      if (abs((lastll-currll)/currll)<EMTol) break()
      if (nem >100) break()
    }
    return(list(ll=currll,start.val=c(currxcoef,currlpoutlier,currtau2,currphi)))    
  }
  
  ll.robustgamma <- function(p){

    xcoef <- p[1:(length(p)-3)]
    lpoutlier <- p[length(p)-2]
    tau2 <- p[length(p)-1]
    phi <- p[length(p)]
    
    poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
    
    lp <- as.vector(x %*% xcoef)
    
    if (!is.null(offset)) lp <- lp+offset
    
    ll1 <- suppressWarnings(dgamma(y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE)+log(1-poutlier))
    
    ll2 <- llrandgammacpp(y, lp, tau2, phi, gh)+log(poutlier)
    
    l <- exp(cbind(ll1,ll2))
    negll <- -sum(log(apply(l,1,sum)))
    if (is.nan(negll)) negll <- NA
    if (!is.finite(negll)) negll <- NA
    return(negll)
  }
  
  if (is.null(starting.values)) {
    maxll <- -Inf
    
    for (i in 1:notrials) {
      try({
        outliers <- rbinom(dim(x)[1],1,0.5)
      thefit <- suppressWarnings(fitonemlreg(y,outliers,x,offset,fixed=NULL))
      if (verbose) print(c(thefit$ll,thefit$start.val))
      if (thefit$ll>maxll) {
        maxll <- thefit$ll
        start.val <- thefit$start.val
      }}
      )
      if (!is.finite(maxll)) stop("No starting values found") 
    }
  } else {
    thefit <- suppressWarnings(fitonemlreg(y,NULL,x,offset,fixed=NULL))
    start.val <- thefit$start.val
  }
  
  thenames <- c(dimnames(x)[[2]],"lpoutlier","tau2", "phi")
  
  names(start.val) <- thenames
  
  parnames(ll.robustgamma) <- names(start.val)
  
  lower.val <- c(rep(-Inf,length(start.val)-3),-Inf,0,0)
  
  names(lower.val) <- names(start.val)

  robustgamma.fit <- mle2(ll.robustgamma,start=start.val,vecpar=TRUE,
                          optimizer="user",optimfun=myoptim,
                          data=list(y=y,x=x,offset=offset),
                            skip.hessian=TRUE,trace=verbose,
                            lower=lower.val)
  
  if (calcHessian) {
    thecoef <- coef(robustgamma.fit)
    ncoef <- length(thecoef)
    if (thecoef[ncoef-1]<0.0001) thecoef[ncoef-1] <- 0.0001
    robustgamma.fit@details$hessian <- hessian(ll.robustgamma,thecoef)
    robustgamma.fit@vcov <- ginv(robustgamma.fit@details$hessian)
  }
  #if (any(is.nan(sqrt(diag(robustgamma.fit@vcov)))))  warning("Error in calculating standard errors.")
    
  xcoef <- matrix(coef(robustgamma.fit)[1:(length(coef(robustgamma.fit))-3)],ncol=1)
  lpoutlier <- coef(robustgamma.fit)[length(coef(robustgamma.fit))-2]
  poutlier <- 1.0/(1+exp(-lpoutlier))
  tau2 <- coef(robustgamma.fit)[length(coef(robustgamma.fit))-1]
  phi <- coef(robustgamma.fit)[length(coef(robustgamma.fit))]
  
  lp <- as.vector(x %*% xcoef)
  
  if (!is.null(offset)) lp <- lp+offset
  
  ll1 <- dgamma(y,shape=1.0/phi,rate=1.0/(phi*exp(lp)),log=TRUE)+log(1-poutlier)  
	ll2 <- llrandgammacpp(y, lp, tau2, phi, gh)+log(poutlier)
  
  l <- exp(cbind(ll1,ll2))
  prop <- l/apply(l,1,sum)
  
  coef.names <- c(dimnames(x)[[2]],"Outlier p.","Tau-sq","Phi")
  return(list(fit=robustgamma.fit,prop=prop,logLik=-robustgamma.fit@min,np=length(coef.names),nobs=dim(x)[1],coef.names=coef.names))  
}
