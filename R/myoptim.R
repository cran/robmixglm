myoptim <- function(par, fn, gr,method,
                          lower, upper,
                          control, hessian,...)
{
  myfn <- function(par) {
    temp <- fn(par)
    if (!is.finite(temp)) temp <- 1.0e100
    return(temp)
  }
  if (length(lower) < length(par)) lower <- rep(lower,length(par))
  if (length(upper) < length(par)) upper <- rep(upper,length(par))
  thenlm <- suppressWarnings(nlminb(par, fn,lower = lower, upper = upper, control= list(eval.max=1000,iter.max=1000)))
  if (thenlm$convergence==0) theoptim <- list(par=thenlm$par,value=thenlm$objective,convergence=thenlm$convergence,message=thenlm$message)
  else theoptim <- suppressWarnings(optim(par, myfn, method="L-BFGS-B",lower=lower,upper=upper,hessian=FALSE,control=list(maxit=1000)))
  if (theoptim$convergence!=0)  {
    ui <- diag(length(par))[is.finite(lower),]
    ci <- lower[is.finite(lower)]
    par[is.finite(lower)] <- ifelse(par[is.finite(lower)]<1.0e-4,1.0e-4,par[is.finite(lower)])
    theoptim <- constrOptim(par, fn, NULL , ui, ci, control=list(maxit=2000))
  }
  # convert any warning to failure
  if (theoptim$convergence!=0) stop(theoptim$message)
  return(theoptim)
}