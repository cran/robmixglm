myoptim <- function(par, fn, gr,method,
                          lower, upper,
                          control, hessian,...)
{
  myfn <- function(par) {
    temp <- fn(par)
    if (!is.finite(temp)) temp <- 1.0e100
    return(temp)
  }
  if (is.null(control$trace)) verbose <- FALSE
  else verbose <- control$trace
  if (verbose) cat("using nlminb","\n")
  if (length(lower) < length(par)) lower <- rep(lower,length(par))
  if (length(upper) < length(par)) upper <- rep(upper,length(par))
  thenlm <- suppressWarnings(nlminb(par, fn,lower = lower, upper = upper, control=control))
  if (thenlm$convergence==0) theoptim <- list(par=thenlm$par,value=thenlm$objective,convergence=thenlm$convergence,message=thenlm$message)
  else {
    if (verbose) cat("using optim","\n")
    control <- list(trace=if(verbose) 1 else 0,maxit=control$iter.max)
    theoptim <- suppressWarnings(optim(par, myfn, method="L-BFGS-B",lower=lower,upper=upper,hessian=FALSE,control=control))
  }
  if (theoptim$convergence!=0)  {
    if (verbose) cat("using constrOptim","\n")
    ui <- diag(length(par))[is.finite(lower),]
    ci <- lower[is.finite(lower)]
    par[is.finite(lower)] <- ifelse(par[is.finite(lower)]<1.0e-4,1.0e-4,par[is.finite(lower)])
    theoptim <- constrOptim(par, fn, NULL , ui, ci, control=control)
  }
  # convert any warning to failure
  if (theoptim$convergence!=0) stop(theoptim$message)
  return(theoptim)
}