negll <- function(p, y, x, offset) {
  lp <- x %*% matrix(p[1:(length(p)-1)],ncol=1)+offset
  theta <- p[length(p)]
  negll <- -sum(dnbinom(y,mu=exp(lp),size=theta,log=TRUE))
  if (negll>1.0e100) negll <- 1.0e100
  if (!is.finite(negll)) stop("Error in calculation log likelihood")
  return(negll)
} 

fitnegbin <- function(y,x,offset) {
  # obtain starting values
  pois.glm <- glm(y~x[,-1], family=poisson, offset=offset)
  initp <- c(coef(pois.glm),1)
  fit <- nlminb(initp, negll, lower = c(rep(-Inf,length(initp)-1),0),   x=x, y=y, offset=offset)
# allow for singular convergence
  if (fit$convergence>1) stop("Convergence failed")
  return(fit)
}
