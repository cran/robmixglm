fitted.robmixglm <- function(object, ...) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  lp <- object$model$X %*% matrix(coef(object),ncol=1) + object$model$offset
  out <- switch (
    object$family,
    gaussian = lp,
    binomial = 1.0/(1.0+exp(-lp)),
    poisson = exp(lp),
    gamma = exp(lp),
    truncpoisson = exp(lp),
    negbinom = exp(lp)
  )
  return(out)
}
