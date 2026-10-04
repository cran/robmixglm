outlierProbsMix <-
  ## Short form for generic function 
  function(object) UseMethod("outlierProbsMix")

outlierProbsMix.robmixglm <- function(object) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  outliers <- object$prop[,2]
  class(outliers) <- "outlierProbsMix"
  return(outliers)
}