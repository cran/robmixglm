outlierProbs <-
  ## Short form for generic function 
  function(object) UseMethod("outlierProbs")

outlierProbs.robmixglm <- function(object) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  
  warning("outlierProbs has been replaced by outlierProbsMix, and will be removed in a future update")
  outliers <- outlierProbsMix.robmixglm(object) 
  return(outliers)
}