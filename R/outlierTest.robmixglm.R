outlierTest <-
  ## Short form for generic function
  function(object, R = 999, cores = max(detectCores(logical = FALSE)-1, 1))
    UseMethod("outlierTest")

outlierTest.robmixglm <- function(object, R = 999, cores = max(detectCores(logical = FALSE)-1, 1)) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
 
  warning("outlierTest has been replaced by outlierTestMix, and will disappear soon") 
  out <- outlierTestMix.robmixglm(object, R, cores) 
  class(out) <- "outlierTestMix"
  return(out)
}
