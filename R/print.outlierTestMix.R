print.outlierTestMix <- function(x, ...) {
  if (!inherits(x, "outlierTestMix"))
    stop("Use only with 'outlierTestMix' objects.\n")
  thestat <- x$t0[1]
  thestat <- ifelse(thestat<0.0,0.0,thestat)
  nullstat <- x$t[,1]
  nullstat <- ifelse(nullstat<0.0,0.0,nullstat)
  pos <- sum(nullstat>thestat)
  test <- (pos+1)/(length(x$t[,1])+1)
  out <- sprintf("p value %1.4f", test)
  cat(out)
  invisible(out)
}

summary.outlierTestMix <- function(object,...) {
  if (!inherits(object, "outlierTestMix"))
    stop("Use only with 'outlierTest' objects.\n")
  class(object) <- "summary.outlierTest"
  return(object)
}

print.summary.outlierTestMix <- function(x,...) {
  if (!inherits(x, "summary.outlierTestMix"))
    stop("Use only with 'summary.outlierTest' objects.\n")
  thestat <- x$t0[1]
  thestat <- ifelse(thestat<0.0,0.0,thestat)
  nullstat <- x$t[,1]
  nullstat <- ifelse(nullstat<0.0,0.0,nullstat)
  pos <- sum(nullstat>thestat)
  test <- (pos+1)/(length(x$t[,1])+1)
  out <- sprintf("p value %1.4f", test)
  cat(out)
  invisible(out)
}
