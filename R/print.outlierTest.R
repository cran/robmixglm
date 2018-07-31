print.outlierTest <- function(x, ...) {
  if (!inherits(x, "outlierTest"))
    stop("Use only with 'outlierTest' objects.\n")
  test <- (x$pos+1)/(x$R+1)
  out <- sprintf("p value %1.4f", test)
  cat(out)
  invisible(out)
}

summary.outlierTest <- function(object,...) {
  if (!inherits(object, "outlierTest"))
    stop("Use only with 'outlierTest' objects.\n")
  class(object) <- "summary.outlierTest"
  return(object)
}

print.summary.outlierTest <- function(x,...) {
  if (!inherits(x, "summary.outlierTest"))
    stop("Use only with 'summary.outlierTest' objects.\n")
  test <- (x$pos+1)/(x$R+1)
  out <- sprintf("p value %1.4f", test)
  cat(out)
  invisible(out)
}
