#' Validate Partial Correlation
#' @param r Partial correlation value
#' @param allow_zero Whether to allow exactly zero (default: TRUE)
#' @param context Context for error messages (optional)
#' @return Validated r value
validate_partial_r <- function(r, allow_zero = TRUE, context = "") {
  if (!is.numeric(r) || any(is.na(r))) {
    stop(paste("Partial correlation must be numeric without NA", context))
  }
  if (any(abs(r) >= 1)) {
    stop(paste("Partial correlation must be between -1 and 1 (exclusive)", context))
  }
  if (!allow_zero && any(r == 0)) {
    stop(paste("Partial correlation cannot be zero for this calculation", context))
  }
  return(r)
}
