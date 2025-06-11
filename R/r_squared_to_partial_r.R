#' Convert R² to Partial Correlation
#' @param r_squared R\00B2 value from regression
#' @param apply_discount Apply framework discount factor
#' @return Partial correlation equivalent
r_squared_to_partial_r <- function(r_squared, apply_discount = FALSE) {
  if (!is.numeric(r_squared) || any(is.na(r_squared))) {
    stop("R^2 must be numeric without NA values")
  }
  if (any(r_squared < 0) || any(r_squared >= 1)) {
    stop("R^2 must be between 0 and 1 (exclusive of 1)")
  }

  r_squared <- apply_discount_factor(r_squared, apply_discount = apply_discount)

  # Simple square root conversion
  r <- sqrt(r_squared)
  return(validate_partial_r(r))
}
