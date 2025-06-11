#' Convert Eta-squared to Partial Correlation
#' @param eta_squared Eta-squared effect size
#' @param apply_discount Apply framework discount factor
#' @return Partial correlation equivalent
eta_squared_to_partial_r <- function(eta_squared, apply_discount = FALSE) {
  if (!is.numeric(eta_squared) || any(is.na(eta_squared))) {
    stop("Eta-squared must be numeric without NA values")
  }
  if (any(eta_squared < 0) || any(eta_squared >= 1)) {
    stop("Eta-squared must be between 0 and 1 (exclusive of 1)")
  }

  eta_squared <- apply_discount_factor(eta_squared, apply_discount = apply_discount)

  # Square root conversion (same as R²)
  r <- sqrt(eta_squared)
  return(validate_partial_r(r))
}
