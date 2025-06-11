#' Convert Cohen's d to Partial Correlation
#' @param d Cohen's d effect size
#' @param apply_discount Apply framework discount factor
#' @return Partial correlation equivalent
cohens_d_to_partial_r <- function(d, apply_discount = FALSE) {
  if (!is.numeric(d) || any(is.na(d))) {
    stop("Cohen's d must be numeric without NA values")
  }
  if (any(abs(d) > 5)) {
    warning("Very large |d| > 5 detected - results may be unreliable")
  }

  d <- apply_discount_factor(d, apply_discount = apply_discount)

  # Standard conversion formula: r = d / sqrt(d^2 + 4)
  r <- d / sqrt(d^2 + 4)
  return(validate_partial_r(r))
}
