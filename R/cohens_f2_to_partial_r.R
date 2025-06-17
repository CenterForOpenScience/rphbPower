#' Convert Cohen's f² to Partial Correlation
#' @param f2 Cohen's f\00B2 effect size
#' @param apply_discount Apply framework discount factor
#' @return Partial correlation equivalent
#' @export
cohens_f2_to_partial_r <- function(f2, apply_discount = FALSE) {
  if (!is.numeric(f2) || any(is.na(f2))) {
    stop("Cohen's f\u00B2 must be numeric without NA values")
  }
  if (any(f2 < 0)) {
    stop("Cohen's f\u00B2 must be non-negative")
  }
  if (any(f2 > 10)) {
    warning("Very large f\u00B2 > 10 detected - results may be unreliable")
  }

  f2 <- apply_discount_factor(f2, apply_discount = apply_discount)

  # Conversion formula: r = sqrt(f² / (1 + f²))
  r <- sqrt(f2 / (1 + f2))
  return(validate_partial_r(r))
}
