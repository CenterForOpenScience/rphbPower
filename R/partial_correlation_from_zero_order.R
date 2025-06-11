#' Estimate Partial Correlation from Zero-Order Correlations
#' @param r_xy Zero-order correlation between X and Y
#' @param r_xz Correlation between X and control variable Z
#' @param r_yz Correlation between Y and control variable Z
#' @return Partial correlation of X and Y controlling for Z
#' @export
partial_correlation_from_zero_order <- function(r_xy, r_xz, r_yz) {
  correlations <- c(r_xy, r_xz, r_yz)
  if (!all(is.numeric(correlations)) || any(is.na(correlations))) {
    stop("All correlations must be numeric without NA values")
  }
  if (any(abs(correlations) >= 1)) {
    stop("All correlations must be between -1 and 1 (exclusive)")
  }

  # Standard partial correlation formula
  numerator <- r_xy - (r_xz * r_yz)
  denominator <- sqrt((1 - r_xz^2) * (1 - r_yz^2))

  if (abs(denominator) < 1e-10) {
    stop("Cannot compute: near-perfect correlation with control variable")
  }

  r_partial <- numerator / denominator
  return(validate_partial_r(r_partial))
}
