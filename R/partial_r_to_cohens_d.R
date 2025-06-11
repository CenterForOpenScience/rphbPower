#' Convert Partial Correlation to Cohen's d
#' @param r Partial correlation
#' @return Cohen's d equivalent
partial_r_to_cohens_d <- function(r) {
  r <- validate_partial_r(r, allow_zero = FALSE, context = "for Cohen's d conversion")

  # Standard conversion formula: d = 2r / sqrt(1 - r^2)
  d <- 2 * r / sqrt(1 - r^2)
  return(d)
}
