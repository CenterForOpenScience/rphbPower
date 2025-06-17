#' Convert Partial Correlation to Cohen's f²
#' @param r Partial correlation
#' @return Cohen's f\00B2 equivalent
#' @export
partial_r_to_cohens_f2 <- function(r) {
  r <- validate_partial_r(r, allow_zero = FALSE, context = "for Cohen's f\u00B2 conversion")

  # Conversion formula: f² = r² / (1 - r²)
  f2 <- r^2 / (1 - r^2)
  return(f2)
}
