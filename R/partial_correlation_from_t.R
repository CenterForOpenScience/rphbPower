#' Estimate Partial Correlation from t-statistic
#' @param t_value t-statistic for coefficient
#' @param df Degrees of freedom
#' @return Partial correlation
#' @export
partial_correlation_from_t <- function(t_value, df) {
  if (!is.numeric(t_value) || !is.numeric(df) || any(is.na(c(t_value, df)))) {
    stop("t-value and df must be numeric without NA values")
  }
  if (any(df <= 0)) {
    stop("Degrees of freedom must be positive")
  }

  # Standard formula: r = t / sqrt(t² + df)
  r_partial <- t_value / sqrt(t_value^2 + df)
  return(validate_partial_r(r_partial))
}
