#' Estimate Partial Correlation from F-statistic
#' @param f_value F-statistic
#' @param df_num Numerator degrees of freedom (default: 1)
#' @param df_den Denominator degrees of freedom
#' @return Partial correlation
partial_correlation_from_f <- function(f_value, df_num = 1, df_den) {
  if (!is.numeric(f_value) || !is.numeric(df_num) || !is.numeric(df_den)) {
    stop("F-statistic parameters must be numeric")
  }
  if (any(f_value < 0)) {
    stop("F-statistic must be non-negative")
  }
  if (any(df_num <= 0) || any(df_den <= 0)) {
    stop("Degrees of freedom must be positive")
  }

  if (df_num == 1) {
    # Convert F to t first, then to r
    t_value <- sqrt(f_value)
    r_partial <- t_value / sqrt(t_value^2 + df_den)
  } else {
    # Multiple df case: convert through f²
    f2 <- (f_value * df_num) / df_den
    r_partial <- sqrt(f2 / (1 + f2))
  }

  return(validate_partial_r(r_partial))
}
