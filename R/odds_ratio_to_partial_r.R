#' Convert Odds Ratio to Partial Correlation
#' @param or The odds ratio from a logistic regression.
#' @param apply_discount Apply framework discount factor (default: FALSE).
#' @return Partial correlation equivalent (r).
#' @export
odds_ratio_to_partial_r <- function(or, apply_discount = FALSE) {
  if (!is.numeric(or) || any(is.na(or)) || any(or <= 0)) {
    stop("Odds ratio (or) must be a positive numeric value.")
  }
  
  or <- apply_discount_factor(or, apply_discount = apply_discount)
  
  # Convert odds ratio to log-odds (beta), then to Cohen's d
  beta <- log(or)
  d <- beta * sqrt(3) / pi
  
  # Convert Cohen's d to partial r using the existing framework function
  r <- cohens_d_to_partial_r(d, apply_discount = FALSE) # Discount already applied
  return(validate_partial_r(r))
}