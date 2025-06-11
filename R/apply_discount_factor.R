#' Apply Framework Discount Factor
#' @param effect_size Input effect size (any metric)
#' @param discount_factor Discount factor (default: 0.75 per framework)
#' @param apply_discount Whether to apply discount (default: TRUE)
#' @return Discounted effect size
#' @export
apply_discount_factor <- function(effect_size, discount_factor = 0.75, apply_discount = TRUE) {
  if (!is.numeric(effect_size) || any(is.na(effect_size))) {
    stop("Effect size must be numeric without NA values")
  }
  if (!is.numeric(discount_factor) || discount_factor <= 0 || discount_factor > 1) {
    stop("Discount factor must be between 0 and 1")
  }

  if (apply_discount) {
    return(effect_size * discount_factor)
  } else {
    return(effect_size)
  }
}
