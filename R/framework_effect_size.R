#' Framework-Integrated Effect Size Input
#' @param effect_size Input effect size
#' @param effect_type Type of effect size ("r", "d", "f2", "r_squared", "eta_squared")
#' @param apply_discount Whether to apply 0.75 discount factor
#' @return Partial correlation with discount applied
#' @export
framework_effect_size <- function(effect_size, effect_type = "r", apply_discount = TRUE) {
  # Apply discount first, then convert
  discounted_effect <- apply_discount_factor(effect_size, apply_discount = apply_discount)

  switch(effect_type,
         "r" = validate_partial_r(discounted_effect),
         "d" = cohens_d_to_partial_r(discounted_effect, apply_discount = FALSE),
         "f2" = cohens_f2_to_partial_r(discounted_effect, apply_discount = FALSE),
         "r_squared" = r_squared_to_partial_r(discounted_effect, apply_discount = FALSE),
         "eta_squared" = eta_squared_to_partial_r(discounted_effect, apply_discount = FALSE),
         "or" = odds_ratio_to_partial_r(discounted_effect, apply_discount = FALSE),
         stop("Unsupported effect_type. Use: r, d, f2, r_squared, eta_squared")
  )
}
