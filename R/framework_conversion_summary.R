#' Framework Effect Size Conversion Summary
#' @param effect_value Input effect size value
#' @param input_type Type of input effect size
#' @param apply_discount Apply framework discount factor
#' @return Named list with all conversions
framework_conversion_summary <- function(effect_value, input_type = "r", apply_discount = TRUE) {
  # Convert to partial correlation first
  r_partial <- framework_effect_size(effect_value, input_type, apply_discount)

  list(
    input_effect = effect_value,
    input_type = input_type,
    discount_applied = apply_discount,
    partial_r = r_partial,
    cohens_d = partial_r_to_cohens_d(r_partial),
    cohens_f2 = partial_r_to_cohens_f2(r_partial),
    r_squared = r_partial^2,
    eta_squared = r_partial^2,
    interpretation = interpret_effect_size(r_partial)
  )
}
