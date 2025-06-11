#' Unified Effect Size Table for Framework
#' @param effect_values Vector of effect sizes
#' @param input_type Type of input effect sizes
#' @param apply_discount Apply framework discount factor
#' @return Data frame with unified conversions
#' @export
unified_effect_size_table <- function(effect_values, input_type = "r", apply_discount = TRUE) {
  if (length(effect_values) == 0) {
    stop("No effect sizes provided")
  }

  # Convert all to partial correlations
  r_values <- sapply(effect_values, function(x) {
    framework_effect_size(x, input_type, apply_discount)
  })

  data.frame(
    Input = effect_values,
    Input_Type = rep(input_type, length(effect_values)),
    Discounted = rep(apply_discount, length(effect_values)),
    Partial_r = round(r_values, 4),
    Cohens_d = round(sapply(r_values, partial_r_to_cohens_d), 4),
    Cohens_f2 = round(sapply(r_values, partial_r_to_cohens_f2), 4),
    R_squared = round(r_values^2, 4),
    Interpretation = sapply(r_values, interpret_effect_size),
    stringsAsFactors = FALSE
  )
}
