#' Effect Size Interpretation (Framework Standard)
#' @param r Partial correlation
#' @param standard Interpretation standard ("cohen", "field")
#' @return Interpretation string
interpret_effect_size <- function(r, standard = "cohen") {
  r <- validate_partial_r(r)
  r_abs <- abs(r)

  if (standard == "cohen") {
    if (r_abs < 0.1) return("Negligible")
    if (r_abs < 0.3) return("Small")
    if (r_abs < 0.5) return("Medium")
    return("Large")
  } else if (standard == "field") {
    if (r_abs < 0.05) return("Negligible")
    if (r_abs < 0.15) return("Small")
    if (r_abs < 0.35) return("Medium")
    if (r_abs < 0.55) return("Large")
    return("Very Large")
  } else {
    stop("Standard must be 'cohen' or 'field'")
  }
}
