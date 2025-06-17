# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Fixed Effects Power Analysis with Framework Integration (v2.2)
#' @param r_partial Within-person partial correlation coefficient (NULL to calculate)
#' @param n_units Number of individuals/units (NULL to calculate)
#' @param n_periods Number of time periods per unit (default = 4)
#' @param power Statistical power (NULL to calculate, default applied when needed = 0.8)
#' @param icc Intraclass correlation coefficient (default = 0.3)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
fixed_effects_power <- function(r_partial = NULL, n_units = NULL, n_periods = 4,
                                power = NULL, icc = 0.3, alpha = 0.05,
                                discount_factor = 0.75, effect_input = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = TRUE)
  }

  if (is.null(power)) {
    power_default_used <- TRUE
    power <- 0.8
  } else {
    power_default_used <- FALSE
  }

  provided_params <- c(
    r_partial = !missing(r_partial) && !is.null(r_partial),
    n_units = !missing(n_units) && !is.null(n_units),
    power = !missing(power) && !is.null(power) && !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n_units, power (or use effect_input)")
  }

  # Validate inputs
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for fixed effects power analysis")
  }
  if (!is.null(n_units) && (n_units < 3 || n_units != round(n_units))) {
    stop("Number of units must be a whole number >= 3.")
  }
  if (!is.null(power) && (power <= 0 || power >= 1)) {
    stop("Power must be between 0 and 1.")
  }
  if (!is.numeric(n_periods) || n_periods < 2 || n_periods != round(n_periods)) {
    stop("Number of periods must be an integer >= 2.")
  }
  if (!is.numeric(icc) || icc < 0 || icc >= 1) {
    stop("ICC must be between 0 and 1.")
  }

  # Perform calculation
  if (!provided_params["power"]) {
    calculated_power <- fixed_effects_power_calculation(r_partial, n_units, n_periods, icc, alpha)
    result <- list(power = calculated_power, calculation_target = "power")
  } else if (!provided_params["n_units"]) {
    calculated_n_units <- fixed_effects_sample_size_calculation(r_partial, power, n_periods, icc, alpha)
    result <- list(n_units = calculated_n_units, calculation_target = "sample_size")
  } else {
    calculated_r_partial <- fixed_effects_effect_size_calculation(n_units, power, n_periods, icc, alpha)
    result <- list(r_partial = calculated_r_partial, calculation_target = "effect_size")
  }

  # Assemble final output object
  final_result <- list(
    analysis_type = paste("fixed_effects", result$calculation_target, sep="_"),
    method = "Fixed Effects Power Analysis (v2.2)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n_units = if(is.null(n_units)) result$n_units else n_units,
    n_periods = n_periods,
    power = if(is.null(power) || power_default_used) result$power else power,
    icc = icc,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = result$calculation_target
  )

  # Add derived values for reporting
  final_result$n_total_observations <- final_result$n_units * final_result$n_periods
  final_result$design_effect <- calculate_fixed_effects_design_effect(final_result$n_periods, final_result$icc)
  final_result$n_effective <- calculate_fixed_effects_effective_n(final_result$n_units, final_result$n_periods, final_result$icc)
  final_result$effect_size_conversions <- framework_conversion_summary(final_result$r_partial, "r", apply_discount = FALSE)
  final_result$interpretation <- interpret_effect_size(final_result$r_partial)
  final_result$panel_details <- list(
    within_person_focus = TRUE,
    time_invariant_controls = "all_unobserved",
    efficiency_vs_cross_sectional = round(1 / final_result$design_effect, 2)
  )

  class(final_result) <- "fixed_effects_power_analysis"
  return(final_result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS - CRITICALLY CORRECTED (v2.2)
# ==============================================================================

#' Calculate Power for Fixed Effects - Corrected (v2.2)
#' @param r_partial Within-person partial correlation
#' @param n_units Number of units
#' @param n_periods Number of periods
#' @param icc Intraclass correlation coefficient (used for reporting, not this calculation)
#' @param alpha Significance level
#' @return Statistical power
#' @export
fixed_effects_power_calculation <- function(r_partial, n_units, n_periods, icc, alpha) {
  # The test of a within-person (time-varying) predictor in a FE model is a t-test.
  # Degrees of freedom are total observations minus number of units (for the fixed effects)
  # minus the number of other predictors. Here, k=1.
  df_den <- (n_units * n_periods) - n_units - 1
  if (df_den <= 0) return(0)

  # The t-statistic for a partial correlation r_partial with the correct df_den.
  t_stat <- r_partial * sqrt(df_den) / sqrt(1 - r_partial^2)
  ncp <- t_stat

  # Critical t-value for a two-tailed test
  t_crit <- stats::qt(1 - alpha / 2, df_den)
  power <- stats::pt(t_crit, df_den, ncp = ncp, lower.tail = FALSE) + stats::pt(-t_crit, df_den, ncp = ncp, lower.tail = TRUE)

  return(pmax(0.001, pmin(0.999, power)))
}

#' Calculate Sample Size for Fixed Effects
#' @param r_partial Within-person partial correlation
#' @param power Target power
#' @param n_periods Number of time periods
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @return Required effective sample size
#' @export
fixed_effects_sample_size_calculation <- function(r_partial, power, n_periods, icc, alpha) {
  power_func <- function(n_units) {
    if ((n_units * n_periods) - n_units - 1 <= 0) return(-1)
    fixed_effects_power_calculation(r_partial, n_units, n_periods, icc, alpha) - power
  }
  min_n <- ceiling(1 + 1/n_periods) + 1 # Min n_units to have df_den > 0
  result <- tryCatch(stats::uniroot(power_func, interval = c(min_n, 100000))$root, error = function(e) NA)
  if(is.na(result)) stop("Could not find sample size.")
  return(ceiling(result))
}

#' Calculate Effect Size for Fixed Effects
#' @param n_units Number of units
#' @param power Target power
#' @param n_periods Number of time periods
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @return Required within-person partial correlation
#' @export
fixed_effects_effect_size_calculation <- function(n_units, power, n_periods, icc, alpha) {
  power_func <- function(r) fixed_effects_power_calculation(r, n_units, n_periods, icc, alpha) - power
  if ((n_units * n_periods) - n_units - 1 <= 0) stop("Insufficient N to calculate effect size.")
  result <- tryCatch(stats::uniroot(power_func, interval = c(0.001, 0.999))$root, error = function(e) NA)
  if(is.na(result)) stop("Could not find effect size.")
  return(result)
}

# ==============================================================================
# HELPER, CONVENIENCE, AND PRINT FUNCTIONS
# ==============================================================================

#' Calculate Design Effect for Fixed Effects (REPORTING ONLY)
#' @param n_periods Number of time periods
#' @param icc Intraclass correlation coefficient
#' @return Design effect multiplier
#' @export
calculate_fixed_effects_design_effect <- function(n_periods, icc) {
  design_effect <- 1 + (n_periods - 1) * icc
  return(design_effect)
}

#' Calculate Effective Sample Size for Fixed Effects (REPORTING ONLY)
#' @param n_units Number of units
#' @param n_periods Number of periods
#' @param icc Intraclass correlation coefficient
#' @return Effective sample size
#' @export
calculate_fixed_effects_effective_n <- function(n_units, n_periods, icc) {
  n_effective <- (n_units * n_periods) / (1 + (n_periods - 1) * icc)
  return(round(n_effective))
}

#' Framework-Integrated Fixed Effects Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n_units Number of units
#' @param n_periods Number of periods
#' @param power Target power
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
#' @export
fixed_effects_framework_power <- function(effect_size, effect_type = "r", n_units = NULL,
                                          n_periods = 4, power = NULL, icc = 0.3,
                                          alpha = 0.05, discount_factor = 0.75) {
  fixed_effects_power(effect_input = effect_size, effect_type = effect_type,
                      n_units = n_units, n_periods = n_periods, power = power,
                      icc = icc, alpha = alpha, discount_factor = discount_factor)
}

#' Quick Fixed Effects Sample Size Calculation
#' @param r_partial Within-person partial correlation
#' @param power Target power (default = 0.8)
#' @param n_periods Number of periods (default = 4)
#' @param icc Intraclass correlation coefficient (default = 0.3)
#' @param alpha Significance level (default = 0.05)
#' @return Required number of units
#' @export
fixed_effects_sample_size <- function(r_partial, power = 0.8, n_periods = 4,
                                      icc = 0.3, alpha = 0.05) {
  result <- fixed_effects_power(r_partial = r_partial, power = power,
                                n_periods = n_periods, icc = icc, alpha = alpha)
  return(result$n_units)
}

#' Quick Fixed Effects Power Calculation
#' @param r_partial Within-person partial correlation
#' @param n_units Number of units
#' @param n_periods Number of periods (default = 4)
#' @param icc Intraclass correlation coefficient (default = 0.3)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
fixed_effects_power_check <- function(r_partial, n_units, n_periods = 4,
                                      icc = 0.3, alpha = 0.05) {
  result <- fixed_effects_power(r_partial = r_partial, n_units = n_units,
                                n_periods = n_periods, icc = icc, alpha = alpha)
  return(result$power)
}

#' Print Method for Fixed Effects Power Analysis
#' @param x Fixed effects power analysis result
#' @param ... Additional arguments
#' @export
print.fixed_effects_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  cat("Within-person partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Number of units:", x$n_units, "\n")
  cat("Number of periods:", x$n_periods, "\n")
  cat("Total observations:", x$n_total_observations, "\n")
  cat("Effective sample size (reporting):", round(x$n_effective), "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")

  cat("\nPanel details:\n")
  cat("  ICC:", round(x$icc, 3), "\n")
  cat("  Design effect:", round(x$design_effect, 2), "\n")
  cat("  Efficiency vs cross-sectional:", x$panel_details$efficiency_vs_cross_sectional, "x\n")

  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat("  Cohen's f\u00B2:", round(conv$cohens_f2, 3), "\n")
    cat("  R\u00B2:", round(conv$r_squared, 3), "\n")
  }

  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("  Note: Uses correct DF and NCP for within-person effects.\n")

  cat("\n")
}
