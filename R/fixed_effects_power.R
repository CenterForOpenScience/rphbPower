# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Fixed Effects Power Analysis with Framework Integration
#' @param r_partial Within-person partial correlation coefficient (NULL to calculate)
#' @param n_units Number of individuals/units (NULL to calculate)
#' @param n_periods Number of time periods per unit (default = 4)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param icc Intraclass correlation coefficient (default = 0.3)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
fixed_effects_power <- function(r_partial = NULL, n_units = NULL, n_periods = 4,
                                power = 0.8, icc = 0.3, alpha = 0.05,
                                discount_factor = 0.75, effect_input = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = TRUE)
  }

  # Use missing() to detect actually provided parameters
  provided_params <- c(
    r_partial = !missing(r_partial) && !is.null(r_partial),
    n_units = !missing(n_units) && !is.null(n_units),
    power = !missing(power) && !is.null(power)
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n_units, power (or use effect_input)")
  }

  # Validate inputs using framework functions
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for fixed effects power analysis")
  }

  if (!is.null(n_units)) {
    if (!is.numeric(n_units) || any(n_units != round(n_units)) || any(n_units < 10)) {
      stop("Number of units must be whole number >= 10")
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(n_periods) || n_periods < 2 || n_periods != round(n_periods)) {
    stop("Number of periods must be integer >= 2")
  }

  if (!is.numeric(icc) || icc < 0 || icc >= 1) {
    stop("ICC must be between 0 and 1")
  }

  # Calculate design effect and effective sample size
  design_effect <- calculate_fixed_effects_design_effect(n_periods, icc)

  if (!is.null(n_units)) {
    n_effective <- calculate_fixed_effects_effective_n(n_units, n_periods, icc)
  } else {
    n_effective <- NULL
  }

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- fixed_effects_power_calculation(r_partial, n_effective, icc, alpha)

    result <- list(
      analysis_type = "fixed_effects_power",
      method = "Fixed Effects Power Analysis",
      r_partial = r_partial,
      n_units = n_units,
      n_periods = n_periods,
      n_total_observations = n_units * n_periods,
      n_effective = n_effective,
      power = calculated_power,
      icc = icc,
      design_effect = design_effect,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else if (is.null(n_units)) {
    # Calculate sample size
    n_effective_needed <- fixed_effects_sample_size_calculation(r_partial, power, icc, alpha)
    calculated_n_units <- calculate_fixed_effects_required_n(n_effective_needed, n_periods, icc)

    result <- list(
      analysis_type = "fixed_effects_sample_size",
      method = "Fixed Effects Sample Size Analysis",
      r_partial = r_partial,
      n_units = calculated_n_units,
      n_periods = n_periods,
      n_total_observations = calculated_n_units * n_periods,
      n_effective = n_effective_needed,
      power = power,
      icc = icc,
      design_effect = design_effect,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- fixed_effects_effect_size_calculation(n_effective, power, icc, alpha)

    result <- list(
      analysis_type = "fixed_effects_effect_size",
      method = "Fixed Effects Effect Size Analysis",
      r_partial = calculated_r_partial,
      n_units = n_units,
      n_periods = n_periods,
      n_total_observations = n_units * n_periods,
      n_effective = n_effective,
      power = power,
      icc = icc,
      design_effect = design_effect,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "effect_size"
    )
  }

  # Add framework integration
  r_for_conversion <- result$r_partial
  result$effect_size_conversions <- framework_conversion_summary(
    r_for_conversion, "r", apply_discount = FALSE
  )
  result$interpretation <- interpret_effect_size(r_for_conversion)

  # Add fixed effects specific information
  result$panel_details <- list(
    within_person_focus = TRUE,
    time_invariant_controls = "all_unobserved",
    efficiency_vs_cross_sectional = round(1 / design_effect, 2)
  )

  class(result) <- "fixed_effects_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Fixed Effects
#' @param r_partial Within-person partial correlation
#' @param n_effective Effective sample size
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @return Statistical power
#' @export
fixed_effects_power_calculation <- function(r_partial, n_effective, icc, alpha) {
  # Degrees of freedom for within-person effect
  df <- n_effective - 2

  if (df <= 0) {
    stop("Insufficient degrees of freedom for analysis")
  }

  # Calculate t-statistic and power
  t_value <- r_partial * sqrt(df / (1 - r_partial^2))
  t_crit <- qt(1 - alpha/2, df)

  # Non-centrality parameter
  ncp <- abs(t_value)

  # Calculate power using non-central t distribution
  power <- 1 - pt(t_crit, df, ncp = ncp) + pt(-t_crit, df, ncp = ncp)

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for Fixed Effects
#' @param r_partial Within-person partial correlation
#' @param power Target power
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @return Required effective sample size
#' @export
fixed_effects_sample_size_calculation <- function(r_partial, power, icc, alpha) {
  # Use iterative approach
  n_min <- 20
  n_max <- 5000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- fixed_effects_power_calculation(r_partial, n_test, icc, alpha)

    if (abs(power_test - power) < tolerance) {
      return(n_test)
    }

    if (power_test < power) {
      n_min <- n_test
    } else {
      n_max <- n_test
    }

    if (n_max - n_min <= 1) break
  }

  return(n_max)
}

#' Calculate Effect Size for Fixed Effects
#' @param n_effective Effective sample size
#' @param power Target power
#' @param icc Intraclass correlation coefficient
#' @param alpha Significance level
#' @return Required within-person partial correlation
#' @export
fixed_effects_effect_size_calculation <- function(n_effective, power, icc, alpha) {
  # Use iterative approach
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- fixed_effects_power_calculation(r_test, n_effective, icc, alpha)

    if (abs(power_test - power) < tolerance) {
      return(r_test)
    }

    if (power_test < power) {
      r_min <- r_test
    } else {
      r_max <- r_test
    }

    if (r_max - r_min <= 0.001) break
  }

  return(r_max)
}

#' Calculate Design Effect for Fixed Effects
#' @param n_periods Number of time periods
#' @param icc Intraclass correlation coefficient
#' @return Design effect multiplier
#' @export
calculate_fixed_effects_design_effect <- function(n_periods, icc) {
  # Standard design effect formula for clustered data
  design_effect <- 1 + (n_periods - 1) * icc
  return(design_effect)
}

#' Calculate Effective Sample Size for Fixed Effects
#' @param n_units Number of units
#' @param n_periods Number of periods
#' @param icc Intraclass correlation coefficient
#' @return Effective sample size
#' @export
calculate_fixed_effects_effective_n <- function(n_units, n_periods, icc) {
  # Total observations
  n_total <- n_units * n_periods

  # Adjust for clustering
  design_effect <- calculate_fixed_effects_design_effect(n_periods, icc)
  n_effective <- n_total / design_effect

  # Additional penalty for fixed effects complexity
  complexity_penalty <- 0.85
  n_effective <- n_effective * complexity_penalty

  return(max(n_effective, n_units * 0.5))
}

#' Calculate Required Number of Units
#' @param n_effective_needed Required effective sample size
#' @param n_periods Number of periods
#' @param icc Intraclass correlation coefficient
#' @return Required number of units
#' @export
calculate_fixed_effects_required_n <- function(n_effective_needed, n_periods, icc) {
  # Reverse the effective sample size calculation
  design_effect <- calculate_fixed_effects_design_effect(n_periods, icc)
  complexity_penalty <- 0.85

  # Calculate required total observations
  n_total_needed <- n_effective_needed * design_effect / complexity_penalty

  # Convert to number of units
  n_units_needed <- ceiling(n_total_needed / n_periods)

  return(n_units_needed)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

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
                                          n_periods = 4, power = 0.8, icc = 0.3,
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

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Fixed Effects Power Analysis
#' @param x Fixed effects power analysis result
#' @param ... Additional arguments
print.fixed_effects_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Within-person partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Number of units:", x$n_units, "\n")
  cat("Number of periods:", x$n_periods, "\n")
  cat("Total observations:", x$n_total_observations, "\n")
  cat("Effective sample size:", round(x$n_effective), "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Panel details
  cat("\nPanel details:\n")
  cat("  ICC:", round(x$icc, 3), "\n")
  cat("  Design effect:", round(x$design_effect, 2), "\n")
  cat("  Efficiency vs cross-sectional:", x$panel_details$efficiency_vs_cross_sectional, "x\n")

  # Framework conversions
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat(paste0("  Cohen's f", cli::symbol$sup_2, ":"), round(conv$cohens_f2, 3), "\n")
    cat(paste0("  R", cli::symbol$sup_2, ":"), round(conv$r_squared, 3), "\n")
  }

  # Analysis details
  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")

  cat("\n")
}
