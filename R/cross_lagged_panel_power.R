# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Cross-Lagged Panel Power Analysis with Framework Integration
#' @param r_partial Cross-lagged partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_waves Number of measurement waves (default = 3)
#' @param stability_coefficient Autoregressive stability coefficient (default = 0.6)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
cross_lagged_panel_power <- function(r_partial = NULL, n = NULL, power = 0.8,
                                     n_waves = 3, stability_coefficient = 0.6,
                                     alpha = 0.05, discount_factor = 0.75,
                                     effect_input = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = TRUE)
  }

  # Use missing() to detect actually provided parameters
  provided_params <- c(
    r_partial = !missing(r_partial) && !is.null(r_partial),
    n = !missing(n) && !is.null(n),
    power = !missing(power) && !is.null(power)
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  # Validate inputs using framework functions
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for cross-lagged panel power analysis")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < 50)) {
      stop("Sample size must be whole number >= 50")
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(n_waves) || n_waves < 2 || n_waves != round(n_waves)) {
    stop("Number of waves must be integer >= 2")
  }

  if (!is.numeric(stability_coefficient) || abs(stability_coefficient) >= 1) {
    stop("Stability coefficient must be between -1 and 1")
  }

  # Calculate effective sample size accounting for longitudinal complexity
  if (!is.null(n)) {
    n_effective <- calculate_cross_lagged_effective_n(n, stability_coefficient, n_waves)
  } else {
    n_effective <- NULL
  }

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- cross_lagged_power_calculation(r_partial, n_effective, stability_coefficient, alpha)

    result <- list(
      analysis_type = "cross_lagged_panel_power",
      method = "Cross-Lagged Panel Power Analysis",
      r_partial = r_partial,
      n = n,
      n_effective = n_effective,
      power = calculated_power,
      n_waves = n_waves,
      stability_coefficient = stability_coefficient,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    n_effective_needed <- cross_lagged_sample_size_calculation(r_partial, power, stability_coefficient, alpha)
    calculated_n <- calculate_cross_lagged_required_n(n_effective_needed, stability_coefficient, n_waves)

    result <- list(
      analysis_type = "cross_lagged_panel_sample_size",
      method = "Cross-Lagged Panel Sample Size Analysis",
      r_partial = r_partial,
      n = calculated_n,
      n_effective = n_effective_needed,
      power = power,
      n_waves = n_waves,
      stability_coefficient = stability_coefficient,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- cross_lagged_effect_size_calculation(n_effective, power, stability_coefficient, alpha)

    result <- list(
      analysis_type = "cross_lagged_panel_effect_size",
      method = "Cross-Lagged Panel Effect Size Analysis",
      r_partial = calculated_r_partial,
      n = n,
      n_effective = n_effective,
      power = power,
      n_waves = n_waves,
      stability_coefficient = stability_coefficient,
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

  # Add cross-lagged specific information
  result$longitudinal_details <- list(
    stability_impact = round(1 - stability_coefficient^2, 3),
    attrition_assumption = 0.15,
    model_complexity = "cross_lagged_panel"
  )

  class(result) <- "cross_lagged_panel_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Cross-Lagged Panel
#' @param r_partial Cross-lagged partial correlation
#' @param n_effective Effective sample size
#' @param stability_coefficient Stability coefficient
#' @param alpha Significance level
#' @return Statistical power
#' @export
cross_lagged_power_calculation <- function(r_partial, n_effective, stability_coefficient, alpha) {
  # Adjust correlation for stability effects
  r_adjusted <- r_partial / sqrt(1 - stability_coefficient^2)

  # Degrees of freedom for cross-lagged path
  df <- n_effective - 3

  if (df <= 0) {
    stop("Insufficient degrees of freedom for analysis")
  }

  # Calculate t-statistic and power
  t_value <- r_adjusted * sqrt(df / (1 - r_adjusted^2))
  t_crit <- qt(1 - alpha/2, df)

  # Non-centrality parameter
  ncp <- abs(t_value)

  # Calculate power using non-central t distribution
  power <- 1 - pt(t_crit, df, ncp = ncp) + pt(-t_crit, df, ncp = ncp)

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for Cross-Lagged Panel
#' @param r_partial Cross-lagged partial correlation
#' @param power Target power
#' @param stability_coefficient Stability coefficient
#' @param alpha Significance level
#' @return Required effective sample size
#' @export
cross_lagged_sample_size_calculation <- function(r_partial, power, stability_coefficient, alpha) {
  # Use iterative approach
  n_min <- 50
  n_max <- 5000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- cross_lagged_power_calculation(r_partial, n_test, stability_coefficient, alpha)

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

#' Calculate Effect Size for Cross-Lagged Panel
#' @param n_effective Effective sample size
#' @param power Target power
#' @param stability_coefficient Stability coefficient
#' @param alpha Significance level
#' @return Required cross-lagged partial correlation
#' @export
cross_lagged_effect_size_calculation <- function(n_effective, power, stability_coefficient, alpha) {
  # Use iterative approach
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- cross_lagged_power_calculation(r_test, n_effective, stability_coefficient, alpha)

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

#' Calculate Effective Sample Size for Cross-Lagged Analysis
#' @param n Nominal sample size
#' @param stability_coefficient Stability coefficient
#' @param n_waves Number of waves
#' @return Effective sample size
#' @export
calculate_cross_lagged_effective_n <- function(n, stability_coefficient, n_waves) {
  # Reductions due to:
  # 1. Partial correlation complexity (controlling for stability)
  # 2. Longitudinal attrition (15% assumed)
  # 3. Model complexity with multiple waves

  stability_penalty <- 1 - (stability_coefficient^2 * 0.3)
  attrition_penalty <- 0.85^(n_waves - 1)
  complexity_penalty <- 0.9

  n_effective <- n * stability_penalty * attrition_penalty * complexity_penalty
  return(max(n_effective, n * 0.5))
}

#' Calculate Required Nominal Sample Size
#' @param n_effective_needed Required effective sample size
#' @param stability_coefficient Stability coefficient
#' @param n_waves Number of waves
#' @return Required nominal sample size
#' @export
calculate_cross_lagged_required_n <- function(n_effective_needed, stability_coefficient, n_waves) {
  # Reverse the effective sample size calculation
  stability_penalty <- 1 - (stability_coefficient^2 * 0.3)
  attrition_penalty <- 0.85^(n_waves - 1)
  complexity_penalty <- 0.9

  total_penalty <- stability_penalty * attrition_penalty * complexity_penalty
  required_n <- ceiling(n_effective_needed / max(total_penalty, 0.5))

  return(required_n)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated Cross-Lagged Panel Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param n_waves Number of waves
#' @param stability_coefficient Stability coefficient
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
#' @export
cross_lagged_panel_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                               power = 0.8, n_waves = 3,
                                               stability_coefficient = 0.6,
                                               alpha = 0.05, discount_factor = 0.75) {
  cross_lagged_panel_power(effect_input = effect_size, effect_type = effect_type,
                           n = n, power = power, n_waves = n_waves,
                           stability_coefficient = stability_coefficient,
                           alpha = alpha, discount_factor = discount_factor)
}

#' Quick Cross-Lagged Panel Sample Size Calculation
#' @param r_partial Cross-lagged partial correlation
#' @param power Target power (default = 0.8)
#' @param n_waves Number of waves (default = 3)
#' @param stability_coefficient Stability coefficient (default = 0.6)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
#' @export
cross_lagged_panel_sample_size <- function(r_partial, power = 0.8, n_waves = 3,
                                           stability_coefficient = 0.6, alpha = 0.05) {
  result <- cross_lagged_panel_power(r_partial = r_partial, power = power,
                                     n_waves = n_waves, stability_coefficient = stability_coefficient,
                                     alpha = alpha)
  return(result$n)
}

#' Quick Cross-Lagged Panel Power Calculation
#' @param r_partial Cross-lagged partial correlation
#' @param n Sample size
#' @param n_waves Number of waves (default = 3)
#' @param stability_coefficient Stability coefficient (default = 0.6)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
cross_lagged_panel_power_check <- function(r_partial, n, n_waves = 3,
                                           stability_coefficient = 0.6, alpha = 0.05) {
  result <- cross_lagged_panel_power(r_partial = r_partial, n = n,
                                     n_waves = n_waves, stability_coefficient = stability_coefficient,
                                     alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Cross-Lagged Panel Power Analysis
#' @param x Cross-lagged panel power analysis result
#' @param ... Additional arguments
#' @export
print.cross_lagged_panel_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Cross-lagged partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Effective sample size:", round(x$n_effective), "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of waves:", x$n_waves, "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Model details
  cat("\nLongitudinal details:\n")
  cat("  Stability coefficient:", round(x$stability_coefficient, 3), "\n")
  cat("  Stability impact:", x$longitudinal_details$stability_impact, "\n")
  cat("  Assumed attrition:", round(x$longitudinal_details$attrition_assumption * 100), "%\n")

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
