# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Repeated Measures Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_timepoints Number of repeated measures (default = 2)
#' @param correlation_between_measures Correlation between measures (default = 0.5)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
repeated_measures_power <- function(r_partial = NULL, n = NULL, power = 0.8,
                                    n_timepoints = 2, correlation_between_measures = 0.5,
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
                                    context = "for repeated measures power analysis")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < 10)) {
      stop("Sample size must be whole number >= 10")
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(n_timepoints) || n_timepoints < 2 || n_timepoints != round(n_timepoints)) {
    stop("Number of timepoints must be integer >= 2")
  }

  if (!is.numeric(correlation_between_measures) || abs(correlation_between_measures) >= 1) {
    stop("Correlation between measures must be between -1 and 1")
  }

  # Calculate design efficiency
  design_efficiency <- calculate_repeated_measures_efficiency(n_timepoints, correlation_between_measures)

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- repeated_measures_power_calculation(r_partial, n, n_timepoints,
                                                            correlation_between_measures, alpha)

    result <- list(
      analysis_type = "repeated_measures_power",
      method = "Repeated Measures Power Analysis",
      r_partial = r_partial,
      n = n,
      power = calculated_power,
      n_timepoints = n_timepoints,
      correlation_between_measures = correlation_between_measures,
      design_efficiency = design_efficiency,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    calculated_n <- repeated_measures_sample_size_calculation(r_partial, power, n_timepoints,
                                                              correlation_between_measures, alpha)

    result <- list(
      analysis_type = "repeated_measures_sample_size",
      method = "Repeated Measures Sample Size Analysis",
      r_partial = r_partial,
      n = calculated_n,
      power = power,
      n_timepoints = n_timepoints,
      correlation_between_measures = correlation_between_measures,
      design_efficiency = design_efficiency,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- repeated_measures_effect_size_calculation(n, power, n_timepoints,
                                                                      correlation_between_measures, alpha)

    result <- list(
      analysis_type = "repeated_measures_effect_size",
      method = "Repeated Measures Effect Size Analysis",
      r_partial = calculated_r_partial,
      n = n,
      power = power,
      n_timepoints = n_timepoints,
      correlation_between_measures = correlation_between_measures,
      design_efficiency = design_efficiency,
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

  class(result) <- "repeated_measures_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Repeated Measures
#' @importFrom stats qf pf
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param n_timepoints Number of timepoints
#' @param correlation_between_measures Correlation between measures
#' @param alpha Significance level
#' @return Statistical power
repeated_measures_power_calculation <- function(r_partial, n, n_timepoints,
                                                correlation_between_measures, alpha) {
  # Convert to F-statistic approach
  f2 <- partial_r_to_cohens_f2(r_partial)

  # Degrees of freedom
  df_num <- n_timepoints - 1
  df_den <- (n - 1) * (n_timepoints - 1)

  if (df_den <= 0) {
    stop("Insufficient degrees of freedom for analysis")
  }

  # Non-centrality parameter
  ncp <- f2 * df_den

  # Critical F-value
  f_crit <- qf(1 - alpha, df_num, df_den)

  # Calculate power using non-central F distribution
  power <- 1 - pf(f_crit, df_num, df_den, ncp = ncp)

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for Repeated Measures
#' @param r_partial Partial correlation
#' @param power Target power
#' @param n_timepoints Number of timepoints
#' @param correlation_between_measures Correlation between measures
#' @param alpha Significance level
#' @return Required sample size
repeated_measures_sample_size_calculation <- function(r_partial, power, n_timepoints,
                                                      correlation_between_measures, alpha) {
  # Use iterative approach
  n_min <- 10
  n_max <- 5000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- repeated_measures_power_calculation(r_partial, n_test, n_timepoints,
                                                      correlation_between_measures, alpha)

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

#' Calculate Effect Size for Repeated Measures
#' @param n Sample size
#' @param power Target power
#' @param n_timepoints Number of timepoints
#' @param correlation_between_measures Correlation between measures
#' @param alpha Significance level
#' @return Required partial correlation
repeated_measures_effect_size_calculation <- function(n, power, n_timepoints,
                                                      correlation_between_measures, alpha) {
  # Use iterative approach
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- repeated_measures_power_calculation(r_test, n, n_timepoints,
                                                      correlation_between_measures, alpha)

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

#' Calculate Design Efficiency for Repeated Measures
#' @param n_timepoints Number of timepoints
#' @param correlation_between_measures Correlation between measures
#' @return Design efficiency multiplier
calculate_repeated_measures_efficiency <- function(n_timepoints, correlation_between_measures) {
  # Formula based on compound symmetry assumption
  efficiency <- 1 + (n_timepoints - 1) * correlation_between_measures

  # Cap efficiency at reasonable bounds
  efficiency <- pmax(0.5, pmin(efficiency, 5.0))

  return(efficiency)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated Repeated Measures Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param n_timepoints Number of timepoints
#' @param correlation_between_measures Correlation between measures
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
repeated_measures_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                              power = 0.8, n_timepoints = 2,
                                              correlation_between_measures = 0.5,
                                              alpha = 0.05, discount_factor = 0.75) {
  repeated_measures_power(effect_input = effect_size, effect_type = effect_type,
                          n = n, power = power, n_timepoints = n_timepoints,
                          correlation_between_measures = correlation_between_measures,
                          alpha = alpha, discount_factor = discount_factor)
}

#' Quick Repeated Measures Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param n_timepoints Number of timepoints (default = 2)
#' @param correlation_between_measures Correlation between measures (default = 0.5)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
repeated_measures_sample_size <- function(r_partial, power = 0.8, n_timepoints = 2,
                                          correlation_between_measures = 0.5, alpha = 0.05) {
  result <- repeated_measures_power(r_partial = r_partial, power = power,
                                    n_timepoints = n_timepoints,
                                    correlation_between_measures = correlation_between_measures,
                                    alpha = alpha)
  return(result$n)
}

#' Quick Repeated Measures Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param n_timepoints Number of timepoints (default = 2)
#' @param correlation_between_measures Correlation between measures (default = 0.5)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
repeated_measures_power_check <- function(r_partial, n, n_timepoints = 2,
                                          correlation_between_measures = 0.5, alpha = 0.05) {
  result <- repeated_measures_power(r_partial = r_partial, n = n,
                                    n_timepoints = n_timepoints,
                                    correlation_between_measures = correlation_between_measures,
                                    alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Repeated Measures Power Analysis
#' @param x Repeated measures power analysis result
#' @param ... Additional arguments
print.repeated_measures_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of timepoints:", x$n_timepoints, "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Design details
  cat("\nDesign details:\n")
  cat("  Correlation between measures:", round(x$correlation_between_measures, 3), "\n")
  cat("  Design efficiency:", round(x$design_efficiency, 2), "x\n")

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
