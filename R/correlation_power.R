# ==============================================================================
# Correlation Power Analysis - Unified Framework Implementation
# ==============================================================================
#
# Power analysis for correlation studies using partial correlations as the
# foundation within the unified regression framework. All analyses integrate
# with core effect size conversion utilities.
#
# Author: Power Analysis Package
# Version: 1.2

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Correlation Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param two_tailed Two-tailed test (default = TRUE)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
correlation_power <- function(r_partial = NULL, n = NULL, power = 0.8,
                              alpha = 0.05, discount_factor = 0.75, two_tailed = TRUE,
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
                                    context = "for correlation power analysis")
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

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1")
  }

  if (!is.logical(two_tailed)) {
    stop("two_tailed must be TRUE or FALSE")
  }

  # Perform power calculation using base R
  if (is.null(power)) {
    # Calculate power
    calculated_power <- correlation_power_calculation(r_partial, n, alpha, two_tailed)

    result <- list(
      analysis_type = "correlation_power",
      method = "Correlation Power Analysis",
      r_partial = r_partial,
      n = n,
      power = calculated_power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    calculated_n <- correlation_sample_size_calculation(r_partial, power, alpha, two_tailed)

    result <- list(
      analysis_type = "correlation_sample_size",
      method = "Correlation Sample Size Analysis",
      r_partial = r_partial,
      n = calculated_n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- correlation_effect_size_calculation(n, power, alpha, two_tailed)

    result <- list(
      analysis_type = "correlation_effect_size",
      method = "Correlation Effect Size Analysis",
      r_partial = calculated_r_partial,
      n = n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      calculation_target = "effect_size"
    )
  }

  # Add framework integration
  r_for_conversion <- result$r_partial
  result$effect_size_conversions <- framework_conversion_summary(
    r_for_conversion, "r", apply_discount = FALSE
  )
  result$interpretation <- interpret_effect_size(r_for_conversion)

  class(result) <- "correlation_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Given Correlation and Sample Size
#' @importFrom stats pt qt
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Statistical power
correlation_power_calculation <- function(r_partial, n, alpha, two_tailed) {
  # Transform correlation to t-statistic
  t_calc <- r_partial * sqrt((n - 2) / (1 - r_partial^2))

  # Critical t-value
  df <- n - 2
  t_crit <- if (two_tailed) {
    qt(1 - alpha/2, df)
  } else {
    qt(1 - alpha, df)
  }

  # Calculate power
  if (two_tailed) {
    power <- 1 - pt(t_crit, df, ncp = abs(t_calc)) + pt(-t_crit, df, ncp = abs(t_calc))
  } else {
    power <- 1 - pt(t_crit, df, ncp = t_calc)
  }

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for Given Correlation and Power
#' @param r_partial Partial correlation
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required sample size
correlation_sample_size_calculation <- function(r_partial, power, alpha, two_tailed) {
  # Use iterative approach to find required sample size
  n_min <- 10
  n_max <- 10000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- correlation_power_calculation(r_partial, n_test, alpha, two_tailed)

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

#' Calculate Effect Size for Given Sample Size and Power
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required partial correlation
correlation_effect_size_calculation <- function(n, power, alpha, two_tailed) {
  # Use iterative approach to find required effect size
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- correlation_power_calculation(r_test, n, alpha, two_tailed)

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

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated Correlation Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
correlation_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                        power = 0.8, alpha = 0.05, discount_factor = 0.75) {
  correlation_power(effect_input = effect_size, effect_type = effect_type,
                    n = n, power = power, alpha = alpha, discount_factor = discount_factor)
}

#' Quick Correlation Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
correlation_sample_size <- function(r_partial, power = 0.8, alpha = 0.05) {
  result <- correlation_power(r_partial = r_partial, power = power, alpha = alpha)
  return(result$n)
}

#' Quick Correlation Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
correlation_power_check <- function(r_partial, n, alpha = 0.05) {
  result <- correlation_power(r_partial = r_partial, n = n, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Correlation Power Analysis
#' @importFrom cli symbol
#' @param x Correlation power analysis result
#' @param ... Additional arguments
print.correlation_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Test type:", if(x$two_tailed) "Two-tailed" else "One-tailed", "\n")

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
