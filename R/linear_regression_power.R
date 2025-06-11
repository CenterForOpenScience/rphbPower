# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Linear Regression Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_predictors Number of predictors in model (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
linear_regression_power <- function(r_partial = NULL, n = NULL, power = 0.8,
                                    n_predictors = 1, alpha = 0.05, discount_factor = 0.75,
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
                                    context = "for linear regression power analysis")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < n_predictors + 10)) {
      stop(paste("Sample size must be whole number >=", n_predictors + 10))
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(n_predictors) || n_predictors < 1 || n_predictors != round(n_predictors)) {
    stop("Number of predictors must be positive integer")
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1")
  }

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- linear_regression_power_calculation(r_partial, n, n_predictors, alpha)

    result <- list(
      analysis_type = "linear_regression_power",
      method = "Linear Regression Power Analysis",
      r_partial = r_partial,
      n = n,
      power = calculated_power,
      n_predictors = n_predictors,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    calculated_n <- linear_regression_sample_size_calculation(r_partial, power, n_predictors, alpha)

    result <- list(
      analysis_type = "linear_regression_sample_size",
      method = "Linear Regression Sample Size Analysis",
      r_partial = r_partial,
      n = calculated_n,
      power = power,
      n_predictors = n_predictors,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- linear_regression_effect_size_calculation(n, power, n_predictors, alpha)

    result <- list(
      analysis_type = "linear_regression_effect_size",
      method = "Linear Regression Effect Size Analysis",
      r_partial = calculated_r_partial,
      n = n,
      power = power,
      n_predictors = n_predictors,
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

  # Add model-specific information
  result$degrees_of_freedom <- list(
    numerator = n_predictors,
    denominator = result$n - n_predictors - 1
  )

  class(result) <- "linear_regression_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Linear Regression
#' @importFrom stats qf pf
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @return Statistical power
linear_regression_power_calculation <- function(r_partial, n, n_predictors, alpha) {
  # Convert partial correlation to F-statistic
  f2 <- partial_r_to_cohens_f2(r_partial)

  # Degrees of freedom
  df_num <- n_predictors
  df_den <- n - n_predictors - 1

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

#' Calculate Sample Size for Linear Regression
#' @param r_partial Partial correlation
#' @param power Target power
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @return Required sample size
linear_regression_sample_size_calculation <- function(r_partial, power, n_predictors, alpha) {
  # Use iterative approach to find required sample size
  n_min <- n_predictors + 10
  n_max <- 10000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- linear_regression_power_calculation(r_partial, n_test, n_predictors, alpha)

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

#' Calculate Effect Size for Linear Regression
#' @param n Sample size
#' @param power Target power
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @return Required partial correlation
linear_regression_effect_size_calculation <- function(n, power, n_predictors, alpha) {
  # Use iterative approach to find required effect size
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- linear_regression_power_calculation(r_test, n, n_predictors, alpha)

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

#' Framework-Integrated Linear Regression Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
linear_regression_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                              power = 0.8, n_predictors = 1, alpha = 0.05,
                                              discount_factor = 0.75) {
  linear_regression_power(effect_input = effect_size, effect_type = effect_type,
                          n = n, power = power, n_predictors = n_predictors,
                          alpha = alpha, discount_factor = discount_factor)
}

#' Quick Linear Regression Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param n_predictors Number of predictors (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
linear_regression_sample_size <- function(r_partial, power = 0.8, n_predictors = 1, alpha = 0.05) {
  result <- linear_regression_power(r_partial = r_partial, power = power,
                                    n_predictors = n_predictors, alpha = alpha)
  return(result$n)
}

#' Quick Linear Regression Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param n_predictors Number of predictors (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
linear_regression_power_check <- function(r_partial, n, n_predictors = 1, alpha = 0.05) {
  result <- linear_regression_power(r_partial = r_partial, n = n,
                                    n_predictors = n_predictors, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Linear Regression Power Analysis
#' @param x Linear regression power analysis result
#' @param ... Additional arguments
print.linear_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of predictors:", x$n_predictors, "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Model details
  cat("\nModel details:\n")
  cat("  Degrees of freedom (num, den):", x$degrees_of_freedom$numerator, ",",
      x$degrees_of_freedom$denominator, "\n")

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
