# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Mediation SEM Power Analysis with Framework Integration
#' @param r_a X -> Mediator path coefficient (NULL to calculate)
#' @param r_b Mediator -> Y path coefficient (controlling for X)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param measurement_error Measurement error proportion (default = 0.1)
#' @param n_indicators_total Total number of indicators across all latent variables (default = 9)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input_a Raw effect size input for path a (alternative to r_a)
#' @param effect_input_b Raw effect size input for path b (alternative to r_b)
#' @param effect_type Type of effect inputs ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
mediation_sem_power <- function(r_a = NULL, r_b = NULL, n = NULL, power = 0.8,
                                measurement_error = 0.1, n_indicators_total = 9,
                                alpha = 0.05, discount_factor = 0.75,
                                effect_input_a = NULL, effect_input_b = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input_a) && !is.null(effect_input_b)) {
    r_a <- framework_effect_size(effect_input_a, effect_type, apply_discount = TRUE)
    r_b <- framework_effect_size(effect_input_b, effect_type, apply_discount = TRUE)
  }

  # Use missing() to detect actually provided parameters
  provided_params <- c(
    r_a = !missing(r_a) && !is.null(r_a),
    r_b = !missing(r_b) && !is.null(r_b),
    n = !missing(n) && !is.null(n),
    power = !missing(power) && !is.null(power)
  )

  # Check that r_a and r_b are provided, and exactly one of n or power
  if (!provided_params["r_a"] || !provided_params["r_b"]) {
    stop("Both r_a and r_b must be provided (or use effect_input_a/b)")
  }

  n_power_provided <- sum(provided_params[c("n", "power")])
  if (n_power_provided != 1) {
    stop("Provide exactly one of: n or power")
  }

  # Validate inputs using framework functions
  if (!is.null(r_a)) {
    r_a <- validate_partial_r(r_a, allow_zero = FALSE, context = "for path a")
  }
  if (!is.null(r_b)) {
    r_b <- validate_partial_r(r_b, allow_zero = FALSE, context = "for path b")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < 100)) {
      stop("Sample size must be whole number >= 100")
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(measurement_error) || measurement_error < 0 || measurement_error >= 0.5) {
    stop("Measurement error must be between 0 and 0.5")
  }

  if (!is.numeric(n_indicators_total) || n_indicators_total < 6 || n_indicators_total != round(n_indicators_total)) {
    stop("Total indicators must be integer >= 6")
  }

  # Calculate SEM-specific adjustments
  sem_complexity_factor <- calculate_sem_complexity_factor(n_indicators_total, measurement_error)

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- mediation_sem_power_calculation(r_a, r_b, n, sem_complexity_factor, alpha)

    result <- list(
      analysis_type = "mediation_sem_power",
      method = "SEM Mediation Power Analysis",
      r_a = r_a,
      r_b = r_b,
      n = n,
      power = calculated_power,
      indirect_effect = r_a * r_b,
      measurement_error = measurement_error,
      n_indicators_total = n_indicators_total,
      sem_complexity_factor = sem_complexity_factor,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else {
    # Calculate sample size
    calculated_n <- mediation_sem_sample_size_calculation(r_a, r_b, power, sem_complexity_factor, alpha)

    result <- list(
      analysis_type = "mediation_sem_sample_size",
      method = "SEM Mediation Sample Size Analysis",
      r_a = r_a,
      r_b = r_b,
      n = calculated_n,
      power = power,
      indirect_effect = r_a * r_b,
      measurement_error = measurement_error,
      n_indicators_total = n_indicators_total,
      sem_complexity_factor = sem_complexity_factor,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )
  }

  # Add framework integration
  result$effect_size_conversions <- list(
    path_a = framework_conversion_summary(r_a, "r", apply_discount = FALSE),
    path_b = framework_conversion_summary(r_b, "r", apply_discount = FALSE),
    indirect_effect = framework_conversion_summary(result$indirect_effect, "r", apply_discount = FALSE)
  )

  result$interpretation <- interpret_effect_size(result$indirect_effect)

  # Add SEM-specific information
  result$sem_details <- list(
    measurement_model_complexity = n_indicators_total,
    estimated_measurement_error = measurement_error,
    model_type = "latent_variable_mediation",
    degrees_of_freedom_penalty = sem_complexity_factor
  )

  class(result) <- "mediation_sem_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for SEM Mediation
#' @param r_a Path a coefficient
#' @param r_b Path b coefficient
#' @param n Sample size
#' @param sem_complexity_factor SEM complexity adjustment
#' @param alpha Significance level
#' @return Statistical power
mediation_sem_power_calculation <- function(r_a, r_b, n, sem_complexity_factor, alpha) {
  # Calculate indirect effect
  indirect_effect <- r_a * r_b

  # Adjust effective sample size for SEM complexity
  n_effective <- n / sem_complexity_factor

  # Use Sobel test statistic adjusted for SEM
  se_indirect <- sqrt((r_b^2 * (1 - r_a^2) + r_a^2 * (1 - r_b^2)) / n_effective)

  # Z-statistic for indirect effect
  z_stat <- abs(indirect_effect) / se_indirect

  # Power calculation (two-tailed)
  z_crit <- qnorm(1 - alpha/2)
  power <- 1 - pnorm(z_crit - z_stat) + pnorm(-z_crit - z_stat)

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for SEM Mediation
#' @param r_a Path a coefficient
#' @param r_b Path b coefficient
#' @param power Target power
#' @param sem_complexity_factor SEM complexity adjustment
#' @param alpha Significance level
#' @return Required sample size
mediation_sem_sample_size_calculation <- function(r_a, r_b, power, sem_complexity_factor, alpha) {
  # Use iterative approach
  n_min <- 100
  n_max <- 2000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- mediation_sem_power_calculation(r_a, r_b, n_test, sem_complexity_factor, alpha)

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

#' Calculate SEM Complexity Factor
#' @param n_indicators_total Total number of indicators
#' @param measurement_error Measurement error proportion
#' @return Complexity adjustment factor
calculate_sem_complexity_factor <- function(n_indicators_total, measurement_error) {
  # Base complexity from number of indicators
  indicator_penalty <- 1 + (n_indicators_total - 6) * 0.05

  # Measurement error penalty
  measurement_penalty <- 1 + measurement_error * 2

  # Combined complexity factor
  complexity_factor <- indicator_penalty * measurement_penalty

  # Cap at reasonable bounds
  return(pmax(1.2, pmin(complexity_factor, 3.0)))
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated SEM Mediation Power Analysis
#' @param effect_size_a Effect size for path a
#' @param effect_size_b Effect size for path b
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param measurement_error Measurement error proportion
#' @param n_indicators_total Total indicators
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
mediation_sem_framework_power <- function(effect_size_a, effect_size_b, effect_type = "r",
                                          n = NULL, power = 0.8, measurement_error = 0.1,
                                          n_indicators_total = 9, alpha = 0.05,
                                          discount_factor = 0.75) {
  mediation_sem_power(effect_input_a = effect_size_a, effect_input_b = effect_size_b,
                      effect_type = effect_type, n = n, power = power,
                      measurement_error = measurement_error,
                      n_indicators_total = n_indicators_total,
                      alpha = alpha, discount_factor = discount_factor)
}

#' Quick SEM Mediation Sample Size Calculation
#' @param r_a Path a coefficient
#' @param r_b Path b coefficient
#' @param power Target power (default = 0.8)
#' @param measurement_error Measurement error (default = 0.1)
#' @param n_indicators_total Total indicators (default = 9)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
mediation_sem_sample_size <- function(r_a, r_b, power = 0.8, measurement_error = 0.1,
                                      n_indicators_total = 9, alpha = 0.05) {
  result <- mediation_sem_power(r_a = r_a, r_b = r_b, power = power,
                                measurement_error = measurement_error,
                                n_indicators_total = n_indicators_total, alpha = alpha)
  return(result$n)
}

#' Quick SEM Mediation Power Calculation
#' @param r_a Path a coefficient
#' @param r_b Path b coefficient
#' @param n Sample size
#' @param measurement_error Measurement error (default = 0.1)
#' @param n_indicators_total Total indicators (default = 9)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
mediation_sem_power_check <- function(r_a, r_b, n, measurement_error = 0.1,
                                      n_indicators_total = 9, alpha = 0.05) {
  result <- mediation_sem_power(r_a = r_a, r_b = r_b, n = n,
                                measurement_error = measurement_error,
                                n_indicators_total = n_indicators_total, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for SEM Mediation Power Analysis
#' @param x SEM mediation power analysis result
#' @param ... Additional arguments
print.mediation_sem_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat(paste0("Path a (X ", symbol$arrow_right, " M):"), round(x$r_a, 4), "\n")
  cat(paste0("Path b (M ", symbol$arrow_right, " Y|X):"), round(x$r_b, 4), "\n")
  cat("Indirect effect:", round(x$indirect_effect, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")

  # SEM details
  cat("\nSEM model details:\n")
  cat("  Total indicators:", x$n_indicators_total, "\n")
  cat("  Measurement error:", round(x$measurement_error, 2), "\n")
  cat("  Complexity factor:", round(x$sem_complexity_factor, 2), "\n")

  # Framework conversions
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions (indirect effect):\n")
    conv <- x$effect_size_conversions$indirect_effect
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
