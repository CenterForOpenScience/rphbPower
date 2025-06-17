# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' SEM Direct Effects Power Analysis (Approximation)
#' @param r_partial Latent path coefficient (standardized beta) as partial correlation (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_predictors Number of predictors of the dependent variable (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param measurement_reliability Average reliability (e.g., Cronbach's alpha) of the constructs (default = 0.85)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
sem_direct_effects_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                     n_predictors = 1, alpha = 0.05,
                                     discount_factor = 0.75, measurement_reliability = 0.85,
                                     effect_input = NULL, effect_type = "r") {

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
    n = !missing(n) && !is.null(n),
    power = !missing(power) && !is.null(power) && !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  # Perform power calculation
  if (!provided_params["power"]) {
    calculated_power <- sem_direct_effects_power_calculation(r_partial, n, n_predictors, alpha, measurement_reliability)
    result <- list(power = calculated_power, calculation_target = "power")
  } else if (!provided_params["n"]) {
    calculated_n <- sem_direct_effects_sample_size_calculation(r_partial, power, n_predictors, alpha, measurement_reliability)
    result <- list(n = calculated_n, calculation_target = "sample_size")
  } else {
    calculated_r_partial <- sem_direct_effects_effect_size_calculation(n, power, n_predictors, alpha, measurement_reliability)
    result <- list(r_partial = calculated_r_partial, calculation_target = "effect_size")
  }

  # Assemble final output object
  final_result <- list(
    analysis_type = paste("sem_direct_effects", result$calculation_target, sep="_"),
    method = "SEM Direct Effects Power Analysis (v2.4, t-test approximation)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n = if(is.null(n)) result$n else n,
    power = if(is.null(power) || power_default_used) result$power else power,
    n_predictors = n_predictors,
    alpha = alpha,
    discount_factor = discount_factor,
    measurement_reliability = measurement_reliability,
    calculation_target = result$calculation_target
  )

  final_result$interpretation <- interpret_effect_size(final_result$r_partial)
  class(final_result) <- "sem_direct_effects_power_analysis"
  return(final_result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS (v2.4)
# ==============================================================================

#' Calculate Power for SEM Direct Effects (Approximation)
#' @param r_partial Partial correlation (standardized path coefficient)
#' @param n Sample size
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @param measurement_reliability Average measurement reliability
#' @return Statistical power
#' @export
sem_direct_effects_power_calculation <- function(r_partial, n, n_predictors, alpha, measurement_reliability) {
  # Attenuate the latent effect size by sqrt of reliability to get the expected observed effect
  r_observed <- r_partial * sqrt(measurement_reliability)

  # Power is calculated as a standard t-test for a regression coefficient.
  df_den <- n - n_predictors - 1
  if (df_den <= 0) return(0)

  ncp <- r_observed * sqrt(df_den) / sqrt(1 - r_observed^2)
  t_crit <- stats::qt(1 - alpha / 2, df_den)
  power <- stats::pt(t_crit, df_den, ncp = ncp, lower.tail = FALSE) + stats::pt(-t_crit, df_den, ncp = ncp, lower.tail = TRUE)
  return(pmax(0.001, pmin(0.999, power)))
}

#' Calculate Sample Size for SEM Direct Effects (Approximation)
#' @param r_partial Partial correlation (standardized path coefficient)
#' @param power Target power
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @param measurement_reliability Average measurement reliability
#' @return Required sample size
#' @export
sem_direct_effects_sample_size_calculation <- function(r_partial, power, n_predictors, alpha, measurement_reliability) {
  power_func <- function(n) {
    sem_direct_effects_power_calculation(r_partial, n, n_predictors, alpha, measurement_reliability) - power
  }
  min_n <- n_predictors + 2
  result <- tryCatch(stats::uniroot(power_func, interval = c(min_n, 100000))$root, error = function(e) NA)
  if(is.na(result)) stop("Could not find sample size. Effect may be too small or power too high.")
  return(ceiling(result))
}

#' Calculate Effect Size for SEM Direct Effects (Approximation)
#' @param n Sample size
#' @param power Target power
#' @param n_predictors Number of predictors
#' @param alpha Significance level
#' @param measurement_reliability Average measurement reliability
#' @return Required partial correlation
#' @export
sem_direct_effects_effect_size_calculation <- function(n, power, n_predictors, alpha, measurement_reliability) {
  power_func <- function(r_latent) {
    sem_direct_effects_power_calculation(r_latent, n, n_predictors, alpha, measurement_reliability) - power
  }
  # Upper bound must be < 1 after attenuation
  upper_bound <- 0.999 / sqrt(measurement_reliability)
  result <- tryCatch(stats::uniroot(power_func, interval = c(0.001, min(0.999, upper_bound)))$root, error = function(e) NA)
  if(is.na(result)) stop("Could not find effect size. Power may be too high for the given N.")
  return(result)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for SEM Direct Effects Power Analysis
#' @param x SEM direct effects power analysis result
#' @param ... Additional arguments
#' @export
print.sem_direct_effects_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  cat("Latent path coefficient (\u03B2):", round(x$r_partial, 4), paste0(" (", x$interpretation, ")"), "\n")
  cat("Required sample size (N):", x$n, "\n")
  cat("Achieved statistical power:", round(x$power, 3), "\n\n")

  cat("Model Parameters:\n")
  cat("  Number of predictors of DV:", x$n_predictors, "\n")
  cat("  Assumed measurement reliability:", round(x$measurement_reliability, 3), "\n")
  cat("  Alpha level:", x$alpha, "\n\n")

  cat("NOTE: This is an approximation. For formal analysis, consider\n")
  cat("      simulation-based power analysis or specialized software (e.g., pwrSEM).\n\n")

  cat("Framework Details:\n")
  cat("  Discount factor applied:", x$discount_factor, "(to initial effect_input)\n")
  cat("  Calculation target:", x$calculation_target, "\n\n")
}
