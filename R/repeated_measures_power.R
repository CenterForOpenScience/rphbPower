# ==============================================================================
# CORE CALCULATION FUNCTIONS (DELEGATED TO pwr)
# ==============================================================================

#' Calculate Power for Repeated Measures
#' @param r_partial Partial correlation coefficient, analogous to sqrt(partial eta-squared) (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param alpha Significance level (default = 0.05)
#' @param correlation_between_measures Average correlation between measures (default = 0.5).
#' @param n_timepoints Number of repeated measures (default = 2). Note: engine uses paired t-test, best for 2 timepoints.
#' @return Power analysis results with framework integration
#' @export
repeated_measures_power_calculation <- function(r_partial, n, alpha, correlation_between_measures, n_timepoints) {
  # Convert r_partial to the equivalent Cohen's d for a paired t-test
  d <- partial_r_to_cohens_d(r_partial)

  # Use pwr.t.test for paired samples as the calculation engine
  power <- pwr::pwr.t.test(n = n, d = d, sig.level = alpha, type = "paired", alternative = "two.sided")$power
  return(power)
}

#' Calculate Sample Size for Repeated Measures
#' @param r_partial Partial correlation coefficient, analogous to sqrt(partial eta-squared) (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param correlation_between_measures Average correlation between measures (default = 0.5).
#' @param n_timepoints Number of repeated measures (default = 2). Note: engine uses paired t-test, best for 2 timepoints.
#' @return Sample Size
#' @export
repeated_measures_sample_size_calculation <- function(r_partial, power, alpha, correlation_between_measures, n_timepoints) {
  d <- partial_r_to_cohens_d(r_partial)
  n <- pwr::pwr.t.test(d = d, power = power, sig.level = alpha, type = "paired", alternative = "two.sided")$n
  return(ceiling(n))
}

#' Calculate Effect Size for Repeated Measures
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param correlation_between_measures Average correlation between measures (default = 0.5).
#' @param n_timepoints Number of repeated measures (default = 2). Note: engine uses paired t-test, best for 2 timepoints.
#' @return Power analysis results with framework integration
#' @export
repeated_measures_effect_size_calculation <- function(n, power, alpha, correlation_between_measures, n_timepoints) {
  # Find the required Cohen's d using the pwr package
  d_target <- pwr::pwr.t.test(n = n, power = power, sig.level = alpha, type = "paired", alternative = "two.sided")$d

  # Convert the resulting Cohen's d back to the framework's r_partial
  r_target <- cohens_d_to_partial_r(d_target)
  return(r_target)
}

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (WRAPPER)
# ==============================================================================
#' Repeated Measures Power Analysis with Framework Integration (v3.0)
#' @param r_partial Partial correlation coefficient, analogous to sqrt(partial eta-squared) (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_timepoints Number of repeated measures (default = 2). Note: engine uses paired t-test, best for 2 timepoints.
#' @param correlation_between_measures Average correlation between measures (default = 0.5).
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d")
#' @return Power analysis results with framework integration
#' @export
repeated_measures_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                    n_timepoints = 2, correlation_between_measures = 0.5,
                                    alpha = 0.05, discount_factor = 0.75,
                                    effect_input = NULL, effect_type = "r") {

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

  if (!provided_params["power"]) {
    result_val <- repeated_measures_power_calculation(r_partial, n, alpha, correlation_between_measures, n_timepoints)
    result <- list(power = result_val, calculation_target = "power")
  } else if (!provided_params["n"]) {
    result_val <- repeated_measures_sample_size_calculation(r_partial, power, alpha, correlation_between_measures, n_timepoints)
    result <- list(n = result_val, calculation_target = "sample_size")
  } else {
    result_val <- repeated_measures_effect_size_calculation(n, power, alpha, correlation_between_measures, n_timepoints)
    result <- list(r_partial = result_val, calculation_target = "effect_size")
  }

  final_result <- list(
    analysis_type = paste("repeated_measures", result$calculation_target, sep="_"),
    method = "Repeated Measures Power Analysis (v3.0, pwr.t.test Engine)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n = if(is.null(n)) result$n else n,
    power = if(is.null(power) || power_default_used) result$power else power,
    n_timepoints = n_timepoints,
    correlation_between_measures = correlation_between_measures,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = result$calculation_target
  )

  r_for_conversion <- final_result$r_partial
  if (exists("partial_r_to_cohens_d")) {
    final_result$effect_size_conversions <- list(
      cohens_d = partial_r_to_cohens_d(r_for_conversion)
    )
  }
  if(exists("interpret_effect_size")) {
    final_result$interpretation <- interpret_effect_size(r_for_conversion)
  }

  class(final_result) <- "repeated_measures_power_analysis"
  return(final_result)
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
#' @export
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
#' @export
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
#' @export
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
#' @export
print.repeated_measures_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0(" (", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of timepoints:", x$n_timepoints, "\n")
  cat("Correlation between measures:", x$correlation_between_measures, "\n")
  cat("Alpha level:", x$alpha, "\n")

  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Equivalent Cohen's d (for paired t-test):", round(conv$cohens_d, 3), "\n")
  }

  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, " (applied to initial effect_input)\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}
