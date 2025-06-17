# ==============================================================================
# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Correlation Power Analysis with Framework Integration (v2.2)
#' @param r_partial Partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param two_tailed Two-tailed test (default = TRUE)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
correlation_power <- function(r_partial = NULL, n = NULL, power = NULL,
                              alpha = 0.05, discount_factor = 0.75, two_tailed = TRUE,
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

  # Validate inputs
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for correlation power analysis")
  }
  if (!is.null(n) && (n < 4 || n != round(n))) {
    stop("Sample size must be a whole number >= 4 for a correlation test.")
  }
  if (!is.null(power) && (power <= 0 || power >= 1)) {
    stop("Power must be between 0 and 1.")
  }
  if (!is.logical(two_tailed)) {
    stop("two_tailed must be TRUE or FALSE")
  }

  # Perform power calculation
  if (!provided_params["power"]) {
    calculated_power <- correlation_power_calculation(r_partial, n, alpha, two_tailed)
    result <- list(power = calculated_power, calculation_target = "power")
  } else if (!provided_params["n"]) {
    calculated_n <- correlation_sample_size_calculation(r_partial, power, alpha, two_tailed)
    result <- list(n = calculated_n, calculation_target = "sample_size")
  } else {
    calculated_r_partial <- correlation_effect_size_calculation(n, power, alpha, two_tailed)
    result <- list(r_partial = calculated_r_partial, calculation_target = "effect_size")
  }

  # Assemble final output object
  final_result <- list(
    analysis_type = paste("correlation", result$calculation_target, sep="_"),
    method = "Correlation Power Analysis (v2.2)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n = if(is.null(n)) result$n else n,
    power = if(is.null(power) || power_default_used) result$power else power,
    alpha = alpha,
    discount_factor = discount_factor,
    two_tailed = two_tailed,
    calculation_target = result$calculation_target
  )

  final_result$effect_size_conversions <- framework_conversion_summary(
    final_result$r_partial, "r", apply_discount = FALSE
  )
  final_result$interpretation <- interpret_effect_size(final_result$r_partial)

  class(final_result) <- "correlation_power_analysis"
  return(final_result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS (REFINED for v2.2)
# ==============================================================================

#' Calculate Power for Given Correlation and Sample Size
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Statistical power
#' @export
correlation_power_calculation <- function(r_partial, n, alpha, two_tailed) {
  df <- n - 2
  if (df <= 0) return(0)

  t_stat <- r_partial * sqrt(df) / sqrt(1 - r_partial^2)
  ncp <- t_stat

  if (two_tailed) {
    t_crit <- stats::qt(1 - alpha / 2, df)
    power <- stats::pt(t_crit, df, ncp = ncp, lower.tail = FALSE) + stats::pt(-t_crit, df, ncp = ncp, lower.tail = TRUE)
  } else {
    t_crit <- stats::qt(1 - alpha, df)
    power <- stats::pt(t_crit, df, ncp = ncp, lower.tail = FALSE)
  }
  return(pmax(0.001, pmin(0.999, power)))
}

#' Calculate Sample Size for Given Correlation and Power
#' @param r_partial Partial correlation
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required sample size
#' @export
correlation_sample_size_calculation <- function(r_partial, power, alpha, two_tailed) {
  power_func <- function(n) correlation_power_calculation(r_partial, n, alpha, two_tailed) - power
  result <- tryCatch(stats::uniroot(power_func, interval = c(4, 100000))$root, error = function(e) NA)
  if (is.na(result)) stop("Could not find a sample size to achieve the desired power. The effect size may be too small.")
  return(ceiling(result))
}

#' Calculate Effect Size for Given Sample Size and Power
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required partial correlation
#' @export
correlation_effect_size_calculation <- function(n, power, alpha, two_tailed) {
  power_func <- function(r) correlation_power_calculation(r, n, alpha, two_tailed) - power
  result <- tryCatch(stats::uniroot(power_func, interval = c(0.001, 0.999))$root, error = function(e) NA)
  if (is.na(result)) stop("Could not find an effect size to achieve the desired power.")
  return(result)
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
#' @export
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
#' @export
correlation_sample_size <- function(r_partial, power = 0.8, alpha = 0.05) {
  result <- correlation_power(r_partial = r_partial, power = power, alpha = alpha)
  return(result$n)
}

#' Quick Correlation Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
correlation_power_check <- function(r_partial, n, alpha = 0.05) {
  result <- correlation_power(r_partial = r_partial, n = n, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Correlation Power Analysis
#' @param x Correlation power analysis result
#' @param ... Additional arguments
#' @export
print.correlation_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Test type:", if(x$two_tailed) "Two-tailed" else "One-tailed", "\n")

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
  cat("\n")
}
