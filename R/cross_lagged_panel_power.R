# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (WRAPPER)
# ==============================================================================
#' Cross-Lagged Panel Power Analysis with Framework Integration (v2.1)
#' @param r_partial Cross-lagged partial correlation coefficient (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param n_waves Number of measurement waves (for user context, default = 3)
#' @param stability_coefficient Autoregressive stability coefficient (for user context, default = 0.6)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d")
#' @return Power analysis results with framework integration
#' @export
cross_lagged_panel_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                     n_waves = 3, stability_coefficient = 0.6,
                                     alpha = 0.05, discount_factor = 0.75,
                                     effect_input = NULL, effect_type = "r") {

  # The power for a single cross-lagged path is equivalent to the power for a
  # single predictor in a multiple regression model where the outcome is
  # predicted by its own lag (the stability path) and the other variable's
  # lag (the cross-lagged path). Thus, n_predictors = 2.
  n_predictors <- 2

  # Delegate the entire calculation to the validated linear regression module.
  # The `n_waves` and `stability_coefficient` arguments are not used in the
  # calculation but are kept in the signature for user context and API consistency.
  result <- linear_regression_power(
    r_partial = r_partial,
    n = n,
    power = power,
    n_predictors = n_predictors,
    alpha = alpha,
    effect_input = effect_input,
    effect_type = effect_type
  )

  # Re-package the result to reflect the cross-lagged panel context
  if (is.null(result$power)) result$power <- if (!is.null(power)) power else 0.8
  if (is.null(result$n)) result$n <- n
  if (is.null(result$r_partial)) result$r_partial <- r_partial

  final_result <- list(
    analysis_type = "cross_lagged_panel",
    method = "Cross-Lagged Panel Power Analysis (v2.1, Regression Engine)",
    r_partial = result$r_partial,
    n = result$n,
    power = result$power,
    n_waves = n_waves,
    stability_coefficient = stability_coefficient,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = names(result)[1] # power, n, or r_partial
  )

  class(final_result) <- "cross_lagged_panel_power_analysis"
  return(final_result)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (RESTORED)
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

  cat("Calculation Target:", x$calculation_target, "\n\n")
  cat("  Cross-lagged partial r :", round(x$r_partial, 4), "\n")
  cat("  Sample size (N)        :", ceiling(x$n), "\n")
  cat("  Statistical Power      :", round(x$power, 3), "\n\n")

  cat("Model Parameters:\n")
  cat("  Number of waves        :", x$n_waves, "(for context)\n")
  cat("  Stability coefficient  :", x$stability_coefficient, "(for context)\n")
  cat("  Alpha level            :", x$alpha, "\n")
  cat("  Note: Power is for one path, controlling for one stability path (2 predictors total).\n\n")

  cat("Framework Details:\n")
  cat("  Discount factor applied  :", x$discount_factor, " (to initial effect_input)\n")
  cat("\n")
}
