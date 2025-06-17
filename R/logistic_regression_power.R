# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (v4.0 - Corrected Engine)
# ==============================================================================

#' Logistic Regression Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient (NULL to calculate).
#' @param n Sample size (NULL to calculate).
#' @param power Statistical power (NULL to calculate, default = 0.8).
#' @param n_predictors Number of predictors. (NOTE: Kept for compatibility, but the current engine only supports a single predictor).
#' @param alpha Significance level (default = 0.05).
#' @param discount_factor Conservative discount factor (default = 0.75).
#' @param mc_reps Kept for compatibility, unused by the new engine.
#' @param search_reps Kept for compatibility, unused by the new engine.
#' @param effect_input Raw effect size input (alternative to r_partial).
#' @param effect_type Type of `effect_input` ("r", "d", "f2", "r_squared").
#' @return A list object containing the power analysis results.
#' @export
logistic_regression_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                      n_predictors = 1, alpha = 0.05, discount_factor = 0.75,
                                      mc_reps = 1000, search_reps = 500, # Compatibility params
                                      effect_input = NULL, effect_type = "r") {

  # --- 1. Parameter Detection & Validation ---
  if (n_predictors > 1) {
    warning("Current engine only supports a single predictor. The 'n_predictors' argument is ignored.")
  }
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = (discount_factor != 1.0))
  }

  # Use the robust parameter detection logic from other validated modules
  if (is.null(power)) {
    power_default_used <- TRUE
    power <- 0.8
  } else {
    power_default_used <- FALSE
  }

  provided_params <- c(
    r_partial = !is.null(r_partial),
    n = !is.null(n),
    power = !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  calculation_target <- names(which(!provided_params))

  # --- 2. Core Calculation ---
  p0 <- 0.5 # Base probability (50/50 outcome)
  p1 <- NULL
  if (!is.null(r_partial)) {
    d <- partial_r_to_cohens_d(r_partial)
    odds_ratio <- exp(d * pi / sqrt(3))
    p1 <- (odds_ratio * p0) / (1 - p0 + odds_ratio * p0)
  }

  # **THE FIX**: Pass NULL to the argument that is the calculation target
  result <- WebPower::wp.logistic(
    n = if (calculation_target == "n") NULL else n,
    p0 = p0,
    p1 = p1,
    power = if (calculation_target == "power") NULL else power,
    alpha = alpha,
    family = "normal"
  )

  # --- 3. Format Output ---
  if (calculation_target == "n") {
    final_result <- list(n = ceiling(result$n), calculation_target = "sample_size")
  } else if (calculation_target == "power") {
    final_result <- list(power = result$power, calculation_target = "power")
  } else { # Solving for r_partial
    logit_p0 <- log(p0 / (1 - p0)); logit_p1 <- log(result$p1 / (1 - result$p1))
    b1 <- logit_p1 - logit_p0; f2 <- (b1^2 * 3) / (pi^2)
    calculated_r <- sqrt(f2 / (1 + f2))
    final_result <- list(r_partial = calculated_r, calculation_target = "effect_size")
  }

  # Assemble the final, consistent output object
  final_output <- list(
    analysis_type = paste("logistic_regression", final_result$calculation_target, sep = "_"),
    method = "Logistic Regression Power Analysis (v4.2, Conforming)",
    r_partial = if (calculation_target != "r_partial") r_partial else final_result$r_partial,
    n = if (calculation_target != "n") n else final_result$n,
    power = if (calculation_target != "power") power else final_result$power,
    n_predictors = 1,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = final_result$calculation_target
  )
  r_for_conversion <- final_output$r_partial
  final_output$effect_size_conversions <- list(
    cohens_d = partial_r_to_cohens_d(r_for_conversion),
    odds_ratio = exp((partial_r_to_cohens_d(r_for_conversion) * pi) / sqrt(3))
  )
  final_output$interpretation <- interpret_effect_size(r_for_conversion)
  class(final_output) <- "logistic_regression_power_analysis"
  return(final_output)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (Unchanged)
# ==============================================================================

#' Quick Logistic Regression Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
#' @export
logistic_regression_sample_size <- function(r_partial, power = 0.8, alpha = 0.05) {
  logistic_regression_power(r_partial = r_partial, power = power, alpha = alpha)$n
}

#' Quick Logistic Regression Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
logistic_regression_power_check <- function(r_partial, n, alpha = 0.05) {
  logistic_regression_power(r_partial = r_partial, n = n, alpha = alpha)$power
}

# ==============================================================================
# PRINT METHOD (Unchanged)
# ==============================================================================

#' Print Method for Logistic Regression Power Analysis
#' @param x Logistic regression power analysis result
#' @param ... Additional arguments
#' @export
print.logistic_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0(" (", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of predictors:", x$n_predictors, "(Note: engine assumes 1)\n")
  cat("Alpha level:", x$alpha, "\n")
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Equivalent Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat("  Equivalent Odds Ratio:", round(conv$odds_ratio, 3), "\n")
  }
  cat("\nFramework details:\n")
  cat("  Discount factor applied via effect_input:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}
