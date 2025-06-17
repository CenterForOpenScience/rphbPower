# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (v2.2 - Final Validated Engine)
# ==============================================================================

#' Linear Regression Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient for the predictor of interest (NULL to calculate).
#' @param n Sample size (NULL to calculate).
#' @param power Statistical power (1 - beta error probability) (NULL to calculate).
#' @param n_predictors Total number of predictors in the model (default = 1).
#' @param alpha The significance level (Type I error probability) (default = 0.05).
#' @param discount_factor The conservative planning discount factor (default = 0.75).
#' @param effect_input Alternative to r_partial for automatic conversion from other effect size metrics.
#' @param effect_type The type of `effect_input` ("r", "d", "f2", "r_squared").
#' @return A list object containing the power analysis results.
#' @export
linear_regression_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                    n_predictors = 1, alpha = 0.05, discount_factor = 0.75,
                                    effect_input = NULL, effect_type = "r") {
  
  # I. Parameter Validation & Setup
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = (discount_factor != 1.0))
  }
  
  target_param <- names(which(sapply(list(r_partial = r_partial, n = n, power = power), is.null)))
  if (length(target_param) != 1) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }
  
  # II. Core Calculation (Delegated to {pwr} package)
  f2 <- if(!is.null(r_partial)) partial_r_to_cohens_f2(r_partial) else NULL
  
  # Numerator degrees of freedom (u) for testing a single coefficient is ALWAYS 1.
  # The total number of predictors (n_predictors) is used for the denominator df (v).
  u_df <- 1
  v_df <- if(!is.null(n)) n - n_predictors - 1 else NULL
  
  power_to_pass <- if(target_param == "power") NULL else power
  
  result <- tryCatch(
    pwr::pwr.f2.test(u = u_df, v = v_df, f2 = f2, sig.level = alpha, power = power_to_pass),
    error = function(e) stop("Could not calculate power. Parameters may be invalid (e.g., n is too small).")
  )
  
  # III. Format Output
  if (target_param == "n") {
    # The result$v is the denominator df (n - k - 1). To get N, we must add k + 1 back.
    calculated_n <- ceiling(result$v + n_predictors + 1)
    final_result <- list(n = calculated_n)
  } else if (target_param == "power") {
    final_result <- list(power = result$power)
  } else { # target_param == "r_partial"
    calculated_r <- cohens_f2_to_partial_r(result$f2)
    final_result <- list(r_partial = calculated_r)
  }
  
  # Assemble the comprehensive final output object
  final_result$calculation_target <- sub("n", "sample_size", target_param)
  
  final_output <- list(
    analysis_type = paste("linear_regression", final_result$calculation_target, sep="_"),
    method = "Linear Regression Power Analysis (v2.2, Final Engine)",
    r_partial = if(target_param != "r_partial") r_partial else final_result$r_partial,
    n = if(target_param != "n") n else final_result$n,
    power = if(target_param != "power") power else final_result$power,
    n_predictors = n_predictors,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = final_result$calculation_target,
    effect_size_conversions = framework_conversion_summary(
      if(target_param != "r_partial") r_partial else final_result$r_partial, "r", apply_discount=FALSE
    )
  )
  final_output$interpretation <- interpret_effect_size(final_output$r_partial)
  final_output$degrees_of_freedom <- list(numerator = u_df, denominator = final_output$n - n_predictors - 1)
  
  class(final_output) <- "linear_regression_power_analysis"
  return(final_output)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Quick Linear Regression Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param n_predictors Number of predictors (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
#' @export
linear_regression_sample_size <- function(r_partial, power = 0.8, n_predictors = 1, alpha = 0.05) {
  linear_regression_power(r_partial = r_partial, power = power, n_predictors = n_predictors, alpha = alpha)$n
}

#' Quick Linear Regression Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param n_predictors Number of predictors (default = 1)
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
linear_regression_power_check <- function(r_partial, n, n_predictors = 1, alpha = 0.05) {
  linear_regression_power(r_partial = r_partial, n = n, n_predictors = n_predictors, alpha = alpha)$power
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Linear Regression Power Analysis
#' @param x Linear regression power analysis result
#' @param ... Additional arguments
#' @export
print.linear_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of predictors:", x$n_predictors, "\n")
  cat("Alpha level:", x$alpha, "\n")
  
  cat("\nModel details:\n")
  cat("  Degrees of freedom (num, den):", x$degrees_of_freedom$numerator, ",",
      x$degrees_of_freedom$denominator, "\n")
  
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat("  Cohen's f\u00B2:", round(conv$cohens_f2, 3), "\n")
    cat("  R\u00B2:", round(conv$r_squared, 3), "\n")
  }
  
  cat("\nFramework details:\n")
  cat("  Discount factor applied via effect_input:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}