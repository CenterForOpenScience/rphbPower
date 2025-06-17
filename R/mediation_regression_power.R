# ==============================================================================
# CORE CALCULATION ENGINE (v4.1 - Corrected)
# ==============================================================================

#' Calculate Power for Mediation Analysis
#' This is the internal engine that performs the power calculation.
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @return Statistical power
#' @export
mediation_power_calculation <- function(r_a, r_b, n, alpha) {
  # This calculation is based on the Aroian test for the significance of an indirect effect.

  # Corrected SE formulas
  se_a_sq <- (1 - r_a^2) / (n - 2)
  se_b_sq <- (1 - r_b^2) / (n - 3)

  # CORRECTED: The third term (se_a_sq * se_b_sq) was missing and has been restored.
  se_ab_sq_val <- r_a^2 * se_b_sq + r_b^2 * se_a_sq + se_a_sq * se_b_sq
  if (is.na(se_ab_sq_val) || se_ab_sq_val <= 0) return(NA)

  z_stat <- (r_a * r_b) / sqrt(se_ab_sq_val)
  z_crit <- stats::qnorm(1 - alpha / 2)
  power <- stats::pnorm(z_stat - z_crit) + stats::pnorm(-z_stat - z_crit)

  return(power)
}

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (v4.1)
# ==============================================================================
#' Mediation Regression Power Analysis with Framework Integration
#'
#' @param r_a Path a (X->M) coefficient as correlation (NULL to calculate).
#' @param r_b Path b (M->Y|X) coefficient as partial correlation (NULL to calculate).
#' @param n Sample size (NULL to calculate).
#' @param power Statistical power (NULL to calculate, default = 0.8).
#' @param alpha Significance level (default = 0.05).
#' @param discount_factor Conservative discount factor applied to input effect sizes (default = 0.75).
#'
#' @return A list object containing the power analysis results.
#' @export
mediation_regression_power <- function(r_a = NULL, r_b = NULL, n = NULL, power = NULL,
                                       alpha = 0.05, discount_factor = 0.75) {

  params <- list(r_a = r_a, r_b = r_b, n = n, power = power)
  null_params <- sapply(params, is.null)

  if (sum(null_params) != 1) {
    stop("Provide values for exactly three of the four parameters: r_a, r_b, n, power.")
  }

  target_param <- names(which(null_params))

  # Apply discount factor to any input effect sizes
  if (!is.null(params$r_a)) params$r_a <- apply_discount_factor(params$r_a, discount_factor)
  if (!is.null(params$r_b)) params$r_b <- apply_discount_factor(params$r_b, discount_factor)

  # Define a target function for uniroot to solve for the missing parameter
  target_function <- function(target_val, param_name, current_params) {
    temp_params <- current_params
    temp_params[[param_name]] <- target_val
    # Power can be NULL here when solving for it, so we need to pass the user-provided value
    power_goal <- if(is.null(current_params$power)) 0.8 else current_params$power
    calculated_power <- mediation_power_calculation(temp_params$r_a, temp_params$r_b, temp_params$n, alpha)
    # uniroot finds where the function crosses zero
    return(calculated_power - power_goal)
  }

  if (target_param == "power") {
    calculated_power <- mediation_power_calculation(params$r_a, params$r_b, params$n, alpha)
    final_result <- list(power = calculated_power, calculation_target = "power")
  } else if (target_param == "n") {
    root <- tryCatch(stats::uniroot(function(x) target_function(x, "n", params), interval=c(5, 100000)),
                     error = function(e) stop("Could not find a sample size. Effect sizes may be too small."))
    final_result <- list(n = ceiling(root$root), calculation_target = "sample_size")
  } else { # Solving for r_a or r_b
    root <- tryCatch(stats::uniroot(function(x) target_function(x, target_param, params), interval=c(0.001, 0.999)),
                     error = function(e) stop("Could not find an effect size. Power may be too high for the given N."))
    final_result <- list(value = root$root, calculation_target = target_param)
  }

  final_output <- list(
    method = "Mediation Regression Power Analysis (v4.1, Analytical Engine)",
    calculation_target = final_result$calculation_target,
    alpha = alpha,
    discount_factor = discount_factor,
    r_a = if(target_param == "r_a") final_result$value else params$r_a,
    r_b = if(target_param == "r_b") final_result$value else params$r_b,
    n = if(target_param == "n") final_result$n else params$n,
    power = if(target_param == "power") final_result$power else params$power
  )
  final_output$indirect_effect <- final_output$r_a * final_output$r_b
  final_output$interpretation <- interpret_mediation_effect(final_output$r_a, final_output$r_b)

  class(final_output) <- "mediation_regression_power_analysis"
  return(final_output)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (Unchanged)
# ==============================================================================

#' Quick Mediation Sample Size Calculation
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param power Target power (default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Framework discount factor
#' @return Required sample size
#' @export
mediation_regression_sample_size <- function(r_a, r_b, power = 0.8, alpha = 0.05, discount_factor = 0.75) {
  mediation_regression_power(r_a = r_a, r_b = r_b, power = power,
                             alpha = alpha, discount_factor = discount_factor)$n
}

#' Quick Mediation Power Calculation
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Framework discount factor
#' @return Statistical power
#' @export
mediation_regression_power_check <- function(r_a, r_b, n, alpha = 0.05, discount_factor = 0.75) {
  mediation_regression_power(r_a = r_a, r_b = r_b, n = n,
                             alpha = alpha, discount_factor = discount_factor)$power
}

# ==============================================================================
# HELPER AND PRINT FUNCTIONS (Unchanged)
# ==============================================================================

#' Interpret Mediation Effect Size
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @return Interpretation string
#' @export
interpret_mediation_effect <- function(r_a, r_b) {
  indirect_effect <- abs(r_a * r_b)
  if (indirect_effect < 0.01) return("Negligible mediation")
  if (indirect_effect < 0.09) return("Small mediation")
  if (indirect_effect < 0.25) return("Medium mediation")
  return("Large mediation")
}

#' Print Method for Mediation Regression Power Analysis
#' @param x Mediation regression power analysis result
#' @param ... Additional arguments
#' @export
print.mediation_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  cat("Calculation Target:", x$calculation_target, "\n")
  cat("Engine: Analytical (Aroian formula for Sobel test)\n")
  cat("Discount Factor Applied to Inputs:", x$discount_factor, "\n\n")

  cat("--- Effective Path Coefficients Used in Calculation ---\n")
  cat("NOTE: These values reflect inputs after the discount factor has been applied.\n\n")

  cat("  A-path (X->M | r_a)      :", round(x$r_a, 4), "\n")
  cat("  B-path (M->Y|X | r_b)    :", round(x$r_b, 4), "\n")
  cat("  Indirect effect (a*b)   :", round(x$indirect_effect, 4), paste0("(", x$interpretation, ")\n\n"))

  cat("--- Results ---\n")
  cat("  Required sample size (N):", ceiling(x$n), "\n")
  cat("  Achieved power          :", round(x$power, 3), "\n")
  cat("  Alpha level             :", x$alpha, "\n\n")
}
