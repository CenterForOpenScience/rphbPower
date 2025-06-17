# ==============================================================================
# CORE CALCULATION ENGINE (v5.1 - Corrected)
# ==============================================================================

#' Calculate Power for SEM Mediation (Analytical Approximation)
#' @param r_a X -> Mediator path coefficient (NULL to calculate)
#' @param r_b Mediator -> Y path coefficient (controlling for X)
#' @param n Sample size (NULL to calculate)
#' @param alpha Significance level (default = 0.05)
#' @param reliability Average reliability of constructs
#' @return Power analysis results with framework integration
#' @export
sem_mediation_power_calculation <- function(r_a, r_b, n, alpha, reliability) {
  if (n < 4) return(0)

  # Attenuate the latent paths to get the expected observed paths
  attenuation_factor <- sqrt(reliability)
  r_a_obs <- r_a * attenuation_factor
  r_b_obs <- r_b * attenuation_factor

  # Corrected SE formulas
  se_a_sq <- (1 - r_a_obs^2) / (n - 2)
  se_b_sq <- (1 - r_b_obs^2) / (n - 3)

  # CORRECTED: The third term (se_a_sq * se_b_sq) was missing and has been restored.
  se_ab_sq_val <- r_a_obs^2 * se_b_sq + r_b_obs^2 * se_a_sq + se_a_sq * se_b_sq
  if (is.na(se_ab_sq_val) || se_ab_sq_val <= 0) return(NA)

  z_ab <- (r_a_obs * r_b_obs) / sqrt(se_ab_sq_val)
  z_crit <- stats::qnorm(1 - alpha / 2)
  power <- stats::pnorm(z_ab - z_crit) + stats::pnorm(-z_ab - z_crit)

  return(pmax(0.001, pmin(0.999, power)))
}

#' SEM Mediation Sample Size Calculation
#' This internal function is used by the main wrapper to solve for sample size.
#' @param r_a X -> Mediator path coefficient (NULL to calculate)
#' @param r_b Mediator -> Y path coefficient (controlling for X)
#' @param power Statistical power
#' @param alpha Significance level (default = 0.05)
#' @param reliability Average reliability of constructs
#' @return Sample size
sem_mediation_sample_size_calculation <- function(r_a, r_b, power, alpha, reliability) {
  power_func <- function(n) sem_mediation_power_calculation(r_a, r_b, n, alpha, reliability) - power
  result <- tryCatch(stats::uniroot(power_func, interval = c(10, 100000))$root, error = function(e) NA)
  if (is.na(result)) stop("Could not find sample size. Effect may be too small after reliability adjustment.")
  return(ceiling(result))
}

#' SEM Mediation Effect Size Calculation
#' This internal function is used by the main wrapper to solve for an effect size.
#' @param target_path NEEDS DESCRIPTION
#' @param other_path_val NEEDS DESCRIPTION
#' @param n Sample Size
#' @param power Statistical power
#' @param alpha Significance level
#' @param reliability Average reliability of constructs
#' @return Effect size
sem_mediation_effect_size_calculation <- function(target_path, other_path_val, n, power, alpha, reliability) {
  power_func <- function(path_val) {
    if (target_path == "a") {
      power_val <- sem_mediation_power_calculation(path_val, other_path_val, n, alpha, reliability)
    } else {
      power_val <- sem_mediation_power_calculation(other_path_val, path_val, n, alpha, reliability)
    }
    return(power_val - power)
  }
  result <- tryCatch(stats::uniroot(power_func, interval = c(0.001, 0.999 / sqrt(reliability)))$root, error = function(e) NA)
  if (is.na(result)) stop("Could not find effect size. Power may be too high for given N and reliability.")
  return(result)
}

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (WRAPPER)
# ==============================================================================
#' Mediation SEM Power Analysis (Analytical Approximation)
#' @param r_a Latent path a (X->M) coefficient (NULL to calculate).
#' @param r_b Latent path b (M->Y|X) coefficient (NULL to calculate).
#' @param n Sample size (NULL to calculate).
#' @param power Statistical power (NULL to calculate, default = 0.8).
#' @param alpha Significance level (default = 0.05).
#' @param discount_factor Conservative discount factor (default = 0.75).
#' @param measurement_reliability Average reliability of constructs (default = 0.9).
#' @return Power analysis results.
#' @export
mediation_sem_power <- function(r_a = NULL, r_b = NULL, n = NULL, power = NULL,
                                alpha = 0.05, discount_factor = 0.75,
                                measurement_reliability = 0.9) {

  args_list <- list(r_a = r_a, r_b = r_b, n = n, power = power)
  null_args <-sapply(args_list, is.null)

  if (sum(null_args) != 1) {
    stop("Provide values for exactly three of the four parameters: r_a, r_b, n, power.")
  }

  target_param <- names(which(null_args))

  # Apply discount factor to any input latent path coefficients
  r_a_effective <- if (!is.null(r_a)) apply_discount_factor(r_a, discount_factor) else NULL
  r_b_effective <- if (!is.null(r_b)) apply_discount_factor(r_b, discount_factor) else NULL

  # Perform calculation based on which parameter is NULL
  if (target_param == "power") {
    result_val <- sem_mediation_power_calculation(r_a_effective, r_b_effective, n, alpha, measurement_reliability)
    result <- list(power = result_val, r_a = r_a_effective, r_b = r_b_effective, n = n)
  } else if (target_param == "n") {
    result_val <- sem_mediation_sample_size_calculation(r_a_effective, r_b_effective, power, alpha, measurement_reliability)
    result <- list(n = result_val, r_a = r_a_effective, r_b = r_b_effective, power = power)
  } else if (target_param == "r_a") {
    result_val <- sem_mediation_effect_size_calculation("a", r_b_effective, n, power, alpha, measurement_reliability)
    result <- list(r_a = result_val, r_b = r_b_effective, n = n, power = power)
  } else { # target_param == "r_b"
    result_val <- sem_mediation_effect_size_calculation("b", r_a_effective, n, power, alpha, measurement_reliability)
    result <- list(r_a = r_a_effective, r_b = result_val, n = n, power = power)
  }

  # Assemble final output, preserving original structure
  final_result <- c(
    list(
      method = "Mediation SEM Power Analysis (v5.1, Corrected Analytical Engine)",
      calculation_target = target_param,
      alpha = alpha,
      discount_factor = discount_factor,
      measurement_reliability = measurement_reliability
    ),
    result
  )
  final_result$indirect_effect <- final_result$r_a * final_result$r_b
  final_result$interpretation <- interpret_mediation_effect(final_result$r_a, final_result$r_b)

  class(final_result) <- "mediation_sem_power_analysis"
  return(final_result)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (Unchanged)
# ==============================================================================

#' Quick Mediation SEM Sample Size Calculation
#'
#' @param r_a Latent path a (X->M) coefficient
#' @param r_b Latent path b (M->Y|X) coefficient
#' @param power Statistical power
#' @param alpha Significance level
#' @param discount_factor Conservative discount factor
#' @param measurement_reliability Average reliability of constructs
#' @return Sample size
mediation_sem_sample_size <- function(r_a, r_b, power = 0.8, alpha = 0.05,
                                      discount_factor = 0.75, measurement_reliability = 0.9) {
  result <- mediation_sem_power(r_a = r_a, r_b = r_b, power = power,
                                alpha = alpha, discount_factor = discount_factor,
                                measurement_reliability = measurement_reliability)
  return(result$n)
}

#' Quick Mediation SEM Power Check
#'
#' @param r_a Latent path a (X->M) coefficient
#' @param r_b Latent path b (M->Y|X) coefficient
#' @param n Sample size
#' @param alpha Significance level
#' @param discount_factor Conservative discount factor
#' @param measurement_reliability Average reliability of constructs
#' @return Power
mediation_sem_power_check <- function(r_a, r_b, n, alpha = 0.05,
                                      discount_factor = 0.75, measurement_reliability = 0.9) {
  result <- mediation_sem_power(r_a = r_a, r_b = r_b, n = n,
                                alpha = alpha, discount_factor = discount_factor,
                                measurement_reliability = measurement_reliability)
  return(result$power)
}


# ==============================================================================
# HELPER AND PRINT FUNCTIONS (Unchanged)
# ==============================================================================

#' Interpret Mediation Effect
#'
#' @param r_a Latent path a (X->M) coefficient
#' @param r_b Latent path b (M->Y|X) coefficient
#'
#' @return Mediation Effect
interpret_mediation_effect <- function(r_a, r_b) {
  indirect_effect <- abs(r_a * r_b)
  if (indirect_effect < 0.01) return("Negligible mediation")
  if (indirect_effect < 0.09) return("Small mediation")
  if (indirect_effect < 0.25) return("Medium mediation")
  return("Large mediation")
}

#' Print Method for SEM Mediation Power Analysis
#' @param x SEM mediation power analysis result
#' @param ... Additional arguments
#' @export
print.mediation_sem_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  cat("Calculation Target:", x$calculation_target, "\n")
  cat("Engine: Analytical Approximation (Aroian test on attenuated paths)\n")
  cat("Discount Factor Applied to Inputs:", x$discount_factor, "\n\n")

  cat("--- Effective Latent Path Coefficients Used in Calculation ---\n")
  cat("NOTE: These values reflect inputs after the discount factor has been applied.\n\n")

  cat("  Latent A-path (X->M | r_a)      :", round(x$r_a, 4), "\n")
  cat("  Latent B-path (M->Y|X | r_b)    :", round(x$r_b, 4), "\n")
  cat("  Latent Indirect effect (a*b)   :", round(x$indirect_effect, 4), paste0("(", x$interpretation, ")\n\n"))

  cat("--- Results ---\n")
  cat("  Required sample size (N)       :", ceiling(x$n), "\n")
  cat("  Achieved power                 :", round(x$power, 3), "\n")
  cat("  Assumed Measurement Reliability:", x$measurement_reliability, "\n")
  cat("  Alpha level                    :", x$alpha, "\n\n")

  cat("NOTE: This is a fast approximation. For final grant proposals or\n")
  cat("      publications, a full Monte Carlo simulation is recommended.\n\n")
}
