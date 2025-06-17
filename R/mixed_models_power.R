# ==============================================================================
# CORE CALCULATION FUNCTIONS (DELEGATED TO pwr and WebPower)
# ==============================================================================

#' Calculate Power for Mixed Models
#' @param r_partial Partial correlation
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param test_level Level being tested
#' @param alpha Significance level
#' @return Statistical power
mixed_models_power_calculation <- function(r_partial, n_groups, n_per_group, icc, test_level, alpha) {
  f2 <- partial_r_to_cohens_f2(r_partial)

  if (test_level == "level1") {
    # Level 1 (within-group) effects are calculated like a cluster-randomized trial
    power <- WebPower::wp.crt2arm(n = n_per_group, J = n_groups, f = sqrt(f2), icc = icc, alpha = alpha)$power
  } else { # test_level == "level2"
    # Level 2 (between-group) effects are like a simple regression on group means
    power <- pwr::pwr.f2.test(u = 1, v = n_groups - 2, f2 = f2, sig.level = alpha)$power
  }
  return(power)
}

#' Calculate Number of Groups for Mixed Models
#' @param r_partial Partial correlation
#' @param power Target power
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param test_level Level being tested
#' @param alpha Significance level
#' @return Required number of groups
#' @export
mixed_models_sample_size_calculation <- function(r_partial, power, n_per_group, icc, test_level, alpha) {
  f2 <- partial_r_to_cohens_f2(r_partial)

  if (test_level == "level1") {
    n_groups <- WebPower::wp.crt2arm(n = n_per_group, f = sqrt(f2), icc = icc, power = power, alpha = alpha)$J
  } else { # test_level == "level2"
    v <- pwr::pwr.f2.test(u = 1, f2 = f2, power = power, sig.level = alpha)$v
    n_groups <- v + 2
  }
  return(ceiling(n_groups))
}

#' Calculate Effect Size for Mixed Models
#' @param n_groups Number of groups
#' @param power Target power
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param test_level Level being tested
#' @param alpha Significance level
#' @return Required partial correlation
#' @export
mixed_models_effect_size_calculation <- function(n_groups, power, n_per_group, icc, test_level, alpha) {
  if (test_level == "level1") {
    # WebPower solves for f, not f-squared
    f_target <- tryCatch(
      WebPower::wp.crt2arm(n = n_per_group, J = n_groups, power = power, icc = icc, alpha = alpha)$f,
      error = function(e) NA
    )
    if(is.na(f_target)) return(NA)
    f2_target <- f_target^2
  } else { # test_level == "level2"
    f2_target <- tryCatch(
      pwr::pwr.f2.test(u = 1, v = n_groups - 2, power = power, sig.level = alpha)$f2,
      error = function(e) NA
    )
    if(is.na(f2_target)) return(NA)
  }

  # Convert f-squared back to the framework's r_partial
  r_target <- sqrt(f2_target / (1 + f2_target))
  return(r_target)
}

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (WRAPPER)
# ==============================================================================
#' Mixed Models Power Analysis - Corrected Formulas
#'
#' Performs power analysis for mixed models, calculating either partial
#' correlation, number of groups, or statistical power based on provided inputs.
#' This version delegates calculations to validated external packages.
#'
#' @param r_partial Partial correlation coefficient (NULL to calculate).
#' @param n_groups Number of groups/clusters (NULL to calculate).
#' @param n_per_group Average observations per group (default = 10).
#' @param power Statistical power (NULL to calculate, default applied when needed = 0.8).
#' @param icc Intraclass correlation coefficient (default = 0.05).
#' @param alpha Significance level (default = 0.05).
#' @param discount_factor Conservative discount factor (default = 0.75).
#' @param test_level Level being tested ("level1" or "level2").
#' @param effect_input Raw effect size input (alternative to r_partial).
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared").
#' @return A list containing power analysis results and framework integrations.
#' @export
mixed_models_power <- function(r_partial = NULL, n_groups = NULL, n_per_group = 10,
                               power = NULL, icc = 0.05, alpha = 0.05,
                               discount_factor = 0.75, test_level = "level1",
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
    n_groups = !missing(n_groups) && !is.null(n_groups),
    power = !missing(power) && !is.null(power) && !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n_groups, power (or use effect_input).")
  }

  if (!provided_params["power"]) {
    result_val <- mixed_models_power_calculation(r_partial, n_groups, n_per_group, icc, test_level, alpha)
    result <- list(power = result_val, calculation_target = "power")
  } else if (!provided_params["n_groups"]) {
    result_val <- mixed_models_sample_size_calculation(r_partial, power, n_per_group, icc, test_level, alpha)
    result <- list(n_groups = result_val, calculation_target = "sample_size")
  } else {
    result_val <- mixed_models_effect_size_calculation(n_groups, power, n_per_group, icc, test_level, alpha)
    result <- list(r_partial = result_val, calculation_target = "effect_size")
  }

  final_result <- list(
    analysis_type = paste("mixed_models", result$calculation_target, test_level, sep="_"),
    method = "Mixed Models Power Analysis (v4.3, pwr/WebPower Engine)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n_groups = if(is.null(n_groups)) result$n_groups else n_groups,
    n_per_group = n_per_group,
    power = if(is.null(power) || power_default_used) result$power else power,
    icc = icc,
    test_level = test_level,
    alpha = alpha,
    discount_factor = discount_factor,
    calculation_target = result$calculation_target
  )

  class(final_result) <- "mixed_models_power_analysis"
  return(final_result)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (RESTORED)
# ==============================================================================

#' Framework-Integrated Mixed Models Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n_groups Number of groups
#' @param power Target power
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param test_level Level being tested
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
#' @export
mixed_models_framework_power <- function(effect_size, effect_type = "r", n_groups = NULL,
                                         power = NULL, n_per_group = 10, icc = 0.05,
                                         test_level = "level1", alpha = 0.05,
                                         discount_factor = 0.75) {
  mixed_models_power(effect_input = effect_size, effect_type = effect_type,
                     n_groups = n_groups, power = power, n_per_group = n_per_group,
                     icc = icc, test_level = test_level, alpha = alpha,
                     discount_factor = discount_factor)
}

#' Quick Mixed Models Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param n_per_group Observations per group (default = 10)
#' @param icc Intraclass correlation (default = 0.05)
#' @param test_level Level being tested (default = "level1")
#' @param alpha Significance level
#' @return Required number of groups
#' @export
mixed_models_sample_size <- function(r_partial, power = 0.8, n_per_group = 10,
                                     icc = 0.05, test_level = "level1", alpha = 0.05) {
  result <- mixed_models_power(r_partial = r_partial, power = power,
                               n_per_group = n_per_group, icc = icc, test_level = test_level,
                               alpha = alpha)
  return(result$n_groups)
}

#' Quick Mixed Models Power Calculation
#' @param r_partial Partial correlation
#' @param n_groups Number of groups
#' @param n_per_group Observations per group (default = 10)
#' @param icc Intraclass correlation (default = 0.05)
#' @param test_level Level being tested (default = "level1")
#' @param alpha Significance level
#' @return Statistical power
#' @export
mixed_models_power_check <- function(r_partial, n_groups, n_per_group = 10,
                                     icc = 0.05, test_level = "level1", alpha = 0.05) {
  result <- mixed_models_power(r_partial = r_partial, n_groups = n_groups,
                               n_per_group = n_per_group, icc = icc, test_level = test_level,
                               alpha = alpha)
  return(result$power)
}


# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Mixed Models Power Analysis
#' @param x Mixed models power analysis result
#' @param ... Additional arguments
#' @export
print.mixed_models_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  cat("Test Level:", x$test_level, "\n")
  cat("Partial correlation:", round(x$r_partial, 4), "\n")
  cat("Number of groups (L2):", x$n_groups, "\n")
  cat("Observations per group (L1):", x$n_per_group, "\n")
  cat("Total sample size:", x$n_groups * x$n_per_group, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Intraclass Correlation (ICC):", round(x$icc, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")

  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, " (applied to initial effect_input)\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}
