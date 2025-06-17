# ==============================================================================
# FILE: mixed_models_power.R
# STATUS: RE-ENGINEERED & DOCUMENTED (v6.1)
# ==============================================================================
#
# v6.1: Adds complete Roxygen documentation for all functions to resolve
#       `devtools::check()` warnings.
#
# ==============================================================================

# Ensure required packages are available for the new engine
if (!require("lme4", quiet = TRUE)) install.packages("lme4")
if (!require("pwr", quiet = TRUE)) install.packages("pwr")

# --- Internal Simulation Engine for Level-1 Effects (Corrected Beta) ---
mixed_model_level1_sim_engine <- function(r_partial, n_groups, n_per_group, icc, alpha, n_sims) {
  significant_results <- 0
  total_n <- n_groups * n_per_group

  # Mathematical Derivation of Beta from r_partial
  if (abs(r_partial) >= 1) stop("r_partial must be between -1 and 1.")
  beta_sq <- (r_partial^2 * (1 - icc)) / (1 - r_partial^2)
  beta <- sqrt(beta_sq) * sign(r_partial)

  cat(paste0("Running Level-1 simulation (", n_sims, " reps)... "))
  pb <- utils::txtProgressBar(min = 0, max = n_sims, style = 3)

  for (i in 1:n_sims) {
    # Data Generation using correct variance components
    sd_intercepts <- sqrt(icc)
    sd_residuals <- sqrt(1 - icc)

    group_id <- rep(1:n_groups, each = n_per_group)
    random_intercepts <- stats::rnorm(n_groups, 0, sd_intercepts)
    x_level1 <- stats::rnorm(total_n, 0, 1)
    error <- stats::rnorm(total_n, 0, sd_residuals)

    y <- random_intercepts[group_id] + beta * x_level1 + error
    sim_data <- data.frame(Y = y, X = x_level1, Group = as.factor(group_id))

    # Model Fitting and P-Value Extraction
    p_value <- tryCatch({
      model <- suppressMessages(lme4::lmer(Y ~ X + (1 | Group), data = sim_data))
      summary(model)$coefficients["X", "Pr(>|t|)"]
    }, error = function(e) { 1.0 })

    if (!is.na(p_value) && p_value < alpha) {
      significant_results <- significant_results + 1
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  cat(" Done.\n")

  return(significant_results / n_sims)
}

#' @title Mixed Models Power Analysis (Final)
#' @description Performs power analysis for mixed-effects (multilevel) models for both Level-1 (within-group) and Level-2 (between-group) effects.
#' @param r_partial Partial correlation for the effect of interest (NULL to calculate).
#' @param n_groups Number of groups (Level-2 units) (NULL to calculate).
#' @param n_per_group Average number of observations per group (Level-1 units).
#' @param power Target statistical power (NULL to calculate).
#' @param icc Intraclass Correlation Coefficient.
#' @param alpha Significance level.
#' @param discount_factor Conservative planning discount factor.
#' @param test_level The level of the predictor: "level1" or "level2".
#' @param n_sims_l1 Number of simulations for Level-1 power estimation.
#' @param effect_input Alternative raw effect size.
#' @param effect_type Type of `effect_input`.
#' @return A list with power analysis results.
#' @export
mixed_models_power <- function(r_partial = NULL, n_groups = NULL, n_per_group = 10,
                               power = NULL, icc = 0.05, alpha = 0.05,
                               discount_factor = 0.75, test_level = "level1",
                               n_sims_l1 = 500, effect_input = NULL, effect_type = "r") {

  # I. Parameter Validation & Setup
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = (discount_factor != 1.0))
  }

  if (is.null(power)) {
    power_default_used <- TRUE
    power_in <- 0.8
  } else {
    power_default_used <- FALSE
    power_in <- power
  }

  provided_params <- c(
    r_partial = !is.null(r_partial),
    n_groups = !is.null(n_groups),
    power = !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n_groups, power (or use effect_input)")
  }
  target_param <- names(which(!provided_params))

  # II. Core Calculation
  if (test_level == "level1") {
    if (target_param != "power") stop("For Level-1 effects, the simulation engine can only solve for power.")

    calculated_power <- mixed_model_level1_sim_engine(r_partial, n_groups, n_per_group, icc, alpha, n_sims_l1)
    result <- list(power = calculated_power, calculation_target = "power")

  } else if (test_level == "level2") {
    f2 <- if(!is.null(r_partial)) partial_r_to_cohens_f2(r_partial) else NULL
    power_to_pass <- if(target_param == "power") NULL else power_in

    pwr_result <- pwr::pwr.f2.test(u = 1, v = if(!is.null(n_groups)) n_groups - 2 else NULL, f2 = f2, sig.level = alpha, power = power_to_pass)

    if (target_param == "n_groups") {
      result <- list(n_groups = ceiling(pwr_result$v + 2), calculation_target = "sample_size")
    } else if (target_param == "power") {
      result <- list(power = pwr_result$power, calculation_target = "power")
    } else {
      result <- list(r_partial = cohens_f2_to_partial_r(pwr_result$f2), calculation_target = "effect_size")
    }
  } else {
    stop("`test_level` must be either 'level1' or 'level2'.")
  }

  # III. Assemble Output
  final_output <- list(
    method = "Mixed Models Power Analysis",
    engine = if(test_level == "level1") "Monte Carlo Simulation" else "Analytical (pwr)",
    calculation_target = result$calculation_target,
    r_partial = if(target_param == "r_partial") result$r_partial else r_partial,
    n_groups = if(target_param == "n_groups") result$n_groups else n_groups,
    n_per_group = n_per_group,
    power = if(target_param == "power") result$power else power_in,
    icc = icc,
    test_level = test_level,
    alpha = alpha,
    discount_factor = discount_factor
  )
  final_output$n_total <- final_output$n_groups * final_output$n_per_group
  class(final_output) <- "mixed_models_power_analysis"
  return(final_output)
}


#' Quick Mixed Models Sample Size
#' @description Convenience function to calculate the number of groups for a Level-2 effect.
#' @param r_partial Partial correlation for the Level-2 predictor.
#' @param power Target power (default = 0.8).
#' @param alpha Significance level (default = 0.05).
#' @param n_per_group Number of observations per group (for context in some calculations).
#' @param icc Intraclass correlation (for context).
#' @param test_level The level of the effect, must be "level2" for this function.
#' @return Required number of groups.
#' @export
mixed_models_sample_size <- function(r_partial, power = 0.8, alpha = 0.05, n_per_group = 10, icc = 0.05, test_level = "level2") {
  if(test_level != "level2") stop("This convenience function only supports solving for N for Level-2 effects.")
  result <- mixed_models_power(r_partial = r_partial, power = power, alpha = alpha, n_per_group = n_per_group, icc = icc, test_level = test_level, n_groups = NULL)
  return(result$n_groups)
}

#' Quick Mixed Models Power Check
#' @description Convenience function to calculate power.
#' @param r_partial Partial correlation for the predictor.
#' @param n_groups Number of groups.
#' @param alpha Significance level (default = 0.05).
#' @param n_per_group Number of observations per group (default = 10).
#' @param icc Intraclass correlation (default = 0.05).
#' @param test_level The level of the effect, "level1" or "level2".
#' @return Achieved statistical power.
#' @export
mixed_models_power_check <- function(r_partial, n_groups, alpha = 0.05, n_per_group = 10, icc = 0.05, test_level = "level1") {
  mixed_models_power(r_partial = r_partial, n_groups = n_groups, alpha = alpha, n_per_group = n_per_group, icc = icc, test_level = test_level, power = NULL)$power
}


#' Print method for mixed models power analysis
#' @param x A `mixed_models_power_analysis` object.
#' @param ... Additional arguments passed to `print`.
#' @export
print.mixed_models_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  cat("Test Level:", x$test_level, "\n")
  cat("Engine:", x$engine, "\n\n")

  cat("Partial correlation:", round(x$r_partial, 4), "\n")
  cat("Number of groups (L2):", x$n_groups, "\n")
  cat("Observations per group (L1):", x$n_per_group, "\n")
  cat("Total sample size:", x$n_total, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Intraclass Correlation (ICC):", round(x$icc, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")

  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}
