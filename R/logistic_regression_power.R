# ==============================================================================
# FILE: logistic_regression_power.R
# STATUS: RE-ENGINEERED WITH SIMULATION ENGINE (v5.2)
# ==============================================================================
#
# v5.2: Corrected the parameter validation logic. The check for the number
#       of provided arguments now correctly accounts for 'effect_input'
#       before throwing an error. This fixes the vignette build failures.
#
# v5.1: Restored the discount_factor functionality.
#
# v5.0: Replaced the simple analytical engine with a full Monte Carlo
#       simulation engine to properly account for n_predictors.
#
# ==============================================================================

#' Internal Simulation Engine for Multiple Logistic Regression (with correlated predictors)
#' @keywords internal
logistic_multi_predictor_sim_engine <- function(r_partial, n, n_predictors, alpha, n_sims, inter_predictor_cor = 0.3) {
  d_target <- partial_r_to_cohens_d(r_partial)
  beta_target <- d_target * pi / sqrt(3)
  significant_results <- 0

  # Handle the n_predictors = 1 edge case first
  if (n_predictors == 1) {
    inter_predictor_cor <- 0
  }

  # Build the correlation matrix for all predictors
  cor_matrix <- matrix(inter_predictor_cor, nrow = n_predictors, ncol = n_predictors)
  diag(cor_matrix) <- 1

  # Cholesky decomposition to create correlated variables
  chol_matrix <- chol(cor_matrix)

  for (i in 1:n_sims) {
    # Generate independent normal variables first
    Z <- matrix(rnorm(n * n_predictors), ncol = n_predictors)
    # Multiply by the Cholesky matrix to induce correlation
    X_correlated <- Z %*% chol_matrix

    # Assume the effect of nuisance predictors is 0 for a standard power calculation
    betas <- c(beta_target, rep(0, n_predictors - 1))
    log_odds <- X_correlated %*% betas
    prob <- 1 / (1 + exp(-log_odds))
    Y <- stats::rbinom(n, 1, prob)

    sim_data <- as.data.frame(cbind(Y, X_correlated))
    model <- suppressWarnings(stats::glm(Y ~ ., data = sim_data, family = "binomial"))

    # Extract p-value for the target predictor (V2, which is the first column of X_correlated)
    p_value <- summary(model)$coefficients[2, "Pr(>|z|)"]

    if (!is.na(p_value) && p_value < alpha) {
      significant_results <- significant_results + 1
    }
  }
  return(significant_results / n_sims)
}


#' Logistic Regression Power Analysis (Simulation-Based)
#' @param r_partial Partial correlation coefficient for the predictor of interest.
#' @param n Sample size.
#' @param power Statistical power.
#' @param n_predictors Total number of predictors in the model.
#' @param alpha The significance level.
#' @param discount_factor The conservative planning discount factor.
#' @param n_sims Number of simulations to run.
#' @param effect_input Alternative to r_partial for automatic conversion.
#' @param effect_type The type of `effect_input`.
#' @param inter_predictor_cor The estimated correlation between predictors (default = 0).
#' @return A list object containing the power analysis results.
#' @importFrom stats rnorm
#' @export
logistic_regression_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                      n_predictors = 1, alpha = 0.05, discount_factor = 0.75,
                                      n_sims = 500, effect_input = NULL, effect_type = "r",
                                      inter_predictor_cor = 0) { # <-- VALUE CHANGED TO 0

  # --- Parameter Validation (remains the same) ---
  effect_provided <- !is.null(r_partial) || !is.null(effect_input)
  provided_args <- c(effect_provided, !is.null(n), !is.null(power))
  if (sum(provided_args) != 2) {
    stop("Provide exactly two of: r_partial (or effect_input), n, power")
  }
  params_list <- list(r_partial = effect_provided, n = !is.null(n), power = !is.null(power))
  target_param <- names(which(!unlist(params_list)))

  # --- Effect Size Setup (remains the same) ---
  r_initial <- NULL
  if (!is.null(effect_input)) {
    r_initial <- framework_effect_size(effect_input, effect_type, apply_discount = FALSE)
  } else if (!is.null(r_partial)) {
    r_initial <- r_partial
  }
  r_for_calc <- if(!is.null(r_initial)) apply_discount_factor(r_initial, discount_factor) else NULL

  # --- Core Calculation (pass new argument to engine) ---
  cat(paste("Running simulation-based power analysis for", target_param, "...\n"))
  if (target_param == "power") {
    calculated_power <- logistic_multi_predictor_sim_engine(r_for_calc, n, n_predictors, alpha, n_sims, inter_predictor_cor)
    result <- list(power = calculated_power)
  } else {
    power_goal <- power
    power_diff_func <- function(val) {
      current_r <- if (target_param == "r_partial") val else r_for_calc
      current_n <- if (target_param == "n") ceiling(val) else n
      current_power <- logistic_multi_predictor_sim_engine(current_r, current_n, n_predictors, alpha, 100, inter_predictor_cor)
      return(current_power - power_goal)
    }
    search_interval <- if(target_param == "n") c(n_predictors + 5, 100000) else c(0.01, 0.99)
    solution <- tryCatch(
      stats::uniroot(power_diff_func, interval = search_interval),
      error = function(e) stop("Could not find a solution. The effect size may be too small or power too high.")
    )
    result <- list(value = solution$root)
  }

  # --- Assemble Output (add new parameter) ---
  final_output <- list(
    method = "Logistic Regression Power Analysis (Simulation Engine)",
    engine_details = paste(n_sims, "simulations per run"),
    calculation_target = target_param,
    r_partial = if (target_param == "r_partial") result$value else r_initial,
    n = if (target_param == "n") ceiling(result$value) else n,
    power = if (target_param == "power") result$power else power,
    n_predictors = n_predictors,
    alpha = alpha,
    discount_factor = discount_factor,
    inter_predictor_cor = inter_predictor_cor # <-- ADD THIS TO OUTPUT
  )
  final_output$interpretation <- interpret_effect_size(final_output$r_partial)
  final_output$effect_size_conversions <- framework_conversion_summary(
    final_output$r_partial, "r", apply_discount=FALSE
  )

  class(final_output) <- "logistic_regression_power_analysis"
  cat("Done.\n")
  return(final_output)
}


#' Quick Logistic Regression Sample Size Calculation
#' @param r_partial Partial correlation for the predictor of interest.
#' @param power Target statistical power (defaults to 0.8).
#' @param n_predictors Total number of predictors in the model.
#' @param ... Other arguments passed to `logistic_regression_power`, such as `effect_input`, `effect_type`, `alpha`, or `inter_predictor_cor`.
#' @return Required sample size (N).
#' @export
logistic_regression_sample_size <- function(r_partial = NULL, power = 0.8, n_predictors = 1, ...) {
  # Capture all arguments into a list
  args <- list(r_partial = r_partial, power = power, n_predictors = n_predictors, ...)
  # Ensure 'n' is the target parameter for this convenience function
  args$n <- NULL
  # Call the main power function with the complete set of arguments
  do.call(logistic_regression_power, args)$n
}


#' Quick Logistic Regression Power Calculation
#' @param r_partial Partial correlation for the predictor of interest.
#' @param n The sample size.
#' @param n_predictors Total number of predictors in the model.
#' @param ... Other arguments passed to `logistic_regression_power`, such as `effect_input`, `effect_type`, `alpha`, or `inter_predictor_cor`.
#' @return Achieved statistical power.
#' @export
logistic_regression_power_check <- function(r_partial = NULL, n = NULL, n_predictors = 1, ...) {
  # Capture all arguments into a list
  args <- list(r_partial = r_partial, n = n, n_predictors = n_predictors, ...)
  # Ensure 'power' is the target parameter for this convenience function
  args$power <- NULL
  # Call the main power function with the complete set of arguments
  do.call(logistic_regression_power, args)$power
}

#' Print Method for Logistic Regression Power Analysis
#' @param x Logistic regression power analysis result
#' @param ... Additional arguments
#' @export
print.logistic_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")
  cat("Input Partial correlation (pre-discount):", round(x$r_partial, 4),
      paste0(" (", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Number of predictors:", x$n_predictors, "\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Engine:", x$engine_details, "\n")

  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat("  Cohen's f\u00B2:", round(conv$cohens_f2, 3), "\n")
  }
  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")
  cat("\n")
}
