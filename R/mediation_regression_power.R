# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Mediation Regression Power Analysis with Framework Integration
#' @param r_a Path coefficient a (X->M) as correlation (NULL to calculate)
#' @param r_b Path coefficient b (M->Y|X) as partial correlation (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param test_type Type of mediation test ("sobel", "joint")
#' @param effect_input_a Raw effect size input for a-path (alternative to r_a)
#' @param effect_input_b Raw effect size input for b-path (alternative to r_b)
#' @param effect_type Type of effect inputs ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
mediation_regression_power <- function(r_a = NULL, r_b = NULL, n = NULL, power = 0.8,
                                       alpha = 0.05, discount_factor = 0.75, test_type = "sobel",
                                       effect_input_a = NULL, effect_input_b = NULL,
                                       effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input_a)) {
    r_a <- framework_effect_size(effect_input_a, effect_type, apply_discount = TRUE)
  }
  if (!is.null(effect_input_b)) {
    r_b <- framework_effect_size(effect_input_b, effect_type, apply_discount = TRUE)
  }

  # Use missing() to detect actually provided parameters
  provided_params <- c(
    r_a = !missing(r_a) && !is.null(r_a),
    r_b = !missing(r_b) && !is.null(r_b),
    n = !missing(n) && !is.null(n),
    power = !missing(power) && !is.null(power)
  )

  if (sum(provided_params) != 3) {
    stop("Provide exactly three of: r_a, r_b, n, power (or use effect_input parameters)")
  }

  # Validate inputs using framework functions
  if (!is.null(r_a)) {
    r_a <- validate_partial_r(r_a, allow_zero = FALSE,
                              context = "for a-path in mediation analysis")
  }

  if (!is.null(r_b)) {
    r_b <- validate_partial_r(r_b, allow_zero = FALSE,
                              context = "for b-path in mediation analysis")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < 50)) {
      stop("Sample size must be whole number >= 50 for mediation analysis")
    }
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1")
  }

  if (!test_type %in% c("sobel", "joint")) {
    stop("test_type must be 'sobel' or 'joint'")
  }

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- mediation_power_calculation(r_a, r_b, n, alpha, test_type)

    result <- list(
      analysis_type = "mediation_power",
      method = "Mediation Regression Power Analysis",
      r_a = r_a,
      r_b = r_b,
      n = n,
      power = calculated_power,
      alpha = alpha,
      discount_factor = discount_factor,
      test_type = test_type,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    calculated_n <- mediation_sample_size_calculation(r_a, r_b, power, alpha, test_type)

    result <- list(
      analysis_type = "mediation_sample_size",
      method = "Mediation Regression Sample Size Analysis",
      r_a = r_a,
      r_b = r_b,
      n = calculated_n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      test_type = test_type,
      calculation_target = "sample_size"
    )

  } else if (is.null(r_a)) {
    # Calculate required a-path
    calculated_r_a <- mediation_effect_size_calculation_a(r_b, n, power, alpha, test_type)

    result <- list(
      analysis_type = "mediation_effect_size_a",
      method = "Mediation Regression A-Path Analysis",
      r_a = calculated_r_a,
      r_b = r_b,
      n = n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      test_type = test_type,
      calculation_target = "effect_size_a"
    )

  } else {
    # Calculate required b-path
    calculated_r_b <- mediation_effect_size_calculation_b(r_a, n, power, alpha, test_type)

    result <- list(
      analysis_type = "mediation_effect_size_b",
      method = "Mediation Regression B-Path Analysis",
      r_a = r_a,
      r_b = calculated_r_b,
      n = n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      test_type = test_type,
      calculation_target = "effect_size_b"
    )
  }

  # Add framework integration
  result$indirect_effect <- result$r_a * result$r_b
  result$a_path_conversions <- framework_conversion_summary(
    result$r_a, "r", apply_discount = FALSE
  )
  result$b_path_conversions <- framework_conversion_summary(
    result$r_b, "r", apply_discount = FALSE
  )
  result$interpretation <- interpret_mediation_effect(result$r_a, result$r_b)

  class(result) <- "mediation_regression_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Mediation Analysis
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @param test_type Type of mediation test
#' @return Statistical power
#' @export
mediation_power_calculation <- function(r_a, r_b, n, alpha, test_type) {
  # Calculate indirect effect
  indirect_effect <- r_a * r_b

  if (test_type == "sobel") {
    # Sobel test standard error
    se_ab <- sqrt(r_a^2 * (1 - r_b^2)/(n - 3) + r_b^2 * (1 - r_a^2)/(n - 2))

    # Z-score for indirect effect
    z_score <- abs(indirect_effect) / se_ab

    # Power calculation
    power <- 1 - pnorm(qnorm(1 - alpha/2) - z_score) + pnorm(-qnorm(1 - alpha/2) - z_score)

  } else if (test_type == "joint") {
    # Joint significance test (test both paths)
    # Power for a-path
    se_a <- sqrt((1 - r_a^2) / (n - 2))
    t_a <- abs(r_a) / se_a
    power_a <- 1 - pt(qt(1 - alpha/2, n - 2), n - 2, ncp = t_a) +
      pt(-qt(1 - alpha/2, n - 2), n - 2, ncp = t_a)

    # Power for b-path
    se_b <- sqrt((1 - r_b^2) / (n - 3))
    t_b <- abs(r_b) / se_b
    power_b <- 1 - pt(qt(1 - alpha/2, n - 3), n - 3, ncp = t_b) +
      pt(-qt(1 - alpha/2, n - 3), n - 3, ncp = t_b)

    # Joint power (both must be significant)
    power <- power_a * power_b
  }

  return(pmax(0, pmin(1, power)))
}

#' Calculate Sample Size for Mediation Analysis
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param power Target power
#' @param alpha Significance level
#' @param test_type Type of mediation test
#' @return Required sample size
#' @export
mediation_sample_size_calculation <- function(r_a, r_b, power, alpha, test_type) {
  # Use iterative approach to find required sample size
  n_min <- 50
  n_max <- 5000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- mediation_power_calculation(r_a, r_b, n_test, alpha, test_type)

    if (abs(power_test - power) < tolerance) {
      return(n_test)
    }

    if (power_test < power) {
      n_min <- n_test
    } else {
      n_max <- n_test
    }

    if (n_max - n_min <= 1) break
  }

  return(n_max)
}

#' Calculate Required A-Path Effect Size
#' @param r_b B-path partial correlation
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param test_type Type of mediation test
#' @return Required a-path correlation
#' @export
mediation_effect_size_calculation_a <- function(r_b, n, power, alpha, test_type) {
  # Use iterative approach to find required a-path
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- mediation_power_calculation(r_test, r_b, n, alpha, test_type)

    if (abs(power_test - power) < tolerance) {
      return(r_test)
    }

    if (power_test < power) {
      r_min <- r_test
    } else {
      r_max <- r_test
    }

    if (r_max - r_min <= 0.001) break
  }

  return(r_max)
}

#' Calculate Required B-Path Effect Size
#' @param r_a A-path correlation
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param test_type Type of mediation test
#' @return Required b-path partial correlation
#' @export
mediation_effect_size_calculation_b <- function(r_a, n, power, alpha, test_type) {
  # Use iterative approach to find required b-path
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- mediation_power_calculation(r_a, r_test, n, alpha, test_type)

    if (abs(power_test - power) < tolerance) {
      return(r_test)
    }

    if (power_test < power) {
      r_min <- r_test
    } else {
      r_max <- r_test
    }

    if (r_max - r_min <= 0.001) break
  }

  return(r_max)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated Mediation Power Analysis
#' @param effect_size_a A-path effect size value
#' @param effect_size_b B-path effect size value
#' @param effect_type Type of effect sizes
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
#' @export
mediation_framework_power <- function(effect_size_a, effect_size_b, effect_type = "r",
                                      n = NULL, power = 0.8, alpha = 0.05,
                                      discount_factor = 0.75) {
  mediation_regression_power(effect_input_a = effect_size_a, effect_input_b = effect_size_b,
                             effect_type = effect_type, n = n, power = power,
                             alpha = alpha, discount_factor = discount_factor)
}

#' Quick Mediation Sample Size Calculation
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param power Target power (default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
#' @export
mediation_sample_size <- function(r_a, r_b, power = 0.8, alpha = 0.05) {
  result <- mediation_regression_power(r_a = r_a, r_b = r_b, power = power, alpha = alpha)
  return(result$n)
}

#' Quick Mediation Power Calculation
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
mediation_power_check <- function(r_a, r_b, n, alpha = 0.05) {
  result <- mediation_regression_power(r_a = r_a, r_b = r_b, n = n, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Interpret Mediation Effect Size
#' @param r_a A-path correlation
#' @param r_b B-path partial correlation
#' @return Interpretation string
#' @export
interpret_mediation_effect <- function(r_a, r_b) {
  indirect_effect <- abs(r_a * r_b)

  if (indirect_effect < 0.01) return("Negligible mediation effect")
  if (indirect_effect < 0.09) return("Small mediation effect")
  if (indirect_effect < 0.25) return("Medium mediation effect")
  return("Large mediation effect")
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Mediation Regression Power Analysis
#' @param x Mediation regression power analysis result
#' @param ... Additional arguments
#' @export
print.mediation_regression_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat(paste0("A-path (X", symbol$arrow_right, "M):"), round(x$r_a, 4), "\n")
  cat(paste0("B-path (M", symbol$arrow_right, "Y|X):"), round(x$r_b, 4), "\n")
  cat("Indirect effect (a x b):", round(x$indirect_effect, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Test type:", x$test_type, "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Framework conversions for paths
  if (!is.null(x$a_path_conversions)) {
    cat("\nA-path conversions:\n")
    conv_a <- x$a_path_conversions
    cat("  Cohen's d:", round(conv_a$cohens_d, 3), "\n")
    cat(paste0("  R", cli::symbol$sup_2, ":"), round(conv_a$r_squared, 3), "\n")
  }

  if (!is.null(x$b_path_conversions)) {
    cat("\nB-path conversions:\n")
    conv_b <- x$b_path_conversions
    cat("  Cohen's d:", round(conv_b$cohens_d, 3), "\n")
    cat(paste0("  R", cli::symbol$sup_2, ":"), round(conv_b$r_squared, 3), "\n")
  }

  # Analysis details
  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")

  cat("\n")
}
