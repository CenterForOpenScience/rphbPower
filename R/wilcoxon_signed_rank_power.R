# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Wilcoxon Signed-Rank Power Analysis with Framework Integration
#' @param r_partial Partial correlation effect size (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param two_tailed Two-tailed test (default = TRUE)
#' @param asymptotic Use asymptotic approximation (default = TRUE)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
wilcoxon_signed_rank_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                       alpha = 0.05, discount_factor = 0.75,
                                       two_tailed = TRUE, asymptotic = TRUE,
                                       effect_input = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = TRUE)
  }

  # Count non-NULL parameters BEFORE applying defaults
  params_provided <- sum(c(!is.null(r_partial), !is.null(n), !is.null(power)))

  if (params_provided != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  # Apply default power value AFTER parameter checking
  if (is.null(power)) power <- 0.8

  if (params_provided != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  # Apply default power value AFTER parameter checking
  if (is.null(power)) power <- 0.8

  # Validate inputs using framework functions
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for Wilcoxon signed-rank power analysis")
  }

  if (!is.null(n)) {
    if (!is.numeric(n) || any(n != round(n)) || any(n < 6)) {
      stop("Sample size must be whole number >= 6 for Wilcoxon signed-rank test")
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

  if (!is.logical(two_tailed)) {
    stop("two_tailed must be TRUE or FALSE")
  }

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- wilcoxon_power_calculation(r_partial, n, alpha, two_tailed, asymptotic)

    result <- list(
      analysis_type = "wilcoxon_power",
      method = "Wilcoxon Signed-Rank Power Analysis",
      r_partial = r_partial,
      n = n,
      power = calculated_power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      asymptotic = asymptotic,
      calculation_target = "power"
    )

  } else if (is.null(n)) {
    # Calculate sample size
    calculated_n <- wilcoxon_sample_size_calculation(r_partial, power, alpha, two_tailed, asymptotic)

    result <- list(
      analysis_type = "wilcoxon_sample_size",
      method = "Wilcoxon Signed-Rank Sample Size Analysis",
      r_partial = r_partial,
      n = calculated_n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      asymptotic = asymptotic,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- wilcoxon_effect_size_calculation(n, power, alpha, two_tailed, asymptotic)

    result <- list(
      analysis_type = "wilcoxon_effect_size",
      method = "Wilcoxon Signed-Rank Effect Size Analysis",
      r_partial = calculated_r_partial,
      n = n,
      power = power,
      alpha = alpha,
      discount_factor = discount_factor,
      two_tailed = two_tailed,
      asymptotic = asymptotic,
      calculation_target = "effect_size"
    )
  }

  # Add framework integration
  r_for_conversion <- result$r_partial
  result$effect_size_conversions <- framework_conversion_summary(
    r_for_conversion, "r", apply_discount = FALSE
  )
  result$interpretation <- interpret_effect_size(r_for_conversion)
  result$are_factor <- calculate_are_factor()
  result$nonparametric_conversions <- wilcoxon_specific_conversions(r_for_conversion)

  class(result) <- "wilcoxon_signed_rank_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Wilcoxon Signed-Rank Test
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @param asymptotic Use asymptotic approximation
#' @return Statistical power
#' @export
wilcoxon_power_calculation <- function(r_partial, n, alpha, two_tailed, asymptotic) {
  # Convert partial correlation to Cohen's d for Wilcoxon calculations
  cohens_d <- partial_r_to_cohens_d(r_partial)

  if (asymptotic) {
    # Use asymptotic relative efficiency approach
    are_factor <- calculate_are_factor()
    effective_n <- n * are_factor

    # Calculate power using normal approximation
    se <- sqrt(1 / (effective_n - 3))
    z_calc <- r_partial / se

    # Critical value
    z_crit <- if (two_tailed) {
      qnorm(1 - alpha/2)
    } else {
      qnorm(1 - alpha)
    }

    # Power calculation
    if (two_tailed) {
      power <- 1 - pnorm(z_crit - abs(z_calc)) + pnorm(-z_crit - abs(z_calc))
    } else {
      power <- 1 - pnorm(z_crit - z_calc)
    }

  } else {
    # Exact calculation using Wilcoxon distribution properties
    power <- exact_wilcoxon_power(cohens_d, n, alpha, two_tailed)
  }

  return(pmax(0.05, pmin(1, power)))
}

#' Calculate Sample Size for Wilcoxon Signed-Rank Test
#' @param r_partial Partial correlation
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @param asymptotic Use asymptotic approximation
#' @return Required sample size
#' @export
wilcoxon_sample_size_calculation <- function(r_partial, power, alpha, two_tailed, asymptotic) {
  # Use iterative approach to find required sample size
  n_min <- 6
  n_max <- 2000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- wilcoxon_power_calculation(r_partial, n_test, alpha, two_tailed, asymptotic)

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

#' Calculate Effect Size for Wilcoxon Signed-Rank Test
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @param asymptotic Use asymptotic approximation
#' @return Required partial correlation
#' @export
wilcoxon_effect_size_calculation <- function(n, power, alpha, two_tailed, asymptotic) {
  # Use iterative approach to find required effect size
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- wilcoxon_power_calculation(r_test, n, alpha, two_tailed, asymptotic)

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
# WILCOXON-SPECIFIC CALCULATIONS
# ==============================================================================

#' Calculate Asymptotic Relative Efficiency Factor
#' @return ARE factor for Wilcoxon vs t-test
#' @export
calculate_are_factor <- function() {
  # Asymptotic relative efficiency of Wilcoxon signed-rank vs t-test
  return(3/pi)  # ≈ 0.955
}

#' Exact Wilcoxon Power Calculation
#' @param cohens_d Cohen's d effect size
#' @param n Sample size
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Statistical power
#' @export
exact_wilcoxon_power <- function(cohens_d, n, alpha, two_tailed) {
  # Convert Cohen's d to probability of favorable outcome
  p_favorable <- pnorm(cohens_d / sqrt(2))

  # Expected value and variance of Wilcoxon statistic under alternative
  expected_w <- n * (n + 1) / 4 * (2 * p_favorable - 1)
  var_w <- n * (n + 1) * (2*n + 1) / 24

  # Normal approximation with continuity correction
  critical_value <- if (two_tailed) {
    qnorm(1 - alpha/2) * sqrt(var_w)
  } else {
    qnorm(1 - alpha) * sqrt(var_w)
  }

  # Calculate power
  z_score <- abs(expected_w) / sqrt(var_w)
  power <- 1 - pnorm(critical_value/sqrt(var_w) - z_score)

  if (two_tailed) {
    power <- 2 * power - 1
    power <- pmax(0, power)
  }

  return(pmax(0.05, pmin(0.99, power)))
}

#' Wilcoxon-Specific Effect Size Conversions
#' @param r_partial Partial correlation
#' @return Named list of Wilcoxon-specific metrics
#' @export
wilcoxon_specific_conversions <- function(r_partial) {
  cohens_d <- partial_r_to_cohens_d(r_partial)

  # Nonparametric effect size metrics
  auc <- pnorm(cohens_d / sqrt(2))
  common_language_es <- auc
  rank_biserial_r <- 2 * auc - 1
  probability_superiority <- auc

  list(
    auc = auc,
    common_language_es = common_language_es,
    rank_biserial_r = rank_biserial_r,
    probability_superiority = probability_superiority
  )
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

#' Framework-Integrated Wilcoxon Power Analysis
#' @param effect_size Effect size value
#' @param effect_type Type of effect size
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param discount_factor Framework discount factor
#' @return Power analysis result
#' @export
wilcoxon_signed_rank_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                                 power = NULL, alpha = 0.05,
                                                 discount_factor = 0.75) {
  # Only pass power if explicitly provided, otherwise let it be calculated
  if (is.null(power)) {
    wilcoxon_signed_rank_power(effect_input = effect_size, effect_type = effect_type,
                               n = n, alpha = alpha, discount_factor = discount_factor)
  } else {
    wilcoxon_signed_rank_power(effect_input = effect_size, effect_type = effect_type,
                               n = n, power = power, alpha = alpha, discount_factor = discount_factor)
  }
}

#' Quick Wilcoxon Sample Size Calculation
#' @param r_partial Partial correlation
#' @param power Target power (default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @return Required sample size
#' @export
wilcoxon_sample_size <- function(r_partial, power = 0.8, alpha = 0.05) {
  result <- wilcoxon_signed_rank_power(r_partial = r_partial, power = power, alpha = alpha)
  return(result$n)
}

#' Quick Wilcoxon Power Calculation
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level (default = 0.05)
#' @return Statistical power
#' @export
wilcoxon_power_check <- function(r_partial, n, alpha = 0.05) {
  result <- wilcoxon_signed_rank_power(r_partial = r_partial, n = n, alpha = alpha)
  return(result$power)
}

# ==============================================================================
# PRINT METHOD
# ==============================================================================

#' Print Method for Wilcoxon Signed-Rank Power Analysis
#' @param x Wilcoxon power analysis result
#' @param ... Additional arguments
#' @export
print.wilcoxon_signed_rank_power_analysis <- function(x, ...) {
  cat("\n", x$method, "\n")
  cat(rep("=", nchar(x$method)), "\n\n")

  # Core results
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Sample size:", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Test type:", if(x$two_tailed) "Two-tailed" else "One-tailed", "\n")
  cat("Calculation method:", if(x$asymptotic) "Asymptotic" else "Exact", "\n")

  # Framework conversions
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat(paste0("  Cohen's f", cli::symbol$sup_2, ":"), round(conv$cohens_f2, 3), "\n")
    cat(paste0("  R", cli::symbol$sup_2, ":"), round(conv$r_squared, 3), "\n")
  }

  # Nonparametric conversions
  if (!is.null(x$nonparametric_conversions)) {
    cat("\nNonparametric conversions:\n")
    np_conv <- x$nonparametric_conversions
    cat("  Area Under Curve:", round(np_conv$auc, 3), "\n")
    cat("  Common Language ES:", round(np_conv$common_language_es, 3), "\n")
    cat("  Rank-biserial r:", round(np_conv$rank_biserial_r, 3), "\n")
    cat("  Probability superiority:", round(np_conv$probability_superiority, 3), "\n")
  }

  # Analysis details
  cat("\nFramework details:\n")
  cat("  ARE factor:", round(x$are_factor, 3), "\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")

  cat("\n")
}
