# ==============================================================================
# CORE CALCULATION FUNCTIONS (DELEGATED TO pwr)
# ==============================================================================

#' Calculate Power for Wilcoxon Signed-Rank Test
#' @param r_partial Partial correlation
#' @param n Sample size
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Statistical power
#' @export
wilcoxon_power_calculation <- function(r_partial, n, alpha, two_tailed) {
  # Convert framework's r_partial to the Cohen's d required by pwr.t.test
  d <- partial_r_to_cohens_d(r_partial)
  alt <- ifelse(two_tailed, "two.sided", "greater")

  # Delegate to the paired t-test power function
  power <- pwr::pwr.t.test(n = n, d = d, sig.level = alpha, type = "paired", alternative = alt)$power
  return(power)
}

#' Calculate Sample Size for Wilcoxon Signed-Rank Test
#' @param r_partial Partial correlation
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required sample size
#' @export
wilcoxon_sample_size_calculation <- function(r_partial, power, alpha, two_tailed) {
  d <- partial_r_to_cohens_d(r_partial)
  alt <- ifelse(two_tailed, "two.sided", "greater")

  # Delegate to the paired t-test power function
  n <- pwr::pwr.t.test(d = d, power = power, sig.level = alpha, type = "paired", alternative = alt)$n
  return(ceiling(n))
}

#' Calculate Effect Size for Wilcoxon Signed-Rank Test
#' @param n Sample size
#' @param power Target power
#' @param alpha Significance level
#' @param two_tailed Two-tailed test flag
#' @return Required partial correlation
#' @export
wilcoxon_effect_size_calculation <- function(n, power, alpha, two_tailed) {
  alt <- ifelse(two_tailed, "two.sided", "greater")

  # Find the required Cohen's d using pwr
  d_target <- pwr::pwr.t.test(n = n, power = power, sig.level = alpha, type = "paired", alternative = alt)$d

  # Convert the result back to the framework's r_partial
  r_target <- cohens_d_to_partial_r(d_target)
  return(r_target)
}

# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION (WRAPPER)
# ==============================================================================

#' Wilcoxon Signed-Rank Power Analysis with Framework Integration (v3.1)
#' @param r_partial Partial correlation effect size (NULL to calculate)
#' @param n Sample size (NULL to calculate)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param two_tailed Two-tailed test (default = TRUE)
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d")
#' @return Power analysis results with framework integration
#' @export
wilcoxon_signed_rank_power <- function(r_partial = NULL, n = NULL, power = NULL,
                                       alpha = 0.05, discount_factor = 0.75,
                                       two_tailed = TRUE,
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
    n = !missing(n) && !is.null(n),
    power = !missing(power) && !is.null(power) && !power_default_used
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n, power (or use effect_input)")
  }

  if (!provided_params["power"]) {
    result_val <- wilcoxon_power_calculation(r_partial, n, alpha, two_tailed)
    result <- list(power = result_val, calculation_target = "power")
  } else if (!provided_params["n"]) {
    result_val <- wilcoxon_sample_size_calculation(r_partial, power, alpha, two_tailed)
    result <- list(n = result_val, calculation_target = "sample_size")
  } else {
    result_val <- wilcoxon_effect_size_calculation(n, power, alpha, two_tailed)
    result <- list(r_partial = result_val, calculation_target = "effect_size")
  }

  final_result <- list(
    analysis_type = paste("wilcoxon_signed_rank", result$calculation_target, sep="_"),
    method = "Wilcoxon Signed-Rank Power Analysis (v3.1, Paired t-test Engine)",
    r_partial = if(is.null(r_partial)) result$r_partial else r_partial,
    n = if(is.null(n)) result$n else n,
    power = if(is.null(power) || power_default_used) result$power else power,
    alpha = alpha,
    discount_factor = discount_factor,
    two_tailed = two_tailed,
    calculation_target = result$calculation_target
  )

  final_result$interpretation <- interpret_effect_size(final_result$r_partial)
  final_result$effect_size_conversions <- list(
    cohens_d = partial_r_to_cohens_d(final_result$r_partial)
  )

  class(final_result) <- "wilcoxon_signed_rank_power_analysis"
  return(final_result)
}

# ==============================================================================
# CONVENIENCE FUNCTIONS (Restored in v3.1)
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
wilcoxon_framework_power <- function(effect_size, effect_type = "r", n = NULL,
                                     power = NULL, alpha = 0.05,
                                     discount_factor = 0.75) {
  wilcoxon_signed_rank_power(effect_input = effect_size, effect_type = effect_type,
                             n = n, power = power, alpha = alpha, discount_factor = discount_factor)
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

  cat("Effect size (r_partial):", round(x$r_partial, 4), paste0(" (", x$interpretation, ")"), "\n")
  cat("Sample size (N of pairs):", x$n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Alpha level:", x$alpha, "\n")
  cat("Test type:", if(x$two_tailed) "Two-tailed" else "One-tailed", "\n\n")

  cat("NOTE: Power is calculated using the parametric equivalent (paired t-test).\n")
  cat("      This provides a standard and robust estimate for study planning.\n\n")

  if (!is.null(x$effect_size_conversions)) {
    cat("Framework Conversions:\n")
    cat("  Equivalent Cohen's d:", round(x$effect_size_conversions$cohens_d, 3), "\n\n")
  }

  cat("Framework Details:\n")
  cat("  Discount factor applied:", x$discount_factor, "(to initial effect_input)\n")
  cat("  Calculation target:", x$calculation_target, "\n\n")
}
