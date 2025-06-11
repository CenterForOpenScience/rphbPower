# ==============================================================================
# MAIN POWER ANALYSIS FUNCTION
# ==============================================================================

#' Mixed Models Power Analysis with Framework Integration
#' @param r_partial Partial correlation coefficient (NULL to calculate)
#' @param n_groups Number of groups/clusters (NULL to calculate)
#' @param n_per_group Average observations per group (default = 10)
#' @param power Statistical power (NULL to calculate, default = 0.8)
#' @param icc Intraclass correlation coefficient (default = 0.05)
#' @param alpha Significance level (default = 0.05)
#' @param discount_factor Conservative discount factor (default = 0.75)
#' @param test_level Level being tested ("level1" or "level2")
#' @param effect_input Raw effect size input (alternative to r_partial)
#' @param effect_type Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")
#' @return Power analysis results with framework integration
#' @export
mixed_models_power <- function(r_partial = NULL, n_groups = NULL, n_per_group = 10,
                               power = 0.8, icc = 0.05, alpha = 0.05,
                               discount_factor = 0.75, test_level = "level1",
                               effect_input = NULL, effect_type = "r") {

  # Parameter detection and validation
  if (!is.null(effect_input)) {
    r_partial <- framework_effect_size(effect_input, effect_type, apply_discount = TRUE)
  }

  # Use missing() to detect actually provided parameters
  provided_params <- c(
    r_partial = !missing(r_partial) && !is.null(r_partial),
    n_groups = !missing(n_groups) && !is.null(n_groups),
    power = !missing(power) && !is.null(power)
  )

  if (sum(provided_params) != 2) {
    stop("Provide exactly two of: r_partial, n_groups, power (or use effect_input)")
  }

  # Validate inputs using framework functions
  if (!is.null(r_partial)) {
    r_partial <- validate_partial_r(r_partial, allow_zero = FALSE,
                                    context = "for mixed models power analysis")
  }

  if (!is.null(n_groups)) {
    if (!is.numeric(n_groups) || any(n_groups != round(n_groups)) || any(n_groups < 3)) {
      stop("Number of groups must be whole number >= 3")
    }
  }

  if (!is.numeric(n_per_group) || n_per_group < 2 || n_per_group != round(n_per_group)) {
    stop("Observations per group must be whole number >= 2")
  }

  if (!is.null(power)) {
    if (!is.numeric(power) || any(power <= 0) || any(power >= 0.999)) {
      stop("Power must be between 0 and 0.999")
    }
  }

  if (!is.numeric(icc) || icc < 0 || icc >= 1) {
    stop("ICC must be between 0 and 1")
  }

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1")
  }

  if (!test_level %in% c("level1", "level2")) {
    stop("test_level must be 'level1' or 'level2'")
  }

  # Calculate design effect and effective sample size
  total_n <- if (!is.null(n_groups)) n_groups * n_per_group else NULL

  # Perform power calculation
  if (is.null(power)) {
    # Calculate power
    calculated_power <- mixed_models_power_calculation(r_partial, n_groups, n_per_group,
                                                       icc, alpha, test_level)

    result <- list(
      analysis_type = "mixed_models_power",
      method = "Mixed Models Power Analysis",
      r_partial = r_partial,
      n_groups = n_groups,
      n_per_group = n_per_group,
      total_n = total_n,
      power = calculated_power,
      icc = icc,
      test_level = test_level,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "power"
    )

  } else if (is.null(n_groups)) {
    # Calculate number of groups
    calculated_n_groups <- mixed_models_sample_size_calculation(r_partial, power, n_per_group,
                                                                icc, alpha, test_level)

    result <- list(
      analysis_type = "mixed_models_sample_size",
      method = "Mixed Models Sample Size Analysis",
      r_partial = r_partial,
      n_groups = calculated_n_groups,
      n_per_group = n_per_group,
      total_n = calculated_n_groups * n_per_group,
      power = power,
      icc = icc,
      test_level = test_level,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "sample_size"
    )

  } else {
    # Calculate effect size
    calculated_r_partial <- mixed_models_effect_size_calculation(n_groups, power, n_per_group,
                                                                 icc, alpha, test_level)

    result <- list(
      analysis_type = "mixed_models_effect_size",
      method = "Mixed Models Effect Size Analysis",
      r_partial = calculated_r_partial,
      n_groups = n_groups,
      n_per_group = n_per_group,
      total_n = total_n,
      power = power,
      icc = icc,
      test_level = test_level,
      alpha = alpha,
      discount_factor = discount_factor,
      calculation_target = "effect_size"
    )
  }

  # Add framework integration
  r_for_conversion <- result$r_partial
  result$effect_size_conversions <- framework_conversion_summary(
    r_for_conversion, "r", apply_discount = FALSE
  )
  result$interpretation <- interpret_effect_size(r_for_conversion)

  # Add mixed models specific information
  result$design_effect <- 1 + (result$n_per_group - 1) * result$icc
  result$effective_n <- effective_sample_size(result$total_n, design_effect = result$design_effect)

  class(result) <- "mixed_models_power_analysis"
  return(result)
}

# ==============================================================================
# CORE CALCULATION FUNCTIONS
# ==============================================================================

#' Calculate Power for Mixed Models
#' @param r_partial Partial correlation
#' @param n_groups Number of groups
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param alpha Significance level
#' @param test_level Level being tested
#' @return Statistical power
#' @export
mixed_models_power_calculation <- function(r_partial, n_groups, n_per_group, icc, alpha, test_level) {
  # Convert partial correlation to F-statistic approach
  f2 <- partial_r_to_cohens_f2(r_partial)

  if (test_level == "level1") {
    # Level-1 (individual-level) effect
    total_n <- n_groups * n_per_group
    design_effect <- 1 + (n_per_group - 1) * icc
    effective_n <- total_n / design_effect

    # Degrees of freedom
    df_num <- 1
    df_den <- effective_n - 2

  } else {
    # Level-2 (group-level) effect
    effective_n <- n_groups
    df_num <- 1
    df_den <- n_groups - 2
  }

  if (df_den <= 0) {
    stop("Insufficient degrees of freedom for analysis")
  }

  # Non-centrality parameter
  ncp <- f2 * df_den

  # Critical F-value
  f_crit <- qf(1 - alpha, df_num, df_den)

  # Calculate power using non-central F distribution
  power <- 1 - pf(f_crit, df_num, df_den, ncp = ncp)

  return(pmax(0, pmin(1, power)))
}

#' Calculate Number of Groups for Mixed Models
#' @param r_partial Partial correlation
#' @param power Target power
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param alpha Significance level
#' @param test_level Level being tested
#' @return Required number of groups
#' @export
mixed_models_sample_size_calculation <- function(r_partial, power, n_per_group, icc, alpha, test_level) {
  # Use iterative approach to find required number of groups
  n_min <- 3
  n_max <- 1000
  tolerance <- 0.001

  for (i in 1:100) {
    n_test <- round((n_min + n_max) / 2)
    power_test <- mixed_models_power_calculation(r_partial, n_test, n_per_group,
                                                 icc, alpha, test_level)

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

#' Calculate Effect Size for Mixed Models
#' @param n_groups Number of groups
#' @param power Target power
#' @param n_per_group Observations per group
#' @param icc Intraclass correlation
#' @param alpha Significance level
#' @param test_level Level being tested
#' @return Required partial correlation
#' @export
mixed_models_effect_size_calculation <- function(n_groups, power, n_per_group, icc, alpha, test_level) {
  # Use iterative approach to find required effect size
  r_min <- 0.01
  r_max <- 0.95
  tolerance <- 0.001

  for (i in 1:100) {
    r_test <- (r_min + r_max) / 2
    power_test <- mixed_models_power_calculation(r_test, n_groups, n_per_group,
                                                 icc, alpha, test_level)

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
                                         power = 0.8, n_per_group = 10, icc = 0.05,
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
#' @return Required number of groups
#' @export
mixed_models_sample_size <- function(r_partial, power = 0.8, n_per_group = 10,
                                     icc = 0.05, test_level = "level1") {
  result <- mixed_models_power(r_partial = r_partial, power = power,
                               n_per_group = n_per_group, icc = icc, test_level = test_level)
  return(result$n_groups)
}

#' Quick Mixed Models Power Calculation
#' @param r_partial Partial correlation
#' @param n_groups Number of groups
#' @param n_per_group Observations per group (default = 10)
#' @param icc Intraclass correlation (default = 0.05)
#' @param test_level Level being tested (default = "level1")
#' @return Statistical power
#' @export
mixed_models_power_check <- function(r_partial, n_groups, n_per_group = 10,
                                     icc = 0.05, test_level = "level1") {
  result <- mixed_models_power(r_partial = r_partial, n_groups = n_groups,
                               n_per_group = n_per_group, icc = icc, test_level = test_level)
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

  # Core results
  cat("Partial correlation:", round(x$r_partial, 4),
      paste0("(", x$interpretation, ")"), "\n")
  cat("Number of groups:", x$n_groups, "\n")
  cat("Observations per group:", x$n_per_group, "\n")
  cat("Total sample size:", x$total_n, "\n")
  cat("Effective sample size:", x$effective_n, "\n")
  cat("Statistical power:", round(x$power, 3), "\n")
  cat("Test level:", x$test_level, "\n")
  cat("Alpha level:", x$alpha, "\n")

  # Design details
  cat("\nDesign characteristics:\n")
  cat("  ICC:", round(x$icc, 3), "\n")
  cat("  Design effect:", round(x$design_effect, 2), "\n")

  # Framework conversions
  if (!is.null(x$effect_size_conversions)) {
    cat("\nFramework conversions:\n")
    conv <- x$effect_size_conversions
    cat("  Cohen's d:", round(conv$cohens_d, 3), "\n")
    cat(paste0("  Cohen's f", cli::symbol$sup_2, ":"), round(conv$cohens_f2, 3), "\n")
    cat(paste0("  R", cli::symbol$sup_2, ":"), round(conv$r_squared, 3), "\n")
  }

  # Analysis details
  cat("\nFramework details:\n")
  cat("  Discount factor:", x$discount_factor, "\n")
  cat("  Calculation target:", x$calculation_target, "\n")

  cat("\n")
}
