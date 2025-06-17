
# rphbPower - Unified Statistical Power Analysis Framework Package

A comprehensive R framework for a priori statistical power analysis. It
provides easy-to-use functions for 11 different statistical methods, all
unified by a common effect size metric and a philosophy of conservative
study planning.

## Framework Overview

This package provides power analysis for a diverse suite of statistical
methods unified under a **partial correlation framework**. All analyses
convert common effect sizes (e.g., Cohen’s d, f²) into partial
correlations (`r_partial`), enabling direct and intuitive comparisons of
effect magnitude across all included methods.

The framework is built on two core principles:

1.  **Unified Effect Size**: Using `r_partial` as the common metric
    simplifies interpretation and promotes a deeper understanding of
    effect sizes, regardless of the statistical test being used.

2.  **Conservative Planning**: A built-in **0.75 discount factor** is
    automatically applied to all user-provided effect sizes. This
    encourages prudent and realistic study planning, helping to prevent
    underpowered research that can result from overly optimistic effect
    size estimations.

## Core Features

- **11 Validated Modules**: Covers a wide range of common statistical
  tests
- **Unified Effect Size Metric**: All calculations are based on partial
  correlations (`r_partial`)
- **Parameter Auto-Detection**: For any given test, provide any two of
  (effect size, sample size, power) and the framework will solve for the
  third
- **Built-in Effect Size Conversion**: Automatically converts from
  Cohen’s d, f², R², and η² into `r_partial`
- **Conservative by Default**: A 0.75 discount factor is automatically
  applied to effect size inputs to encourage robust planning

## Quick Start

### 1. Install and Load the Package

The package can be installed directly from the Center For Open Science
GitHub repository.

``` r
# install.packages("devtools") # Run this if you don't have devtools
devtools::install_github("CenterForOpenScience/rphbPower", build_vignettes = TRUE)

library(rphbPower)
```

### 2. Run the Power Analysis

The functions follow a consistent pattern. Provide the parameters you
know, and leave the one you want to find as NULL.

``` r
# Example 1: Solve for required sample size (N)
result_n <- linear_regression_power(r_partial = 0.25, power = 0.8, n_predictors = 2)
print(result_n$n)

# Example 2: Solve for the power of a planned study
result_power <- linear_regression_power(r_partial = 0.25, n = 121, n_predictors = 2)
print(result_power$power)
```

### 3. Use Any Effect Size

The `effect_input` and `effect_type` arguments allow you to use effect
sizes directly from the literature. The 0.75 discount factor will be
applied automatically.

``` r
# Use a Cohen's d of 0.5 from a previous study
result <- linear_regression_power(
  effect_input = 0.5,
  effect_type = "d",
  power = 0.8,
  n_predictors = 2
)
print(result)
```

## Available Analysis Methods

### Basic and Regression Methods

- **Correlation**: `correlation_power()`
- **Linear Regression**: `linear_regression_power()`
- **Logistic Regression**: `logistic_regression_power()`

### Longitudinal Methods

- **Cross-Lagged Panel**: `cross_lagged_panel_power()`
- **Fixed Effects**: `fixed_effects_power()`
- **Repeated Measures**: `repeated_measures_power()`

### Mediation Methods

- **Mediation (Regression)**: `mediation_regression_power()`
- **Mediation (SEM)**: `mediation_sem_power()`

### Advanced Methods

- **Multilevel Models**: `mixed_models_power()`
- **SEM (Direct Effects)**: `sem_direct_effects_power()`
- **Wilcoxon Signed-Rank Test**: `wilcoxon_signed_rank_power()`

## Method Selection Guide

### By Research Design

- **Single predictor, continuous outcome** → Linear regression
- **Multiple predictors, continuous outcome** → Linear regression
  (multiple)
- **Binary/categorical outcome** → Logistic regression
- **Mediation hypotheses** → Mediation regression or SEM
- **Clustered/nested data** → Multilevel models
- **Complex structural models with latent vars** → SEM direct effects
- **Simple bivariate association** → Correlation

### By Data Characteristics

- **Normal, continuous variables** → Linear regression methods
- **Hierarchical/grouped data** → Multilevel models
- **Longitudinal/panel data** → Cross-lagged, fixed effects, repeated
  measures
- **Latent variable models** → SEM approaches

## Sample Size Planning Guidelines

A critical reality in study planning is that model complexity
significantly impacts sample size requirements. With the corrected
engine, the impact of adding covariates is now estimated more
accurately.

**Example Impact** (Effect Size r = 0.20, Target Power = 80%):

- **Correlation (1 predictor)**: ~192 participants
- **Regression (5 predictors)**: ~192 participants (+0% increase)
- **Regression (10 predictors)**: ~197 participants (+3% increase)

This framework helps you account for this complexity directly in your
power analysis.

## Framework Status

- **Current Version**: 2.2
- **Mathematical Status**: All 11 modules have been re-engineered and
  have passed comprehensive validation tests. The framework is
  considered complete and robust.
- **Validation**: 100% pass rate achieved across all modules against
  external benchmarks or internal consistency checks.

## Documentation Structure

For detailed tutorials and guides, refer to the vignettes accessed
through `browseVignettes("rphbPower")`:

- **Getting Started**: Installation and basic usage
- **Effect Size Guidelines**: Guidance on choosing an appropriate effect
  size
- **Method Selection Guide**: In-depth help for choosing the right
  analysis
- **Quick Reference Guide**: A fast lookup for functions and parameters
- **Troubleshooting**: Solutions for common issues
