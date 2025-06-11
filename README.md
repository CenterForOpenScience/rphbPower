
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rphbPower - Unified Power Analysis Package

<!-- badges: start -->
<!-- badges: end -->

Comprehensive power analysis functions for regression-based statistical
methods using partial correlations as the standardized effect size
metric with integrated conservative planning.

## Framework Overview

This package provides power analysis for diverse statistical methods
unified under a **partial correlation framework**. All analyses convert
effect sizes to partial correlations, enabling direct comparison across
correlation, regression, mediation, SEM, multilevel, and longitudinal
approaches.

**Core Features:** - **Unified Effect Size Metric**: Partial
correlations across all analysis types - **Conservative Planning**:
Built-in 0.75 discount factor for realistic sample sizes  
- **Auto-Detection**: Provide any 2 of (effect size, sample size,
power) - calculates the third - **Framework Integration**: Conversion
between Cohen’s d, f², R², eta-squared - **Mathematical Precision**:
Accurate F-distribution calculations for precise power analysis

## Quick Start

### Basic Power Analysis Pattern

``` r
# Provide any two parameters - calculates the third
result <- linear_regression_power(r_partial = 0.25, n = 122, n_predictors = 1)          # → power
result <- linear_regression_power(r_partial = 0.25, power = 0.8, n_predictors = 1)      # → sample size  
result <- linear_regression_power(n = 192, power = 0.8, n_predictors = 1)               # → effect size

print(result)  # Comprehensive output with framework conversions
```

### Framework Effect Size Integration

``` r
# Direct conversion from any effect size type using effect_input
result <- linear_regression_power(
  effect_input = 0.4,          # Cohen's d from literature
  effect_type = "d",           # Automatic conversion to partial r
  power = 0.8,
  n_predictors = 3
)

# Multiple effect size types supported
linear_regression_power(effect_input = 0.15, effect_type = "f2", power = 0.8, n_predictors = 2)
linear_regression_power(effect_input = 0.09, effect_type = "r_squared", power = 0.8, n_predictors = 2)
linear_regression_power(effect_input = 0.25, effect_type = "r", power = 0.8, n_predictors = 2)
```

## Available Analysis Methods

| Analysis Type | Location | Key Function | Effect Size Input |
|----|----|----|----|
| **Correlation** | `05_methods/5.1_correlation/` | `correlation_power()` | r (zero-order correlation) |
| **Linear Regression** | `05_methods/5.2_regression/` | `linear_regression_power()` | r_partial |
| **Logistic Regression** | `05_methods/5.2_regression/` | `logistic_regression_power()` | r_partial or OR |
| **Cross-Lagged Panel** | `05_methods/5.3_longitudinal/` | `cross_lagged_panel_power()` | r_partial |
| **Fixed Effects** | `05_methods/5.3_longitudinal/` | `fixed_effects_power()` | r_partial |
| **Repeated Measures** | `05_methods/5.3_longitudinal/` | `repeated_measures_power()` | r_partial |
| **Mediation (Regression)** | `05_methods/5.4_mediation/` | `mediation_regression_power()` | r_a, r_b paths |
| **Mediation (SEM)** | `05_methods/5.4_mediation/` | `mediation_sem_power()` | r_a, r_b paths |
| **Multilevel Models** | `05_methods/5.5_multilevel/` | `mixed_models_power()` | r_partial |
| **SEM Direct Effects** | `05_methods/5.6_sem/` | `sem_direct_effects_power()` | r_partial |
| **Nonparametric** | `05_methods/5.7_nonparametric/` | `wilcoxon_signed_rank_power()` | r_partial |

## Method Selection Guide

### By Research Design

- **Single predictor, continuous outcome** → Linear regression
- **Multiple predictors, continuous outcome** → Linear regression
  (multiple)
- **Binary/categorical outcome** → Logistic regression
- **Mediation hypotheses** → Mediation regression or SEM
- **Clustered/nested data** → Multilevel models
- **Complex structural models** → SEM direct effects
- **Simple bivariate association** → Correlation

### By Data Characteristics

- **Normal, continuous variables** → Linear regression methods
- **Hierarchical/grouped data** → Multilevel models
- **Longitudinal/panel data** → Cross-lagged, fixed effects, repeated
  measures
- **Latent variable models** → SEM approaches

## Core Framework Functions

### Effect Size Conversions

``` r
# Direct effect size conversion using power analysis functions
result <- linear_regression_power(effect_input = 0.4, effect_type = "d", power = 0.8, n_predictors = 2)
r_partial <- result$r_partial  # Converted partial correlation

# Available effect types: "r", "d", "f2", "r_squared", "eta_squared"
load_analysis_method("correlation")
correlation_power(effect_input = 0.15, effect_type = "f2", power = 0.8)

load_analysis_method("mediation_regression")
mediation_regression_power(effect_input_a = 0.09, effect_type = "r_squared", r_b = 0.3, power = 0.8)
```

### Working with Results

``` r
# All power analysis functions return comprehensive results
result <- linear_regression_power(r_partial = 0.25, power = 0.8, n_predictors = 3)

# Access key components
print(result$n)                           # Required sample size
print(result$effect_size_conversions)     # Cohen's d, f², R²
print(result$interpretation)              # Effect size interpretation
print(result$discount_factor)             # Conservative planning factor
```

## Framework Quality Assurance

### Mathematical Foundation

The framework uses accurate F-distribution calculations with proper
non-centrality parameter computation (`ncp = f² × df_denominator`)
ensuring precise power analysis across all regression-based methods.

### Conservative Planning

- **Automatic Discount Factor**: 0.75 applied to all effect sizes for
  realistic planning
- **Cross-Method Validation**: Partial correlations enable direct
  comparison across analyses
- **Literature Integration**: Systematic conversion from diverse effect
  size metrics

### Sample Size Planning Guidelines

**Critical Planning Reality**: Model complexity significantly impacts
sample size requirements. Always specify actual predictor count rather
than using simple correlation estimates.

**Example Impact (r = 0.20, 80% power)**: - Single predictor model: ~192
participants - Five predictor model: ~320 participants (+67%) - Ten
predictor model: ~410 participants (+113%)

**Framework Recommendations**: 1. Use framework auto-detection for
precise calculations 2. Apply conservative discount factors (default
behavior) 3. Account for model complexity in planning 4. Validate effect
sizes through framework conversions

## Documentation Structure

### Core Documentation

- `01_README.md` - Package overview (this file)
- `02_setup.R` - Working setup system
- `03_Getting Started.md` - Tutorial guide
- `04_core/` - Framework utilities and mathematical functions

### Method Documentation

- `05_methods/` - Complete analysis methods (11 total)
- Each method includes: power analysis functions, examples, and
  documentation
- All methods integrate with unified framework architecture

### Comprehensive Guides

- `06_documentation/` - Advanced guides and reference materials
- Assumption validation, effect size guidelines, troubleshooting
- Quick reference with sample size estimates

## Framework Validation

✅ **Mathematical Accuracy**: All F-distribution calculations verified
and validated

✅ **Framework Integration**: Unified partial correlation approach
across 11 methods

✅ **Conservative Planning**: Automatic discount factor system for
realistic estimates

✅ **Cross-Method Consistency**: Direct comparison enabled through
standardized effect sizes

✅ **Documentation Quality**: All sample size calculations verified

## Getting Help

### Quick Reference

Use `06_documentation/quick_reference_guide.md` for fast lookup of
functions, parameters, and sample size estimates.

### Method Selection

Use `06_documentation/method_selection_guide.md` for guidance on
choosing appropriate analysis approaches.

### Troubleshooting

Use `06_documentation/troubleshooting.md` for common issues and
framework-specific solutions.

## Version Information

**Current Version**: 1.2

**Mathematical Status**: F-distribution calculations verified and
validated (v1.2)

**Documentation Status**: Sample size calculations verified

**Framework Status**: All 11 methods mathematically verified

The unified framework provides sophisticated power analysis with simple,
consistent patterns across diverse statistical methods while maintaining
mathematical precision and conservative planning principles.
