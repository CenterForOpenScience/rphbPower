% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_power.R
\name{correlation_power}
\alias{correlation_power}
\title{Correlation Power Analysis with Framework Integration (v2.2)}
\usage{
correlation_power(
  r_partial = NULL,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  discount_factor = 0.75,
  two_tailed = TRUE,
  effect_input = NULL,
  effect_type = "r"
)
}
\arguments{
\item{r_partial}{Partial correlation coefficient (NULL to calculate)}

\item{n}{Sample size (NULL to calculate)}

\item{power}{Statistical power (NULL to calculate, default = 0.8)}

\item{alpha}{Significance level (default = 0.05)}

\item{discount_factor}{Conservative discount factor (default = 0.75)}

\item{two_tailed}{Two-tailed test (default = TRUE)}

\item{effect_input}{Raw effect size input (alternative to r_partial)}

\item{effect_type}{Type of effect_input ("r", "d", "f2", "r_squared", "eta_squared")}
}
\value{
Power analysis results with framework integration
}
\description{
Correlation Power Analysis with Framework Integration (v2.2)
}
