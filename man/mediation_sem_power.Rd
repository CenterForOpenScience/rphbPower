% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediation_sem_power.R
\name{mediation_sem_power}
\alias{mediation_sem_power}
\title{Mediation SEM Power Analysis (Analytical Approximation)}
\usage{
mediation_sem_power(
  r_a = NULL,
  r_b = NULL,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  discount_factor = 0.75,
  measurement_reliability = 0.9
)
}
\arguments{
\item{r_a}{Latent path a (X->M) coefficient (NULL to calculate).}

\item{r_b}{Latent path b (M->Y|X) coefficient (NULL to calculate).}

\item{n}{Sample size (NULL to calculate).}

\item{power}{Statistical power (NULL to calculate, default = 0.8).}

\item{alpha}{Significance level (default = 0.05).}

\item{discount_factor}{Conservative discount factor (default = 0.75).}

\item{measurement_reliability}{Average reliability of constructs (default = 0.9).}
}
\value{
Power analysis results.
}
\description{
Mediation SEM Power Analysis (Analytical Approximation)
}
