#' Rigorously Corrected Function to Generate Mediation Data
#'
#' This function generates simulated data (X, M, Y) for a pure mediation model (X -> M -> Y).
#' It is constructed to ensure the *population* parameters precisely match the inputs:
#' - `cor(X, M) = r_a`
#' - `partial_cor(M, Y, controlling_for=X) = r_b`
#' - `partial_cor(X, Y, controlling_for=M) = 0` (pure mediation, c_prime = 0)
#' This is achieved by deriving the necessary zero-order correlation matrix.
#'
#' @param n Sample size
#' @param r_a Standardized path a coefficient (X -> M), i.e., cor(X,M)
#' @param r_b Standardized path b coefficient (M -> Y | X), i.e., partial_cor(M,Y,X)
#' @return Data frame with X, M, Y variables (standardized)
#' @export
generate_mediation_data_corrected <- function(n, r_a, r_b) {
  # Validate inputs
  r_a <- validate_partial_r(r_a, allow_zero = TRUE, context = "for r_a")
  r_b <- validate_partial_r(r_b, allow_zero = TRUE, context = "for r_b")

  # Mathematical Derivation of the Zero-Order Correlation Matrix
  # We need to find the zero-order correlations (r_xm, r_xy, r_my) that satisfy the
  # population conditions for pure mediation (c_prime = 0).
  # 1. r_xm = r_a (by definition)
  # 2. From c_prime=0 (partial cor(X,Y|M)=0), we know r_xy = r_xm * r_my.
  #    So, r_xy = r_a * r_my.
  # 3. From partial cor(M,Y|X)=r_b, we solve for r_my:
  #    The formula is: (r_my - r_mx * r_xy) / sqrt((1-r_mx^2)(1-r_xy^2)) = r_b
  #    Substituting (1) and (2) gives:
  #    (r_my - r_a * (r_a * r_my)) / sqrt((1-r_a^2)(1-(r_a*r_my)^2)) = r_b
  #    Solving this equation for r_my yields:
  #    r_my = r_b / sqrt(1 - r_a^2 * (1 - r_b^2))

  # Calculate the required zero-order correlations
  r_xm <- r_a

  # Handle edge case where denominator might be zero or negative
  denom_r_my <- 1 - r_a^2 * (1 - r_b^2)
  if (denom_r_my <= 0) {
    stop(paste("Invalid correlation combination for r_a=", r_a, "and r_b=", r_b))
  }
  r_my <- r_b / sqrt(denom_r_my)

  r_xy <- r_a * r_my

  # Construct the population correlation matrix (Sigma)
  sigma <- matrix(c(
    1.0,  r_xm, r_xy,
    r_xm, 1.0,  r_my,
    r_xy, r_my, 1.0
  ), nrow = 3, byrow = TRUE)

  # Validate the final matrix before generation
  if (any(abs(c(r_xm, r_xy, r_my)) > 1.0)) {
    # This can happen in extreme cases if floating point errors push a value just over 1
    # It indicates an impossible correlation structure.
    stop(paste("Cannot generate data: derived correlation matrix is invalid for r_a=", r_a, ", r_b=", r_b))
  }

  # Generate multivariate normal data using Cholesky decomposition
  L <- tryCatch(chol(sigma), error = function(e) {
    stop(paste("Cannot generate data: correlation matrix is not positive-definite for r_a=", r_a, ", r_b=", r_b, ". Error:", e$message))
  })

  Z <- matrix(stats::rnorm(n * 3), ncol = 3)
  data_sim <- Z %*% L

  # Name and scale columns for robustness
  colnames(data_sim) <- c("X", "M", "Y")
  data.frame(
    X = scale(data_sim[,1])[,1],
    M = scale(data_sim[,2])[,1],
    Y = scale(data_sim[,3])[,1]
  )
}
