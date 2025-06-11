#' Sample Size Adjustment for Design Effects
#' @param n_nominal Nominal sample size
#' @param design_effect Design effect multiplier (default: 1.0)
#' @param icc Intracluster correlation (optional)
#' @param cluster_size Average cluster size (optional)
#' @return Effective sample size
#' @export
effective_sample_size <- function(n_nominal, design_effect = 1.0, icc = NULL, cluster_size = NULL) {
  if (!is.numeric(n_nominal) || any(n_nominal <= 0)) {
    stop("Nominal sample size must be positive")
  }
  if (!is.numeric(design_effect) || any(design_effect < 1)) {
    stop("Design effect must be >= 1")
  }

  if (!is.null(icc) && !is.null(cluster_size)) {
    if (!is.numeric(icc) || !is.numeric(cluster_size)) {
      stop("ICC and cluster size must be numeric")
    }
    if (any(icc < 0) || any(icc >= 1)) {
      stop("ICC must be between 0 and 1")
    }
    if (any(cluster_size <= 1)) {
      stop("Cluster size must be > 1")
    }

    # Calculate design effect from ICC and cluster size
    design_effect <- 1 + (cluster_size - 1) * icc
  }

  n_effective <- n_nominal / design_effect
  return(round(n_effective))
}
