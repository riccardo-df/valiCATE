#' Correction Term for CWATE
#'
#' Computes the correction term C_hat_i for each observation, which enters the CWATE influence function and accounts for the
#' estimation error in the weight function.
#'
#' @param Gamma Numeric vector of AIPW pseudo-outcomes.
#' @param tau_hat Numeric vector of CATE predictions.
#' @param wt_info List from \code{\link{compute_weights}} with elements \code{omega}, \code{v}, and \code{nabla}.
#' @param spec Weight specification from \code{\link{get_weight_spec}}.
#'
#' @return
#' Numeric vector of correction terms C_hat_i.
#'
#' @details
#' For CDF-based weights (AUTOC, AUC-HVL, QINI), the correction term involves a sum over all observations that is computed
#' in O(n log n) via sorting and cumulative sums, rather than the naive O(n^2) double loop.\cr
#'
#' For BLP, the correction simplifies to an O(n) vectorized expression.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{correction_ncwate}}, \code{\link{estimate_cwate}}
#'
#' @keywords internal
correction_cwate <- function(Gamma, tau_hat, wt_info, spec) {
  n <- length(Gamma)

  if (!spec$cdf_based) {
    return(-mean(Gamma) * (tau_hat - mean(tau_hat)))
  }

  ## CDF-based weights: O(n log n) computation.
  a <- Gamma * wt_info$nabla
  ord <- order(tau_hat)

  # Prefix sum of a in sorted order.
  a_sorted <- a[ord]
  cumsum_a <- cumsum(a_sorted)

  # For each observation i, sum_{j: tau_hat_j <= tau_hat_i} a_j.
  rnk <- rank(tau_hat, ties.method = "max")
  prefix_sum_i <- cumsum_a[rnk]

  # Constant term: sum_j a_j * F_hat_j.
  const <- sum(a * wt_info$v)

  ## Output.
  return((prefix_sum_i - const) / n)
}


#' Correction Term for NCWATE
#'
#' Computes the correction term C^gamma_hat_i for each observation, which enters the NCWATE influence function and accounts for
#' the estimation error in the weight function.
#'
#' @param Gamma Numeric vector of AIPW pseudo-outcomes.
#' @param tau_hat Numeric vector of CATE predictions.
#' @param gamma_hat Scalar, the NCWATE point estimate.
#' @param wt_info List from \code{\link{compute_weights}} with elements \code{omega}, \code{v}, and \code{nabla}.
#' @param spec Weight specification from \code{\link{get_weight_spec}}.
#'
#' @return
#' Numeric vector of correction terms C^gamma_hat_i.
#'
#' @details
#' For CDF-based weights (AUTOC, AUC-HVL, QINI), the correction term involves a sum over all observations that is computed
#' in O(n log n) via sorting and cumulative sums, rather than the naive O(n^2) double loop.\cr
#'
#' For BLP, the correction simplifies to an O(n) vectorized expression.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{correction_cwate}}, \code{\link{estimate_ncwate}}
#'
#' @keywords internal
correction_ncwate <- function(Gamma, tau_hat, gamma_hat, wt_info, spec) {
  n <- length(Gamma)
  resid <- Gamma - gamma_hat * tau_hat

  if (!spec$cdf_based) {
    ## BLP: C^gamma_hat_i = -(tau_hat_i - tau_bar) * mean(Gamma - gamma_hat * tau_hat).
    return(-mean(resid) * (tau_hat - mean(tau_hat)))
  }

  ## CDF-based weights: O(n log n) computation.
  a <- resid * wt_info$nabla
  ord <- order(tau_hat)

  a_sorted <- a[ord]
  cumsum_a <- cumsum(a_sorted)

  rnk <- rank(tau_hat, ties.method = "max")
  prefix_sum_i <- cumsum_a[rnk]

  const <- sum(a * wt_info$v)

  ## Output.
  return((prefix_sum_i - const) / n)
}
