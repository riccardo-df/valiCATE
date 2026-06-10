#' CWATE Estimation
#'
#' Estimates the Centered-Weighted Average Treatment Effect (CWATE) for a given weight function.
#'
#' @param Gamma Numeric vector of AIPW pseudo-outcomes.
#' @param tau_hat Numeric vector of CATE predictions.
#' @param weight_name Character string, one of \code{"AUTOC"}, \code{"AUC-HVL"}, \code{"BLP"}, \code{"QINI"}.
#' @param alpha Significance level for the confidence interval.
#'
#' @return
#' A list with elements \code{estimate}, \code{se}, \code{ci} (named vector with \code{lower} and \code{upper}),
#' \code{test_stat}, and \code{p_value}.
#'
#' @details
#' The CWATE is estimated as theta_hat = mean(omega_hat * Gamma), where omega_hat denotes the weight function evaluated at the
#' CATE predictions and Gamma denotes the AIPW pseudo-outcomes.\cr
#'
#' A one-sided test of H0: theta <= 0 is reported. Rejection signals that the predicted CATEs capture genuine heterogeneity
#' in the correct direction.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{estimate_ncwate}}, \code{\link{valiCATE}}
#'
#' @keywords internal
estimate_cwate <- function(Gamma, tau_hat, weight_name, alpha = 0.05) {
  ## 1.) Compute weights.
  n <- length(Gamma)
  spec <- get_weight_spec(weight_name)
  wt_info <- compute_weights(tau_hat, spec)

  ## 2.) Point estimate: theta_hat = mean(omega * Gamma).
  theta_hat <- mean(wt_info$omega * Gamma)

  ## 3.) Correction term.
  C_hat <- correction_cwate(Gamma, tau_hat, wt_info, spec)

  ## 4.) Influence function: phi_hat_i = omega_i * Gamma_i - theta_hat + C_hat_i.
  phi_hat <- wt_info$omega * Gamma - theta_hat + C_hat

  ## 5.) Variance and standard error.
  V_hat <- mean(phi_hat^2)
  se <- sqrt(V_hat / n)

  ## Output.
  return(assemble_inference(theta_hat, se, alpha, type = "cwate"))
}
