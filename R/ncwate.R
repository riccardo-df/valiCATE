#' NCWATE Estimation
#'
#' Estimates the Normalized Centered-Weighted Average Treatment Effect (NCWATE) for a given weight function.
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
#' The NCWATE is estimated as gamma_hat = theta_hat / D_hat, where theta_hat is the CWATE estimate and
#' D_hat = mean(omega_hat * tau_hat) is the denominator.\cr
#'
#' A two-sided test of H0: gamma = 1 is reported. Non-rejection signals that the predicted CATEs recover the true CATEs.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{estimate_cwate}}, \code{\link{valiCATE}}
#'
#' @keywords internal
estimate_ncwate <- function(Gamma, tau_hat, weight_name, alpha = 0.05) {
  ## 1.) Compute weights.
  n <- length(Gamma)
  spec <- get_weight_spec(weight_name)
  wt_info <- compute_weights(tau_hat, spec)

  ## 2.) CWATE point estimate.
  theta_hat <- mean(wt_info$omega * Gamma)

  ## 3.) Denominator: D_hat = mean(omega * tau_hat).
  D_hat <- mean(wt_info$omega * tau_hat)

  ## 4.) NCWATE point estimate: gamma_hat = theta_hat / D_hat.
  gamma_hat <- theta_hat / D_hat

  ## 5.) Correction term for NCWATE.
  C_gamma_hat <- correction_ncwate(Gamma, tau_hat, gamma_hat, wt_info, spec)

  ## 6.) Influence function: psi_hat_i = (1/D_hat) * (omega_i * (Gamma_i - gamma_hat * tau_hat_i) + C_gamma_hat_i).
  psi_hat <- (1 / D_hat) * (wt_info$omega * (Gamma - gamma_hat * tau_hat) + C_gamma_hat)

  ## 7.) Variance and standard error.
  V_hat <- mean(psi_hat^2)
  se <- sqrt(V_hat / n)

  ## Output.
  return(assemble_inference(gamma_hat, se, alpha, type = "ncwate"))
}
