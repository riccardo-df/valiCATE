#' Weight Function Specifications
#'
#' Returns a list describing the weight function for a given weight name.
#'
#' @param weight_name Character string, one of \code{"AUTOC"}, \code{"AUC-HVL"}, \code{"BLP"}, \code{"QINI"}.
#'
#' @return
#' A list with the following elements:\cr
#'
#' \code{name}: Character name of the weight.\cr
#' \code{h}: Function h(t, v) defining omega(x) = h(tau_hat(x), v(x)).\cr
#' \code{nabla_v_h}: Function computing the partial derivative of h with respect to v.\cr
#' \code{cdf_based}: Logical. If \code{TRUE}, v(x) = F(tau_hat(x)) and m_x(t) = I(t <= tau_hat(x)).
#'   If \code{FALSE} (BLP), v(x) = E[tau_hat] and m_x(t) = t.
#'
#' @details
#' The four weight functions are defined as follows:\cr
#'
#' AUTOC: h(t, v) = -log(1 - v) - 1, with nabla_v_h = 1 / (1 - v).\cr
#' AUC-HVL: h(t, v) = log(v / (1 - v)), with nabla_v_h = 1 / (v * (1 - v)).\cr
#' BLP: h(t, v) = t - v, with nabla_v_h = -1.\cr
#' QINI: h(t, v) = v - 1/2, with nabla_v_h = 1.\cr
#'
#' For CDF-based weights (AUTOC, AUC-HVL, QINI), v is computed using a safe empirical CDF that avoids 0 and 1 values.
#' For BLP, v is the sample mean of the CATE predictions.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{compute_weights}}, \code{\link{valiCATE}}
#'
#' @keywords internal
get_weight_spec <- function(weight_name) {
  ## Output.
  switch(weight_name,
    "AUTOC" = list(
      name = "AUTOC",
      h = function(t, v) -log(1 - v) - 1,
      nabla_v_h = function(t, v) 1 / (1 - v),
      cdf_based = TRUE
    ),
    "AUC-HVL" = list(
      name = "AUC-HVL",
      h = function(t, v) log(v / (1 - v)),
      nabla_v_h = function(t, v) 1 / (v * (1 - v)),
      cdf_based = TRUE
    ),
    "BLP" = list(
      name = "BLP",
      h = function(t, v) t - v,
      nabla_v_h = function(t, v) -1,
      cdf_based = FALSE
    ),
    "QINI" = list(
      name = "QINI",
      h = function(t, v) v - 0.5,
      nabla_v_h = function(t, v) 1,
      cdf_based = TRUE
    ),
    stop(paste0("Invalid 'weight_name'. Unknown weight: '", weight_name, "'."), call. = FALSE)
  )
}


#' Compute Weight Values
#'
#' Given CATE predictions and a weight specification, computes the weight values omega_hat(X_i), the values v_hat(X_i),
#' and the gradient values nabla_v_h(X_i) for each observation.
#'
#' @param tau_hat Numeric vector of CATE predictions.
#' @param spec Weight specification list from \code{\link{get_weight_spec}}.
#'
#' @return
#' A list with the following elements:\cr
#'
#' \code{omega}: Numeric vector of weight values omega_hat(X_i).\cr
#' \code{v}: Numeric vector of centering values v_hat(X_i).\cr
#' \code{nabla}: Numeric vector of gradient values nabla_v_h(X_i).
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{get_weight_spec}}, \code{\link{valiCATE}}
#'
#' @keywords internal
compute_weights <- function(tau_hat, spec) {
  ## 1.) Compute v_hat.
  n <- length(tau_hat)

  if (spec$cdf_based) {
    v <- ecdf_safe(tau_hat)
  } else {
    v <- rep(mean(tau_hat), n)
  }

  ## 2.) Compute omega_hat and nabla_v_h.
  omega <- spec$h(tau_hat, v)
  nabla <- rep(spec$nabla_v_h(tau_hat, v), length.out = n)

  ## Output.
  return(list("omega" = omega, "v" = v, "nabla" = nabla))
}
