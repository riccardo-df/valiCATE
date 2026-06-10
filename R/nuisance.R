#' AIPW Pseudo-Outcomes
#'
#' Computes AIPW pseudo-outcomes from nuisance function estimates.
#'
#' @param Y Outcome vector.
#' @param D Treatment indicator vector, binary 0/1 with 1 for treated.
#' @param pscore Propensity score estimates of P(D = 1 | X).
#' @param mu0 Conditional mean estimates of E[Y | D = 0, X].
#' @param mu1 Conditional mean estimates of E[Y | D = 1, X].
#'
#' @return
#' Numeric vector of AIPW pseudo-outcomes.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{estimate_nuisances}}, \code{\link{valiCATE}}
#'
#' @keywords internal
compute_aipw_scores <- function(Y, D, pscore, mu0, mu1) {
  ## Output.
  return(mu1 - mu0 + D / pscore * (Y - mu1) - (1 - D) / (1 - pscore) * (Y - mu0))
}


#' Nuisance Estimation via Cross-Fitting
#'
#' Estimates propensity scores and conditional means via K-fold cross-fitting
#' using \code{\link[grf]{regression_forest}}, and construct AIPW pseudo-outcomes.
#'
#' @param Y Outcome vector.
#' @param D Treatment indicator vector, binary 0/1 with 1 for treated.
#' @param X Covariate matrix or data frame, without intercept.
#' @param n_folds Number of cross-fitting folds.
#'
#' @return
#' A list with the following elements:\cr
#'
#' \code{scores}: Numeric vector of AIPW pseudo-outcomes.\cr
#' \code{pscore}: Numeric vector of propensity score estimates.\cr
#' \code{mu0}: Numeric vector of conditional mean estimates for control units.\cr
#' \code{mu1}: Numeric vector of conditional mean estimates for treated units.
#'
#' @details
#' Nuisance functions are estimated using honest \code{\link[grf]{regression_forest}}s. Propensity scores are clipped to
#' [0.025, 0.975] to avoid extreme inverse probability weights.\cr
#'
#' @import grf
#' @importFrom stats predict
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{compute_aipw_scores}}, \code{\link{valiCATE}}
#'
#' @keywords internal
estimate_nuisances <- function(Y, D, X, n_folds = 5) {
  ## 0.) Store useful quantities.
  n <- length(Y)
  X <- as.matrix(X)

  ## 1.) Create fold assignments.
  fold_ids <- sample(rep(seq_len(n_folds), length.out = n))

  ## 2.) Initialize output vectors.
  pscore_hat <- numeric(n)
  mu0_hat <- numeric(n)
  mu1_hat <- numeric(n)

  ## 3.) Cross-fitting.
  for (k in seq_len(n_folds)) {
    test_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    X_train <- X[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]
    Y_train <- Y[train_idx]
    D_train <- D[train_idx]

    # 3a.) Propensity score: E[D|X].
    pscore_forest <- grf::regression_forest(X_train, D_train)
    pscore_hat[test_idx] <- stats::predict(pscore_forest, X_test)$predictions

    # 3b.) mu0: E[Y|D=0,X] using control units.
    ctrl_idx <- which(D_train == 0)
    mu0_forest <- grf::regression_forest(X_train[ctrl_idx, , drop = FALSE], Y_train[ctrl_idx])
    mu0_hat[test_idx] <- stats::predict(mu0_forest, X_test)$predictions

    # 3c.) mu1: E[Y|D=1,X] using treated units.
    trt_idx <- which(D_train == 1)
    mu1_forest <- grf::regression_forest(X_train[trt_idx, , drop = FALSE], Y_train[trt_idx])
    mu1_hat[test_idx] <- stats::predict(mu1_forest, X_test)$predictions
  }

  ## 4.) Clip propensity scores to avoid extreme weights.
  pscore_hat <- pmax(pmin(pscore_hat, 0.975), 0.025)

  ## 5.) Compute AIPW scores.
  scores <- compute_aipw_scores(Y, D, pscore_hat, mu0_hat, mu1_hat)

  ## Output.
  return(list("scores" = scores, "pscore" = pscore_hat, "mu0" = mu0_hat, "mu1" = mu1_hat))
}
