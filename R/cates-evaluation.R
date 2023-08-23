#' CATEs Evaluation
#'
#' Evaluates the quality of CATEs estimates by implementing several Generic Machine Learning methodologies.
#'
#' @param y Observed outcomes.
#' @param D Treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param cates Estimated CATEs. CATEs must be estimated using only the training sample.
#' @param is_train Logical vector denoting which observations belong to the training sample.
#' @param pscore Estimated propensity scores. They must be estimated using only the training sample. If not provided by the user, they are estimated internally via an honest regression forest.
#' @param strategy String indicating the desired identification strategy. One of "wr", "ht", "aipw".
#' @param denoise String denoting whether to include optional constructed covariates to the regressions. One of "none", "cddf1", "cddf2", "mck1", "mck2", "mck3". Useful to reduce the variance of the estimation.
#'
#' @return
#' Return.
#'
#' @examples
#' ## Generate data.
#' set.seed(1986)
#'
#' n <- 3000
#' k <- 3
#'
#' @md
#' @details
#' Currently, the routine supports only single-splits.\cr
#' Regression functions are estimated via honest \code{\link[grf]{regression_forest}}s.\cr
#' If strategy is set to \code{"aipw"}, doubly-robust scores are estimated internally using the validation sample via 5-fold cross fitting and using honest regression forests
#' (see the \code{\link[aggTrees]{dr_scores}} function).\cr
#' If not provided, \code{cates} are estimated internally using a causal forest estimator and the training sample.\cr
#' If not provided, \code{pscore} are estimated internally using a regression forest.\cr
#'
#' @import grf
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
evalue_cates <- function(y, D, X, cates, is_train, pscore = NULL, strategy = "wr", denoise = "none") {
  ## 0.) Handling inputs and checks.
  if (is.logical(D)) D <- as.numeric(D)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("Invalid 'X'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.logical(is_train)) stop("Invalid 'is_train'. This must be a logical vector.", call. = FALSE)
  if (!(strategy %in% c("wr", "ht", "aipw"))) stop("Invalid 'strategy'. This must be one of 'wr', 'ht', 'aipw'.", call. = FALSE)
  if (!(denoise %in% c("none", "cddf1", "cddf2", "mck1", "mck2", "mck3"))) stop("Invalid 'strategy'. This must be one of 'none', 'cddf1', 'cddf2', 'mck1', 'mck2', 'mck3'.", call. = FALSE)
  if (strategy == "aipw" & denoise != "none") stop("If 'aipw' is the selected strategy, then 'denoise' must be 'none'.", call. = FALSE)
  if (strategy != "ht" & denoise %in% c("mck2", "mck3")) stop("'denoise' can be set to 'mck2' or 'mck3' only if 'strategy' is 'ht'.", call. = FALSE)

  # mu_val <- mu0_val <- mu1_val <- NULL

  ## 1.) Split the sample as indicated by the user. If necessary, estimate cates and propensity score using the training sample and doubly-robust scores using the validation sample.
  train_idx <- which(is_train)
  val_idx <- which(!is_train)

  y_tr <- y[train_idx]
  D_tr <- D[train_idx]
  X_tr <- X[train_idx, ]
  cates_tr <- cates[train_idx]

  y_val <- y[val_idx]
  D_val <- D[val_idx]
  X_val <- X[val_idx, ]
  cates_val <- cates[val_idx]

  if (is.null(pscore)) {
    pscore_forest <- grf::regression_forest(X_tr, D_tr)
    pscore_val <- stats::predict(pscore_forest, X_val)$predictions
  }

  if (strategy == "aipw") scores <- aggTrees::dr_scores(y_val, D_val, X_val) else scores <- NULL

  ## 2.) Estimate regression functions for the optional covariates in training sample.
  if (denoise != "none") {
    # Estimate functions that we always need when including denoise terms.
    mu0_forest <- grf::regression_forest(X_tr[D_tr == 0, ], y_tr[D_tr == 0])
    mu0_val <- stats::predict(mu0_forest, X_val)$predictions

    # We need to estimate \mu only under one (strategy, denoise) combination. Similarly, we need to estimate \mu_1 only under two (strategy, denoise) combinations.
    if (strategy == "wr" & denoise == "mck1") {
      mu_forest <- grf::regression_forest(X_tr, y_tr)
      mu_val <- stats::predict(mu_forest, X_val)$predictions
    } else if (strategy == "ht" & denoise %in% c("mck2", "mck3")) {
      mu1_forest <- grf::regression_forest(X_tr[D_tr == 1, ], y_tr[D_tr == 1])
      mu1_val <- stats::predict(mu1_forest, X_val)$predictions
    }
  }

  ## 3.) Estimate BLP and GATES in validation sample.
  blp_results <- blp_regression(y_val, D_val, cates_val, strategy, denoise, pscore_val, mu_val, mu0_val, mu1_val, scores)

  ## Output.

}
