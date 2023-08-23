#' CATEs Evaluation
#'
#' Evaluates the quality of CATEs estimates by implementing several Generic Machine Learning methodologies.
#'
#' @param y Observed outcomes.
#' @param D Treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param cates Estimated CATEs. CATEs must be estimated using only the training sample.
#' @param is_train Logical vector denoting which observations belong to the training sample.
#' @param pscore Propensity scores. If unknown, they must be estimated using only the training sample. If not provided by the user, they are estimated internally via an honest \code{\link[grf]{regression_forest}}.
#' @param verbose Logical, set to FALSE to prevent the function from printing the progresses.
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
evalue_cates <- function(y, D, X, cates, is_train, pscore = NULL, verbose = TRUE) {
  ## 0.) Handling inputs and checks.
  if (is.logical(D)) D <- as.numeric(D)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("Invalid 'X'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.logical(is_train)) stop("Invalid 'is_train'. This must be a logical vector.", call. = FALSE)

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

  ## 2.) Estimate nuisance functions using the training sample and doubly-robust scores using the validation sample and cross-fitting.
  if (verbose) cat("Estimating nuisance functions and AIPW scores; \n")
  mu_forest <- grf::regression_forest(X_tr, y_tr)
  mu0_forest <- grf::regression_forest(X_tr[D_tr == 0, ], y_tr[D_tr == 0])
  mu1_forest <- grf::regression_forest(X_tr[D_tr == 1, ], y_tr[D_tr == 1])
  
  mu_val <- stats::predict(mu_forest, X_val)$predictions
  mu0_val <- stats::predict(mu0_forest, X_val)$predictions
  mu1_val <- stats::predict(mu1_forest, X_val)$predictions
  
  scores <- aggTrees::dr_scores(y_val, D_val, X_val)
    
  ## 3.) Estimate BLP and GATES in validation sample.
  if (verbose) cat("BLP estimation; \n")
  blp_results <- blp_estimation(y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores)

  ## Output.
  out <- list("BLP" = blp_results)
  class_out <- "evalue_cates"
  return(out)
}
