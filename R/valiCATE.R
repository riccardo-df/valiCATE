#' CATEs Validation
#'
#' Validates machine learning predictions of conditional average treatment effects (CATEs) by estimating
#' Centered-Weighted Average Treatment Effects (CWATEs) and their normalized counterparts (NCWATEs).
#'
#' @param Y Outcome vector (validation sample).
#' @param D Treatment indicator vector, binary 0/1 with 1 for treated (validation sample).
#' @param X Covariate matrix or data frame, without intercept (validation sample).
#' @param cates Named list of CATE prediction vectors on the validation sample produced by different models.
#'   Each element must be a numeric vector of the same length as \code{Y}. CATE models must be estimated using only the training
#'   sample for valid inference.
#' @param weights Character vector controlling which weight functions to use. Admitted values are \code{"AUTOC"}, \code{"AUC-HVL"},
#'   \code{"BLP"}, and \code{"QINI"}.
#' @param scores Optional, pre-computed AIPW pseudo-outcomes. If not provided by the user, scores are estimated internally via
#'   cross-fitting. Useful to save computational time if scores have already been estimated.
#' @param n_folds Optional, number of cross-fitting folds for nuisance estimation. Default is 5. Ignored if \code{scores} is provided.
#' @param alpha Optional, significance level for confidence intervals and hypothesis tests. Default is 0.05.
#' @param verbose Optional, set to \code{FALSE} to prevent the function from printing the progresses.
#'
#' @return
#' Object of class \code{valiCATE}.
#' 
#' @examples
#' \donttest{## Generate data.
#' set.seed(1986)
#'
#' n <- 1000
#' k <- 2
#'
#' X <- matrix(rnorm(n * k), ncol = k)
#' colnames(X) <- paste0("x", seq_len(k))
#' D <- rbinom(n, size = 1, prob = 0.5)
#' mu0 <- 0.5 * X[, 1]
#' mu1 <- 0.5 * X[, 1] + X[, 2]
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n, sd = 0.5)
#'
#' ## Split into training and validation samples.
#' train_idx <- sample(1:n, n / 2)
#' val_idx <- setdiff(1:n, train_idx)
#'
#' ## Estimate CATEs on the training sample, predict on the validation sample.
#' library(grf)
#' cf <- causal_forest(X[train_idx, ], Y[train_idx], D[train_idx])
#' cates <- predict(cf, X[val_idx, ])$predictions
#'
#' ## Validate using the validation sample.
#' result <- valiCATE(Y[val_idx], D[val_idx], X[val_idx, ],
#'                    cates = list("causal_forest" = cates))
#'
#' ## We can also compare multiple models.
#' cf2 <- causal_forest(X[train_idx, ], Y[train_idx], D[train_idx],
#'                      num.trees = 500)
#' cates2 <- predict(cf2, X[val_idx, ])$predictions
#'
#' result_multi <- valiCATE(Y[val_idx], D[val_idx], X[val_idx, ],
#'                          cates = list("cf_default" = cates,
#'                                       "cf_500" = cates2))
#'
#' ## We have compatibility with generic S3-methods.
#' summary(result)
#' summary(result, latex = TRUE)}
#'
#' @details
#' The user must provide observations on the outcomes, the treatment status, and the covariates of units in the validation sample
#' using the first three arguments. The user must also provide CATE predictions on the validation sample as a named list
#' storing predictions produced by different models. Be careful, CATE models must be estimated using only the training sample to 
#' achieve valid inference.\cr
#'
#' Estimation is based on AIPW pseudo-outcomes, which require nuisance function estimates (propensity score, conditional mean of the
#' outcome for treated and control units). By default, nuisance functions are estimated internally via honest
#' \code{\link[grf]{regression_forest}}s using K-fold cross-fitting in the validation sample. Alternatively, the user can supply
#' pre-computed AIPW scores via the \code{scores} argument; in this case, scores should be cross-fitted in the validation sample
#' to achieve valid inference.\cr
#'
#' For the CWATE, a one-sided test of H0: theta <= 0 is reported. Rejection signals that the predicted CATEs capture genuine
#' heterogeneity in the correct direction. For the NCWATE, a two-sided test of H0: gamma = 1 is reported. Non-rejection signals
#' that the predicted CATEs recover the true CATEs.
#'
#' @import grf
#' @importFrom stats predict var qnorm pnorm
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{summary.valiCATE}}, \code{\link{print.valiCATE}}, \code{\link{plot.valiCATE}}
#'
#' @export
valiCATE <- function(Y, D, X, cates,
                     weights = c("AUTOC", "AUC-HVL", "BLP", "QINI"),
                     scores = NULL,
                     n_folds = 5,
                     alpha = 0.05,
                     verbose = TRUE) {
  ## 0.) Handling inputs and checks.
  cl <- match.call()
  validate_inputs(Y, D, X, cates, weights, scores, n_folds, alpha, verbose)
  n <- length(Y)

  ## 1.) Estimate nuisance functions and AIPW scores.
  if (is.null(scores)) {
    if (verbose) cat("Estimating nuisance functions via", n_folds, "-fold cross-fitting; \n")
    nuis <- estimate_nuisances(Y, D, X, n_folds)
    scores <- nuis$scores
  } else {
    if (verbose) cat("Using pre-computed AIPW scores; \n")
  }

  ## 2.) Estimate CWATEs and NCWATEs.
  model_names <- names(cates)
  results <- list()

  for (model_nm in model_names) {
    if (verbose) cat("CWATE/NCWATE estimation for model:", model_nm, "; \n")
    tau_hat <- cates[[model_nm]]
    results[[model_nm]] <- list()

    for (wt in weights) {
      cwate_res <- estimate_cwate(scores, tau_hat, wt, alpha)
      ncwate_res <- estimate_ncwate(scores, tau_hat, wt, alpha)
      results[[model_nm]][[wt]] <- list("cwate" = cwate_res, "ncwate" = ncwate_res)
    }
  }

  ## 3.) Construct valiCATE object.
  if (verbose) cat("Output. \n\n")

  out <- list("results" = results, "scores" = scores, "cates" = cates,
              "weights" = weights, "alpha" = alpha, "n" = n, "call" = cl)
  class(out) <- "valiCATE"

  ## Output.
  return(out)
}
