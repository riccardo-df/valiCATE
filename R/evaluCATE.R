#' CATEs Evaluation
#'
#' Evaluates the quality of CATE predictions by estimating the best linear predictor (BLP) of the actual CATEs using the estimated
#' CATEs, the sorted group average treatment effects (GATES), and the rank-weighted average treatment effects (RATEs) induced by
#' the estimated CATEs.
#'
#' @param Y_tr Observed outcomes for the training sample.
#' @param Y_val Observed outcomes for the validation sample.
#' @param D_tr Treatment indicator for the training sample.
#' @param D_val Treatment indicator for the validation sample.
#' @param X_tr Covariate matrix for the training sample (no intercept).
#' @param X_val Covariate matrix for the validation sample (no intercept).
#' @param cates_val Named list storing CATE predictions on the validation sample produced by different models. Models must be estimated using only the training sample.
#' @param strategies Character vector controlling the identification and estimation strategies to implement for BLP and GATES. Admitted values are \code{"WR"}, \code{"HT"}, and \code{"AIPW"}.
#' @param denoising Character vector controlling if and which additional covariates to include in the regressions for BLP and GATES to reduce the variance of the estimation. Admitted values are \code{"none"}, \code{"cddf1"}, \code{"cddf2"}, \code{"mck1"}, \code{"mck2"}, and \code{"mck3"}.
#' @param pscore_val Propensity score predictions on the validation sample. Must be produced by a model estimated using only the training sample (unless the propensity score is known, in which case we provide the true values).
#' @param mu_val Conditional mean predictions on the validation sample. Must be produced by a model estimated using only the training sample.
#' @param mu0_val Control units' conditional mean predictions on the validation sample. Must be produced by a model estimated using only the training sample.
#' @param mu1_val Treated units' conditional mean predictions on the validation sample. Must be produced by a model estimated using only the training sample.
#' @param n_groups Number of groups to be formed for the GATES analysis.
#' @param beneficial Logical, whether the treatment is beneficial to units. If \code{TRUE}, units are ranked according to decreasing values of \code{cates_val} to estimate the RATEs, otherwise they are ranked according to increasing values of \code{cates_val}.
#' @param n_boot Number of bootstrap replications to estimate the standard error of the RATE estimates.
#' @param verbose Logical, set to FALSE to prevent the function from printing the progresses.
#'
#' @return
#' An \code{evaluCATE} object.
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
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#' 
#' ## Sample split.
#' train_idx <- sample(c(TRUE, FALSE), length(Y), replace = TRUE)
#' 
#' X_tr <- X[train_idx, ]
#' X_val <- X[!train_idx, ]
#' 
#' D_tr <- D[train_idx]
#' D_val <- D[!train_idx]
#' 
#' Y_tr <- Y[train_idx]
#' Y_val <- Y[!train_idx]
#' 
#' ## CATEs estimation. Models are estimated with training sample.
#' # T-learner.
#' library(grf)
#' 
#' forest_treated <- regression_forest(X_tr[D_tr == 1, ], Y_tr[D_tr == 1])
#' forest_control <- regression_forest(X_tr[D_tr == 0, ], Y_tr[D_tr == 0])
#' cates_val_t <- predict(forest_treated, X_val)$predictions - 
#'                predict(forest_control, X_val)$predictions 
#' 
#' # Grf.
#' forest_grf <- causal_forest(X_tr, Y_tr, D_tr) 
#' cates_val_grf <- predict(forest_grf, X_val)$predictions 
#' 
#' ## CATEs evaluation. 
#' # Use all strategies with no denoising.
#' strategies <- c("WR", "HT")
#' denoising <- "none"
#' 
#' # We know true pscore.
#' pscore_val <- rep(0.5, length(Y_val))
#' 
#' # Construct CATEs list.
#' cates_val <- list("T-learner" = cates_val_t,
#'                   "grf" = cates_val_grf)
#' 
#' # Call main function.
#' evaluation <- evaluCATE(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, 
#'                         strategies = strategies, denoising = denoising, 
#'                         pscore_val = pscore_val)
#' 
#' ## Generic S3 methods.
#' summary(evaluation, target = "BLP")
#' summary(evaluation, target = "BLP", latex = TRUE)
#'
#' summary(evaluation, target = "GATES")
#' summary(evaluation, target = "GATES", latex = TRUE)
#'
#' summary(evaluation, target = "RATE")
#' summary(evaluation, target = "RATE", latex = TRUE)
#'
#' plot(evaluation, target = "GATES")
#' plot(evaluation, target = "TOC")
#' plot(evaluation, target = "RATE")}
#'
#' @md
#' @details
#' The user must provide observations on the outcomes, the treatment status, and the covariates of units in the training and validation samples separately
#' using the first six arguments. The user must also provide a named list storing CATE predictions on the validation sample produced from different models (the Example section below shows how to construct such a list).
#' Be careful, models must be estimated using only the training sample to achieve valid inference.\cr
#' 
#' The \code{\link{evaluCATE}} function allows the implementation of three different strategies for BLP and GATES identification and estimation: a) Weighted Residuals (WR), Horwitz-Thompson (HT), and
#' Augmented Inverse-Probability Weighting (AIPW). The user can choose their preferred strategies by controlling the \code{strategies} argument. This has no impact on RATEs estimation. GATES are also 
#' always estimated using a nonparametric approach.\cr
#' 
#' Most of the BLP and GATES estimation strategies involve fitting a suitable linear model. For each model, there exist various sets of constructed covariates that one can add to reduce the variance of 
#' the estimation. The user can choose whether to add these additional covariates by controlling the \code{denoising} argument (check the online 
#' \href{https://riccardo-df.github.io/evaluCATE/articles/denoising.html}{denoising vignette} for details). This has no impact on RATEs estimation and on the results from the nonparametric GATES estimation strategy. 
#' 
#' The constructed covariates depend on particular nuisance functions, e.g., propensity score and conditional mean of the outcome. The user can supply predictions on the validation sample of these functions 
#' by using the  optional arguments \code{pscore_val}, \code{mu_val}, \code{mu0_val}, and \code{mu1_val}. Be careful, as these predictions must be produced by models estimated using only the training sample. If not 
#' provided by the user, these functions are estimated internally  via honest \code{\link[grf]{regression_forest}}s using only the training sample. \cr
#' 
#' For the linear models, standard errors are estimated using the Eicker-Huber-White estimator. Under our careful sample splitting procedure, these standard errors can then used to test various 
#' hypotheses of effect heterogeneity. For the GATES, we focus on three distinct hypotheses: whether all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, 
#' and whether the differences in the GATES across all pairs of groups are zero (for the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure and report the median of 
#' the adjusted p-values). The nonparametric  approach tests only the first of these hypotheses. Check the \href{https://riccardo-df.github.io/evaluCATE/articles/hypotheses-testing.html}{hypotheses testing vignette} 
#' for details.\cr
#' 
#' To estimate the BLP and GATES using the AIPW strategy, doubly-robust scores are estimated internally using the validation sample via 5-fold cross fitting and honest regression forests (see the 
#' \code{\link[aggTrees]{dr_scores}} function for details). The same doubly-robust scores are also used to estimate the RATEs.\cr
#' 
#' Groups for the GATES analysis are constructed by cutting the distribution of \code{cates_val} into \code{n_groups} quantiles. If this leads to one or more groups composed of only treated or only control units, 
#' the function raises an error. Possible solutions include: a) change your original training-validation sample split; b) increase the fraction of observations allocated to the validation sample.\cr
#' 
#' The \code{\link{evaluCATE}} function estimates two different RATEs: AUTOC and QINI coefficients. Sample-averaging estimators are employed. Standard errors are estimated by the standard deviation of the bootstrap 
#' estimates obtained using the half-sample bootstrap.\cr
#' 
#' Check the online \href{https://riccardo-df.github.io/evaluCATE/articles/evaluCATE-short-tutorial.html}{short tutorial} for a guided usage of this function.\cr 
#'
#' @import grf
#' @importFrom stats predict
#' @importFrom stats var
#'
#' @author Riccardo Di Francesco
#'
#' @export
evaluCATE <- function(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, 
                      strategies = c("WR", "HT", "AIPW"), denoising = c("none", "cddf1", "cddf2", "mck1", "mck2", "mck3"),
                      pscore_val = NULL, mu_val = NULL, mu0_val = NULL, mu1_val = NULL,
                      n_groups = 5, beneficial = TRUE, n_boot = 200, verbose = TRUE) {
  ## 0.) Handling inputs and checks.
  if (is.logical(D_tr)) D_tr <- as.numeric(D_tr)
  if (is.logical(D_val)) D_val <- as.numeric(D_val)
  if (any(!(D_tr %in% c(0, 1)))) stop("Invalid 'D_tr'. Only binary treatments are allowed.", call. = FALSE)
  if (any(!(D_val %in% c(0, 1)))) stop("Invalid 'D_val'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X_tr) & !is.data.frame(X_tr)) stop("Invalid 'X_tr'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.matrix(X_val) & !is.data.frame(X_val)) stop("Invalid 'X_val'. This must be either a matrix or a data frame.", call. = FALSE)
  if (any(sapply(cates_val, stats::var) == 0)) stop("No variation in at least one element of 'cates_val'.", call. = FALSE)
  if (any(!(strategies %in% c("WR", "HT", "AIPW")))) stop("Invalid 'strategies'. Must be one of 'WR', 'HT', 'AIPW'.", call. = FALSE)
  if (any(!(denoising %in% c("none", "cddf1", "cddf2", "mck1", "mck2", "mck3")))) stop("Invalid 'denoising'. Must be one of 'none', 'cddf1', 'cddf2', 'mck1', 'mck2', 'mck3'.", call. = FALSE)
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)
  if (n_boot <= 1 | n_boot %% 1 != 0) stop("Invalid 'n_boot'. This must be an integer greater than 1.", call. = FALSE)
  if (!is.logical((beneficial))) stop("Invalid 'beneficial'. This must be either TRUE or FALSE.", call. = FALSE)
  if (!is.logical((verbose))) stop("Invalid 'verbose'. This must be either TRUE or FALSE.", call. = FALSE)
  
  ## 1.) Estimate necessary nuisance functions using the training sample. Propensity score always required.
  if (verbose) cat("Estimating nuisance functions; \n")
  
  scores_val <- aggTrees::dr_scores(Y_val, D_val, X_val)
  
  if (is.null(pscore_val)) {
    pscore_forest <- grf::regression_forest(X_tr, D_tr)
    pscore_val <- stats::predict(pscore_forest, X_val)$predictions
  }
  
  if (any(denoising %in% c("cddf1", "cddf2", "mck2", "mck3"))) {
    if (is.null(mu0_val)) {
      mu0_forest <- grf::regression_forest(as.matrix(X_tr[D_tr == 0, ], ncol = dim(X_tr)[2]), Y_tr[D_tr == 0])
      mu0_val <- stats::predict(mu0_forest, X_val)$predictions
    }
  } 
  
  if (any(denoising %in% c("mck2", "mck3"))) {
    if (is.null(mu1_val)) {
      mu1_forest <- grf::regression_forest(as.matrix(X_tr[D_tr == 1, ], ncol = dim(X_tr)[2]), Y_tr[D_tr == 1])
      mu1_val <- stats::predict(mu1_forest, X_val)$predictions
    }
  } 
  
  if (any(denoising == "mck1")) {
    if (is.null(mu_val)) {
      mu_forest <- grf::regression_forest(X_tr, Y_tr)
      mu_val <- stats::predict(mu_forest, X_val)$predictions
    }
  }
  
  ## 3.) Estimate BLP, GATES, and RATE in validation sample.
  if (verbose) cat("BLP estimation; \n")
  blp_results <- lapply(cates_val, function(x) { blp_estimation(Y_val, D_val, x, pscore_val, mu_val, mu0_val, mu1_val, scores_val, strategies, denoising) })
    
  if (verbose) cat("GATES estimation; \n")
  gates_results <- lapply(cates_val, function(x) { gates_estimation(Y_val, D_val, x, pscore_val, mu_val, mu0_val, mu1_val, scores_val, n_groups, strategies, denoising) })
    
  if (verbose) cat("RATE estimation; \n")
  rate_results <- lapply(cates_val, function(x) { rate_estimation(x, scores_val, beneficial, n_boot) }) 
    
  ## 4.) Output.
  if (verbose) cat("Output. \n\n")
  out <- list("BLP" = blp_results, "GATES" = gates_results, "RATE" = rate_results)
  class(out) <- "evaluCATE"
  return(out)
}
