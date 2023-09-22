#' CATEs Evaluation
#'
#' Evaluates the quality of CATEs estimates by estimating the best linear predictor (BLP) of the actual CATEs using the estimated
#' CATEs, the sorted group average treatment effects (GATES), and the rank-weighted average treatment effect (RATE) induced by
#' the estimated CATEs.
#'
#' @param Y Observed outcomes.
#' @param D Treatment indicator.
#' @param X Covariate matrix (no intercept).
#' @param cates Estimated CATEs. Must be estimated using only the training sample.
#' @param is_train Logical vector denoting which observations belong to the training sample.
#' @param pscore Propensity scores. If unknown, must be estimated using only the training sample. 
#' @param mu Estimated regression function. Must be estimated using only the training sample.
#' @param mu0 Estimated regression function for control units. Must be estimated using only the training sample.
#' @param mu1 Estimated regression function for treated units. Must be estimated using only the training sample.
#' @param n_groups Number of groups to be formed for the GATES analysis.
#' @param beneficial Logical, whether the treatment is beneficial to units. If \code{TRUE}, units are ranked according to decreasing values of \code{cates} to estimate the RATEs, otherwise they are ranked according to increasing values of \code{cates}.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates <- predict(forest, X)$predictions # We predict on the whole sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore <- rep(0.5, length(Y))
#' evaluation <- evalu_cates(Y, D, X, cates, train_idx, pscore = pscore)
#' 
#' ## Generic S3 methods.
#' summary(evaluation, target = "BLP")
#' summary(evaluation, target = "BLP", latex = TRUE)
#' 
#' summary(evaluation, target = "GATES")
#' summary(evaluation, target = "GATES", latex = TRUE)
#'
#' plot(evaluation, target = "GATES")
#' plot(evaluation, target = "TOC")}
#'
#' @md
#' @details
#' To estimate BLP, GATES, and RATEs, the user must provide observations on the outcomes, the treatment status, and the covariates of units in the whole sample, as well 
#' as their estimated CATEs. Be careful, as the CATEs must be estimated only with part of the sample, which we call the training sample (see the example section below).\cr
#' 
#' To let the function know which observations were used for the CATEs estimation, the user must also provide a logical vector with the \code{TRUE}s denoting observations in the 
#' training sample. This way, \code{\link{evalu_cates}} knows which observations to use to post-process the CATEs estimates.\cr
#' 
#' \code{\link{evalu_cates}} implements a number of strategies to estimate the BLP and the GATES. Most of them involve fitting a suitable linear model. The linear models differ according to the
#' different identification strategies. Furthermore, for each strategy, there exist various sets of constructed covariates that one can add to reduce the variance of the estimation. \code{\link{evalue_cates}}
#' fits and returns all these possible models. GATES are also estimated using a nonparametric approach. Check the online 
#' \href{https://riccardo-df.github.io/evaluCATE/articles/evalu-cates-short-tutorial.html}{short tutorial} for details.\cr 
#' 
#' For the linear models, standard errors are estimated using the Eicker-Huber-White estimator. These standard errors are then used to test three distinct hypotheses of effect heterogeneity: whether
#' all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, and whether the differences in the GATES across all pairs of groups are zero.
#' For the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure and report the median of the adjusted p-values. The nonparametric approach tests only the first
#' of these hypotheses. Check the \href{https://riccardo-df.github.io/evaluCATE/articles/hypotheses-testing.html}{hypotheses testing vignette} for details.\cr
#' 
#' Some of the linear models involve covariates that depend on particular nuisance functions, e.g., propensity score and conditional mean of the outcome 
#' (check the online \href{https://riccardo-df.github.io/evaluCATE/articles/denoising.html}{denoising vignette} for details about these covariates). The user can supply estimates of these functions by using the 
#' optional arguments \code{pscore}, \code{mu}, \code{mu0}, and \code{mu1}. Be careful, as these must be obtained using only the training sample. If not provided by the user, these functions are estimated internally 
#' via honest \code{\link[grf]{regression_forest}}s using only the training sample. \cr
#' 
#' To estimate the BLP and GATES using the AIPW strategy, doubly-robust scores are estimated internally using the validation sample via 5-fold cross fitting and honest regression forests (see the 
#' \code{\link[aggTrees]{dr_scores}} function for details). The same doubly-robust scores are also used to estimate the RATEs.\cr
#' 
#' Groups are constructed by cutting the distribution of \code{cates} into \code{n_groups} quantiles. If this leads to one or more groups composed of only treated or only control units, the function raises an error.\cr
#' 
#' Two different RATEs are estimated: AUTOC and QINI coefficient. Sample-averaging estimators are employed. Standard errors are estimated by the standard deviation of the bootstrap estimates obtained using the half-sample bootstrap. 
#'
#' @import grf
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
evalu_cates <- function(Y, D, X, cates, is_train,
                        pscore = NULL, mu = NULL, mu0 = NULL, mu1 = NULL,
                        n_groups = 5, beneficial = TRUE, n_boot = 200, verbose = TRUE) {
  ## 0.) Handling inputs and checks.
  if (is.logical(D)) D <- as.numeric(D)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("Invalid 'X'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.logical(is_train)) stop("Invalid 'is_train'. This must be a logical vector.", call. = FALSE)
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)
  if (n_boot <= 1 | n_boot %% 1 != 0) stop("Invalid 'n_boot'. This must be an integer greater than 1.", call. = FALSE)
  if (!is.logical((beneficial))) stop("Invalid 'beneficial'. This must be either TRUE or FALSE.", call. = FALSE)
  if (!is.logical((verbose))) stop("Invalid 'verbose'. This must be either TRUE or FALSE.", call. = FALSE)
  
  ## 1.) Split the sample as indicated by the user.  and doubly-robust scores using the validation sample.
  train_idx <- which(is_train)
  val_idx <- which(!is_train)
  
  Y_tr <- Y[train_idx]
  D_tr <- D[train_idx]
  X_tr <- matrix(X[train_idx, ], ncol = dim(X)[2])
  cates_tr <- cates[train_idx]
  
  Y_val <- Y[val_idx]
  D_val <- D[val_idx]
  X_val <- matrix(X[val_idx, ], ncol = dim(X)[2])
  cates_val <- cates[val_idx]
  
  ## 2.) If necessary, estimate nuisance functions using the training sample. Then,estimated AIPW scores in validation sample via cross-fitting.
  if (verbose) cat("Estimating nuisance functions and AIPW scores; \n")
  
  if (is.null(pscore)) {
    pscore_forest <- grf::regression_forest(X_tr, D_tr)
    pscore <- stats::predict(pscore_forest, X)$predictions
  }
  
  if (is.null(mu)) {
    mu_forest <- grf::regression_forest(X_tr, Y_tr)
    mu <- stats::predict(mu_forest, X)$predictions
  }
  
  if (is.null(mu0)) {
    mu0_forest <- grf::regression_forest(matrix(X_tr[D_tr == 0, ], ncol = dim(X_tr)[2]), Y_tr[D_tr == 0])
    mu0 <- stats::predict(mu0_forest, X)$predictions
  }
  
  if (is.null(mu1)) {
    mu1_forest <- grf::regression_forest(matrix(X_tr[D_tr == 1, ], ncol = dim(X)[2]), Y_tr[D_tr == 1])
    mu1 <- stats::predict(mu1_forest, X)$predictions
  }
  
  pscore_tr <- pscore[train_idx]
  pscore_val <- pscore[val_idx]
  
  mu_tr <- mu[train_idx]
  mu_val <- mu[val_idx]
  
  mu0_tr <- mu0[train_idx]
  mu0_val <- mu0[val_idx]
  
  mu1_tr <- mu1[train_idx]
  mu1_val <- mu1[val_idx]
  
  scores_val <- aggTrees::dr_scores(Y_val, D_val, X_val)
  
  ## 3.) Estimate BLP, GATES, and RATE in validation sample.
  if (verbose) cat("BLP estimation; \n")
  blp_results <- blp_estimation(Y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores_val)
  
  if (verbose) cat("GATES estimation; \n")
  gates_results <- gates_estimation(Y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores_val, n_groups)
  
  if (verbose) cat("RATE estimation; \n")
  rate_results <- rate_estimation(cates_val, scores_val, beneficial, n_boot)
  
  ## Output.
  if (verbose) cat("Output. \n\n")
  out <- list("BLP" = blp_results, "GATES" = gates_results, "RATE" = rate_results)
  class(out) <- "evaluCATE"
  return(out)
}
