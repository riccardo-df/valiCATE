#' CATEs Evaluation
#'
#' Evaluates the quality of CATEs estimates by estimating the best linear predictor (BLP) of the actual CATEs using the estimated
#' CATEs, the sorted group average treatment effects (GATES), and the rank-weighted average treatment effect (RATE) induced by
#' the estimated CATEs.
#'
#' @param Y_tr Observed outcomes for the training sample.
#' @param Y_val Observed outcomes for the validation sample.
#' @param D_tr Treatment indicator for the training sample.
#' @param D_val Treatment indicator for the validation sample.
#' @param X_tr Covariate matrix for the training sample (no intercept).
#' @param X_val Covariate matrix for the validation sample (no intercept).
#' @param cates_val CATE predictions on the validation sample. Must be produced by a model estimated using only the training sample.
#' @param strategies Character vector controlling the identification strategies to implement for BLP and GATES. Admitted values are \code{"WR"} (weighted residuals), \code{"HT"} (Horwitz-Thompson), and \code{"AIPW"} (augmented inverse-probability weighting).
#' @param denoising Character vector controlling if and which additional covariates to include in the regressions to reduce the variance of the estimation. Admitted values are \code{"none"}, \code{"cddf1"}, \code{"cddf2"}, \code{"mck1"}, \code{"mck2"}, and \code{"mck3"}.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates_val <- predict(forest, X_val)$predictions # We predict on the validation sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. Do not use denoising.
#' strategies <- c("WR", "HT", "AIPW")
#' denoising <- "none"
#' 
#' pscore_val <- rep(0.5, length(Y_val))
#' 
#' evaluation <- evaluCATE(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, 
#'                         strategies = strategies, denoising = denoising, pscore_val = pscore_val)
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
#' To estimate BLP, GATES, and RATEs, the user must provide observations on the outcomes, the treatment status, and the covariates of units in the training and validation samples separately.
#' Additionally, the user must provide CATE predictions on the validation sample obtained from a model estimated using only the training sample.\cr
#' 
#' \code{\link{evaluCATE}} implements a number of strategies to estimate the BLP and the GATES. Most of them involve fitting a suitable linear model. The linear models differ according to the
#' different identification strategies. Furthermore, for each strategy, there exist various sets of constructed covariates that one can add to reduce the variance of the estimation. \code{\link{evaluCATE}}
#' fits and returns all these possible models. GATES are also estimated using a nonparametric approach. Check the online 
#' \href{https://riccardo-df.github.io/evaluCATE/articles/evaluCATE-short-tutorial.html}{short tutorial} for details.\cr 
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
#' @importFrom stats predict
#' @importFrom stats var
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
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
  if (stats::var(cates_val) == 0) stop("No variation in 'cates_val'.", call. = FALSE)
  if (any(!(strategies %in% c("WR", "HT", "AIPW")))) stop("Invalid 'strategies'. Must be one of 'none', 'cddf1', 'cddf2', 'mck1', 'mck2', or 'mck3'.", call. = FALSE)
  if (any(!(denoising %in% c("none", "cddf1", "cddf2", "mck1", "mck2", "mck3")))) stop("Invalid 'strategies'. Must be one of 'WR', 'HT', or 'AIPW'.", call. = FALSE)
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)
  if (n_boot <= 1 | n_boot %% 1 != 0) stop("Invalid 'n_boot'. This must be an integer greater than 1.", call. = FALSE)
  if (!is.logical((beneficial))) stop("Invalid 'beneficial'. This must be either TRUE or FALSE.", call. = FALSE)
  if (!is.logical((verbose))) stop("Invalid 'verbose'. This must be either TRUE or FALSE.", call. = FALSE)
  
  ## 1.) Estimate necessary nuisance functions using the training sample. Propensity score always required.
  if (verbose) cat("Estimating nuisance functions; \n")
  
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
  
  ## 2.) Estimate AIPW scores in validation sample via cross-fitting if necessary. 
  if (any(strategies == "AIPW")) {
    if (verbose) cat("Estimating AIPW scores; \n")
    scores_val <- aggTrees::dr_scores(Y_val, D_val, X_val)
  } else {
    scores_val <- NULL
  }
  
  ## 3.) Estimate BLP, GATES, and RATE in validation sample.
  if (verbose) cat("BLP estimation; \n")
  blp_results <- blp_estimation(Y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores_val, strategies, denoising)
  
  if (verbose) cat("GATES estimation; \n")
  gates_results <- gates_estimation(Y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores_val, n_groups, strategies, denoising)
  
  if (verbose) cat("RATE estimation; \n")
  rate_results <- rate_estimation(cates_val, scores_val, beneficial, n_boot)
  
  ## 4.) Output.
  if (verbose) cat("Output. \n\n")
  out <- list("BLP" = blp_results, "GATES" = gates_results, "RATE" = rate_results)
  class(out) <- "evaluCATE"
  return(out)
}
