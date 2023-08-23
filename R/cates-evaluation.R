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
#' @param n_groups Number of groups to be formed.
#' @param verbose Logical, set to FALSE to prevent the function from printing the progresses.
#'
#' @return
#' Return.
#'
#' @examples
#' ## Generate data.
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
#' y <- mu0 + D * (mu1 - mu0) + rnorm(n)
#' 
#' ## Sample split.
#' train_idx <- sample(c(TRUE, FALSE), length(y), replace = TRUE)
#' 
#' X_tr <- X[train_idx, ]
#' X_val <- X[!train_idx, ]
#' 
#' D_tr <- D[train_idx]
#' D_val <- D[!train_idx]
#' 
#' y_tr <- y[train_idx]
#' y_val <- y[!train_idx]
#' 
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, y_tr, D_tr) # We use only the training sample.
#' cates <- predict(forest, X)$predictions # We predict on the whole sample.
#' 
#' ## CATEs evaluation.  
#' pscore <- rep(0.5, length(y_val))
#' n_groups <- 5
#' 
#' evaluation <- evalue_cates(y, D, X, cates, train_idx, pscore, n_groups)
#' 
#' ## Compare results for a given model.
#' blp_model <- evaluation$BLP$aipw
#' gates_model <- evaluation$GATES$aipw
#' 
#' # True ATE vs estimated ATE.
#' cat("True ATE      : ", round(mean(mu1 - mu0), 3), " 
#' Estimated ATE : ", round(blp_model$coefficients["beta1"], 3), " [", round(blp_model$conf.low["beta1"], 3), ", ", round(blp_model$conf.high["beta1"], 3), "]", sep = "")
#' 
#' # True "quality" of estimated CATEs vs estimated quality.
#' # (We can do this because we know that, by DGP, we have heterogeneous effects.)
#' cat("True quality      : ", cor(mu1[!train_idx] - mu0[!train_idx], cates[!train_idx]), " 
#' Estimated quality : ", round(blp_model$coefficients["beta2"], 3), " [", round(blp_model$conf.low["beta2"], 3), ", ", round(blp_model$conf.high["beta2"], 3), "]", sep = "") 
#' 
#' # True GATES with estimated GATES.
#' cuts <- seq(0, 1, length = n_groups+1)[-c(1, n_groups+1)]
#' group_indicators <- GenericML::quantile_group(cates[!train_idx], cutoffs = cuts)
#' colnames(group_indicators) <- paste0(1:n_groups)
#' true_gates <- apply(group_indicators, 2, function(x) {mean(cates[!train_idx][x])})
#' 
#' library(ggplot2)
#' 
#' plot_dta <- data.frame("group" = 1:n_groups, "true_gate" = true_gates, 
#'                        "estimated_gate" = gates_model$coefficients[1:n_groups], 
#'                        "se" = gates_model$std.error[1:n_groups])
#' 
#' ggplot(plot_dta, aes(x = group, y = true_gate)) +
#'   geom_point(aes(color = "True")) +
#'   geom_point(aes(y = estimated_gate, color = "Estimated")) +
#'   geom_errorbar(aes(x = group, ymin = estimated_gate - 1.96 * se, ymax = estimated_gate + 1.96 * se), color = "black") +
#'   xlab("Group") + ylab("GATES") + 
#'   scale_color_manual(name = "", breaks = c("True", "Estimated"), values = c("True" = "tomato", "Estimated" = "dodgerblue")) +
#'   theme_bw() + 
#'   theme(legend.position = c(0.2, 0.85))
#'
#' @md
#' @details
#' \code{\link{evalue_cates}} targets the estimation of the best linear predictor (BLP) of the actual CATEs using the estimated CATEs and of the sorted group average treatment effects (GATES).
#' To this end, the user must provide observations on the outcomes, the treatment status, and the covariates of units in the whole sample, as well as their estimated CATEs. Be careful,
#' as the CATEs must be estimated only with part of the sample, which we call the training sample (see the example section below).\cr
#' 
#' To let the function know which observations were used for the CATEs estimation, the user must also provide a logical vector with the \code{TRUE}s denoting observations in the 
#' training sample. This way, \code{\link{evalue_cates}} knows which observations to use to post-process the CATEs estimates.\cr
#' 
#' \code{\link{evalue_cates}} implements a number of strategies to estimate the BLP and the GATES. Most of them involve fitting a suitable linear model. The linear models differ according to the
#' different strategies. Furthermore, for each strategy, there exist various sets of constructed covariates that one can add to reduce the variance of the estimation. \code{\link{evalue_cates}}
#' fits and returns all these possible models. Check the online \href{https://riccardo-df.github.io/evalueCATE/articles/evalue-cates-short-tutorial.html}{short tutorial} for details.\cr 
#' 
#' Some of these models involve covariates that depend on particular nuisance functions. These functions are estimated internally via honest \code{\link[grf]{regression_forest}}s.
#' Check the online \href{https://riccardo-df.github.io/evalueCATE/articles/denoising.html}{denoising vignette} for details about these covariates.\cr
#' 
#' For the linear models, standard errors are estimated using the Eicker-Huber-White estimator.\cr
#' 
#' To estimate the BLP and GATES using the AIPW strategy, doubly-robust scores are estimated internally using the validation sample via 5-fold cross fitting and honest regression 
#' forests (see the \code{\link[aggTrees]{dr_scores}} function).\cr
#' 
#' The estimated GATES are sorted to enforce monotonicity.
#'
#' @import grf
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
evalue_cates <- function(y, D, X, cates, is_train, pscore = NULL, n_groups = 5, verbose = TRUE) {
  ## 0.) Handling inputs and checks.
  if (is.logical(D)) D <- as.numeric(D)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)
  if (!is.matrix(X) & !is.data.frame(X)) stop("Invalid 'X'. This must be either a matrix or a data frame.", call. = FALSE)
  if (!is.logical(is_train)) stop("Invalid 'is_train'. This must be a logical vector.", call. = FALSE)
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)

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
    pscore <- stats::predict(pscore_forest, X)$predictions
  }
  
  pscore_tr <- pscore[train_idx]
  pscore_val <- pscore[!train_idx]

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
  
  if (verbose) cat("GATES estimation; \n")
  gates_results <- gates_estimation(y_val, D_val, cates_val, pscore_val, mu_val, mu0_val, mu1_val, scores, n_groups)

  ## Output.
  if (verbose) cat("Output. \n\n")
  out <- list("BLP" = blp_results, "GATES" = gates_results)
  class_out <- "evalue_cates"
  return(out)
}
