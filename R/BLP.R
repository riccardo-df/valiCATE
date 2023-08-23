#' BLP Estimation
#'
#' Estimates the best linear predictor of the actual CATEs using the estimated CATEs.
#'
#' @param y Observed outcomes.
#' @param D Treatment indicator.
#' @param cates Estimated CATEs. CATEs must be estimated with different observations than those in \code{y} and \code{D}.
#' @param pscore Propensity scores. If unknown, they must be estimated using different observations than those in \code{y} and \code{D}. 
#' @param mu Estimated regression function. It must be estimated with different observations than those in \code{y} and \code{D}. 
#' @param mu0 Estimated regression function for control units. It must be estimated with different observations than those in \code{y} and \code{D}.
#' @param mu1 Estimated regression function for treated units. It must be estimated with different observations than those in \code{y} and \code{D}. 
#' @param scores Estimated doubly-robust scores. They must be estimated via K-fold cross-fitting. 
#'
#' @return
#' A list of fitted models.
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
#' ## CATEs and nuisance functions estimation.
#' library(grf)
#' 
#' cates_forest <- causal_forest(X_tr, y_tr, D_tr) # We use only the training sample.
#' mu_forest <- regression_forest(X_tr, y_tr)
#' mu0_forest <- regression_forest(X_tr[D_tr == 0, ], y_tr[D_tr == 0])
#' mu1_forest <- regression_forest(X_tr[D_tr == 1, ], y_tr[D_tr == 1])
#' 
#' cates_val <- predict(cates_forest, X_val)$predictions # We predict on the validation sample.
#' mu_val <- predict(mu_forest, X_val)$predictions
#' mu0_val <- predict(mu0_forest, X_val)$predictions
#' mu1_val <- predict(mu1_forest, X_val)$predictions
#' 
#' library(aggTrees)
#' scores_val <- dr_scores(y_val, D_val, X_val)
#' 
#' ## BLP estimation. Here, we know true pscores. Otherwise, estimate them in training sample.
#' pscore_val <- rep(0.5, length(y_val)) 
#' blp_results <- blp_estimation(y_val, D_val, cates_val, 
#'                               pscore_val, mu_val, mu0_val, mu1_val, scores_val)
#'
#' ## Compare true ATE vs estimated ATE (i.e., \hat{\beta}_1).
#' cat("True ATE      : ", round(mean(mu1 - mu0), 3), " 
#' Estimated ATE 
#'     - wr_none :  ", round(blp_results$wr_none$coefficients["beta1"], 3), "
#'     - wr_cddf1: ", round(blp_results$wr_cddf1$coefficients["beta1"], 3), "
#'     - wr_cddf2: ", round(blp_results$wr_cddf2$coefficients["beta1"], 3), " 
#'     - wr_mck1 : ", round(blp_results$wr_mck1$coefficients["beta1"], 3), " 
#'     - ht_none : ", round(blp_results$ht_none$coefficients["beta1"], 3), " 
#'     - ht_cddf1: ", round(blp_results$ht_cddf1$coefficients["beta1"], 3), " 
#'     - ht_cddf2: ", round(blp_results$ht_cddf2$coefficients["beta1"], 3), " 
#'     - ht_mck1 : ", round(blp_results$ht_mck1$coefficients["beta1"], 3), " 
#'     - ht_mck2 : ", round(blp_results$ht_mck2$coefficients["beta1"], 3), " 
#'     - ht_mck3 : ", round(blp_results$ht_mck3$coefficients["beta1"], 3), " 
#'     - aipw    : ", round(blp_results$aipw$coefficients["beta1"], 3), sep = "")
#' 
#' ## Compare true "quality" of CATEs estimated vs estimated quality (i.e., \hat{\beta}_2).
#' ## We can do this because we know that, by DGP, we have heterogeneous effects.
#' cat("True quality      : ", cor(mu1[!train_idx] - mu0[!train_idx], cates_val), " 
#' Estimated quality 
#'     - wr_none     : ", round(blp_results$wr_none$coefficients["beta2"], 3), "
#'     - wr_cddf1    : ", round(blp_results$wr_cddf1$coefficients["beta2"], 3), "
#'     - wr_cddf2    : ", round(blp_results$wr_cddf2$coefficients["beta2"], 3), " 
#'     - wr_mck1     : ", round(blp_results$wr_mck1$coefficients["beta2"], 3), " 
#'     - ht_none     : ", round(blp_results$ht_none$coefficients["beta2"], 3), " 
#'     - ht_cddf1    : ", round(blp_results$ht_cddf1$coefficients["beta2"], 3), " 
#'     - ht_cddf2    : ", round(blp_results$ht_cddf2$coefficients["beta2"], 3), " 
#'     - ht_mck1     : ", round(blp_results$ht_mck1$coefficients["beta2"], 3), " 
#'     - ht_mck2     : ", round(blp_results$ht_mck2$coefficients["beta2"], 3), " 
#'     - ht_mck3     : ", round(blp_results$ht_mck3$coefficients["beta2"], 3), " 
#'     - aipw        : ", round(blp_results$aipw$coefficients["beta2"], 3), sep = "")
#'
#' @md
#' @details
#' \code{\link{blp_estimation}} estimates the best linear predictor (BLP) of the actual CATEs using the estimated CATEs. To this end, the user must provide observations on the outcomes and the treatment status of units in 
#' the validation sample, as well as their estimated cates and nuisance functions. These estimates must be obtained by using only observations from the training sample (see the example section below).\cr
#' 
#' The BLP is estimated using three different strategies, all involving fitting suitable linear models. Check the \href{PUT LINK HERE}{online vignette} for details.\cr
#' 
#' Standard errors are estimated using the Eicker-Huber-White estimator.
#'
#' @import estimatr stats
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
blp_estimation <- function(y, D, cates, pscore, mu, mu0, mu1, scores) {
  ## 1.) Construct covariates 
  wr_weights <- (pscore * (1 - pscore))^(-1)
  D_residual <- D - pscore
  demeaned_cates <- cates - mean(cates)
  interaction_D_cates <- D_residual * demeaned_cates
  interaction_pscore_cates <- pscore * cates
  H <- D_residual * wr_weights
  HY <- H * y
  Hmu0 <- H * mu0
  Hpscore <- H * pscore
  Hinteraction_pscore_cates <- Hpscore * cates
  new_mck_covariate <- H * (1 - pscore) * cates
  Hmu0_pscore <- H * mu0 * pscore
  Hmu1_pscore <- H * mu1 * (1 - pscore)
  Hmu0_pscore_mu1_pscore <- Hmu0_pscore + Hmu1_pscore
  
  ## 2.) Fit linear models via OLS.
  wr_none_model <- estimatr::lm_robust(y ~ 0 + ., data = data.frame("y" = y, "beta1" = D_residual, "beta2" = interaction_D_cates), weights = wr_weights, se_type = "HC1") 
  wr_cddf1_model <- estimatr::lm_robust(y ~ 0 + ., data = data.frame("y" = y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu0" = mu0), weights = wr_weights, se_type = "HC1") 
  wr_cddf2_model <- estimatr::lm_robust(y ~ 0 + ., data = data.frame("y" = y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu0" = mu0, "constant" = rep(1, length(y)), "pscore" = pscore, "pscore.tauhat" = interaction_pscore_cates), weights = wr_weights, se_type = "HC1") 
  wr_mck1_model <- estimatr::lm_robust(y ~ 0 + ., data = data.frame("y" = y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu" = mu), weights = wr_weights, se_type = "HC1") 
  ht_none_model <- estimatr::lm_robust(Hy ~ 0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates), se_type = "HC1") 
  ht_cddf1_model <- estimatr::lm_robust(Hy ~ 0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates, "Hmu0" = Hmu0), se_type = "HC1") 
  ht_cddf2_model <- estimatr::lm_robust(Hy ~ 0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates, "Hmu0" = Hmu0, "Hpscore" = Hpscore, "Hpscore_tauhat" = Hinteraction_pscore_cates), se_type = "HC1") 
  ht_mck1_model <- estimatr::lm_robust(Hy ~ 0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates, "Hmu0.pscore" = Hmu0_pscore, "H1_pscore.tauhat" = new_mck_covariate), se_type = "HC1") 
  ht_mck2_model <- estimatr::lm_robust(Hy ~0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates, "Hpscore" = Hpscore, "Hmu0.pscore" = Hmu0_pscore, "Hmu1.1_pscore" = Hmu1_pscore), se_type = "HC1") 
  ht_mck3_model <- estimatr::lm_robust(Hy ~ 0 + ., data = data.frame("Hy" = HY, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates, "Hpscore" = Hpscore, "Hmu0.pscore+Hmu1.1_pscore" = Hmu0_pscore_mu1_pscore), se_type = "HC1") 
  aipw_model <- estimatr::lm_robust(aipw ~ 0 + ., data = data.frame("aipw" = scores, "beta1" = rep(1, length(y)), "beta2" = demeaned_cates), se_type = "HC1") 
  
  ## 3.) Output.
  out <- list("wr_none" = wr_none_model, "wr_cddf1" = wr_cddf1_model, "wr_cddf2" = wr_cddf2_model, "wr_mck1" = wr_mck1_model, 
              "ht_none" = ht_none_model, "ht_cddf1" = ht_cddf1_model, "ht_cddf2" = ht_cddf2_model, "ht_mck1" = ht_mck1_model, "ht_mck2" = ht_mck2_model, "ht_mck3" = ht_mck3_model,
              "aipw" = aipw_model)
  return(out)
}
