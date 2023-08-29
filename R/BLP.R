#' BLP Estimation
#'
#' Estimates the best linear predictor (BLP) of the actual CATEs using the estimated CATEs.
#'
#' @param Y Observed outcomes.
#' @param D Treatment indicator.
#' @param cates Estimated CATEs. Must be estimated with different observations than those in \code{Y} and \code{D}.
#' @param pscore Propensity scores. If unknown, they must be estimated using different observations than those in \code{Y} and \code{D}. 
#' @param mu Estimated regression function. Must be estimated with different observations than those in \code{Y} and \code{D}. 
#' @param mu0 Estimated regression function for control units. Must be estimated with different observations than those in \code{Y} and \code{D}.
#' @param mu1 Estimated regression function for treated units. Must be estimated with different observations than those in \code{Y} and \code{D}. 
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting with the same observations as in \code{Y} and \code{D}. 
#'
#' @return
#' A list of fitted models as \code{\link[estimatr]{lm_robust}} objects.
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
#' ## CATEs and nuisance functions estimation.
#' ## We use only the training sample for estimation.
#' ## We predict on the validation sample.
#' library(grf)
#' 
#' cates_forest <- causal_forest(X_tr, Y_tr, D_tr) 
#' mu_forest <- regression_forest(X_tr, Y_tr)
#' mu0_forest <- regression_forest(X_tr[D_tr == 0, ], Y_tr[D_tr == 0])
#' mu1_forest <- regression_forest(X_tr[D_tr == 1, ], Y_tr[D_tr == 1])
#' 
#' cates_val <- predict(cates_forest, X_val)$predictions #
#' mu_val <- predict(mu_forest, X_val)$predictions
#' mu0_val <- predict(mu0_forest, X_val)$predictions
#' mu1_val <- predict(mu1_forest, X_val)$predictions
#' 
#' ## AIPW scores estimation.
#' ## Cross-fitting on the validation sample.
#' library(aggTrees)
#' scores_val <- dr_scores(Y_val, D_val, X_val)
#' 
#' ## BLP estimation. 
#' ## Here, we know true pscores. Otherwise, estimate them in training sample.
#' pscore_val <- rep(0.5, length(Y_val)) 
#' blp_results <- blp_estimation(Y_val, D_val, cates_val, 
#'                               pscore_val, mu_val, mu0_val, mu1_val, 
#'                               scores_val)
#'
#' ## Compare true ATE vs estimated ATE.
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
#' ## Compare true "quality" of estimated CATEs vs estimated quality.
#' ## (We can do this because we know that, by DGP, we have heterogeneous effects.)
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
#' @details
#' \code{\link{blp_estimation}} estimates the BLP of the actual CATEs using the estimated CATEs. To this end, the user must provide observations on the outcomes and the treatment status of units in 
#' the validation sample, as well as their estimated cates and nuisance functions. These estimates must be obtained by using only observations from the training sample (see the example section below).\cr
#' 
#' The BLP is estimated using three different strategies, all involving fitting suitable linear models. Check the \href{https://riccardo-df.github.io/evaluCATE/articles/evalue-cates-short-tutorial.html}{online vignette} 
#' for details.\cr
#' 
#' Standard errors are estimated using the Eicker-Huber-White estimator.
#'
#' @import estimatr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
blp_estimation <- function(Y, D, cates, pscore, mu, mu0, mu1, scores) {
  ## 1.) Construct covariates 
  wr_weights <- (pscore * (1 - pscore))^(-1)
  D_residual <- D - pscore
  demeaned_cates <- cates - mean(cates)
  interaction_D_cates <- D_residual * demeaned_cates
  interaction_pscore_cates <- pscore * cates
  H <- D_residual * wr_weights
  HY <- H * Y
  Hmu0 <- H * mu0
  Hpscore <- H * pscore
  Hinteraction_pscore_cates <- Hpscore * cates
  new_mck_covariate <- H * (1 - pscore) * cates
  Hmu0_pscore <- H * mu0 * pscore
  Hmu1_pscore <- H * mu1 * (1 - pscore)
  Hmu0_pscore_mu1_pscore <- Hmu0_pscore + Hmu1_pscore
  
  ## 2.) Define specifications.
  wr_none_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates)
  wr_cddf1_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu0" = mu0)
  wr_cddf2_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "constant" = rep(1, length(Y)), "mu0" = mu0, "pscore" = pscore, "pscore.tauhat" = interaction_pscore_cates)
  wr_mck1_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu" = mu)
  
  ht_none_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates)
  ht_cddf1_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0)
  ht_cddf2_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.pscore" = Hpscore, "H.pscore.tauhat" = Hinteraction_pscore_cates)
  ht_mck1_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.1_pscore.tauhat" = new_mck_covariate)
  ht_mck2_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "H.mu1.1_pscore" = Hmu1_pscore)
  ht_mck3_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
  
  aipw_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates)
  
  ## 3.) Fit linear models.
  wr_none_model <- estimatr::lm_robust(Y ~ 0 + ., wr_none_dta, weights = wr_weights, se_type = "HC1") 
  wr_cddf1_model <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf1_dta, weights = wr_weights, se_type = "HC1") 
  wr_cddf2_model <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf2_dta, weights = wr_weights, se_type = "HC1") 
  wr_mck1_model <- estimatr::lm_robust(Y ~ 0 + ., wr_mck1_dta, weights = wr_weights, se_type = "HC1") 
  
  ht_none_model <- estimatr::lm_robust(HY ~ 0 + ., ht_none_dta, se_type = "HC1") 
  ht_cddf1_model <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf1_dta, se_type = "HC1") 
  ht_cddf2_model <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf2_dta, se_type = "HC1") 
  ht_mck1_model <- estimatr::lm_robust(HY ~ 0 + ., ht_mck1_dta, se_type = "HC1") 
  ht_mck2_model <- estimatr::lm_robust(HY ~0 + ., ht_mck2_dta, se_type = "HC1") 
  ht_mck3_model <- estimatr::lm_robust(HY ~ 0 + ., ht_mck3_dta, se_type = "HC1") 
  
  aipw_model <- estimatr::lm_robust(aipw ~ 0 + ., aipw_dta, se_type = "HC1") 
  
  ## 4.) Output.
  out <- list("wr_none" = wr_none_model, "wr_cddf1" = wr_cddf1_model, "wr_cddf2" = wr_cddf2_model, "wr_mck1" = wr_mck1_model, 
              "ht_none" = ht_none_model, "ht_cddf1" = ht_cddf1_model, "ht_cddf2" = ht_cddf2_model, "ht_mck1" = ht_mck1_model, "ht_mck2" = ht_mck2_model, "ht_mck3" = ht_mck3_model,
              "aipw" = aipw_model)
  return(out)
}
