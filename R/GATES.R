#' GATES Estimation
#'
#' Estimates the group average treatment effects, with the groups formed by cutting the distribution of the estimated CATEs into K quantiles.
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
#' We rearrange the estimate GATES to obey the monotonicity property.
#'
#' @import estimatr stats
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
blp_regression <- function(y, D, cates, pscore, mu, mu0, mu1, scores) {
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
  wr_none_model <- estimatr::lm_robust(y ~ 0 + D_residual + interaction_D_cates, weights = wr_weights, se_type = "HC1") 
  wr_cddf1_model <- estimatr::lm_robust(y ~ 0 + D_residual + interaction_D_cates + mu0 , weights = wr_weights, se_type = "HC1") 
  wr_cddf2_model <- estimatr::lm_robust(y ~ D_residual + interaction_D_cates + mu0 + interaction_pscore_cates, weights = wr_weights, se_type = "HC1") 
  wr_mck1_model <- estimatr::lm_robust(y_test ~ 0 + D_residual + interaction_D_cates + mu, weights = wr_weights, se_type = "HC1") 
  ht_none_model <- estimatr::lm_robust(HY ~ demeaned_cates, se_type = "HC1") 
  ht_cddf1_model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0, se_type = "HC1") 
  ht_cddf2_model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0 + Hpscore + Hinteraction_pscore_cates, se_type = "HC1") 
  ht_mck1_model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0 + new_mck_covariate, se_type = "HC1") 
  ht_mck2_model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0_pscore + Hmu1_pscore + Hpscore, se_type = "HC1") 
  ht_mck3_model <- estimatr::lm_robust(HY ~ demeaned_cates + Hpscore + Hmu0_pscore_mu1_pscore, se_type = "HC1") 
  aipw_model <- estimatr::lm_robust(scores ~ demeaned_cates, se_type = "HC1") 
  
  ## 3.) Output.
  out <- list("wr_none" = wr_none_model, "wr_cddf1" = wr_cddf1_model, "wr_cddf2" = wr_cddf2_model, "wr_mck1" = wr_mck1_model, 
              "ht_none" = ht_none_model, "ht_cddf1" = ht_cddf1_model, "ht_cddf2" = ht_cddf2_model, "ht_mck1" = ht_mck1_model, "ht_mck2" = ht_mck2_model, "ht_mck3" = ht_mck3_model,
              "aipw" = aipw_model)
  return(out)
}
