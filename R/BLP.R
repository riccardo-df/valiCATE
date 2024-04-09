#' BLP Estimation
#'
#' Estimates the best linear predictor (BLP) of the actual CATEs using the estimated CATEs.
#'
#' @param Y Observed outcomes.
#' @param D Treatment indicator.
#' @param cates CATE predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}.
#' @param pscore Propensity scores predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D} (unless the propensity score is known, in which case we provide the true values).
#' @param mu Conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param mu0 Control units' conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param mu1 Treated units' conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting using the same observations as in \code{Y} and \code{D}. 
#' @param strategies Character vector controlling the identification strategies to implement for BLP and GATES. Admitted values are \code{"WR"} (weighted residuals), \code{"HT"} (Horwitz-Thompson), and \code{"AIPW"} (augmented inverse-probability weighting).
#' @param denoising Character vector controlling if and which additional covariates to include in the regressions to reduce the variance of the estimation. Admitted values are \code{"none"}, \code{"cddf1"}, \code{"cddf2"}, \code{"mck1"}, \code{"mck2"}, and \code{"mck3"}.
#' 
#' @return
#' A list of fitted models as \code{\link[estimatr]{lm_robust}} objects.
#'
#' @details
#' To estimate the BLP of the actual CATEs using the estimated CATEs, the user must provide observations on the outcomes and the treatment status of units in 
#' the validation sample, as well as their estimated cates and nuisance functions. Be careful, as these estimates must be obtained using only observations from the training sample (see the example section below).
#' Additionally, the user must provide doubly-robust scores estimated in the validation sample using K-fold cross fitting.\cr
#' 
#' The BLP is estimated using three different strategies, all involving fitting suitable linear models. For each of these strategis, different model specifications are considered that differ in additional and
#' optional covariates that can be included in the regressions to reduce the estimation variance. Check the online \href{https://riccardo-df.github.io/evaluCATE/articles/evaluCATE-short-tutorial.html}{short tutorial}
#' and \href{https://riccardo-df.github.io/evaluCATE/articles/denoising.html}{denoising vignette} for details.\cr
#' 
#' Standard errors are estimated using the Eicker-Huber-White estimator.
#'
#' @import estimatr
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{gates_estimation}}, \code{\link{toc_estimation}}, \code{\link{rate_estimation}}
#'
#' @keywords internal
blp_estimation <- function(Y, D, cates, pscore, mu, mu0, mu1, scores, strategies, denoising) {
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
  
  ## 2.) Fit linear models. Save them in a list for flexible output.
  model_list <- list()
  counter <- 1
  
  if (any(strategies == "WR")) {
    if (any(denoising == "none")) {
      wr_none_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_none_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_none" 
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf1")) {
      wr_cddf1_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu0" = mu0)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf1_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_cddf1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf2")) {
      wr_cddf2_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "constant" = rep(1, length(Y)), "mu0" = mu0, "pscore" = pscore, "pscore.tauhat" = interaction_pscore_cates)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf2_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_cddf2"  
      counter <- counter + 1
    }  
    
    if (any(denoising == "mck1")) {
      wr_mck1_dta <- data.frame("Y" = Y, "beta1" = D_residual, "beta2" = interaction_D_cates, "mu" = mu)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_mck1_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_mck1"  
      counter <- counter + 1
    } 
  }
  
  if (any(strategies == "HT")) {
    if (any(denoising == "none")) {
      ht_none_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_none_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_none"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf1")) {
      ht_cddf1_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_cddf1"  
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf2")) {
      ht_cddf2_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.pscore" = Hpscore, "H.pscore.tauhat" = Hinteraction_pscore_cates)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_cddf2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck1")) {
      ht_mck1_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.1_pscore.tauhat" = new_mck_covariate)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck2")) {
      ht_mck2_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "H.mu1.1_pscore" = Hmu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck3")) {
      ht_mck3_dta <- data.frame("HY" = HY, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck3_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck3"  
      counter <- counter + 1
    } 
  }
  
  if (any(strategies == "AIPW")) {
    if (any(denoising == "none")) {
      aipw_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_none"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf1")) {
      aipw_cddf1_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_cddf1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_cddf1"  
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf2")) {
      aipw_cddf2_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.pscore" = Hpscore, "H.pscore.tauhat" = Hinteraction_pscore_cates)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_cddf2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_cddf2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck1")) {
      aipw_mck1_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.mu0" = Hmu0, "H.1_pscore.tauhat" = new_mck_covariate)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck2")) {
      aipw_mck2_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "H.mu1.1_pscore" = Hmu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck3")) {
      aipw_mck3_dta <- data.frame("aipw" = scores, "beta1" = rep(1, length(Y)), "beta2" = demeaned_cates, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck3_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck3"  
      counter <- counter + 1
    } 
  } 

  ## 3.) Output.
  return(model_list)
}
