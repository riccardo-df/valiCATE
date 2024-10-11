#' Testing for Effect Heterogeneity using the estimated GATES
#'
#' Tests for effect heterogeneity using the estimated parametric model for the sorted group average treatment effects (GATES).
#'
#' @param model Estimated parametric model, as one of those returned by \code{\link{gates_estimation}}.
#'
#' @return
#' A list with the p-values for the three hypotheses tested.
#'
#' @details
#' \code{model} must consist of a \code{lm_robust} object where the coefficients identifying the GATES must be called \code{"group1"}, \code{"group2"}, and so on.\cr
#' 
#' Three distinct hypotheses of effect heterogeneity are tested: whether all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, 
#' and whether the differences in the GATES across all pairs of groups are zero. For the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure 
#' and report the median of the adjusted p-values.
#'
#' @import estimatr GenericML evalITR
#' @importFrom car linearHypothesis
#' @importFrom stats p.adjust
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{gates_estimation}}
#'
#' @keywords internal
test_parametric_gates <- function(model) {
  ## 0.) Extract number of groups.
  regressors <- names(model$coefficients)
  n_groups <- sum(grepl("group", regressors))
  
  ## 1.) Test equality of coefficients.
  all_equal <- car::linearHypothesis(model, paste0("group1 = group", 2:n_groups), test = "F")$`Pr(>F)`[2]
  
  ## 2.) Test if most and least affected groups differ.
  top_group <- paste0("group", n_groups)
  bottom_group <- "group1"
  
  largest_difference <- car::linearHypothesis(model, paste(top_group, "-", bottom_group, "= 0"), test = "F")$`Pr(>F)`[2]
  
  ## 3.) Test all pairwise differences. Adjust p-values using Holm's procedure.
  p_values <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(p_values) <- paste0("group", seq_len(n_groups))
  colnames(p_values) <- paste0("group", seq_len(n_groups))
  
  for (i in seq_len(n_groups-1)) {
    for (j in seq_len(n_groups)[-c(1:i)]) {
      group_i <- paste0("group", i)
      group_j <- paste0("group", j)
      
      p_values[j, i] <- car::linearHypothesis(model, paste(group_i, "=", group_j), test = "F")$`Pr(>F)`[2]
    }
  }
  
  p_values_vec <- c(p_values) # First column, then second column, then third column ...
  p_values_holm_vec <- stats::p.adjust(p_values_vec, method = "holm")
  p_values_holm <- matrix(p_values_holm_vec, nrow = n_groups, ncol = n_groups)
  
  pairwise_differences <- p_values_holm
  
  ## 4.) Output.
  return(list("all_equal" = all_equal, "largest_difference" = largest_difference, "pairwise_differences" = pairwise_differences))
}


#' GATES Estimation
#'
#' Estimates the sorted group average treatment effects (GATES), with the groups formed by cutting the distribution of the estimated CATEs into K quantiles.
#'
#' @param Y Observed outcomes.
#' @param D Treatment indicator.
#' @param cates CATE predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}.
#' @param pscore Propensity scores predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D} (unless the propensity score is known, in which case we provide the true values).
#' @param mu Conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param mu0 Control units' conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param mu1 Treated units' conditional mean predictions. Must be produced by a model estimated using different observations than those in \code{Y} and \code{D}
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting using the same observations as in \code{Y} and \code{D}. 
#' @param n_groups Number of groups to be formed.
#' @param strategies Character vector controlling the identification strategies to implement for BLP and GATES. Admitted values are \code{"WR"} (weighted residuals), \code{"HT"} (Horwitz-Thompson), and \code{"AIPW"} (augmented inverse-probability weighting).
#' @param denoising Character vector controlling if and which additional covariates to include in the regressions to reduce the variance of the estimation. Admitted values are \code{"none"}, \code{"cddf1"}, \code{"cddf2"}, \code{"mck1"}, \code{"mck2"}, and \code{"mck3"}.
#' 
#' @return
#' A list of fitted models as \code{\link[estimatr]{lm_robust}} objects and a data frame with point estimates and standard errors for the nonparametric estimator.
#'
#' @details
#' To estimate the GATES, the user must provide observations on the outcomes and the treatment status of units in 
#' the validation sample, as well as their estimated cates and nuisance functions. Be careful, as these estimates must be obtained using only observations from the training sample 
#' (see the example section below). Additionally, the user must provide doubly-robust scores estimated in the validation sample using K-fold cross fitting.\cr
#' 
#' Groups are constructed by cutting the distribution of \code{cates} into \code{n_groups} quantiles. If this leads to one or more groups composed of only treated or only control units, the function raises an error.\cr
#' 
#' The GATES are estimated using four different strategies: three involving fitting suitable linear models, and one nonparametric approach. Check the online 
#' \href{https://riccardo-df.github.io/valiCATE/articles/valiCATE-short-tutorial.html}{short tutorial} for details.\cr
#' 
#' For the linear models, standard errors are estimated using the Eicker-Huber-White estimator. These standard errors are then used to test three distinct hypotheses of effect heterogeneity: whether
#' all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, and whether the differences in the GATES across all pairs of groups are zero.
#' For the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure and report the median of the adjusted p-values. The nonparametric approach tests only the first
#' of these hypotheses. Check \href{https://riccardo-df.github.io/valiCATE/articles/hypotheses-testing.html}{hypotheses testing vignette} for details.\cr
#' 
#' Each strategy based on linear models supports different model specifications 
#' that differ in additional and optional covariates that can be included in the regressions to reduce the estimation variance.
#' Check \href{https://riccardo-df.github.io/valiCATE/articles/denoising.html}{denoising vignette} for details.
#'
#' @import estimatr GenericML evalITR 
#' @importFrom stats median
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{blp_estimation}}, \code{\link{toc_estimation}}, \code{\link{rate_estimation}}, \code{\link{test_parametric_gates}}
#'
#' @keywords internal
gates_estimation <- function(Y, D, cates, pscore, mu, mu0, mu1, scores, n_groups = 5, strategies, denoising) {
  ## 0.) Handling inputs and checks.
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)

  ## 1.) Generate groups by cutting the distribution of the CATEs. If we have homogeneous groups (i.e., only treated or only control units), raise an error.  
  cuts <- seq(0, 1, length = n_groups+1)[-c(1, n_groups+1)]
  group_indicators <- GenericML::quantile_group(cates, cutoffs = cuts)
  class(group_indicators) <- "numeric"
  colnames(group_indicators) <- paste0("group", 1:dim(group_indicators)[2])
  
  out_condition <- FALSE
  for (g in seq_len(dim(group_indicators)[2])) {
    if (sum(D * group_indicators[, g]) %in% c(0, sum(group_indicators[, g]))) out_condition <- TRUE 
  }
  
  if (out_condition) stop("We have one or more homogeneous groups. Please try a different 'n_groups' or a different sample split.", call. = FALSE)
  
  ## 2.) Construct covariates.
  wr_weights <- (pscore * (1 - pscore))^(-1) 
  D_residual <- D - pscore
  D_residual_interaction <- D_residual * group_indicators
  colnames(D_residual_interaction) <- paste0("group", 1:dim(D_residual_interaction)[2])
  pscore_interaction <- pscore * group_indicators
  colnames(pscore_interaction) <- paste0("pscore", 1:dim(D_residual_interaction)[2])
  H <- D_residual * wr_weights
  HY <- H * Y
  Hmu0 <- H * mu0
  Hpscore_interaction <- H * pscore_interaction
  colnames(Hpscore_interaction) <- paste0("H.pscore", 1:dim(group_indicators)[2])
  new_mck_covariate <- H * (1 - pscore) * cates
  Hpscore <- H * pscore
  Hmu0_pscore <- H * mu0 * pscore
  Hmu1_pscore <- H * mu1 * (1 - pscore)
  Hmu0_pscore_mu1_pscore <- Hmu0_pscore + Hmu1_pscore
  
  ## 3.) Fit linear models. Save them in a list for flexible output.
  model_list <- list()
  counter <- 1
  
  if (any(strategies == "WR")) {
    if (any(denoising == "none")) {
      wr_none_dta <- data.frame("Y" = Y, D_residual_interaction)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_none_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_none" 
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf1")) {
      wr_cddf1_dta <- data.frame("Y" = Y, D_residual_interaction, mu0)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf1_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_cddf1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf2")) {
      wr_cddf2_dta <- data.frame("Y" = Y, D_residual_interaction, mu0, pscore_interaction)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_cddf2_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_cddf2"  
      counter <- counter + 1
    }  
    
    if (any(denoising == "mck1")) {
      wr_mck1_dta <- data.frame("Y" = Y, D_residual_interaction, mu)
      model_list[[counter]] <- estimatr::lm_robust(Y ~ 0 + ., wr_mck1_dta, weights = wr_weights, se_type = "HC1") 
      names(model_list)[[counter]] <- "wr_mck1"  
      counter <- counter + 1
    } 
  }
  
  if (any(strategies == "HT")) {
    if (any(denoising == "none")) {
      ht_none_dta <- data.frame("HY" = HY, group_indicators)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_none_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_none"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf1")) {
      ht_cddf1_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_cddf1"  
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf2")) {
      ht_cddf2_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0, Hpscore_interaction)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_cddf2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_cddf2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck1")) {
      ht_mck1_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0, "H.1-pscore.tauhat" = new_mck_covariate)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck2")) {
      ht_mck2_dta <- data.frame("HY" = HY, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "H.mu1.1_pscore" = Hmu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck3")) {
      ht_mck3_dta <- data.frame("HY" = HY, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(HY ~ 0 + ., ht_mck3_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "ht_mck3"  
      counter <- counter + 1
    } 
  }
  
  if (any(strategies == "AIPW")) {
    if (any(denoising == "none")) {
      aipw_dta <- data.frame("aipw" = scores, group_indicators)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_none"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "cddf1")) {
      aipw_cddf1_dta <- data.frame("aipw" = scores, group_indicators, "H.mu0" = Hmu0)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_cddf1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_cddf1"  
      counter <- counter + 1
    }
    
    if (any(denoising == "cddf2")) {
      aipw_cddf2_dta <- data.frame("aipw" = scores, group_indicators, "H.mu0" = Hmu0, Hpscore_interaction)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_cddf2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_cddf2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck1")) {
      aipw_mck1_dta <- data.frame("aipw" = scores, group_indicators, "H.mu0" = Hmu0, "H.1-pscore.tauhat" = new_mck_covariate)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck1_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck1"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck2")) {
      aipw_mck2_dta <- data.frame("aipw" = scores, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "H.mu1.1_pscore" = Hmu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck2_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck2"  
      counter <- counter + 1
    } 
    
    if (any(denoising == "mck3")) {
      aipw_mck3_dta <- data.frame("aipw" = scores, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
      model_list[[counter]] <- estimatr::lm_robust(aipw ~ 0 + ., aipw_mck3_dta, se_type = "HC1") 
      names(model_list)[[counter]] <- "aipw_mck3"  
      counter <- counter + 1
    } 
  } 
  
  ## 4.) Nonparametric estimator.
  imai_li <- evalITR::GATE(D, cates, Y, n_groups)
  imai_li_results <- data.frame("group" = 1:n_groups, "GATE" = imai_li$gate, "SE" = imai_li$sd)
  
  ## 5.) Hypotheses testing.
  parametric_models <- model_list
  new_parametric_models <- list()
  counter2 <- 1
  
  for (model in parametric_models) {
    temp_model <- model
    tests <- test_parametric_gates(temp_model)
    
    temp_model$p.value.gates.all.equal <- tests$all_equal
    temp_model$p.value.gates.largest.difference <- tests$largest_difference
    temp_model$p.value.gates.pairwise.differences.median <- stats::median(tests$pairwise_differences, na.rm = TRUE)
    
    new_parametric_models[[counter2]] <- temp_model 
    names(new_parametric_models)[[counter2]] <-  names(parametric_models)[[counter2]] 
    counter2 <- counter2 + 1
  }
  
  imai_li_results$p.value.gates.all.equal <- as.numeric(evalITR::het.test(D, cates, Y, n_groups)$pval)
  imai_li_results$p.value.gates.largest.difference <- NA
  imai_li_results$p.value.gates.pairwise.differences.median <- NA
  
  ## 7.) Output.
  out <- append(new_parametric_models, list(imai_li_results))
  names(out)[[counter2]] <- c("imai_li")
  return(out)
}
