#' Testing for Effect Heterogeneity using the estimated GATES
#'
#' Tests for effect heterogeneity using the estimated parametric model for the sorted group average treatment effects (GATES).
#'
#' @param model Estimated parametric model, as one of those returned by \code{\link{gates_estimation}}.
#'
#' @return
#' A list with the p-values for the three hypotheses tested.
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
#' cates_val <- predict(cates_forest, X_val)$predictions 
#' mu_val <- predict(mu_forest, X_val)$predictions
#' mu0_val <- predict(mu0_forest, X_val)$predictions
#' mu1_val <- predict(mu1_forest, X_val)$predictions
#' 
#' ## AIPW scores estimation.
#' ## Cross-fitting on the validation sample.
#' library(aggTrees)
#' scores_val <- dr_scores(Y_val, D_val, X_val)
#' 
#' ## GATEs estimation. Use default of five groups.
#' pscore_val <- rep(0.5, length(Y_val)) # We know true pscores.
#' gates_results <- gates_estimation(Y_val, D_val, cates_val, 
#'                                   pscore_val, mu_val, mu0_val, mu1_val, 
#'                                   scores_val)}
#'
#' @details
#' \code{model} must consist of a \code{lm_robust} object where the coefficients identifying the GATES must be called \code{"group1"}, \code{"group2"}, and so on.\cr
#' 
#' Three distinct hypotheses of effect heterogeneity are tested: whether all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, 
#' and whether the differences in the GATES across all pairs of groups are zero. For the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure 
#' and report the median of the adjusted p-values.
#'
#' @import estimatr GenericML evalITR car
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{gates_estimation}}
#'
#' @export
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
#' @param cates Estimated CATEs. CATEs must be estimated with different observations than those in \code{y} and \code{D}.
#' @param cates Estimated CATEs. Must be estimated with different observations than those in \code{Y} and \code{D}.
#' @param pscore Propensity scores. If unknown, must be estimated using different observations than those in \code{Y} and \code{D}. 
#' @param mu Estimated regression function. Must be estimated with different observations than those in \code{Y} and \code{D}. 
#' @param mu0 Estimated regression function for control units. Must be estimated with different observations than those in \code{Y} and \code{D}.
#' @param mu1 Estimated regression function for treated units. Must be estimated with different observations than those in \code{Y} and \code{D}. 
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting with the same observations as in \code{Y} and \code{D}.
#' @param n_groups Number of groups to be formed.
#'
#' @return
#' A list of fitted models as \code{\link[estimatr]{lm_robust}} objects and a data frame with point estimates and standard errors for the nonparametric estimator.
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
#' cates_val <- predict(cates_forest, X_val)$predictions 
#' mu_val <- predict(mu_forest, X_val)$predictions
#' mu0_val <- predict(mu0_forest, X_val)$predictions
#' mu1_val <- predict(mu1_forest, X_val)$predictions
#' 
#' ## AIPW scores estimation.
#' ## Cross-fitting on the validation sample.
#' library(aggTrees)
#' scores_val <- dr_scores(Y_val, D_val, X_val)
#' 
#' ## GATEs estimation. Use default of five groups.
#' pscore_val <- rep(0.5, length(Y_val)) # We know true pscores.
#' gates_results <- gates_estimation(Y_val, D_val, cates_val, 
#'                                   pscore_val, mu_val, mu0_val, mu1_val, 
#'                                   scores_val)}
#'
#' @details
#' To estimate the GATES, the user must provide observations on the outcomes and the treatment status of units in 
#' the validation sample, as well as their estimated cates and nuisance functions. Be careful, as these estimates must be obtained using only observations from the training sample 
#' (see the example section below). Additionally, the user must provide doubly-robust scores estimated in the validation sample using K-fold cross fitting.\cr
#' 
#' Groups are constructed by cutting the distribution of \code{cates} into \code{n_groups} quantiles. If this leads to one or more groups composed of only treated or only control units, the function raises an error.\cr
#' 
#' The GATES are estimated using four different strategies: three involving fitting suitable linear models, and one nonparametric approach. Check the online 
#' \href{https://riccardo-df.github.io/evaluCATE/articles/evalu-cates-short-tutorial.html}{short tutorial} for details.\cr
#' 
#' For the linear models, standard errors are estimated using the Eicker-Huber-White estimator. These standard errors are then used to test three distinct hypotheses of effect heterogeneity: whether
#' all GATES are equal to each other, whether the largest and the smallest GATES are different from each other, and whether the differences in the GATES across all pairs of groups are zero.
#' For the last test, we adjust p-values to account for multiple hypotheses testing using Holm's procedure and report the median of the adjusted p-values. The nonparametric approach tests only the first
#' of these hypotheses. Check \href{https://riccardo-df.github.io/evaluCATE/articles/hypotheses-testing.html}{hypotheses testing vignette} for details.\cr
#' 
#' Each strategy based on linear models supports different model specifications 
#' that differ in additional and optional covariates that can be included in the regressions to reduce the estimation variance.
#' Check \href{https://riccardo-df.github.io/evaluCATE/articles/denoising.html}{denoising vignette} for details.
#'
#' @import estimatr GenericML evalITR
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{blp_estimation}}, \code{\link{toc_estimation}}, \code{\link{rate_estimation}}, \code{\link{test_parametric_gates}}
#'
#' @export
gates_estimation <- function(Y, D, cates, pscore, mu, mu0, mu1, scores, n_groups = 5) {
  ## 0.) Handling inputs and checks.
  if (n_groups <= 1 | n_groups %% 1 != 0) stop("Invalid 'n_groups'. This must be an integer greater than 1.", call. = FALSE)

  ## 1.) Generate groups by cutting the distribution of the CATEs. If we have homogeneous groups (i.e., only treated or only control units), raise an error.  
  cuts <- seq(0, 1, length = n_groups+1)[-c(1, n_groups+1)]
  group_indicators <- GenericML::quantile_group(cates, cutoffs = cuts)
  class(group_indicators) <- "numeric"
  colnames(group_indicators) <- paste0("group", 1:dim(group_indicators)[2])
  
  out_condition <- FALSE
  for (g in seq_len(dim(group_indicators)[2])) {
    if (sum(D[group_indicators[, g] == 0]) %in% c(0, dim(group_indicators)[1])) out_condition <- TRUE 
  }
  
  if (out_condition) stop("We have one or more homogeneous groups. Please try a different 'k' or a different sample split.", call. = FALSE)
  
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
  
  ## 3.) Define specifications.
  wr_none_dta <- data.frame("Y" = Y, D_residual_interaction)
  wr_cddf1_dta <- data.frame("Y" = Y, D_residual_interaction, mu0)
  wr_cddf2_dta <- data.frame("Y" = Y, D_residual_interaction, mu0, pscore_interaction)
  wr_mck1_dta <- data.frame("Y" = Y, D_residual_interaction, mu)
  
  ht_none_dta <- data.frame("HY" = HY, group_indicators)
  ht_cddf1_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0)
  ht_cddf2_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0, Hpscore_interaction)
  ht_mck1_dta <- data.frame("HY" = HY, group_indicators, "H.mu0" = Hmu0, "H.1-pscore.tauhat" = new_mck_covariate)
  ht_mck2_dta <- data.frame("HY" = HY, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore" = Hmu0_pscore, "h.mu1.1_pscore" = Hmu1_pscore)
  ht_mck3_dta <- data.frame("HY" = HY, group_indicators, "H.pscore" = Hpscore, "H.mu0.pscore+H.mu1.1_pscore" = Hmu0_pscore_mu1_pscore)
  
  aipw_dta <- data.frame("aipw" = scores, group_indicators)
  
  ## 4.) Fit linear models.
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
  
  ## 5.) Nonparametric estimator.
  imai_li <- evalITR::GATE(D, cates, Y, n_groups)
  imai_li_results <- data.frame("group" = 1:n_groups, "GATE" = imai_li$gate, "SE" = imai_li$sd)
  
  ## 6.) Hypotheses testing.
  parametric_models <- list(wr_none_model, wr_cddf1_model, wr_cddf2_model, wr_mck1_model, 
                 ht_none_model, ht_cddf1_model, ht_cddf2_model, ht_mck1_model, ht_mck2_model, ht_mck3_model,
                 aipw_model)
  
  new_parametric_models <- list()
  counter <- 1
  
  for (model in parametric_models) {
    temp_model <- model
    tests <- test_parametric_gates(temp_model)
    
    temp_model$p.value.gates.all.equal <- tests$all_equal
    temp_model$p.value.gates.largest.difference <- tests$largest_difference
    temp_model$p.value.gates.pairwise.differences.median <- median(tests$pairwise_differences, na.rm = TRUE)
    
    new_parametric_models[[counter]] <- temp_model 
    counter <- counter + 1
  }
  
  imai_li_results$p.value.gates.all.equal <- as.numeric(evalITR::het.test(D, cates, Y, n_groups)$pval)
  imai_li_results$p.value.gates.largest.difference <- NA
  imai_li_results$p.value.gates.pairwise.differences.median <- NA
  
  ## 7.) Output.
  out <- append(new_parametric_models, list(imai_li_results))
  names(out) <- c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1",
                  "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3",
                  "aipw",
                  "imai_li")
  return(out)
}
