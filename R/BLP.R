#' BLP Estimation Strategies
#'
#' Estimates the best linear predictor of the actual CATEs using the estimated CATEs.
#'
#' @param y Observed outcomes.
#' @param D Treatment indicator.
#' @param cates Estimated CATEs. CATEs must be estimated with different observations than those in \code{y} and \code{D}.
#' @param strategy String indicating the desired identification strategy. One of "wr", "ht", "aipw".
#' @param denoise String denoting whether to include optional constructed covariates to the regression. One of "none", "cddf1", "cddf2", "mck1", "mck2", "mck3". Useful to reduce the variance of the estimation.
#' @param pscore Estimated propensity scores. They must be estimated with different observations than those in \code{y} and \code{D}. Necessary only under certain combinations of \code{strategy} and \code{denoise}.
#' @param mu Estimated regression function. It must be estimated with different observations than those in \code{y} and \code{D}. Necessary only under certain combinations of \code{strategy} and \code{denoise}.
#' @param mu0 Estimated regression function for control units. It must be estimated with different observations than those in \code{y} and \code{D}. Necessary only under certain combinations of \code{strategy} and \code{denoise}.
#' @param mu1 Estimated regression function for treated. It must be estimated with different observations than those in \code{y} and \code{D}. Necessary only under certain combinations of \code{strategy} and \code{denoise}.
#' @param scores Estimated doubly-robust scores. They must be estimated with different observations than those in \code{y} and \code{D}. Necessary only if \code{strategy} is set to \code{"aipw"}.
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
#' Currently, the routine supports only single-splits.
#'
#' @import estimatr aggTrees stats
#'
#' @author Riccardo Di Francesco
#'
#' @seealso Other functions
#'
#' @export
blp_regression <- function(y, D, cates, strategy = "wr", denoise = "none", pscore = NULL, mu  = NULL, mu0 = NULL, mu1 = NULL, scores = NULL) {
  ## 0.) Handling inputs and checks.
  if (!(strategy %in% c("wr", "ht", "aipw"))) stop("Invalid 'strategy'. This must be one of 'wr', 'ht', 'aipw'.", call. = FALSE)
  if (!(denoise %in% c("none", "cddf1", "cddf2", "mck1", "mck2", "mck3"))) stop("Invalid 'strategy'. This must be one of 'none', 'cddf1', 'cddf2', 'mck1', 'mck2', 'mck3'.", call. = FALSE)
  if (!is.null(denoise)) {
    if (is.null(mu0)) stop("When 'denoise' is not 'none', you must provide 'mu0'.", call. = FALSE)
    if (denoise != "cddf1" & is.null(pscore)) stop("When 'denoise' is not 'none' or 'cddf1', you must provide 'pscore'.", call. = FALSE)
    if (strategy == "wr" & denoise == "mck1" & is.null(mu0)) stop("When 'strategy' is 'wr' and 'denoise' is 'mck1', you must provide 'mu0'.", call. = FALSE)
    if (strategy == "ht" & denoise %in% c("mck2", "mck3") & is.null(mu1)) stop("When 'strategy' is 'ht' and 'denoise' is 'mck2' or 'mck3', you must provide 'mu1'.", call. = FALSE)
  }
  if (strategy == "aipw" & !is.null(denoise)) stop("Strategy 'aipw' supports only 'denoise' set to 'none'.", call. = FALSE)
  if (strategy == "aipw" & is.null(scores)) stop("When strategy is 'aipw', you need to supply the estimated 'scores'.", call. = FALSE)

  ## 1.) Depending on the combination (strategy, denoise), construct covariates and fit linear model via OLS.
  if (strategy == "wt") {
    wr_weights <- (pscore * (1 - pscore))^(-1)
    D_residual <- D - pscore
    demeaned_cates <- cates - mean(cates)
    interaction_D_cates <- D_residual * demeaned_cates

    if (denoise == "none") {
      model <- estimatr::lm_robust(y ~ 0 + D_residual + interaction_D_cates, weights = wr_weights, se_type = "HC1")
    } else if (denoise == "cddf1") {
      model <- estimatr::lm_robust(y ~ 0 + D_residual + interaction_D_cates + mu0 , weights = wr_weights, se_type = "HC1")
    } else if (denoise == "cddf2") {
      interaction_pscore_cates <- pscore * cates

      if (length(unique(pscore)) == 1) { # If pscore is constant, drop it.
        model <- estimatr::lm_robust(y ~ D_residual + interaction_D_cates + mu0 + interaction_pscore_cates, weights = wr_weights, se_type = "HC1")
      } else {
        model <- estimatr::lm_robust(y ~ D_residual + interaction_D_cates + mu0 + pscore + interaction_pscore_cates, weights = wr_weights, se_type = "HC1")
      }
    } else if (denoise == "mck1") {
      model <- estimatr::lm_robust(y_test ~ 0 + D_residual + interaction_D_cates + mu, weights = wr_weights, se_type = "HC1")
    }
  } else if (strategy == "ht") {
    wr_weights <- (pscore * (1 - pscore))^(-1)
    D_residual <- D - pscore
    H <- D_residual * wr_weights
    HY <- H * y
    demeaned_cates <- cates - mean(cates)

    if (denoise == "none") {
      model <- estimatr::lm_robust(HY ~ demeaned_cates, se_type = "HC1")
    } else if (denoise == "cddf1") {
      Hmu0 <- H * mu0

      model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0, se_type = "HC1")
    } else if (denoise == "cddf2") {
      Hmu0 <- H * mu0
      Hpscore <- H * pscore
      Hinteraction_pscore_cates <- Hpscore * cates

      model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0 + Hpscore + Hinteraction_pscore_cates, se_type = "HC1")
    } else if (denoise == "mck1") {
      Hmu0 <- H * mu0
      new_mck_covariate <- H * (1 - pscore) * cates

      model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0 + new_mck_covariate, se_type = "HC1")
    } else if (denoise == "mck2") {
      Hpscore <- H * pscore
      Hmu0_pscore <- H * mu0 * pscore
      Hmu1_pscore <- H * mu1 * (1 - pscore)

      model <- estimatr::lm_robust(HY ~ demeaned_cates + Hmu0_pscore + Hmu1_pscore + Hpscore, se_type = "HC1")
    } else if (denoise == "mck3") {
      Hpscore <- H * pscore
      Hmu0_pscore <- H * mu0 * pscore
      Hmu1_pscore <- H * mu1 * (1 - pscore)
      Hmu0_pscore_mu1_pscore <- Hmu0_pscore + Hmu1_pscore

      model <- estimatr::lm_robust(HY ~ demeaned_cates + Hpscore + Hmu0_pscore_mu1_pscore, se_type = "HC1")
    }
  } else if (strategy == "aipw") {
    demeaned_cates <- cates - mean(cates)

    model <- estimatr::lm_robust(scores ~ demeaned_cates, se_type = "HC1")
  }

  ## 2.) Output.
  out <- model
  return(out)
}
