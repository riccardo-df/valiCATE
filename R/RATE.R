#' TOC Estimation
#'
#' Estimates the targeting operator characteristic (TOC) induced by the estimated CATEs.
#'
#' @param cates Estimated CATEs. Must be out-of-sample predictions on the validation sample.
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting using the validation sample. 
#' @param beneficial Logical, whether the treatment is beneficial to units. If \code{TRUE}, units are ranked according to decreasing values of \code{cates}, otherwise they are ranked according to increasing values of \code{cates}.
#'
#' @return
#' A vector of estimated TOCs sorted according to the increasing or decreasing values of \code{cates}.
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
#' cates_val <- predict(cates_forest, X_val)$predictions 
#' 
#' ## AIPW scores estimation.
#' ## Cross-fitting on the validation sample.
#' library(aggTrees)
#' scores_val <- dr_scores(Y_val, D_val, X_val)
#' 
#' ## TOC estimation. 
#' toc_results <- toc_estimation(cates_val, scores_val)
#'
#' @details
#' To estimate the TOC induced by the estimated CATEs, the user must provide the estimated CATEs and doubly-robust scores. Be careful, as the CATEs must be estimated only with 
#' using the training sample, while the doubly-robust scores must be estimated in the validation sample using K-fold cross fitting (see the example section below).\cr
#' 
#' The TOC is estimated using a sample-averaging estimator. Check the \href{https://riccardo-df.github.io/evaluCATE/articles/evalue-cates-short-tutorial.html}{online vignette} 
#' for details.
#'
#' @author Riccardo Di Francesco
#'
#' @seealso \code{\link{blp_estimation}}, \code{\link{gates_estimation}}, \code{\link{rate_estimation}}
#'
#' @export
toc_estimation <- function(cates, scores, beneficial = TRUE) {
  sorted_scores <- scores[order(cates, decreasing = beneficial)]
  return(cumsum(sorted_scores) / cumsum(rep(1, length(scores))) - mean(sorted_scores))
}


#' RATE Estimation
#'
#' Estimates the rank-weighted average treatment effects (RATEs) induced by the estimated CATEs.
#'
#' @param cates Estimated CATEs. Must be out-of-sample predictions on the validation sample.
#' @param scores Estimated doubly-robust scores. Must be estimated via K-fold cross-fitting using a validation sample. 
#' @param n_boot Number of bootstrap replications to estimate standard errors.
#' @param beneficial Logical, whether the treatment is beneficial to units. If \code{TRUE}, units are ranked according to decreasing values of \code{cates}, otherwise they are ranked according to increasing values of \code{cates}.
#' 
#' @return
#' A list storing a data frame for the TOC results, a list with two data frames for the AUTOC and QINI results, and a
#' data frame for the bootstrap results.
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
#' cates_val <- predict(cates_forest, X_val)$predictions 
#' 
#' ## AIPW scores estimation.
#' ## Cross-fitting on the validation sample.
#' library(aggTrees)
#' scores_val <- dr_scores(Y_val, D_val, X_val)
#' 
#' ## RATE estimation. 
#' rate_results <- rate_estimation(cates_val, scores_val)
#'
#' @details
#' To estimate the RATEs induced the estimated CATEs, the user must provide the estimated CATEs and doubly-robust scores. Be careful, as the CATEs must be estimated only with 
#' using the training sample, while the doubly-robust scores must be estimated in the validation sample using K-fold cross fitting (see the example section below).\cr
#' 
#' \code{\link{rate_estimation}} first calls the \code{\link{toc_estimation}} function to estimate the TOCs, and then computes two RATEs, the AUTOC and the QINI coefficient, 
#' as weighted averages of the TOCs. Check the \href{https://riccardo-df.github.io/evaluCATE/articles/evalue-cates-short-tutorial.html}{online vignette} 
#' for details.\cr
#' 
#' Standard errors are estimated using the standard deviation of the bootstrap estimates obtained using the half-sample bootstrap.
#'
#' @author Riccardo Di Francesco
#' 
#' @importFrom dplyr bind_rows
#'
#' @seealso \code{\link{blp_estimation}}, \code{\link{gates_estimation}}, \code{\link{toc_estimation}}
#'
#' @export
rate_estimation <- function(cates, scores, beneficial = TRUE, n_boot = 200) {
  ## 1.) Estimate the TOC for a grid of u.
  toc_results <- toc_estimation(cates, scores)
  
  ## 2.) Estimate the AUTOC and the QINI coefficient.
  n <- length(cates)
  
  autoc <- mean(toc_results)
  qini <- 1 / n * sum(cumsum(rep(1, n)) / sum(rep(1, n)) * toc_results)
  
  ## 3.) Half-sample bootstrap to obtain standard errors.
  boot_results <- list()

  for (b in seq_len(n_boot)) {
    boot_idx <- sample(1:n, floor(n/2), replace = FALSE) 
    
    boot_cates <- cates[boot_idx]
    boot_scores <- scores[boot_idx]
    
    boot_toc_results <- toc_estimation(boot_cates, boot_scores)
    
    boot_autoc <- mean(boot_toc_results)
    boot_qini <- 1 / floor(n/2) * sum(cumsum(rep(1, floor(n/2))) / sum(rep(1, floor(n/2))) * boot_toc_results)
    
    boot_results[[b]] <- data.frame("AUTOC" = boot_autoc, "QINI" = boot_qini)
  }
  
  boot_results <- dplyr::bind_rows(boot_results)
  
  se <- apply(boot_results, 2, sd)
  
  ## 4.) Output.
  out <- list("toc_results" = toc_results,
              "rate_results" = list("AUTOC" = data.frame("rate" = autoc, "se" = se["AUTOC"]), "QINI" = data.frame("rate" = qini, "se" = se["QINI"])),
              "boot_results" = boot_results)
  
  return(out)
}
