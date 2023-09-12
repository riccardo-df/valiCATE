#' Summary Method for evaluCATE Objects
#'
#' Summarizes an \code{evaluCATE} object.
#'
#' @param object An \code{evaluCATE} object.
#' @param latex String to denote whether and which table to print in LATEX code. Must be one of \code{"none"}, \code{"BLP"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates <- predict(forest, X)$predictions # We predict on the whole sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore <- rep(0.5, length(Y))
#' evaluation <- evalue_cates(Y, D, X, cates, train_idx, pscore = pscore)
#' 
#' ## Summary.
#' summary(evaluation)
#' summary(evaluation, latex = "BLP")
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#' 
#' @seealso \code{\link{evalue_cates}}
#' 
#' @import stats
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.evaluCATE <- function(object, latex = "none", ...) {
  if (!(latex %in% c("none", "BLP"))) stop("Invalid 'latex'. This must be one of 'none', 'BLP'.", call. = FALSE)
  
  n_models <- length(object$BLP)
  
  if (latex == "none") {
    max_chars <- max(sapply(names(object$BLP), nchar))
    
    cat("CATEs evaluation results \n\n")
    
    cat("Estimated ATE + 95% confidence intervals: \n")
    
    for (i in seq_len(n_models)) {
      model <- object$BLP[[i]]
      model_name <- names(object$BLP)[i]
      model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
      estimated_ate <- round(stats::coef(model)[names(stats::coef(model)) == "beta1"], 2)
      lower_ci <- round(model$conf.low["beta1"], 3)
      upper_ci <- round(model$conf.high["beta1"], 3)
      cat(model_name_space, ": ", estimated_ate, " [", lower_ci, ", ", upper_ci, "] \n", sep = "")
    }
    
    cat("\n")
    
    cat("Estimated HET + 95% confidence intervals: \n")
    
    for (i in seq_len(n_models)) {
      model <- object$BLP[[i]]
      model_name <- names(object$BLP)[i]
      model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
      estimated_ate <- round(stats::coef(model)[names(stats::coef(model)) == "beta2"], 2)
      lower_ci <- round(model$conf.low["beta2"], 3)
      upper_ci <- round(model$conf.high["beta2"], 3)
      cat(model_name_space, ": ", estimated_ate, " [", lower_ci, ", ", upper_ci, "] \n", sep = "")
    }
    
    cat("\n")
    
    cat("RATEs results + 95% confidence intervals: \n")
    
    autoc <- round(object$RATE$rate_results$AUTOC$rate, 2)
    autoc_lower_ci <- round(object$RATE$rate_results$AUTOC$rate - 1.96 * object$RATE$rate_results$AUTOC$se, 3)
    autoc_upper_ci <- round(object$RATE$rate_results$AUTOC$rate + 1.96 * object$RATE$rate_results$AUTOC$se, 3)
    
    qini <- round(object$RATE$rate_results$QINI$rate, 2)
    qini_lower_ci <- round(object$RATE$rate_results$QINI$rate - 1.96 * object$RATE$rate_results$QINI$se, 3)
    qini_upper_ci <- round(object$RATE$rate_results$QINI$rate + 1.96 * object$RATE$rate_results$QINI$se, 3)
    
    
    cat("AUTOC: ", autoc, " [", autoc_lower_ci, ", ", autoc_upper_ci, "]
QINI:  ", qini, " [", qini_lower_ci, ", ", qini_upper_ci, "] \n", sep = "")
    
  } else if (latex == "BLP") {
    estimated_ates <- sapply(object$BLP, function(x) { round(stats::coef(x)[names(stats::coef(x)) == "beta1"], 2) })
    estimated_ates_lower_ci <- sapply(object$BLP, function(x) { round(x$conf.low["beta1"], 3) })
    estimated_ates_upper_ci <- sapply(object$BLP, function(x) { round(x$conf.high["beta1"], 3) })
    
    estimated_hets <- sapply(object$BLP, function(x) { round(stats::coef(x)[names(stats::coef(x)) == "beta2"], 2) })
    estimated_hets_lower_ci <- sapply(object$BLP, function(x) { round(x$conf.low["beta2"], 3) })
    estimated_hets_upper_ci <- sapply(object$BLP, function(x) { round(x$conf.high["beta2"], 3) })
    
    cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
        \\centering
        \\begin{adjustbox}{width = 1\\textwidth}
        \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", n_models), "}
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        & ", paste0(rename_latex(names(object$BLP[-n_models])), " & "), rename_latex(names(object$BLP[n_models])), " \\\\
        \\addlinespace[2pt]
        \\hline \\\\[-1.8ex] \n\n", sep = "")
    
    cat("        ATE ($\\beta_1$) & ", paste0(estimated_ates[-n_models], " & "), estimated_ates[n_models], " \\\\ \n", sep = "")
    cat("                         &  ", paste0("[", estimated_ates_lower_ci[-n_models], ", ", estimated_ates_upper_ci[-n_models],  "] & "), paste0("[", estimated_ates_lower_ci[n_models], ", ", estimated_ates_upper_ci[n_models],  "]"), " \\\\ \n", sep = "")
    
    cat("        HET ($\\beta_2$) & ", paste0(estimated_hets[-n_models], " & "), estimated_hets[n_models], " \\\\ \n", sep = "")
    cat("                         &  ", paste0("[", estimated_hets_lower_ci[-n_models], ", ", estimated_hets_upper_ci[-n_models],  "] & "), paste0("[", estimated_hets_lower_ci[n_models], ", ", estimated_hets_upper_ci[n_models],  "]"), " \\\\ \n", sep = "")
    
    cat("\n        \\addlinespace[3pt]
        \\\\[-1.8ex]\\hline
        \\hline \\\\[-1.8ex]
        \\end{tabular}
        \\end{adjustbox}
        \\caption{BLP results. $95\\%$ confidence intervals are displayed in brackets under each point estimate.}
        \\label{table_blp_results}
    \\end{table}
\\endgroup")
  } 
}


#' Print Method for evaluCATE Objects
#'
#' Prints an \code{evaluCATE} object.
#'
#' @param x An \code{evaluCATE} object.
#' @param latex String to denote whether and which table to print in LATEX code. Must be one of \code{"none"}, \code{"BLP"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prints an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates <- predict(forest, X)$predictions # We predict on the whole sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore <- rep(0.5, length(Y))
#' evaluation <- evalue_cates(Y, D, X, cates, train_idx, pscore = pscore)
#' 
#' ## Print.
#' print(evaluation)
#' print(evaluation, latex = "BLP")
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.evaluCATE <- function(x, latex = "none", ...) {
  summary.evaluCATE(x, latex, ...)
}


#' Plot Method for evaluCATE Objects
#'
#' Plots an \code{evaluCATE} object.
#'
#' @param x An \code{evaluCATE} object.
#' @param target String controlling which plot to display. Must be either \code{"GATES"} or \code{"RATE"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Plots an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates <- predict(forest, X)$predictions # We predict on the whole sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore <- rep(0.5, length(Y))
#' evaluation <- evalue_cates(Y, D, X, cates, train_idx, pscore = pscore)
#' 
#' ## Plot.
#' plot(evaluation, target = "GATES")
#' plot(evaluation, target = "RATE")
#'
#' @details
#' If \code{target == "GATES"}, the estimated GATES and 95% confidence intervals are displayed.\cr
#' 
#' If \code{target == "RATE"}, the estimated TOCs for the considered threshold values are displayed.\cr
#'
#' @import dplyr ggplot2 ggsci stats
#'
#' @author Riccardo Di Francesco
#'
#' @export
plot.evaluCATE <- function(x, target = "GATES", ...) {
  ## Checks.
  if (!(target %in% c("GATES", "RATE"))) stop("Invalid 'target'. This must be either 'GATES' or 'RATE'", call. = FALSE)
  
  group <- NULL
  estimator <- NULL
  estimated_gate <- NULL
  se <- NULL
  
  u <- NULL
  TOC <- NULL
  
  ## Plot.
  if (target == "GATES") {
    ## Extract estimated GATES and standard errors.
    gates_list <- lapply(x$GATES[1:length(x$GATES)-1], function(y) {y$coefficients[grep("gamma", names(y$coefficients))]}) %>%
      append(list(c(x$GATES$imai_li$GATE)))
    names(gates_list)[length(gates_list)] <- "imai_li"
    
    gates_se_list <- lapply(x$GATES[1:length(x$GATES)-1], function(y) {summary(y)$coefficients[grep("gamma", names(y$coefficients)), 2]}) %>%
      append(list(c(x$GATES$imai_li$SE)))
    names(gates_list)[length(gates_list)] <- "imai_li"
    
    ## Arrange plot data and call ggplot.
    n_groups <- length(stats::coef(x$GATES$wr_none))
    plot_dta <- data.frame("group" = rep(1:n_groups, length(gates_list)), 
                           "estimator" = factor(rep(names(gates_list), each = n_groups)),
                           "estimated_gate" = unlist(gates_list), 
                           "se" = unlist(gates_se_list))
    
    ggplot2::ggplot(plot_dta, ggplot2::aes(x = group, y = estimated_gate, colour = estimator)) +
      ggplot2::geom_point() + 
      ggplot2::geom_errorbar(aes(x = group, ymin = estimated_gate - 1.96 * se, ymax = estimated_gate + 1.96 * se)) +
      # ggsci::scale_fill_rickandmorty() + 
      ggplot2::xlab("Group") +ggplot2:: ylab("GATES") + ggplot2::labs(colour = "Estimator") +
      ggplot2::theme_bw() 
  } else if (target == "RATE") {
    ggplot2::ggplot(x$RATE$toc_results, ggplot2::aes(x = u, y = TOC)) +
      ggplot2::geom_line(color = "tomato", linewidth = 1.5) + 
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + 
      ggplot2::xlab("Fraction treated") +ggplot2:: ylab("TOC") +
      ggplot2::theme_bw() 
  }
}
