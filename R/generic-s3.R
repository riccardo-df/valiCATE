#' Summary Method for evaluCATE Objects
#'
#' Summarizes an \code{evaluCATE} object.
#'
#' @param object An \code{evaluCATE} object.
#' @param target String controlling which plot to display. Must be either \code{"BLP"} or \code{"GATES"}.
#' @param latex Logical, whether to print LATEX code for a table. Different tables are produced according to the \code{target} argument.
#' @param which_models Character vector with the names of the models to report. Admitted values are \code{"wr_none"}, \code{"wr_cddf1"}, \code{"wr_cddf2"}, \code{"wr_mck1"}, \code{"ht_none"}, \code{"ht_cddf1"}, \code{"ht_cddf2"}, \code{"ht_mck1"}, \code{"ht_mck2"}, \code{"ht_mck3"}, \code{"aipw"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates_val <- predict(forest, X_val)$predictions # We predict on the validation sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore_val <- rep(0.5, length(Y_val))
#' evaluation <- evaluCATE(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, pscore_val = pscore_val)
#' 
#' ## Summary.
#' summary(evaluation, target = "BLP")
#' summary(evaluation, target = "BLP", latex = TRUE)
#' 
#' summary(evaluation, target = "GATES")
#' summary(evaluation, target = "GATES", latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#' 
#' @seealso \code{\link{evaluCATE}}
#' 
#' @importFrom stats coef
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.evaluCATE <- function(object, target = "BLP", 
                              latex = FALSE, which_models = c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1", "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3", "aipw"), ...) {
  if (!(target %in% c("BLP", "GATES"))) stop("Invalid 'target'. This must be either 'GATES' or 'TOC'", call. = FALSE)
  if (!is.logical(latex)) stop("Invalid 'latex'. This must be either TRUE or FALSE.", call. = FALSE)
  if (any(!(which_models %in% c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1", "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3", "aipw")))) stop("Invalid 'which_models'. Check the documentation for admitted values.", call. = FALSE)
  if (target == "BLP" & sum(!(which_models %in% names(object$BLP))) > 0) stop("Invalid 'which_models'. This must be a character vector with one or more names from 'names(object$target)'.", call. = FALSE)
  if (target == "GATES" & sum(!(which_models %in% names(object$GATES))) > 0) stop("Invalid 'which_models'. This must be a character vector with one or more names from 'names(object$target)'.", call. = FALSE)
  
  n_models_blp <- length(object$BLP)
  n_models_gates <- length(object$GATES)
  n_groups <- length(object$GATES$wr_none$coefficients)
  
  if (!latex) {
    if (target == "BLP") {
      max_chars <- max(sapply(names(object$BLP), nchar))
      
      cat("BLP and RATEs results \n\n")
      
      cat("Estimated ATE + 95% confidence intervals: \n")
      
      for (i in seq_len(n_models_blp)) {
        model <- object$BLP[[i]]
        model_name <- names(object$BLP)[i]
        
        if (!(model_name %in% which_models)) next
        
        model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
        estimated_ate <- format(round(stats::coef(model)[names(stats::coef(model)) == "beta1"], 2), nsmall = 2)
        lower_ci <- format(round(model$conf.low["beta1"], 3), nsmall = 3)
        upper_ci <- format(round(model$conf.high["beta1"], 3), nsmall = 3)
        cat(model_name_space, if (estimated_ate < 0) ": " else ":  ", estimated_ate, if (lower_ci < 0) " [" else " [ ", lower_ci, if (upper_ci < 0) ", " else ",  ", upper_ci, "] \n", sep = "")
      }
      
      cat("\n")
      
      cat("Estimated HET + 95% confidence intervals: \n")
      
      for (i in seq_len(n_models_blp)) {
        model <- object$BLP[[i]]
        model_name <- names(object$BLP)[i]
        
        if (!(model_name %in% which_models)) next
        
        model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
        estimated_het <- format(round(stats::coef(model)[names(stats::coef(model)) == "beta2"], 2), nsmall = 2)
        lower_ci <- format(round(model$conf.low["beta2"], 3), nsmall = 3)
        upper_ci <- format(round(model$conf.high["beta2"], 3), nsmall = 3)
        cat(model_name_space, if (estimated_het < 0) ": " else ":  ", estimated_het, if (lower_ci < 0) " [" else " [ ", lower_ci, if (upper_ci < 0) ", " else ",  ", upper_ci, "] \n", sep = "")
      }
      
      cat("\n")
      
      cat("RATEs results + 95% confidence intervals: \n")
      
      autoc <- format(round(object$RATE$rate_results$AUTOC$rate, 2), nsmall = 2)
      autoc_lower_ci <- format(round(object$RATE$rate_results$AUTOC$rate - 1.96 * object$RATE$rate_results$AUTOC$se, 3), nsmall = 3)
      autoc_upper_ci <- format(round(object$RATE$rate_results$AUTOC$rate + 1.96 * object$RATE$rate_results$AUTOC$se, 3), nsmall = 3)
      
      qini <- format(round(object$RATE$rate_results$QINI$rate, 2), nsmall = 2)
      qini_lower_ci <- format(round(object$RATE$rate_results$QINI$rate - 1.96 * object$RATE$rate_results$QINI$se, 3), nsmall = 3)
      qini_upper_ci <- format(round(object$RATE$rate_results$QINI$rate + 1.96 * object$RATE$rate_results$QINI$se, 3), nsmall = 3)
      
      
      cat(if (autoc < 0) "AUTOC: " else "AUTOC:  ", autoc, if (autoc_lower_ci < 0) " [" else "[ ", autoc_lower_ci, if (autoc_upper_ci < 0) ", " else ",  ", autoc_upper_ci, "] \n",
          if (qini < 0) "QINI:  " else "QINI:   ", qini, if (qini_lower_ci < 0) " [" else "[ ", qini_lower_ci, if (qini_upper_ci < 0) ", " else ",  ", qini_upper_ci, "] \n", sep = "")
    } else if (target == "GATES") {
      max_chars <- max(sapply(names(object$GATES), nchar))
      
      cat("GATES results \n\n")
      
      cat("Estimated GATES + 95% confidence intervals: \n")
      
      for (i in seq_len(n_models_gates)) {
        model <- object$GATES[[i]]
        model_name <- names(object$GATES)[i]
        
        if (!(model_name %in% which_models)) next
        
        model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
        
        if (model_name != "imai_li") {
          estimated_gates <- format(round(stats::coef(model)[grepl("group", names(stats::coef(model)))], 2), nsmall = 2)
          lower_ci <- format(round(model$conf.low[grepl("group", names(stats::coef(model)))], 3), nsmall = 3)
          upper_ci <- format(round(model$conf.high[grepl("group", names(stats::coef(model)))], 3), nsmall = 3)
          
          cat(model_name_space, ": \n", sep = "") 
          
          for (k in seq_len(n_groups)) {
            cat(if (estimated_gates[k] < 0) "    Group " else "     Group ", k, " ", estimated_gates[k], " [", lower_ci[k], ", ", upper_ci[k], "] \n", sep = "")
          }
          
        } else if (model_name == "imai_li") {
          estimated_gates <- format(round(model$GATE, 2), nsmall = 2)
          estimated_ses <- format(round(model$SE, 2), nsmall = 3)
          
          lower_ci <- format(round(as.numeric(estimated_gates) - 1.96 * as.numeric(estimated_ses), 3), nsmall = 3)
          upper_ci <- format(round(as.numeric(estimated_gates) + 1.96 * as.numeric(estimated_ses), 3), nsmall = 3)
          
          cat(if (estimated_gates[k] < 0) "    Group " else "     Group ", k, " ", estimated_gates[k], " [", lower_ci[k], ", ", upper_ci[k], "] \n", sep = "")          
        }
      }
      
      cat("\n")
      
      cat("Hypotheses testing results (p-values): \n")
      
      for (i in seq_len(n_models_gates)) {
        model <- object$GATES[[i]]
        model_name <- names(object$GATES)[i]
        
        if (!(model_name %in% which_models)) next
        
        model_name_space <- paste0(model_name, strrep(" ", max_chars - nchar(model_name)))
        
        if (model_name != "imai_li") {
          gates_all_equal <- format(round(model$p.value.gates.all.equal, 3), nsmall = 3)
          gates_largest_diff <- format(round(model$p.value.gates.largest.difference, 3), nsmall = 3)

          cat(model_name_space, ": 
    GATES_1 = GATES_2 = ... = GATES_K : ",  gates_all_equal, "
    GATES_K = GATES_1                 : ",  gates_largest_diff, "\n", sep = "") 
        } else if (model_name == "imai_li") {
          gates_all_equal <- format(round(model$p.value.gates.all.equal, 3), nsmall = 3)[1]
          gates_largest_diff <- format(round(model$p.value.gates.largest.difference, 3), nsmall = 3)[1]
          
          cat(model_name_space, ": 
    GATES_1 = GATES_2 = ... = GATES_K : ",  gates_all_equal, "
    GATES_K = GATES_1                 : ",  gates_largest_diff, "\n", sep = "")          
        }
      }
    }
  } else if (latex) {
    if (target == "BLP") {
      estimated_ates <- sapply(object$BLP, function(x) { format(round(stats::coef(x)[names(stats::coef(x)) == "beta1"], 2), nsmall = 2) })
      estimated_ates_lower_ci <- sapply(object$BLP, function(x) { format(round(x$conf.low["beta1"], 3), nsmall = 3) })
      estimated_ates_upper_ci <- sapply(object$BLP, function(x) { format(round(x$conf.high["beta1"], 3), nsmall = 3) })
      
      estimated_hets <- sapply(object$BLP, function(x) { format(round(stats::coef(x)[names(stats::coef(x)) == "beta2"], 2), nsmall = 2) })
      estimated_hets_lower_ci <- sapply(object$BLP, function(x) { format(round(x$conf.low["beta2"], 3), nsmall = 3) })
      estimated_hets_upper_ci <- sapply(object$BLP, function(x) { format(round(x$conf.high["beta2"], 3), nsmall = 3) })
      
      names(estimated_ates) <- names(estimated_ates_lower_ci) <- names (estimated_ates_upper_ci) <- names(estimated_hets) <- names(estimated_hets_lower_ci) <- names(estimated_hets_upper_ci) <- unlist(strsplit(names(estimated_ates), ".beta1"))
      
      estimated_ates <- estimated_ates[which_models]
      estimated_ates_lower_ci <- estimated_ates_lower_ci[which_models]
      estimated_ates_upper_ci <- estimated_ates_upper_ci[which_models]

      estimated_hets <- estimated_hets[which_models]
      estimated_hets_lower_ci <- estimated_hets_lower_ci[which_models]
      estimated_hets_upper_ci <- estimated_hets_upper_ci[which_models]   
      
      cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
      \\centering
      \\begin{adjustbox}{width = 1\\textwidth}
      \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", n_models_blp), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      & ATE ($\\beta_1$) & HET ($\\beta_2$)  \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")
      
      for (i in seq_len(length(estimated_ates))) {
        cat("      ", rename_latex(names(estimated_ates)[i]), " & ", estimated_ates[i], " & ", estimated_hets[i], " \\\\ \n", sep = "")
        cat("                &  ", paste0("[", estimated_ates_lower_ci[i], ", ", estimated_ates_upper_ci[i],  "] & "), paste0("[", estimated_hets_lower_ci[i], ", ", estimated_hets_upper_ci[i],  "]"), " \\\\ \n", sep = "")
        
      } ## ADD RATES
      
      cat("\n      \\addlinespace[3pt]
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      \\end{tabular}
      \\end{adjustbox}
      \\caption{BLP results. $95\\%$ confidence intervals are displayed in brackets under each point estimate.}
      \\label{table_blp_results}
    \\end{table}
\\endgroup")
    } else if (target == "GATES") {
      gates_list <- lapply(object$GATES[1:length(object$GATES)-1], function(y) { format(round(y$coefficients[grep("group", names(y$coefficients))], 2), nsmall = 2) }) %>%
        append(list(c(format(round(object$GATES$imai_li$GATE, 2), nsmall = 2))))
      names(gates_list)[length(gates_list)] <- "imai_li"
      
      gates_se_list <- lapply(object$GATES[1:length(object$GATES)-1], function(y) { format(round(summary(y)$coefficients[grep("group", names(y$coefficients)), 2], 3), nsmall = 3) }) %>%
        append(list(c(format(round(object$GATES$imai_li$SE, 3), nsmall = 3))))
      names(gates_se_list)[length(gates_se_list)] <- "imai_li"
      
      gates_cil_list <- mapply(function(point, se) { format(round(as.numeric(point) - 1.96 * as.numeric(se), 3), nsmall = 3) }, point = gates_list, se = gates_se_list, SIMPLIFY = FALSE)
      gates_ciu_list <- mapply(function(point, se) { format(round(as.numeric(point) + 1.96 * as.numeric(se), 3), nsmall = 3) }, point = gates_list, se = gates_se_list, SIMPLIFY = FALSE)
      
      if (sum(which_models != "all") > 0) {
        gates_list <- gates_list[which_models]
        gates_se_list <- gates_list[which_models]
        gates_cil_list <- gates_list[which_models]
        gates_ciu_list <- gates_list[which_models]
      }
      
      cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
      \\centering
      \\begin{adjustbox}{width = 1\\textwidth}
      \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", n_groups), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      & ", paste0("\\textit{Group ", seq_len(n_groups - 1), "} & "), paste0("\\textit{Group ", n_groups, "}"), " \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")
      
      for (i in seq_len(length(gates_list))) {
        cat("      ", rename_latex(names(gates_list)[i]), " & ", paste0(gates_list[[i]][1:n_groups-1], " & "), paste0(gates_list[[i]][n_groups], "\\\\ \n", sep = ""))
        cat("                &  ", paste0("[", gates_cil_list[[i]][1:n_groups-1], ", ", gates_ciu_list[[i]][1:n_groups-1],  "] & "), paste0("[", gates_cil_list[[i]][n_groups], ", ", gates_ciu_list[[i]][n_groups],  "]"), " \\\\ \n", sep = "")
      }
      
      cat("\n      \\addlinespace[3pt]
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      \\end{tabular}
      \\end{adjustbox}
      \\caption{GATES results. $95\\%$ confidence intervals are displayed in brackets under each point estimate.}
      \\label{table_blp_results}
    \\end{table}
\\endgroup")
    }
  } 
}


#' Print Method for evaluCATE Objects
#'
#' Prints an \code{evaluCATE} object.
#'
#' @param x An \code{evaluCATE} object.
#' @param target String controlling which plot to display. Must be either \code{"BLP"} or \code{"GATES"}.
#' @param latex Logical, whether to print LATEX code for a table. Different tables are produced according to the \code{target} argument.
#' @param which_models Character vector with the names of the models to report. Admitted values are \code{"wr_none"}, \code{"wr_cddf1"}, \code{"wr_cddf2"}, \code{"wr_mck1"}, \code{"ht_none"}, \code{"ht_cddf1"}, \code{"ht_cddf2"}, \code{"ht_mck1"}, \code{"ht_mck2"}, \code{"ht_mck3"}, \code{"aipw"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prints an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates_val <- predict(forest, X_val)$predictions # We predict on the validation sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore_val <- rep(0.5, length(Y_val))
#' evaluation <- evaluCATE(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, pscore_val = pscore_val)
#' 
#' ## Print.
#' print(evaluation, target = "BLP")
#' print(evaluation, target = "BLP", latex = TRUE)
#' 
#' print(evaluation, target = "GATES")
#' print(evaluation, target = "GATES", latex = TRUE)}
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.evaluCATE <- function(x, target = "BLP", 
                            latex = FALSE, which_models = c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1", "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3", "aipw"), ...) {
  summary.evaluCATE(x, target, latex, which_models, ...)
}


#' Plot Method for evaluCATE Objects
#'
#' Plots an \code{evaluCATE} object.
#'
#' @param x An \code{evaluCATE} object.
#' @param target String controlling which plot to display. Must be either \code{"GATES"} or \code{"TOC"}.
#' @param which_models Character vector with the names of the models to report. Admitted values are \code{"wr_none"}, \code{"wr_cddf1"}, \code{"wr_cddf2"}, \code{"wr_mck1"}, \code{"ht_none"}, \code{"ht_cddf1"}, \code{"ht_cddf2"}, \code{"ht_mck1"}, \code{"ht_mck2"}, \code{"ht_mck3"}, \code{"aipw"}. Ignored if \code{target == "TOC"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Plots an \code{evaluCATE} object.
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
#' ## CATEs estimation.
#' library(grf)
#' 
#' forest <- causal_forest(X_tr, Y_tr, D_tr) # We use only the training sample.
#' cates_val <- predict(forest, X_val)$predictions # We predict on the validation sample.
#' 
#' ## CATEs evaluation. Estimate all nuisances internally. 
#' pscore_val <- rep(0.5, length(Y_val))
#' evaluation <- evaluCATE(Y_tr, Y_val, D_tr, D_val, X_tr, X_val, cates_val, pscore_val = pscore_val)
#' 
#' ## Plot.
#' plot(evaluation, target = "GATES")
#' plot(evaluation, target = "TOC")}
#'
#' @import dplyr ggplot2 ggsci
#' @importFrom stats coef
#'
#' @author Riccardo Di Francesco
#'
#' @export
plot.evaluCATE <- function(x, target = "GATES", which_models = c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1", "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3", "aipw"), ...) {
  ## Checks.
  if (!(target %in% c("GATES", "TOC"))) stop("Invalid 'target'. This must be either 'GATES' or 'TOC'", call. = FALSE)
  if (any(!(which_models %in% c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1", "ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3", "aipw")))) stop("Invalid 'which_models'. Check the documentation for admitted values.", call. = FALSE)
  
  group <- NULL
  estimator <- NULL
  estimated_gate <- NULL
  se <- NULL
  strategy <- NULL
  denoise <- NULL
  
  u <- NULL
  TOC <- NULL
  
  ## Plot.
  if (target == "GATES") {
    ## Extract estimated GATES and standard errors.
    gates_list <- lapply(x$GATES[1:length(x$GATES)-1], function(y) {y$coefficients[grep("group", names(y$coefficients))]}) %>%
      append(list(c(x$GATES$imai_li$GATE)))
    names(gates_list)[length(gates_list)] <- "imai_li"
    
    gates_se_list <- lapply(x$GATES[1:length(x$GATES)-1], function(y) {summary(y)$coefficients[grep("group", names(y$coefficients)), 2]}) %>%
      append(list(c(x$GATES$imai_li$SE)))
    names(gates_se_list)[length(gates_se_list)] <- "imai_li"
    
    ## Arrange plot data and call ggplot.
    n_groups <- length(stats::coef(x$GATES$wr_none))
    plot_dta <- data.frame("group" = rep(1:n_groups, length(gates_list)), 
                           "estimator" = rep(names(gates_list), each = n_groups),
                           "estimated_gate" = unlist(gates_list), 
                           "se" = unlist(gates_se_list)) %>%
      dplyr::filter(estimator %in% which_models)
    
    plot_dta$strategy <- NA
    plot_dta$strategy[plot_dta$estimator %in% c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1")] <- "WR"
    plot_dta$strategy[plot_dta$estimator %in% c("ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3")] <- "HT"
    plot_dta$strategy[plot_dta$estimator == "aipw"] <- "AIPW"
    
    plot_dta$denoise <- NA
    plot_dta$denoise[plot_dta$estimator %in% c("wr_none", "ht_none", "aipw")] <- "None"
    plot_dta$denoise[plot_dta$estimator %in% c("wr_cddf1", "ht_cddf1")] <- "cddf1"
    plot_dta$denoise[plot_dta$estimator %in% c("wr_cddf2", "ht_cddf2")] <- "cddf2"
    plot_dta$denoise[plot_dta$estimator %in% c("wr_mck1", "ht_mck1")] <- "mck1"
    plot_dta$denoise[plot_dta$estimator %in% c("ht_mck2")] <- "mck2"
    plot_dta$denoise[plot_dta$estimator %in% c("ht_mck3")] <- "mck3"
    
    facet_condition <- length(unique(factor(plot_dta$denoise, levels = c("None", "cddf1", "cddf2", "mck1", "mck2", "mck3")))) > 1
    
    ggplot2::ggplot(plot_dta, ggplot2::aes(x = group, y = estimated_gate, colour = factor(estimator, levels = which_models))) +
      ggplot2::geom_point() + 
      ggplot2::geom_errorbar(aes(x = group, ymin = estimated_gate - 1.96 * se, ymax = estimated_gate + 1.96 * se)) +
      ggplot2::facet_grid(cols = vars(factor(strategy, levels = c("WR", "HT", "AIPW"))), rows = if (facet_condition) vars(factor(denoise, levels = c("None", "cddf1", "cddf2", "mck1", "mck2", "mck3"))) else NULL, scales = "fixed") +
      ggplot2::xlab("Group") + ggplot2::ylab("GATES") + ggplot2::labs(colour = "Estimator") +
      ggplot2::theme_bw() + 
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"), legend.position = "none")
  } else if (target == "TOC") {
    plot_dta <- data.frame("unit" = (1:length(x$RATE$toc_results) / length(x$RATE$toc_results)), "TOC" = x$RATE$toc_results)
    ggplot2::ggplot(plot_dta, ggplot2::aes(x = unit, y = TOC)) +
      ggplot2::geom_line(color = "tomato", linewidth = 1.5) + 
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + 
      ggplot2::xlab("Fraction treated") +ggplot2:: ylab("TOC") +
      ggplot2::theme_bw() 
  }
}
