#' Summary Method for evaluCATE Objects
#'
#' Summarizes an \code{\link{evaluCATE}} object.
#'
#' @param object An \code{\link{evaluCATE}} object.
#' @param target String controlling which parameters we are interested in. Must be one of \code{"BLP"}, \code{"GATES"}, \code{"RATE"}.
#' @param latex Logical, whether to print LATEX code for a table. Different tables are produced according to the \code{target} argument.
#' @param which_models Character vector controlling which results to display. Admitted values are those stored in \code{names(object$GATES[[1]])}. Ignored if \code{target == "RATE"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Summarizes an \code{\link{evaluCATE}} object.
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}, \code{multirow}.
#' 
#' @seealso \code{\link{evaluCATE}}
#' 
#' @import cli
#' @importFrom stats coef
#' @importFrom stringr str_sub
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
summary.evaluCATE <- function(object, target = "BLP", 
                              latex = FALSE, which_models = names(object$GATES[[1]]), ...) {
  if (!(target %in% c("BLP", "GATES", "RATE"))) stop("Invalid 'target'. This must be one of 'BLP', 'GATES', 'RATE'", call. = FALSE)
  if (!is.logical(latex)) stop("Invalid 'latex'. This must be either TRUE or FALSE.", call. = FALSE)
  if (any(!(which_models %in% names(object$GATES[[1]])))) stop("Invalid 'which_models'. Admitted values are those stored in 'names(object$GATES)'.", call. = FALSE)

  which_models_blp <- setdiff(which_models, "imai_li")
  n_models_blp <- length(which_models_blp)
  n_models_gates <- length(which_models)
  n_groups <- length(object$GATES[[1]]$wr_none$coefficients)
  
  cate_models <- names(object$BLP)
  n_cate_models <- length(cate_models)
  
  if (!latex) {
    if (target == "BLP") {
      limits_table <- 6
      cli::cli_h1("BLP RESULTS")
      
      for (j in seq_len(length(cate_models))) {
        cli::cli_h2(cate_models[j])
        
        cat(paste0(strrep(" ", limits_table), "MODEL", paste0(strrep(" ", limits_table))), "|", paste0(strrep(" ", limits_table), "ATE", paste0(strrep(" ", limits_table))), "|", paste0(strrep(" ", limits_table), "HET", paste0(strrep(" ", limits_table))), "\n")
        cat(strrep("-", limits_table * 2 + nchar("MODEL")), "|", strrep("-", limits_table * 2 + nchar("ATE")), "|", strrep("-", limits_table * 2 + nchar("HET")), "| \n")
        
        for (i in seq_len(n_models_blp)) {
          model <- object$BLP[[j]][[i]]
          model_name <- names(object$BLP[[j]])[i]
          
          if (!(model_name %in% which_models_blp)) next
        
          estimated_ate <- format(round(stats::coef(model)[names(stats::coef(model)) == "beta1"], 2), nsmall = 2)
          ate_lower_ci <- format(round(model$conf.low["beta1"], 3), nsmall = 3)
          ate_upper_ci <- format(round(model$conf.high["beta1"], 3), nsmall = 3)
          
          estimated_het <- format(round(stats::coef(model)[names(stats::coef(model)) == "beta2"], 2), nsmall = 2)
          het_lower_ci <- format(round(model$conf.low["beta2"], 3), nsmall = 3)
          het_upper_ci <- format(round(model$conf.high["beta2"], 3), nsmall = 3)
          
          ate_ci_print <- paste0("[", ate_lower_ci, ", ", ate_upper_ci, "]")
          het_ci_print <- paste0("[", het_lower_ci, ", ", het_upper_ci, "]")
          
          fill_model <- limits_table * 2 + (nchar("MODEL") - nchar(model_name))
          fill_ate <- limits_table * 2 + (nchar("ATE") - nchar(estimated_ate))
          fill_het <- limits_table * 2 + (nchar("HET") - nchar(estimated_het))
          fill_ate_ci <- limits_table * 2 + (nchar("ATE") - nchar(ate_ci_print))
          fill_het_ci <- limits_table * 2 + (nchar("HET") - nchar(het_ci_print))
          
          if (fill_model %% 2 == 0) {
            fill_model_left <- fill_model / 2 
            fill_model_right <- fill_model / 2 
          } else {
            fill_model_left <- (fill_model + 1) / 2 
            fill_model_right <- (fill_model - 1) / 2 
          }
          
          if (fill_ate %% 2 == 0) {
            fill_ate_left <- fill_ate / 2 
            fill_ate_right <- fill_ate / 2 
          } else {
            fill_ate_left <- (fill_ate + 1) / 2 
            fill_ate_right <- (fill_ate - 1) / 2 
          }
          
          if (fill_het %% 2 == 0) {
            fill_het_left <- fill_het / 2 
            fill_het_right <- fill_het / 2 
          } else {
            fill_het_left <- (fill_het + 1) / 2 
            fill_het_right <- (fill_het - 1) / 2 
          }
          
          if (fill_ate_ci %% 2 == 0) {
            fill_ate_ci_left <- fill_ate_ci / 2 
            fill_ate_ci_right <- fill_ate_ci / 2 
          } else {
            fill_ate_ci_left <- (fill_ate_ci + 1) / 2 
            fill_ate_ci_right <- (fill_ate_ci - 1) / 2 
          }
          
          if (fill_het_ci %% 2 == 0) {
            fill_het_ci_left <- fill_het_ci / 2 
            fill_het_ci_right <- fill_het_ci / 2 
          } else {
            fill_het_ci_left <- (fill_het_ci + 1) / 2 
            fill_het_ci_right <- (fill_het_ci - 1) / 2 
          }
          
          cat(paste0(strrep(" ", fill_model_left), model_name, paste0(strrep(" ", fill_model_right))), "|", paste0(strrep(" ", fill_ate_left), estimated_ate, paste0(strrep(" ", fill_ate_right))), "|", paste0(strrep(" ", fill_het_left), estimated_het, paste0(strrep(" ", fill_het_right))), "| \n")
          cat(strrep(" ", limits_table * 2 + nchar("MODEL")), "|", paste0(strrep(" ", fill_ate_ci_left), ate_ci_print, paste0(strrep(" ", fill_ate_ci_right))), "|", paste0(strrep(" ", fill_het_ci_left), het_ci_print, paste0(strrep(" ", fill_het_ci_right))), "| \n")    
        }
      }
      
      cli::cli_h1("")
    } else if (target == "GATES") {
      limits_table <- 6
      # cli::cli_h1("GATES RESULTS")
      # 
      # for (j in seq_len(length(cate_models))) {
      #   cli::cli_h2(cate_models[j])
      # 
      #   group_string <- paste0(strrep(" ", limits_table), "GROUP", 1:n_groups, paste0(strrep(" ", limits_table)), "|")
      # 
      #   cat(paste0(strrep(" ", limits_table), "MODEL", paste0(strrep(" ", limits_table))), "|", group_string, "\n")
      #   cat(strrep("-", limits_table * 2 + nchar("MODEL ")), "|", strrep(paste0(strrep("-", limits_table * 2 + nchar("groupk ")), "|"), n_groups), "\n", sep = "")
      #   
      #   for (i in seq_len(n_models_gates)) {
      #     model <- object$GATES[[j]][[i]]
      #     model_name <- names(object$GATES[[j]])[i]
      #     
      #     if (!(model_name %in% which_models)) next
      #     
      #     if (model_name != "imai_li") {
      #       estimated_gates <- format(round(stats::coef(model)[grepl("group", names(stats::coef(model)))], 2), nsmall = 2)
      #       lower_ci <- format(round(model$conf.low[grepl("group", names(stats::coef(model)))], 3), nsmall = 3)
      #       upper_ci <- format(round(model$conf.high[grepl("group", names(stats::coef(model)))], 3), nsmall = 3)
      #     } else if (model_name == "imai_li") {
      #       estimated_gates <- format(round(model$GATE, 2), nsmall = 2)
      #       estimated_ses <- format(round(model$SE, 2), nsmall = 3)
      #       
      #       lower_ci <- format(round(as.numeric(estimated_gates) - 1.96 * as.numeric(estimated_ses), 3), nsmall = 3)
      #       upper_ci <- format(round(as.numeric(estimated_gates) + 1.96 * as.numeric(estimated_ses), 3), nsmall = 3)
      #     }
      #     
      #     cis_print <- paste0("[", lower_ci, ", ", upper_ci, "]")
      # 
      #     fill_model <- limits_table * 2 + (nchar("MODEL") - nchar(model_name))
      #     fill_gates <- (limits_table * 2 + (nchar("groupk") - nchar(estimated_gates)))[1]
      #     fill_cis <- (limits_table * 2 + (nchar("groupk") - nchar(cis_print)))[1]
      #     
      #     if (fill_model %% 2 == 0) {
      #       fill_model_left <- fill_model / 2 
      #       fill_model_right <- fill_model / 2 
      #     } else {
      #       fill_model_left <- (fill_model + 1) / 2 
      #       fill_model_right <- (fill_model - 1) / 2 
      #     }
      #     
      #     if (fill_gates %% 2 == 0) {
      #       fill_gates_left <- fill_gates / 2 
      #       fill_gates_right <- fill_gates / 2 
      #     } else {
      #       fill_gates_left <- (fill_gates + 1) / 2 
      #       fill_gates_right <- (fill_gates - 1) / 2 
      #     }
      #     
      #     if (fill_cis %% 2 == 0) {
      #       fill_cis_left <- fill_cis / 2 
      #       fill_cis_right <- fill_cis / 2 
      #     } else {
      #       fill_cis_left <- (fill_cis + 1) / 2 
      #       fill_cis_right <- (fill_cis - 1) / 2 
      #     }
      #     
      #     gates_string <- paste0(strrep(" ", fill_gates_left), estimated_gates, paste0(strrep(" ", fill_gates_right)), "|")
      #     cis_string <- paste0(strrep(" ", fill_cis_left), cis_print, paste0(strrep(" ", fill_cis_right)), "|")
      #     
      #     cat(paste0(strrep(" ", fill_model_left), model_name, paste0(strrep(" ", fill_model_right))), "|", gates_string, "\n")
      #     cat(strrep(" ", limits_table * 2 + nchar("MODEL")), "|", cis_string, " \n")    
      #   }
      # }
      # 
      # cli::cli_h1("")
      # 
      # cat("\n")
      # 
      cli::cli_h1("HYPOTHESIS TESTING RESULTS (p-values)")

      for (j in seq_len(length(cate_models))) {
        cli::cli_h2(cate_models[j])

        cat(paste0(strrep(" ", limits_table), "MODEL", paste0(strrep(" ", limits_table))), "|", paste0(strrep(" ", limits_table), "GATES_1 = GATES_2 = ... = GATES_K", paste0(strrep(" ", limits_table))), "|", paste0(strrep(" ", limits_table), "GATES_K = GATES_1", paste0(strrep(" ", limits_table))), "\n")
        cat(strrep("-", limits_table * 2 + nchar("MODEL")), "|", strrep("-", limits_table * 2 + nchar("GATES_1 = GATES_2 = ... = GATES_K")), "|", strrep("-", limits_table * 2 + nchar("GATES_K = GATES_1")), "| \n")
        
        for (i in seq_len(n_models_gates)) {
          model <- object$GATES[[j]][[i]]
          model_name <- names(object$GATES[[j]])[i]
          
          if (!(model_name %in% which_models)) next
          
          if (model_name != "imai_li") {
            gates_all_equal <- format(round(model$p.value.gates.all.equal, 3), nsmall = 3)
            gates_largest_diff <- format(round(model$p.value.gates.largest.difference, 3), nsmall = 3)
          } else if (model_name == "imai_li") {
            gates_all_equal <- format(round(model$p.value.gates.all.equal, 3), nsmall = 3)[1]
            gates_largest_diff <- format(round(model$p.value.gates.largest.difference, 3), nsmall = 3)[1]
          }
          
          fill_model <- limits_table * 2 + (nchar("MODEL") - nchar(model_name))
          fill_gates_all_equal <- limits_table * 2 + (nchar("GATES_1 = GATES_2 = ... = GATES_K") - nchar(gates_all_equal))
          fill_gates_largest_diff <- limits_table * 2 + (nchar("GATES_K = GATES_1") - nchar(gates_largest_diff))
          
          if (fill_model %% 2 == 0) {
            fill_model_left <- fill_model / 2 
            fill_model_right <- fill_model / 2 
          } else {
            fill_model_left <- (fill_model + 1) / 2 
            fill_model_right <- (fill_model - 1) / 2 
          }
          
          if (fill_gates_all_equal %% 2 == 0) {
            fill_gates_all_equal_left <- fill_gates_all_equal / 2 
            fill_gates_all_equal_right <- fill_gates_all_equal / 2 
          } else {
            fill_gates_all_equal_left <- (fill_gates_all_equal + 1) / 2 
            fill_gates_all_equal_right <- (fill_gates_all_equal - 1) / 2 
          }
          
          if (fill_gates_largest_diff %% 2 == 0) {
            fill_gates_largest_diff_left <- fill_gates_largest_diff / 2 
            fill_gates_largest_diff_right <- fill_gates_largest_diff / 2 
          } else {
            fill_gates_largest_diff_left <- (fill_gates_largest_diff + 1) / 2 
            fill_gates_largest_diff_right <- (fill_gates_largest_diff - 1) / 2 
          }
          
          cat(paste0(strrep(" ", fill_model_left), model_name, paste0(strrep(" ", fill_model_right))), "|", paste0(strrep(" ", fill_gates_all_equal_left), gates_all_equal, paste0(strrep(" ", fill_gates_all_equal_right))), "|", paste0(strrep(" ", fill_gates_largest_diff_left), gates_largest_diff, paste0(strrep(" ", fill_gates_largest_diff_right))), "| \n")
        } 
      }
      
      cli::cli_h1("")
    } else if (target == "RATE") {
      limits_table <- 6
      cli::cli_h1("RATEs RESULTS")
      
      for (j in seq_len(length(cate_models))) {
        cli::cli_h2(cate_models[j])
        
        cat(paste0(strrep(" ", limits_table), "WEIGHT", paste0(strrep(" ", limits_table))), "|", paste0(strrep(" ", limits_table), "RATE", paste0(strrep(" ", limits_table))), "|", "\n")
        cat(strrep("-", limits_table * 2 + nchar("WEIGHT")), "|", strrep("-", limits_table * 2 + nchar("RATE")), "|", " \n")
      
        autoc <- format(round(object$RATE[[cate_models[j]]]$rate_results$AUTOC$rate, 2), nsmall = 2)
        autoc_lower_ci <- format(round(object$RATE[[cate_models[j]]]$rate_results$AUTOC$rate - 1.96 * object$RATE[[cate_models[j]]]$rate_results$AUTOC$se, 3), nsmall = 3)
        autoc_upper_ci <- format(round(object$RATE[[cate_models[j]]]$rate_results$AUTOC$rate + 1.96 * object$RATE[[cate_models[j]]]$rate_results$AUTOC$se, 3), nsmall = 3)
      
        qini <- format(round(object$RATE[[cate_models[j]]]$rate_results$QINI$rate, 2), nsmall = 2)
        qini_lower_ci <- format(round(object$RATE[[cate_models[j]]]$rate_results$QINI$rate - 1.96 * object$RATE[[cate_models[j]]]$rate_results$QINI$se, 3), nsmall = 3)
        qini_upper_ci <- format(round(object$RATE[[cate_models[j]]]$rate_results$QINI$rate + 1.96 * object$RATE[[cate_models[j]]]$rate_results$QINI$se, 3), nsmall = 3)
      
        autoc_ci_print <- paste0("[", autoc_lower_ci, ", ", autoc_upper_ci, "]")
        qini_ci_print <- paste0("[", qini_lower_ci, ", ", qini_upper_ci, "]")
        
        fill_autoc <- limits_table * 2 + (nchar("RATE") - nchar(autoc))
        fill_qini <- limits_table * 2 + (nchar("RATE") - nchar(qini))
        fill_autoc_ci <- limits_table * 2 + (nchar("RATE") - nchar(autoc_ci_print))
        fill_qini_ci <- limits_table * 2 + (nchar("RATE") - nchar(qini_ci_print))
        
        if (fill_autoc %% 2 == 0) {
          fill_autoc_left <- fill_autoc / 2 
          fill_autoc_right <- fill_autoc / 2 
        } else {
          fill_autoc_left <- (fill_autoc + 1) / 2 
          fill_autoc_right <- (fill_autoc - 1) / 2 
        }
        
        if (fill_qini %% 2 == 0) {
          fill_qini_left <- fill_qini / 2 
          fill_qini_right <- fill_qini / 2 
        } else {
          fill_qini_left <- (fill_qini + 1) / 2 
          fill_qini_right <- (fill_qini - 1) / 2 
        }
        
        if (fill_autoc_ci %% 2 == 0) {
          fill_autoc_ci_left <- fill_autoc_ci / 2 
          fill_autoc_ci_right <- fill_autoc_ci / 2 
        } else {
          fill_autoc_ci_left <- (fill_autoc_ci + 1) / 2 
          fill_autoc_ci_right <- (fill_autoc_ci - 1) / 2 
        }
        
        if (fill_qini_ci %% 2 == 0) {
          fill_qini_ci_left <- fill_qini_ci / 2 
          fill_qini_ci_right <- fill_qini_ci / 2 
        } else {
          fill_qini_ci_left <- (fill_qini_ci + 1) / 2 
          fill_qini_ci_right <- (fill_qini_ci - 1) / 2 
        }
        
        cat("      ",   "autoc", "     ", "|", paste0(strrep(" ", fill_autoc_left), autoc, paste0(strrep(" ", fill_autoc_right))), "| \n")
        cat(strrep(" ", limits_table * 2 + nchar("WEIGHT")), "|", paste0(strrep(" ", fill_autoc_ci_left), autoc_ci_print, paste0(strrep(" ", fill_autoc_ci_right))), "| \n")    
        cat("      ",   "qini ", "     ", "|", paste0(strrep(" ", fill_qini_left), qini, paste0(strrep(" ", fill_qini_right))), "| \n")
        cat(strrep(" ", limits_table * 2 + nchar("WEIGHT")), "|", paste0(strrep(" ", fill_qini_ci_left), qini_ci_print, paste0(strrep(" ", fill_qini_ci_right))), "| \n")
      }
    }
  } else if (latex) {
    if (target == "BLP") {
      ates <- rep(0, n_cate_models * n_models_blp)
      ates_lower_ci <- rep(0, n_cate_models * n_models_blp)
      ates_upper_ci <- rep(0, n_cate_models * n_models_blp)
      
      hets <- rep(0, n_cate_models * n_models_blp)
      hets_lower_ci <- rep(0, n_cate_models * n_models_blp)
      hets_upper_ci <- rep(0, n_cate_models * n_models_blp)
      
      for (i in seq_len(n_cate_models)) {
        estimated_ates <- sapply(object$BLP[[i]], function(x) { format(round(stats::coef(x)[names(stats::coef(x)) == "beta1"], 2), nsmall = 2) })
        estimated_ates_lower_ci <- sapply(object$BLP[[i]], function(x) { format(round(x$conf.low["beta1"], 3), nsmall = 3) })
        estimated_ates_upper_ci <- sapply(object$BLP[[i]], function(x) { format(round(x$conf.high["beta1"], 3), nsmall = 3) })
        
        estimated_hets <- sapply(object$BLP[[i]], function(x) { format(round(stats::coef(x)[names(stats::coef(x)) == "beta2"], 2), nsmall = 2) })
        estimated_hets_lower_ci <- sapply(object$BLP[[i]], function(x) { format(round(x$conf.low["beta2"], 3), nsmall = 3) })
        estimated_hets_upper_ci <- sapply(object$BLP[[i]], function(x) { format(round(x$conf.high["beta2"], 3), nsmall = 3) })
        
        names(estimated_ates) <- names(estimated_ates_lower_ci) <- names (estimated_ates_upper_ci) <- names(estimated_hets) <- names(estimated_hets_lower_ci) <- names(estimated_hets_upper_ci) <- unlist(strsplit(names(estimated_ates), ".beta1"))
        
        which_models_blp <- setdiff(which_models, "imai_li")
        
        estimated_ates <- estimated_ates[which_models_blp]
        estimated_ates_lower_ci <- estimated_ates_lower_ci[which_models_blp]
        estimated_ates_upper_ci <- estimated_ates_upper_ci[which_models_blp]
        
        estimated_hets <- estimated_hets[which_models_blp]
        estimated_hets_lower_ci <- estimated_hets_lower_ci[which_models_blp]
        estimated_hets_upper_ci <- estimated_hets_upper_ci[which_models_blp] 
        
        ates[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_ates
        ates_lower_ci[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_ates_lower_ci
        ates_upper_ci[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_ates_upper_ci
        
        hets[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_hets
        hets_lower_ci[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_hets_lower_ci
        hets_upper_ci[seq_len(n_models_blp) + (n_models_blp * (i - 1))] <- estimated_hets_upper_ci
      }
      
      model_names_print <- rep(NA, length(which_models_blp))
      model_names_print[which_models_blp == "wr_none"] <- "WR"
      model_names_print[which_models_blp == "wr_cddf1"] <- "WR\\_{cddf1}"
      model_names_print[which_models_blp == "wr_cddf2"] <- "WR\\_{cddf2}"
      model_names_print[which_models_blp == "wr_mck1"] <- "WR\\_{mck1}"
      model_names_print[which_models_blp == "ht_none"] <- "HT"
      model_names_print[which_models_blp == "ht_cddf1"] <- "HT\\_{cddf1}"
      model_names_print[which_models_blp == "ht_cddf2"] <- "HT\\_{cddf2}"
      model_names_print[which_models_blp == "ht_mck1"] <- "HT\\_{mck1}"
      model_names_print[which_models_blp == "ht_mck2"] <- "HT\\_{mck2}"
      model_names_print[which_models_blp == "ht_mck3"] <- "HT\\_{mck3}"
      model_names_print[which_models_blp == "aipw"] <- "AIPW"
      
      cmid_points_start <- seq(2, n_models_blp * n_cate_models, by = n_models_blp)
      cmid_points_end <- seq(2+n_models_blp-1, n_models_blp * n_cate_models+1, by = n_models_blp)
      
      cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
      \\centering
      \\begin{adjustbox}{width = 1\\textwidth}
      \\begin{tabular}{@{\\extracolsep{5pt}}l", rep(" c", n_models_blp * n_cate_models), "}
      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      & ", stringr::str_sub(paste0("\\multicolumn{", n_models_blp, "}{c}{\\textit{", cate_models, "}} & ", collapse = " "), end = -3), " \\\\ ", paste0("\\cmidrule{", cmid_points_start, "-", cmid_points_end, "} "), "
      & ", stringr::str_sub(strrep(stringr::str_sub(paste0("\\textit{", model_names_print, "} & ", collapse = " "), end = -2), n_cate_models), end = -2), " \\\\
      \\addlinespace[2pt]
      \\hline \\\\[-1.8ex] \n\n", sep = "")
      
      cat("      \\multirow{2}{6em}{ATE ($\\hat{\\beta}_1$)} & ", stringr::str_sub(paste(paste0(ates, " & "), collapse = " "), end = -3), "\\\\ \n", sep = "")
      cat("      & ", stringr::str_sub(paste0("[", ates_lower_ci, ", ", ates_upper_ci,  "] & ", collapse = " "), end = -3), "\\\\ \n", sep = "")
      
      cat("      \\multirow{2}{6em}{HET ($\\hat{\\beta}_2$)} & ", stringr::str_sub(paste(paste0(hets, " & "), collapse = " "), end = -3), "\\\\ \n", sep = "")
      cat("      & ", stringr::str_sub(paste0("[", hets_lower_ci, ", ", hets_upper_ci,  "] & ", collapse = " "), end = -3), "\\\\ \n", sep = "")
      
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
      cat("NOT IMPLEMENTED, AS WE FEEL ANY TABLE WOULD BE MESSY. \n")
      cat("WE APPRECIATE ANY SUGGESTIONS ON NICE FORMATS YOU WOULD LIKE TO HAVE.")
      cat("MAY WE SUGGEST TO USE THE PRINT METHOD TO PRODUCE A NICE PLOT? \n")
    } else if (target == "RATE") {
      cat("YET TO BE IMPLEMENTED, SORRY :( \n")
      cat("WE APPRECIATE ANY SUGGESTIONS ON THE FORMAT.")
    }
  } 
}


#' Print Method for evaluCATE Objects
#'
#' Prints an \code{\link{evaluCATE}} object.
#'
#' @param x An \code{\link{evaluCATE}} object.
#' @param target String controlling which parameters we are interested in. Must be one of \code{"BLP"}, \code{"GATES"}, \code{"RATE"}.
#' @param latex Logical, whether to print LATEX code for a table. Different tables are produced according to the \code{target} argument.
#' @param which_models Character vector controlling which results to display. Admitted values are those stored in \code{names(object$GATES)}. Ignored if \code{target == "RATE"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' Prints an \code{\link{evaluCATE}} object.
#' 
#' @details 
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}, \code{multirow}.
#' 
#' @author Riccardo Di Francesco
#' 
#' @export
print.evaluCATE <- function(x, target = "BLP", 
                            latex = FALSE, which_models = names(x$GATES[[1]]), ...) {
  summary.evaluCATE(x, target, latex, which_models, ...)
}


#' Plot Method for evaluCATE Objects
#'
#' Plots an \code{\link{evaluCATE}} object.
#'
#' @param x An \code{\link{evaluCATE}} object.
#' @param target String controlling which plot to display. Must be one of \code{"GATES"}, \code{"TOC"}, or \code{"RATE"}.
#' @param which_models Character vector controlling which results to display. Admitted values are those stored in \code{names(x$GATES[[1]])}. Ignored if \code{target != "GATES"}.
#' @param gates_hline Logical, whether to display an horizontal line at zero in the GATES plot. Ignored if \code{target != "GATES"}.
#' @param toc_smoother Integer number, controls the amount of smoothing in the \code{TOC} and \code{RATE} plots. Smoothing is achieved by plotting the TOCs only for every other \code{toc_smoother} unit (in order of treatment benefit). Set to 1 to plot the TOC of all units, to 2 to plot the TOC of every other unit, ecc. Ignored if \code{target == "GATES"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Plots an \code{\link{evaluCATE}} object.
#' 
#' @details
#' Check the online \href{https://riccardo-df.github.io/evaluCATE/articles/more-on-plotting.html}{plotting vignette} for details.
#' 
#'
#' @import dplyr ggplot2 ggsci
#' @importFrom stats coef
#' @importFrom gridExtra grid.arrange
#'
#' @author Riccardo Di Francesco
#'
#' @export
plot.evaluCATE <- function(x, target = "GATES", which_models = names(x$GATES[[1]]), 
                           gates_hline = TRUE, toc_smoother = 1, ...) {
  ## Checks.
  if (!(target %in% c("GATES", "TOC", "RATE"))) stop("Invalid 'target'. This must be one of 'GATES', 'TOC', or 'RATE'.", call. = FALSE)
  if (any(!(which_models %in% names(x$GATES[[1]])))) stop("Invalid 'which_models'. Check the documentation for admitted values.", call. = FALSE)
  if (!is.logical(gates_hline)) stop("Invalid 'gates_hline'. This must be either 'TRUE' or 'FALSE'.", call. = FALSE)
  if (!is.numeric(toc_smoother) & toc_smoother < 1 & toc_smoother %% 1 != 0) stop("Invalid 'toc_smoother'. This must be an integer greater than or equal to one.", call. = FALSE)
  
  group <- NULL
  estimator <- NULL
  estimated_gate <- NULL
  se <- NULL
  strategy <- NULL
  denoise <- NULL
  value <- NULL
  variable <- NULL
  model <- NULL
  CATE_estimator <- NULL
  
  u <- NULL
  TOC <- NULL
  
  cate_models <- names(x$BLP)
  n_cate_models <- length(cate_models)
  
  ## Plot.
  if (target == "GATES") {
    ## Extract information from main call.
    n_models_gates <- length(which_models)
    n_groups <- length(x$GATES[[1]]$wr_none$coefficients)
    
    ## Extract estimated GATES and standard errors.
    list_plot_dta <- list()
    
    for (j in seq_len(n_cate_models)) {
      gates_list <- lapply(x$GATES[[j]][1:length(x$GATES[[j]])-1], function(y) {y$coefficients[grep("group", names(y$coefficients))]}) %>%
        append(list(c(x$GATES[[j]]$imai_li$GATE)))
      names(gates_list)[length(gates_list)] <- "imai_li"
      
      gates_se_list <- lapply(x$GATES[[j]][1:length(x$GATES[[j]])-1], function(y) {summary(y)$coefficients[grep("group", names(y$coefficients)), 2]}) %>%
        append(list(c(x$GATES[[j]]$imai_li$SE)))
      names(gates_se_list)[length(gates_se_list)] <- "imai_li"
      
      list_plot_dta[[j]] <- data.frame("group" = rep(1:n_groups, length(gates_list)), 
                                       "CATE_estimator" = cate_models[j],
                                       "model" = rep(names(gates_list), each = n_groups),
                                       "estimated_gate" = unlist(gates_list),
                                       "se" = unlist(gates_se_list)) %>%
        dplyr::filter(model %in% which_models)
    }
    
    ## Arrange plot data and call ggplot.
    plot_dta <- dplyr::bind_rows(list_plot_dta)
    
    plot_dta$strategy <- NA
    plot_dta$strategy[plot_dta$model %in% c("wr_none", "wr_cddf1", "wr_cddf2", "wr_mck1")] <- "WR"
    plot_dta$strategy[plot_dta$model %in% c("ht_none", "ht_cddf1", "ht_cddf2", "ht_mck1", "ht_mck2", "ht_mck3")] <- "HT"
    plot_dta$strategy[plot_dta$model == "aipw"] <- "AIPW"
    plot_dta$strategy[plot_dta$model == "imai_li"] <- "Nonparametric"
    
    plot_dta$denoise <- NA
    plot_dta$denoise[plot_dta$model %in% c("wr_none", "ht_none", "aipw", "imai_li")] <- "None"
    plot_dta$denoise[plot_dta$model %in% c("wr_cddf1", "ht_cddf1")] <- "cddf1"
    plot_dta$denoise[plot_dta$model %in% c("wr_cddf2", "ht_cddf2")] <- "cddf2"
    plot_dta$denoise[plot_dta$model %in% c("wr_mck1", "ht_mck1")] <- "mck1"
    plot_dta$denoise[plot_dta$model %in% c("ht_mck2")] <- "mck2"
    plot_dta$denoise[plot_dta$model %in% c("ht_mck3")] <- "mck3"
    
    y_axis_top <- max(plot_dta$estimated_gate + 1.96 * plot_dta$se) + 0.01 * mean(plot_dta$estimated_gate)
    y_axis_bottom <- min(plot_dta$estimated_gate - 1.96 * plot_dta$se) - 0.01 * mean(plot_dta$estimated_gate)
    
    facet_condition <- length(unique(factor(plot_dta$denoise, levels = c("None", "cddf1", "cddf2", "mck1", "mck2", "mck3")))) > 1
    
    plot_list <- list()
    
    for (j in seq_len(n_cate_models)) {
      plot_list[[j]] <- plot_dta %>%
        dplyr::filter(CATE_estimator == cate_models[j]) %>%
        ggplot2::ggplot(ggplot2::aes(x = group, y = estimated_gate, colour = factor(model, levels = which_models))) +
        ggplot2::geom_point() + 
        ggplot2::geom_hline(yintercept = if (gates_hline) 0 else NULL, linetype = "dashed") + 
        ggplot2::geom_errorbar(aes(x = group, ymin = estimated_gate - 1.96 * se, ymax = estimated_gate + 1.96 * se)) +
        ggplot2::facet_wrap(vars(CATE_estimator)) +
        ggplot2::facet_grid(cols = vars(factor(strategy, levels = c("WR", "HT", "AIPW", "Nonparametric"))), rows = if (facet_condition) vars(factor(denoise, levels = c("None", "cddf1", "cddf2", "mck1", "mck2", "mck3"))) else NULL, scales = "fixed") +
        ggplot2::xlab("") + ggplot2::ylab("GATES") + ggplot2::labs(colour = "Estimator") + ggplot2::ggtitle(cate_models[j]) +
        ggplot2::ylim(y_axis_bottom, y_axis_top) +
        ggplot2::theme_bw() + 
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"), legend.position = "none",
                       plot.title = ggplot2::element_text(size = 15, face = "italic", colour = "black", hjust = 0.5))
    }
    
    gridExtra::grid.arrange(grobs = plot_list) 
  } else if (target == "TOC") {
    list_plot_dta <- list()
    
    for (j in seq_len(n_cate_models)) {
      list_plot_dta[[j]] <- data.frame("unit" = (1:length(x$RATE[[j]]$toc_results) / length(x$RATE[[j]]$toc_results)), 
                                       "TOC" = x$RATE[[j]]$toc_results,
                                       "CATE_estimator" = cate_models[j])
    }
    
    plot_dta <- dplyr::bind_rows(list_plot_dta)
    
    if (toc_smoother == 1) {
      plot_dta %>%
        ggplot2::ggplot(ggplot2::aes(x = unit, y = TOC)) +
        ggplot2::geom_line(color = "tomato", linewidth = 1.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_grid(cols = vars(factor(CATE_estimator, levels = cate_models))) +
        ggplot2::xlab("Fraction treated") +ggplot2:: ylab("TOC") +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"))
    } else {
      plot_dta %>%
        dplyr::slice(which(row_number() %% toc_smoother == 1)) %>%
        ggplot2::ggplot(ggplot2::aes(x = unit, y = TOC)) +
        ggplot2::geom_line(color = "tomato", linewidth = 1.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_grid(cols = vars(factor(CATE_estimator, levels = cate_models))) +
        ggplot2::xlab("Fraction treated") +ggplot2:: ylab("TOC") +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"))
    }
  } else if (target == "RATE") {
    list_plot_dta <- list()
    
    for (j in seq_len(n_cate_models)) {
      list_plot_dta[[j]] <- data.frame("unit" = (1:length(x$RATE[[j]]$toc_results) / length(x$RATE[[j]]$toc_results)), 
                                       "TOC" = x$RATE[[j]]$toc_results,
                                       "CATE_estimator" = cate_models[j])
    }
    
    plot_dta <- dplyr::bind_rows(list_plot_dta)
    
    if (toc_smoother == 1) {
      plot_dta %>%
        dplyr::mutate(TOCu = TOC * unit) %>% 
        reshape2::melt(id.vars = c("unit", "CATE_estimator"), measure.vars = c("TOC", "TOCu")) %>%
        ggplot2::ggplot(ggplot2::aes(x = unit, y = value, group = variable, colour = variable)) +
        ggplot2::geom_line(linewidth = 1.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_color_hue(labels = c(expression("TOC"), expression("u" %*% "TOC"))) +
        ggplot2::facet_grid(cols = vars(factor(CATE_estimator, levels = cate_models))) +
        ggplot2::xlab("Fraction treated") + ggplot2:: ylab("") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = c(0.9, 0.9), legend.title = ggplot2::element_blank(), legend.direction = "vertical", legend.text = element_text(size = 7),
                       strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"))
    } else {
      plot_dta %>%
        dplyr::mutate(TOCu = TOC * unit) %>% 
        dplyr::slice(which(row_number() %% toc_smoother == 1)) %>%
        reshape2::melt(id.vars = c("unit", "CATE_estimator"), measure.vars = c("TOC", "TOCu")) %>%
        ggplot2::ggplot(ggplot2::aes(x = unit, y = value, group = variable, colour = variable)) +
        ggplot2::geom_line(linewidth = 1.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_color_hue(labels = c(expression("TOC"), expression("u" %*% "TOC"))) +
        ggplot2::facet_grid(cols = vars(factor(CATE_estimator, levels = cate_models))) +
        ggplot2::xlab("Fraction treated") + ggplot2:: ylab("") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = c(0.9, 0.9), legend.title = ggplot2::element_blank(), legend.direction = "vertical", legend.text = element_text(size = 7),
                       strip.text.x = ggplot2::element_text(size = 10, face = "italic"), strip.text.y = ggplot2::element_text(size = 10, face = "italic"))
    }
  }
}
