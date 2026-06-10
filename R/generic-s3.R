#' Summary Method for valiCATE Objects
#'
#' Summarizes an \code{\link{valiCATE}} object.
#'
#' @param object An \code{\link{valiCATE}} object.
#' @param latex Logical, whether to print LATEX code for a table. If \code{TRUE}, a ready-to-compile LATEX table with CWATE and
#'   NCWATE point estimates, confidence intervals, and p-values is printed to the console.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Summarizes an \code{\link{valiCATE}} object.
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
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n, sd = 0.5)
#'
#' ## Split into training and validation samples.
#' train_idx <- sample(1:n, n / 2)
#' val_idx <- setdiff(1:n, train_idx)
#'
#' ## Estimate CATEs on the training sample, predict on the validation sample.
#' library(grf)
#' cf <- causal_forest(X[train_idx, ], Y[train_idx], D[train_idx])
#' cates <- predict(cf, X[val_idx, ])$predictions
#'
#' ## Validate using the validation sample.
#' result <- valiCATE(Y[val_idx], D[val_idx], X[val_idx, ],
#'                    cates = list("causal_forest" = cates))
#'
#' ## Summarize.
#' summary(result)
#' summary(result, latex = TRUE)}
#'
#' @details
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{valiCATE}}
#'
#' @export
summary.valiCATE <- function(object, latex = FALSE, ...) {
  ## Handling inputs and checks.
  if (!is.logical(latex)) stop("Invalid 'latex'. This must be either TRUE or FALSE.", call. = FALSE)

  if (!latex) {
    summary_console(object)
  } else {
    summary_latex(object)
  }

  invisible(object)
}


#' Print Method for valiCATE Objects
#'
#' Prints an \code{\link{valiCATE}} object.
#'
#' @param x An \code{\link{valiCATE}} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' Prints an \code{\link{valiCATE}} object.
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{valiCATE}}
#'
#' @export
print.valiCATE <- function(x, ...) {
  summary.valiCATE(x, ...)
}


#' Plot Method for valiCATE Objects
#'
#' Produces a coefficient plot of CWATE or NCWATE estimates with confidence intervals.
#'
#' @param x An \code{\link{valiCATE}} object.
#' @param type Character, either \code{"cwate"} or \code{"ncwate"}. Controls which estimand is plotted.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' A \code{ggplot} object, returned invisibly.
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
#' Y <- mu0 + D * (mu1 - mu0) + rnorm(n, sd = 0.5)
#'
#' ## Split into training and validation samples.
#' train_idx <- sample(1:n, n / 2)
#' val_idx <- setdiff(1:n, train_idx)
#'
#' ## Estimate CATEs on the training sample, predict on the validation sample.
#' library(grf)
#' cf <- causal_forest(X[train_idx, ], Y[train_idx], D[train_idx])
#' cates <- predict(cf, X[val_idx, ])$predictions
#'
#' ## Validate using the validation sample.
#' result <- valiCATE(Y[val_idx], D[val_idx], X[val_idx, ],
#'                    cates = list("causal_forest" = cates))
#'
#' ## Plot.
#' plot(result)
#' plot(result, type = "ncwate")}
#'
#' @details
#' This method requires the \code{ggplot2} package. The plot displays point estimates as dots and
#' confidence intervals as vertical whiskers, faceted by model. A horizontal reference line is drawn
#' at 0 for the CWATE (null: no heterogeneity) and at 1 for the NCWATE (null: recovery).
#'
#' @author Riccardo Di Francesco
#'
#' @references
#' \itemize{
#'   \item Di Francesco, R., & Knaus, M. C. (2025). Validating ML Predictions of Heterogeneous Treatment Effects via CWATE.
#' }
#'
#' @seealso \code{\link{valiCATE}}, \code{\link{summary.valiCATE}}
#'
#' @export
plot.valiCATE <- function(x, type = "cwate", ...) {
  ## 1.) Handling inputs and checks.
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it with install.packages('ggplot2').", call. = FALSE)
  }

  valid_types <- c("cwate", "ncwate")
  if (!(type %in% valid_types)) {
    stop("Invalid 'type'. This must be either \"cwate\" or \"ncwate\".", call. = FALSE)
  }

  ## 2.) Extract results into a data frame.
  model_names <- names(x$results)
  rows <- list()
  for (model_nm in model_names) {
    for (wt in x$weights) {
      res <- x$results[[model_nm]][[wt]][[type]]
      rows[[length(rows) + 1]] <- data.frame(
        model = model_nm,
        weight = wt,
        estimate = res$estimate,
        ci_lower = res$ci["lower"],
        ci_upper = res$ci["upper"],
        stringsAsFactors = FALSE
      )
    }
  }

  plot_df <- do.call(rbind, rows)
  rownames(plot_df) <- NULL

  ## 3.) Preserve ordering.
  plot_df$model <- factor(plot_df$model, levels = model_names)
  plot_df$weight <- factor(plot_df$weight, levels = x$weights)

  ## 4.) Color palette.
  palette <- c("AUTOC" = "#2166AC", "AUC-HVL" = "#B2182B",
               "BLP" = "#1B7837", "QINI" = "#E08214")
  colors <- palette[x$weights]

  ## 5.) Reference line (0 for CWATE, 1 for NCWATE).
  ref_value <- if (type == "cwate") 0 else 1

  ## 6.) Build plot.
  estimate <- ci_lower <- ci_upper <- weight <- model <- NULL
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = weight, y = estimate, color = weight)) +
    ggplot2::geom_hline(yintercept = ref_value, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                           width = 0.15, linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::facet_wrap(~ model, nrow = 1) +
    ggplot2::labs(x = NULL, y = toupper(type)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      axis.text.x = ggplot2::element_text(color = "black", angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4)
    )

  ## Output.
  print(p)
  invisible(p)
}


#' Console Summary
#'
#' Displays CWATE and NCWATE results in the console.
#'
#' @param object An \code{\link{valiCATE}} object.
#'
#' @return
#' Console output.
#'
#' @keywords internal
summary_console <- function(object) {
  alpha <- object$alpha
  ci_level <- paste0(round((1 - alpha) * 100), "%")

  cate_models <- names(object$results)

  for (model_nm in cate_models) {
    cat(strrep("=", 70), "\n")
    cat("Model:", model_nm, "\n")
    cat(strrep("=", 70), "\n\n")

    ## CWATE table.
    cat("CWATE (H0: theta <= 0, one-sided test)\n")
    cat(strrep("-", 70), "\n")
    cat(sprintf("  %-10s %10s %10s %16s %10s\n", "Weight", "Estimate", "SE", paste0(ci_level, " CI"), "p-value"))
    cat(strrep("-", 70), "\n")

    for (wt in object$weights) {
      res <- object$results[[model_nm]][[wt]]$cwate
      ci_str <- sprintf("[%.4f, %.4f]", res$ci["lower"], res$ci["upper"])
      cat(sprintf("  %-10s %10.4f %10.4f %16s %10.4f\n", wt, res$estimate, res$se, ci_str, res$p_value))
    }

    cat("\n")

    ## NCWATE table.
    cat("NCWATE (H0: gamma = 1, two-sided test)\n")
    cat(strrep("-", 70), "\n")
    cat(sprintf("  %-10s %10s %10s %16s %10s\n", "Weight", "Estimate", "SE", paste0(ci_level, " CI"), "p-value"))
    cat(strrep("-", 70), "\n")

    for (wt in object$weights) {
      res <- object$results[[model_nm]][[wt]]$ncwate
      ci_str <- sprintf("[%.4f, %.4f]", res$ci["lower"], res$ci["upper"])
      cat(sprintf("  %-10s %10.4f %10.4f %16s %10.4f\n", wt, res$estimate, res$se, ci_str, res$p_value))
    }

    cat("\n")
  }
}


#' LaTeX Summary
#'
#' Produces LATEX table code for CWATE and NCWATE results.
#'
#' @param object An \code{\link{valiCATE}} object.
#'
#' @return
#' LATEX code printed to console.
#'
#' @details
#' Compilation of the LATEX code requires the following packages: \code{booktabs}, \code{float}, \code{adjustbox}.
#'
#' @keywords internal
summary_latex <- function(object) {
  alpha <- object$alpha
  ci_level <- paste0(round((1 - alpha) * 100), "\\%")
  model_names <- names(object$results)
  n_models <- length(model_names)

  ## Escape underscores in model names for LATEX.
  model_names_tex <- gsub("_", "\\\\_", model_names)

  cat("\\begingroup
    \\setlength{\\tabcolsep}{8pt}
    \\renewcommand{\\arraystretch}{1.1}
    \\begin{table}[H]
      \\centering
      \\caption{CWATE and NCWATE results.}
      \\vspace{-0.3cm}
      \\label{table_cwate_ncwate}
      \\begin{adjustbox}{width = 1\\textwidth}\n")

  col_spec <- paste0("l", paste(rep(" c c", n_models), collapse = ""))
  cat("      \\begin{tabular}{@{\\extracolsep{5pt}}", col_spec, "}\n")
  cat("      \\\\[-1.8ex]\\hline\n")
  cat("      \\hline \\\\[-1.8ex]\n")

  ## Header row 1: model names spanning 2 columns each.
  header1 <- paste0("\\multicolumn{2}{c}{\\textit{", model_names_tex, "}}", collapse = " & ")
  cmidrule <- paste0(sapply(seq_len(n_models), function(i) {
    start <- 2 * (i - 1) + 2
    paste0("\\cmidrule{", start, "-", start + 1, "}")
  }), collapse = " ")
  cat("      & ", header1, " \\\\ ", cmidrule, "\n")

  ## Header row 2: CWATE / NCWATE for each model.
  header2 <- paste(rep("CWATE & NCWATE", n_models), collapse = " & ")
  cat("      & ", header2, " \\\\\n")
  cat("      \\addlinespace[2pt]\n")
  cat("      \\hline \\\\[-1.8ex] \n\n")

  ## Data rows.
  for (wt in object$weights) {
    wt_tex <- gsub("-", "--", wt)

    ## Point estimates.
    vals <- sapply(model_names, function(m) {
      cwate_est <- sprintf("%.4f", object$results[[m]][[wt]]$cwate$estimate)
      ncwate_est <- sprintf("%.4f", object$results[[m]][[wt]]$ncwate$estimate)
      paste(cwate_est, "&", ncwate_est)
    })
    cat("      ", wt_tex, " & ", paste(vals, collapse = " & "), " \\\\\n")

    ## Confidence intervals.
    cis <- sapply(model_names, function(m) {
      cwate_ci <- object$results[[m]][[wt]]$cwate$ci
      ncwate_ci <- object$results[[m]][[wt]]$ncwate$ci
      cwate_str <- sprintf("[%.3f, %.3f]", cwate_ci["lower"], cwate_ci["upper"])
      ncwate_str <- sprintf("[%.3f, %.3f]", ncwate_ci["lower"], ncwate_ci["upper"])
      paste(cwate_str, "&", ncwate_str)
    })
    cat("      & ", paste(cis, collapse = " & "), " \\\\\n")

    ## p-values.
    pvals <- sapply(model_names, function(m) {
      cwate_p <- sprintf("(%.4f)", object$results[[m]][[wt]]$cwate$p_value)
      ncwate_p <- sprintf("(%.4f)", object$results[[m]][[wt]]$ncwate$p_value)
      paste(cwate_p, "&", ncwate_p)
    })
    cat("      & ", paste(pvals, collapse = " & "), " \\\\\n")
    cat("      \\addlinespace[3pt]\n")
  }

  cat("\n      \\\\[-1.8ex]\\hline
      \\hline \\\\[-1.8ex]
      \\end{tabular}
      \\end{adjustbox}
      \\begin{minipage}{1\\textwidth}
      \\scriptsize
      \\renewcommand{\\baselineskip}{11pt}
      \\vspace{0.05cm}
      \\textit{Notes.} ", ci_level,
      " confidence intervals in brackets, $p$-values in parentheses. CWATE: $\\mathcal{H}_0\\!: \\theta \\leq 0$ (one-sided). NCWATE: $\\mathcal{H}_0\\!: \\gamma = 1$ (two-sided).
      \\end{minipage}
    \\end{table}
\\endgroup\n", sep = "")
}
