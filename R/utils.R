#' Input Validation
#'
#' Validates inputs to the main \code{\link{valiCATE}} function.
#'
#' @param Y Outcome vector.
#' @param D Treatment indicator vector (binary).
#' @param X Covariate matrix or data frame (no intercept).
#' @param cates Named list of CATE predictions.
#' @param weights Character vector of weight names.
#' @param scores Pre-computed AIPW pseudo-outcomes, or \code{NULL}.
#' @param n_folds Number of cross-fitting folds.
#' @param alpha Significance level.
#' @param verbose Logical flag for printing progress.
#'
#' @return
#' Invisibly returns \code{TRUE} if all checks pass.
#'
#' @keywords internal
validate_inputs <- function(Y, D, X, cates, weights, scores,
                            n_folds, alpha, verbose) {
  ## Store useful quantities.
  n <- length(Y)

  ## Run checks.
  # Y.
  if (!is.numeric(Y)) stop("Invalid 'Y'. This must be a numeric vector.", call. = FALSE)

  # D.
  if (!is.numeric(D)) stop("Invalid 'D'. This must be a numeric vector.", call. = FALSE)
  if (length(D) != n) stop("Invalid 'D'. This must have the same length as 'Y'.", call. = FALSE)
  if (any(!(D %in% c(0, 1)))) stop("Invalid 'D'. Only binary treatments are allowed.", call. = FALSE)

  # X.
  if (!is.matrix(X) & !is.data.frame(X)) stop("Invalid 'X'. This must be either a matrix or a data frame.", call. = FALSE)
  if (nrow(X) != n) stop("Invalid 'X'. This must have the same number of rows as the length of 'Y'.", call. = FALSE)

  # cates.
  if (!is.list(cates)) stop("Invalid 'cates'. This must be a named list.", call. = FALSE)
  if (is.null(names(cates)) | any(names(cates) == "")) stop("Invalid 'cates'. This must be a named list with all elements named.", call. = FALSE)
  for (nm in names(cates)) {
    if (!is.numeric(cates[[nm]]) | length(cates[[nm]]) != n) stop(paste0("Invalid 'cates'. Element '", nm, "' must be a numeric vector of length ", n, "."), call. = FALSE)
    if (stats::var(cates[[nm]]) == 0) stop(paste0("No variation in element '", nm, "' of 'cates'."), call. = FALSE)
  }

  # weights.
  valid_weights <- c("AUTOC", "AUC-HVL", "BLP", "QINI")
  if (any(!(weights %in% valid_weights))) stop(paste0("Invalid 'weights'. Admitted values are: ", paste(valid_weights, collapse = ", "), "."), call. = FALSE)

  # scores.
  if (!is.null(scores)) if (!is.numeric(scores) | length(scores) != n) stop("Invalid 'scores'. This must be a numeric vector of the same length as 'Y'.", call. = FALSE)

  # n_folds.
  if (!is.numeric(n_folds) | length(n_folds) != 1 | n_folds < 2 | n_folds != round(n_folds)) stop("Invalid 'n_folds'. This must be an integer greater than or equal to 2.", call. = FALSE)

  # alpha.
  if (!is.numeric(alpha) | length(alpha) != 1 | alpha <= 0 | alpha >= 1) stop("Invalid 'alpha'. This must be a number strictly between 0 and 1.", call. = FALSE)

  # verbose.
  if (!is.logical(verbose) | length(verbose) != 1) stop("Invalid 'verbose'. This must be a single logical value.", call. = FALSE)

  ## Output.
  invisible(TRUE)
}


#' Safe Empirical CDF
#'
#' Computes the empirical CDF using the Hazen plotting position \code{(rank(x) - 0.5) / n}, clipped to \code{[eps, 1 - eps]}
#' to avoid boundary singularities in weight functions like AUTOC and AUC-HVL.
#'
#' @param x Numeric vector.
#'
#' @return
#' Numeric vector of empirical CDF values in (0, 1).
#'
#' @keywords internal
ecdf_safe <- function(x) {
  n <- length(x)
  eps <- 1e-3

  ecdf <- (rank(x, ties.method = "average") - 0.5) / n

  ## Output.
  return(pmin(pmax(ecdf, eps), 1 - eps))
}
