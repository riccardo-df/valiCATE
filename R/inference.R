#' Confidence Interval Construction
#'
#' Builds a two-sided confidence interval from an estimate, standard error, and significance level.
#'
#' @param estimate Point estimate.
#' @param se Standard error.
#' @param alpha Significance level.
#'
#' @return
#' Named numeric vector of length 2 with elements \code{lower} and \code{upper}.
#'
#' @keywords internal
make_ci <- function(estimate, se, alpha) {
  z <- stats::qnorm(1 - alpha / 2)

  ## Output.
  return(c("lower" = estimate - z * se, "upper" = estimate + z * se))
}


#' One-Sided Test for CWATE
#'
#' Tests H0: theta <= 0 against H1: theta > 0.
#'
#' @param estimate Point estimate theta_hat.
#' @param se Standard error.
#'
#' @return
#' A list with elements \code{test_stat} (z-statistic) and \code{p_value}.
#'
#' @keywords internal
test_cwate <- function(estimate, se) {
  z <- estimate / se
  p <- 1 - stats::pnorm(z)

  ## Output.
  return(list("test_stat" = z, "p_value" = p))
}


#' Two-Sided Test for NCWATE
#'
#' Tests H0: gamma = 1 against H1: gamma != 1.
#'
#' @param estimate Point estimate gamma_hat.
#' @param se Standard error.
#'
#' @return
#' A list with elements \code{test_stat} (z-statistic) and \code{p_value}.
#'
#' @keywords internal
test_ncwate <- function(estimate, se) {
  z <- (estimate - 1) / se
  p <- 2 * (1 - stats::pnorm(abs(z)))

  ## Output.
  return(list("test_stat" = z, "p_value" = p))
}


#' Assemble Inference Results
#'
#' Packages point estimate, standard error, confidence interval, and hypothesis test results into a single list.
#'
#' @param estimate Point estimate.
#' @param se Standard error.
#' @param alpha Significance level.
#' @param type Character, either \code{"cwate"} for a one-sided test or \code{"ncwate"} for a two-sided test.
#'
#' @return
#' A list with elements \code{estimate}, \code{se}, \code{ci} (named vector with \code{lower} and \code{upper}),
#' \code{test_stat}, and \code{p_value}.
#'
#' @keywords internal
assemble_inference <- function(estimate, se, alpha, type = c("cwate", "ncwate")) {
  type <- match.arg(type)

  ## 1.) Confidence interval.
  ci <- make_ci(estimate, se, alpha)

  ## 2.) Hypothesis test.
  if (type == "cwate") {
    tst <- test_cwate(estimate, se)
  } else {
    tst <- test_ncwate(estimate, se)
  }

  ## Output.
  return(list("estimate" = estimate, "se" = se, "ci" = ci, "test_stat" = tst$test_stat, "p_value" = tst$p_value))
}
