test_that("valiCATE works with pre-computed scores", {
  set.seed(42)
  n <- 200
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  tau <- X[, 2]
  Y <- 0.5 * X[, 1] + D * tau + rnorm(n)
  tau_hat <- tau + rnorm(n, sd = 0.3)
  scores <- tau + rnorm(n, sd = 0.5)

  res <- valiCATE(Y, D, X, cates = list("model1" = tau_hat),
                  scores = scores, verbose = FALSE)

  expect_s3_class(res, "valiCATE")
  expect_equal(res$n, n)
  expect_equal(res$weights, c("AUTOC", "AUC-HVL", "BLP", "QINI"))
  expect_true("model1" %in% names(res$results))

  for (wt in res$weights) {
    expect_true(is.finite(res$results$model1[[wt]]$cwate$estimate))
    expect_true(is.finite(res$results$model1[[wt]]$ncwate$estimate))
  }
})

test_that("valiCATE rejects a bare numeric vector for cates", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  tau_hat <- rnorm(n)
  scores <- rnorm(n)

  expect_error(
    valiCATE(Y, D, X, cates = tau_hat, scores = scores, verbose = FALSE),
    "named list"
  )
})

test_that("valiCATE works with multiple models", {
  set.seed(42)
  n <- 200
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  cates_list <- list("model_A" = rnorm(n), "model_B" = rnorm(n) + 0.5)
  res <- valiCATE(Y, D, X, cates = cates_list, scores = scores, verbose = FALSE)

  expect_equal(length(res$results), 2)
  expect_true("model_A" %in% names(res$results))
  expect_true("model_B" %in% names(res$results))
})

test_that("valiCATE works with subset of weights", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores,
                  weights = c("BLP", "QINI"), verbose = FALSE)

  expect_equal(res$weights, c("BLP", "QINI"))
  expect_equal(length(res$results$model), 2)
})

test_that("print.valiCATE runs without error", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores, verbose = FALSE)
  expect_output(print(res), "CWATE")
})

test_that("summary.valiCATE runs without error", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores, verbose = FALSE)
  expect_output(summary(res), "CWATE")
  expect_output(summary(res), "NCWATE")
})

test_that("plot.valiCATE returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores, verbose = FALSE)
  p <- plot(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot.valiCATE works with type = 'ncwate' and multiple models", {
  skip_if_not_installed("ggplot2")
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("m1" = rnorm(n), "m2" = rnorm(n)),
                  scores = scores, verbose = FALSE)
  p <- plot(res, type = "ncwate")
  expect_s3_class(p, "ggplot")
})

test_that("plot.valiCATE rejects invalid type", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores, verbose = FALSE)
  expect_error(plot(res, type = "invalid"), "Invalid 'type'")
})

test_that("summary.valiCATE latex output runs without error", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  scores <- rnorm(n)

  res <- valiCATE(Y, D, X, cates = list("model" = rnorm(n)), scores = scores, verbose = FALSE)
  expect_output(summary(res, latex = TRUE), "begin\\{table\\}")
})
