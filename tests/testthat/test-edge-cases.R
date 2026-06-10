test_that("valiCATE rejects non-binary D", {
  X <- matrix(rnorm(20), ncol = 2)
  colnames(X) <- c("x1", "x2")
  expect_error(
    valiCATE(rnorm(10), c(0, 1, 2, 0, 1, 0, 1, 0, 1, 0), X, cates = list("m" = rnorm(10)), scores = rnorm(10)),
    "binary"
  )
})

test_that("valiCATE rejects invalid weights", {
  X <- matrix(rnorm(20), ncol = 2)
  colnames(X) <- c("x1", "x2")
  expect_error(
    valiCATE(rnorm(10), rbinom(10, 1, 0.5), X, cates = list("m" = rnorm(10)), scores = rnorm(10), weights = "INVALID"),
    "Invalid 'weights'"
  )
})

test_that("valiCATE rejects zero-variance cates", {
  X <- matrix(rnorm(20), ncol = 2)
  colnames(X) <- c("x1", "x2")
  expect_error(
    valiCATE(rnorm(10), rbinom(10, 1, 0.5), X, cates = list("m" = rep(1, 10)), scores = rnorm(10)),
    "No variation"
  )
})

test_that("valiCATE rejects mismatched lengths", {
  X <- matrix(rnorm(20), ncol = 2)
  colnames(X) <- c("x1", "x2")
  expect_error(
    valiCATE(rnorm(10), rbinom(10, 1, 0.5), X, cates = list("m" = rnorm(5)), scores = rnorm(10)),
    "length"
  )
})

test_that("valiCATE rejects invalid alpha", {
  X <- matrix(rnorm(20), ncol = 2)
  colnames(X) <- c("x1", "x2")
  expect_error(
    valiCATE(rnorm(10), rbinom(10, 1, 0.5), X, cates = list("m" = rnorm(10)), scores = rnorm(10), alpha = 1.5),
    "alpha"
  )
})

test_that("inference functions work correctly", {
  ci <- valiCATE:::make_ci(1.0, 0.5, 0.05)
  expect_equal(unname(ci["lower"]), 1.0 - 1.96 * 0.5, tolerance = 0.01)
  expect_equal(unname(ci["upper"]), 1.0 + 1.96 * 0.5, tolerance = 0.01)

  ## One-sided test: positive estimate should have small p-value.
  tst <- valiCATE:::test_cwate(3.0, 1.0)
  expect_true(tst$p_value < 0.01)

  ## Two-sided test: estimate near null should have large p-value.
  tst2 <- valiCATE:::test_ncwate(1.0, 0.5)
  expect_true(tst2$p_value > 0.5)
})
