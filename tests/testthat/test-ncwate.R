test_that("NCWATE returns correct structure", {
  set.seed(42)
  n <- 200
  tau_hat <- rnorm(n) + 1
  Gamma <- tau_hat + rnorm(n, sd = 0.5)

  res <- valiCATE:::estimate_ncwate(Gamma, tau_hat, "BLP")
  expect_true(is.list(res))
  expect_named(res, c("estimate", "se", "ci", "test_stat", "p_value"))
  expect_true(res$se > 0)
})

test_that("NCWATE is near 1 when tau_hat = tau", {
  set.seed(42)
  n <- 2000
  tau_hat <- rnorm(n, mean = 1)
  Gamma <- tau_hat + rnorm(n, sd = 0.3)

  res <- valiCATE:::estimate_ncwate(Gamma, tau_hat, "BLP")
  expect_true(abs(res$estimate - 1) < 0.3)
})

test_that("NCWATE works for all weight types", {
  set.seed(42)
  n <- 200
  tau_hat <- rnorm(n) + 1
  Gamma <- tau_hat + rnorm(n, sd = 0.5)

  for (wt in c("AUTOC", "AUC-HVL", "BLP", "QINI")) {
    res <- valiCATE:::estimate_ncwate(Gamma, tau_hat, wt)
    expect_true(is.finite(res$estimate))
    expect_true(is.finite(res$se))
    expect_true(res$se > 0)
  }
})
