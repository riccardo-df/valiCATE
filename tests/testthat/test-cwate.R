test_that("CWATE returns correct structure", {
  set.seed(42)
  n <- 200
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n)

  res <- valiCATE:::estimate_cwate(Gamma, tau_hat, "BLP")
  expect_true(is.list(res))
  expect_named(res, c("estimate", "se", "ci", "test_stat", "p_value"))
  expect_true(res$se > 0)
  expect_equal(length(res$ci), 2)
})

test_that("CWATE with null DGP has p-value > 0.05 on average", {
  set.seed(100)
  pvals <- numeric(100)
  for (r in seq_len(100)) {
    n <- 300
    Gamma <- rnorm(n)
    tau_hat <- rnorm(n)
    res <- valiCATE:::estimate_cwate(Gamma, tau_hat, "BLP")
    pvals[r] <- res$p_value
  }
  ## Under the null, rejection rate should be near alpha = 0.05.
  reject_rate <- mean(pvals < 0.05)
  expect_true(reject_rate < 0.15)
})

test_that("CWATE detects positive signal", {
  set.seed(42)
  n <- 1000
  tau_hat <- rnorm(n)
  tau_true <- tau_hat + 0.1 * rnorm(n)
  Gamma <- tau_true + rnorm(n, sd = 0.5)

  res <- valiCATE:::estimate_cwate(Gamma, tau_hat, "BLP")
  expect_true(res$estimate > 0)
})

test_that("CWATE works for all weight types", {
  set.seed(42)
  n <- 200
  Gamma <- rnorm(n, mean = 1)
  tau_hat <- rnorm(n) + 1

  for (wt in c("AUTOC", "AUC-HVL", "BLP", "QINI")) {
    res <- valiCATE:::estimate_cwate(Gamma, tau_hat, wt)
    expect_true(is.finite(res$estimate))
    expect_true(is.finite(res$se))
    expect_true(res$se > 0)
  }
})
