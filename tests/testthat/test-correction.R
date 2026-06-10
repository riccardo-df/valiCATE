test_that("BLP correction is mean-zero", {
  set.seed(42)
  n <- 500
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n)
  spec <- valiCATE:::get_weight_spec("BLP")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)
  C_hat <- valiCATE:::correction_cwate(Gamma, tau_hat, wt_info, spec)
  expect_true(abs(mean(C_hat)) < 1e-10)
})

test_that("CDF-based correction matches naive O(n^2) for QINI", {
  set.seed(123)
  n <- 50
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n)
  spec <- valiCATE:::get_weight_spec("QINI")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)

  ## Fast version.
  C_fast <- valiCATE:::correction_cwate(Gamma, tau_hat, wt_info, spec)

  ## Naive O(n^2) loop.
  C_naive <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      C_naive[i] <- C_naive[i] +
        Gamma[j] * wt_info$nabla[j] * (as.numeric(tau_hat[i] >= tau_hat[j]) - wt_info$v[j])
    }
    C_naive[i] <- C_naive[i] / n
  }

  expect_equal(C_fast, C_naive, tolerance = 1e-10)
})

test_that("CDF-based correction matches naive O(n^2) for AUTOC", {
  set.seed(456)
  n <- 50
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n)
  spec <- valiCATE:::get_weight_spec("AUTOC")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)

  C_fast <- valiCATE:::correction_cwate(Gamma, tau_hat, wt_info, spec)

  C_naive <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      C_naive[i] <- C_naive[i] +
        Gamma[j] * wt_info$nabla[j] * (as.numeric(tau_hat[i] >= tau_hat[j]) - wt_info$v[j])
    }
    C_naive[i] <- C_naive[i] / n
  }

  expect_equal(C_fast, C_naive, tolerance = 1e-10)
})

test_that("NCWATE correction matches naive O(n^2) for QINI", {
  set.seed(789)
  n <- 50
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n) + 1
  gamma_hat <- 0.8
  spec <- valiCATE:::get_weight_spec("QINI")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)

  C_fast <- valiCATE:::correction_ncwate(Gamma, tau_hat, gamma_hat, wt_info, spec)

  resid <- Gamma - gamma_hat * tau_hat
  C_naive <- numeric(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      C_naive[i] <- C_naive[i] +
        resid[j] * wt_info$nabla[j] * (as.numeric(tau_hat[i] >= tau_hat[j]) - wt_info$v[j])
    }
    C_naive[i] <- C_naive[i] / n
  }

  expect_equal(C_fast, C_naive, tolerance = 1e-10)
})

test_that("BLP NCWATE correction is mean-zero", {
  set.seed(42)
  n <- 500
  Gamma <- rnorm(n)
  tau_hat <- rnorm(n) + 1
  gamma_hat <- 1.0
  spec <- valiCATE:::get_weight_spec("BLP")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)
  C_gamma <- valiCATE:::correction_ncwate(Gamma, tau_hat, gamma_hat, wt_info, spec)
  expect_true(abs(mean(C_gamma)) < 1e-10)
})
