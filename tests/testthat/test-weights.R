test_that("get_weight_spec returns correct structures for all weights", {
  for (wt in c("AUTOC", "AUC-HVL", "BLP", "QINI")) {
    spec <- valiCATE:::get_weight_spec(wt)
    expect_true(is.list(spec))
    expect_equal(spec$name, wt)
    expect_true(is.function(spec$h))
    expect_true(is.function(spec$nabla_v_h))
    expect_true(is.logical(spec$cdf_based))
  }
})

test_that("get_weight_spec errors on unknown weight", {
  expect_error(valiCATE:::get_weight_spec("INVALID"), "Unknown weight")
})

test_that("AUTOC weight values are correct", {
  spec <- valiCATE:::get_weight_spec("AUTOC")
  expect_equal(spec$h(1, 0.5), -log(0.5) - 1)
  expect_equal(spec$nabla_v_h(1, 0.5), 2)
  expect_true(spec$cdf_based)
})

test_that("AUC-HVL weight values are correct", {
  spec <- valiCATE:::get_weight_spec("AUC-HVL")
  expect_equal(spec$h(1, 0.5), 0)
  expect_equal(spec$nabla_v_h(1, 0.5), 4)
  expect_true(spec$cdf_based)
})

test_that("BLP weight values are correct", {
  spec <- valiCATE:::get_weight_spec("BLP")
  expect_equal(spec$h(3, 2), 1)
  expect_equal(spec$nabla_v_h(1, 1), -1)
  expect_false(spec$cdf_based)
})

test_that("QINI weight values are correct", {
  spec <- valiCATE:::get_weight_spec("QINI")
  expect_equal(spec$h(1, 0.5), 0)
  expect_equal(spec$nabla_v_h(1, 1), 1)
  expect_true(spec$cdf_based)
})

test_that("compute_weights produces centered weights for CDF-based", {
  set.seed(42)
  tau_hat <- rnorm(200)
  for (wt in c("AUTOC", "AUC-HVL", "QINI")) {
    spec <- valiCATE:::get_weight_spec(wt)
    wt_info <- valiCATE:::compute_weights(tau_hat, spec)
    expect_equal(length(wt_info$omega), 200)
    expect_equal(length(wt_info$v), 200)
    expect_equal(length(wt_info$nabla), 200)
  }
})

test_that("BLP weights are tau_hat - mean(tau_hat)", {
  set.seed(42)
  tau_hat <- rnorm(100)
  spec <- valiCATE:::get_weight_spec("BLP")
  wt_info <- valiCATE:::compute_weights(tau_hat, spec)
  expect_equal(wt_info$omega, tau_hat - mean(tau_hat))
  expect_true(abs(mean(wt_info$omega)) < 1e-12)
})

test_that("ecdf_safe returns values in (0,1)", {
  x <- c(1, 2, 3, 4, 5)
  v <- valiCATE:::ecdf_safe(x)
  expect_true(all(v > 0 & v < 1))
  expect_equal(v, (1:5 - 0.5) / 5)
})
