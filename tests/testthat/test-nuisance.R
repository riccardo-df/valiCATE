test_that("compute_aipw_scores has correct formula", {
  Y <- c(1, 0, 1, 0)
  D <- c(1, 0, 1, 0)
  pscore <- c(0.5, 0.5, 0.5, 0.5)
  mu0 <- c(0.2, 0.2, 0.2, 0.2)
  mu1 <- c(0.8, 0.8, 0.8, 0.8)

  scores <- valiCATE:::compute_aipw_scores(Y, D, pscore, mu0, mu1)

  ## Manual calculation:
  ## i=1: D=1: 0.8-0.2 + 1/0.5*(1-0.8) - 0 = 0.6 + 0.4 = 1.0
  ## i=2: D=0: 0.8-0.2 + 0 - 1/0.5*(0-0.2) = 0.6 + 0.4 = 1.0
  ## i=3: D=1: 0.8-0.2 + 1/0.5*(1-0.8) - 0 = 0.6 + 0.4 = 1.0
  ## i=4: D=0: 0.8-0.2 + 0 - 1/0.5*(0-0.2) = 0.6 + 0.4 = 1.0
  expect_equal(scores, c(1.0, 1.0, 1.0, 1.0))
})

test_that("compute_aipw_scores is unbiased for E[tau|X]", {
  set.seed(42)
  n <- 5000
  X <- rnorm(n)
  pscore <- rep(0.5, n)
  D <- rbinom(n, 1, pscore)
  mu0 <- 0.5 * X
  mu1 <- 0.5 * X + X
  Y <- mu0 + D * (mu1 - mu0) + rnorm(n, sd = 0.1)

  scores <- valiCATE:::compute_aipw_scores(Y, D, pscore, mu0, mu1)
  ## E[Gamma | X] = tau(X) = X, so mean(Gamma) should be near mean(X) = 0.
  expect_true(abs(mean(scores) - mean(X)) < 0.1)
})
