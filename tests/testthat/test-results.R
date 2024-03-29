test_that("pagfl results", {
  set.seed(1)
  alpha_0 <- matrix(c(-1, -5, .3, 1.5, -1, 1.5), nc = 2)
  sim <- sim_DGP(N = 50, n_periods = 500, p = 2, n_groups = 3, alpha_0 = alpha_0)
  y <- sim$y
  X <- sim$X
  estim <- pagfl(y = y, X = X, n_periods = 500, lambda = 5)
  expect_equal(estim$K_hat, 3)
  expect_equal(estim$groups_hat, sim$groups, ignore_attr = TRUE)
  expect_equal(estim$alpha_hat, alpha_0, ignore_attr = TRUE, tolerance = 0.05)
})

test_that("tv_pagfl results", {
  set.seed(1)
  sim <- sim_tv_DGP(N = 50, n_periods = 500, DGP = 1, sd_error = .5)
  y <- sim$y
  X <- sim$X
  estim <- tv_pagfl(y = y, X = X, n_periods = 500, lambda = 5)
  expect_equal(estim$K_hat, 3)
  expect_equal(estim$groups_hat, sim$groups, ignore_attr = TRUE)
  for (k in 1:3) {
    sim$alpha[, , k] <- sim$alpha[, , k] - mean(sim$alpha[, , k])
    estim$alpha_hat[, , k] <- estim$alpha_hat[, , k] - mean(estim$alpha_hat[, , k])
  }
  expect_equal(estim$alpha_hat, sim$alpha, ignore_attr = TRUE, tolerance = 0.5)
})
