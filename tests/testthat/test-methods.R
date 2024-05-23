test_that("S3 pagfl", {
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  colnames(X) <- c("a", "b")
  estim <- pagfl(y ~ X, n_periods = 150, lambda = 5)

  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(20, 2))
  expect_equal(colnames(coef_res), colnames(X))
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  expect_snapshot(estim,  cran = FALSE)
})


test_that("S3 tv_pagfl", {
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  estim <- tv_pagfl(y ~ X, n_periods = 100, lambda = 7)

  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(100, 2, 10))
  expect_equal(colnames(coef_res)[-1], "X1")
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  expect_snapshot(estim,  cran = FALSE)
})
