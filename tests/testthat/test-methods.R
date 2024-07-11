test_that("S3 pagfl", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  colnames(X) <- c("a", "b")
  data <- data.frame(y = y, X)
  estim <- pagfl(y ~ a + b, data = data, n_periods = 150, lambda = 5)
  # print
  expect_snapshot_output(estim, cran = FALSE)
  # coef
  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(20, 2))
  expect_equal(colnames(coef_res), colnames(X))
  # resid
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  # fitted
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  # summary
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  # formula
  formula_estim <- formula(estim)
  expect_type(formula_estim, "language")
  expect_equal(all.vars(formula_estim), c("y", colnames(X)))
  # df.residual
  expect_equal(df.residual(estim), length(y) - length(y) / 150 - ncol(X) * estim$groups$n_groups)
  # summary with only one group and unbalanced data
  data$i_index <- rep(1:(length(y) / 150), each = 150)
  data$t_index <- rep(1:150, length(y) / 150)
  data <- data[-c(1, 200), ]
  estim <- pagfl(y ~ a + b, data = data, lambda = 1e3, index = c("i_index", "t_index"))
  expect_snapshot(summary(estim))
})

test_that("S3 grouped_plm", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  groups_0 <- sim$groups
  colnames(X) <- c("a", "b")
  data <- data.frame(y = y, X)
  estim <- grouped_plm(y ~ a + b, data = data, groups = groups_0, n_periods = 150)
  # print
  expect_snapshot_output(estim, cran = FALSE)
  # coef
  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(20, 2))
  expect_equal(colnames(coef_res), colnames(X))
  # resid
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  # fitted
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  # summary
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  # formula
  formula_estim <- formula(estim)
  expect_type(formula_estim, "language")
  expect_equal(all.vars(formula_estim), c("y", colnames(X)))
  # df.residual
  expect_equal(df.residual(estim), length(y) - length(y) / 150 - ncol(X) * estim$groups$n_groups)
  # summary with only one group and unbalanced data
  data$i_index <- rep(1:(length(y) / 150), each = 150)
  data$t_index <- rep(1:150, length(y) / 150)
  data <- data[-c(1, 200), ]
  estim <- grouped_plm(y ~ a + b, data = data, groups = rep(1, length(groups_0)), index = c("i_index", "t_index"))
  expect_snapshot(summary(estim))
})

test_that("S3 tv_pagfl", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  data <- data.frame(y = y, X1 = X)
  data$a <- stats::rnorm(length(y))
  estim <- tv_pagfl(y ~ X1, data = data, n_periods = 100, lambda = 7)
  estim_const <- tv_pagfl(y ~ X1 + a, data = data, n_periods = 100, lambda = 7, const_coef = "a")
  # print
  expect_snapshot(estim,  cran = FALSE)
  # coef
  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(100, 2, 10))
  expect_equal(colnames(coef_res)[-1], "X1")
  coef_res_const <- coef(estim_const)
  # resid
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  # fitted
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  # summary
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  # formula
  formula_estim <- formula(estim)
  expect_type(formula_estim, "language")
  expect_equal(all.vars(formula_estim), c("y", "X1"))
  # df.residual
  expect_equal(df.residual(estim), length(y) - length(y) / 100 - 2 * (estim$args$M + estim$args$d + 1) * estim$groups$n_groups)
})

test_that("S3 grouped_tv_plm", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  groups_0 <- sim$groups
  data <- data.frame(y = y, X1 = X)
  data$a <- stats::rnorm(length(y))
  estim <- grouped_tv_plm(y ~ X1, data = data, groups = groups_0, n_periods = 100)
  estim_const <- grouped_tv_plm(y ~ X1 + a, data = data, groups = groups_0, n_periods = 100, const_coef = "a")
  # print
  expect_snapshot(estim,  cran = FALSE)
  # coef
  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(100, 2, 10))
  expect_equal(colnames(coef_res)[-1], "X1")
  coef_res_const <- coef(estim_const)
  # resid
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(length(y), 3))
  # fitted
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(length(y), 3))
  # summary
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  # formula
  formula_estim <- formula(estim)
  expect_type(formula_estim, "language")
  expect_equal(all.vars(formula_estim), c("y", "X1"))
  # df.residual
  expect_equal(df.residual(estim), length(y) - length(y) / 100 - 2 * (estim$args$M + estim$args$d + 1) * estim$groups$n_groups)
})

test_that("S3 tv_pagfl const coef unbalanced", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  df <- data.frame(y = y, X)
  set.seed(1)
  df$a <- rnorm(length(y))
  df$i_index <- rep(1:(length(y) / 100), each = 100)
  df$t_index <- rep(1:100, length(y) / 100)
  rm_index <- as.logical(rbinom(n = length(y), prob = 0.7, size = 1))
  df <- df[rm_index, ]
  estim <- tv_pagfl(y ~ X + a, const_coef = "a", data = df, index = c("i_index", "t_index"), lambda = 25)
  # coef
  coef_res <- coef(estim)
  expect_equal(dim(coef_res), c(100, 3, 10))
  expect_equal(colnames(coef_res), c("(Intercept)", "X", "a"))
  # resid
  resid_res <- resid(estim)
  expect_equal(dim(resid_res), c(nrow(df), 3))
  # fitted
  fitted_res <- fitted(estim)
  expect_equal(dim(fitted_res), c(nrow(df), 3))
  # summary
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  expect_snapshot(estim,  cran = FALSE)
  # formula
  formula_estim <- formula(estim)
  expect_type(formula_estim, "language")
  expect_equal(all.vars(formula_estim), c("y", "X", "a"))
  # df.residual
  expect_equal(df.residual(estim), nrow(df) - max(df$i_index) - (1 + 2 * (estim$args$M + estim$args$d + 1)) * estim$groups$n_groups)
})

test_that("S3 grouped_tv_plm const coef unbalanced summary", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  groups_0 <- sim$groups
  df <- data.frame(y = y, X)
  set.seed(1)
  df$a <- rnorm(length(y))
  df$i_index <- rep(1:(length(y) / 100), each = 100)
  df$t_index <- rep(1:100, length(y) / 100)
  rm_index <- as.logical(rbinom(n = length(y), prob = 0.7, size = 1))
  df <- df[rm_index, ]
  # vanilla
  estim <- grouped_tv_plm(y ~ X + a, const_coef = "a", data = df, groups = groups_0, index = c("i_index", "t_index"))
  summary_estim <- summary(estim)
  expect_snapshot_output(summary_estim, cran = FALSE)
  expect_snapshot(estim,  cran = FALSE)
  # only one group
  estim <- grouped_tv_plm(y ~ ., data = df, groups = rep(1, length(groups_0)), index = c("i_index", "t_index"))
  expect_no_error(summary(estim))
  expect_no_error(coef(estim))
})
