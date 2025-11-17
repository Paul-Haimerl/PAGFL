check_pagfl_pls <- function(estim, groups_0, alpha_0) {
  expect_equal(estim$groups$groups, groups_0, ignore_attr = TRUE)
  expect_equal(estim$groups$n_groups, max(groups_0), ignore_attr = TRUE)
  expect_equal(estim$coefficients, alpha_0, ignore_attr = TRUE, tolerance = 0.05)
  expect_equal(round(estim$IC$IC, 6), 1.045192)
  expect_equal(round(estim$IC$msr, 6), 0.983799)
  resid_0 <- readRDS(test_path("fixtures", "pagfl_pls_resid.rds"))
  expect_equal(estim$residuals, resid_0)
}

check_pagfl_pgmm <- function(estim, groups_0, alpha_0) {
  expect_equal(estim$groups$groups, groups_0, ignore_attr = TRUE)
  expect_equal(estim$groups$n_groups, max(groups_0), ignore_attr = TRUE)
  expect_equal(estim$coefficients, alpha_0, ignore_attr = TRUE, tolerance = 0.05)
  expect_equal(round(estim$IC$IC, 6), 2.025255)
  expect_equal(round(estim$IC$msr, 6), 1.963861)
  resid_0 <- readRDS(test_path("fixtures", "pagfl_pgmm_resid.rds"))
  expect_equal(estim$residuals, resid_0)
}

check_pagfl_output <- function(estim, X, i_index = NULL, t_index = NULL, oracle = FALSE) {
  alpha_hat <- estim$coefficients
  expect_equal(colnames(alpha_hat), colnames(X))
  expect_equal(rownames(alpha_hat), paste("Group", 1:estim$groups$n_groups))
  if (!oracle) {
    expect_length(estim$convergence, 2)
    expect_length(estim$IC, 3)
    expect_length(estim$args, 10)
  } else {
    expect_length(estim$IC, 2)
    expect_length(estim$args, 5)
  }
  if (!is.null(i_index)) {
    all.equal(estim$args$labs$i, i_index)
    all.equal(estim$args$labs$t, t_index)
    expect_equal(names(estim$groups$groups), unique(i_index))
  } else {
    expect_equal(names(estim$groups$groups), as.character(1:20))
  }
}

check_fuse_time <- function(estim, groups_0, alpha_0) {
  expect_equal(estim$groups$groups, groups_0, ignore_attr = TRUE)
  expect_equal(estim$groups$n_groups, max(groups_0), ignore_attr = TRUE)
  expect_equal(estim$coefficients$tv[20:80, 2, ], alpha_0[20:80, 2, ], ignore_attr = TRUE, tolerance = 0.1)
  expect_equal(round(estim$IC$IC, 6), 0.313991)
  expect_equal(round(estim$IC$msr, 6), 1.053228)
  resid_0 <- readRDS(test_path("fixtures", "fuse_time_resid.rds"))
  expect_equal(estim$residuals, resid_0)
}


check_fuse_time_output <- function(estim, X, i_index = NULL, t_index = NULL, oracle = FALSE) {
  alpha_hat <- estim$coefficients$tv
  expect_equal(colnames(alpha_hat)[-1], colnames(X))
  expect_equal(dimnames(alpha_hat)[[3]], paste("Group", 1:estim$groups$n_groups))
  if (!oracle) {
    expect_length(estim$convergence, 2)
    expect_length(estim$IC, 3)
    expect_length(estim$args, 10)
  } else {
    expect_length(estim$IC, 2)
    expect_length(estim$args, 5)
  }

  if (!is.null(i_index)) {
    all.equal(estim$args$labs$i, i_index)
    all.equal(estim$args$labs$t, t_index)
    expect_equal(names(estim$groups$groups), as.character(unique(i_index)))
  } else {
    expect_equal(names(estim$groups$groups), as.character(1:10))
  }
}

check_pagfl_sim <- function(sim, N, n_periods, p, q = NULL, alpha_0 = NULL, K) {
  expect_length(sim$y, N * n_periods)
  expect_equal(dim(sim$X), c(N * n_periods, p))
  expect_equal(max(sim$groups), K)
  expect_length(sim$groups, N)
  if (is.null(alpha_0)) {
    expect_equal(dim(sim$alpha), c(K, p))
  } else {
    expect_equal(sim$alpha, alpha_0, ignore_attr = TRUE)
  }
  if (is.null(q)) {
    expect_null(sim$Z)
  } else {
    expect_equal(dim(sim$Z), c(N * n_periods, q))
  }
}

check_fuse_time_sim <- function(sim, N, n_periods, p, K) {
  expect_length(sim$y, N * n_periods)
  expect_equal(dim(sim$X), c(N * n_periods, p))
  expect_equal(max(sim$groups), K)
  expect_length(sim$groups, N)
  expect_equal(dim(sim$alpha), c(n_periods, p, K))
  expect_equal(dim(sim$beta), c(n_periods, p, N))
}
