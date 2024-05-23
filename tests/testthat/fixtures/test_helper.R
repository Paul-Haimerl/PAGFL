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

check_pagfl_output <- function(estim, i_index = NULL, t_index = NULL) {
  alpha_hat <- estim$coefficients
  expect_equal(colnames(alpha_hat), colnames(X))
  expect_equal(rownames(alpha_hat), paste("Group", 1:estim$groups$n_groups))
  expect_length(estim$convergence, 2)
  expect_length(estim$IC, 3)
  expect_length(estim$args, 9)

  if (!is.null(i_index)) {
    all.equal(estim$args$labs$i, i_index)
    all.equal(estim$args$labs$t, t_index)
    expect_equal(names(estim$groups$groups), unique(i_index))
  } else {
    expect_equal(names(estim$groups$groups), as.character(1:20))
  }
}

check_tv_pagfl <- function(estim, groups_0, alpha_0) {
  expect_equal(estim$groups$groups, groups_0, ignore_attr = TRUE)
  expect_equal(estim$groups$n_groups, max(groups_0), ignore_attr = TRUE)
  expect_equal(estim$coefficients$tv[20:80, 2, ], alpha_0[20:80, 2, ], ignore_attr = TRUE, tolerance = 0.1)
  expect_equal(round(estim$IC$IC, 6), 1.511957)
  expect_equal(round(estim$IC$msr, 6), 1.053228)
  resid_0 <- readRDS(test_path("fixtures", "tv_pagfl_resid.rds"))
  expect_equal(estim$residuals, resid_0)
}

check_pagfl_output <- function(estim, X, i_index = NULL, t_index = NULL) {
  alpha_hat <- estim$coefficients
  expect_equal(colnames(alpha_hat), colnames(X))
  expect_equal(rownames(alpha_hat), paste("Group", 1:estim$groups$n_groups))
  expect_length(estim$convergence, 2)
  expect_length(estim$IC, 3)
  expect_length(estim$args, 9)

  if (!is.null(i_index)) {
    all.equal(estim$args$labs$i, i_index)
    all.equal(estim$args$labs$t, t_index)
    expect_equal(names(estim$groups$groups), unique(i_index))
  } else {
    expect_equal(names(estim$groups$groups), as.character(1:20))
  }
}


check_tv_pagfl_output <- function(estim, X, i_index = NULL, t_index = NULL) {
  alpha_hat <- estim$coefficients$tv
  expect_equal(colnames(alpha_hat)[-1], colnames(X))
  expect_equal(dimnames(alpha_hat)[[3]], paste("Group", 1:estim$groups$n_groups))
  expect_length(estim$convergence, 2)
  expect_length(estim$IC, 3)
  expect_length(estim$args, 9)

  if (!is.null(i_index)) {
    all.equal(estim$args$labs$i, i_index)
    all.equal(estim$args$labs$t, t_index)
    expect_equal(names(estim$groups$groups), as.character(unique(i_index)))
  } else {
    expect_equal(names(estim$groups$groups), as.character(1:10))
  }
}
