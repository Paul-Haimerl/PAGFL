test_that("Incorrect n_periods and index", {
  N <- 5
  n_periods <- 10
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, DGP = 2)
  y <- sim$y
  X <- sim$X[, -1]
  # Wrong number of time periods
  expect_error(pagfl(y = y, X = X, n_periods = 11, lambda = 1, verbose = F))
  expect_error(tv_pagfl(y = y, X = X, n_periods = 11, lambda = 10, verbose = F))
  i_index <- rep(1:N, each = n_periods)
  t_index <- rep(1:n_periods, N)
  y <- cbind(y, i_index = i_index, t_index = t_index)
  X <- cbind(X, i_index = i_index, t_index = t_index)
  # Wrong index variables
  expect_error(pagfl(y = y, X = X, index = c("i_index", "a"), lambda = 1, verbose = F))
  expect_error(tv_pagfl(y = y, X = X, index = c("i_index", "a"), lambda = 10, verbose = F))
  # Indices do not match
  y <- cbind(sim$y, i_index = rep((N + 1):(2 * N), each = n_periods), t_index = t_index)
  expect_error(pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 1, verbose = F))
  expect_error(tv_pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 10, verbose = F))
})

test_that("Row and colnames of the output", {
  N <- 5
  n_periods <- 10
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, DGP = 2)
  y <- sim$y
  X <- as.matrix(sim$X[, -1])
  i_index <- rep(1:N, each = n_periods)
  t_index <- rep(1:n_periods, N)
  y <- cbind(y, i_index = i_index, t_index = t_index)
  colnames(X) <- "a"
  X <- cbind(X, i_index = i_index, t_index = t_index)
  estim <- tv_pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 15, verbose = F)
  estim_const <- pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 1, verbose = F)
  expect_equal(colnames(estim$alpha_hat), colnames(X)[1])
  expect_equal(colnames(estim_const$alpha_hat), colnames(X)[1])
  expect_equal(rownames(estim$alpha_hat), as.character(unique(t_index)))
  expect_equal(names(estim$groups_hat), as.character(unique(i_index)))
  expect_equal(names(estim_const$groups_hat), as.character(unique(i_index)))
})

test_that("Output dimensions with missings", {
  N <- 5
  n_periods <- 20
  delete_index <- as.logical(rbinom(n = N * n_periods, prob = 0.8, size = 1))
  # Ensure the first and last period of at least one individual is not discarded
  delete_index[1] <- TRUE
  delete_index[n_periods] <- TRUE
  # Ensure at that the first period of all remaining individuals is omitted
  delete_index[1 + 1:(N - 1) * n_periods] <- FALSE
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, DGP = 1)
  y <- sim$y
  X <- sim$X
  i_index <- rep(1:N, each = n_periods)
  t_index <- rep(1:n_periods, N)
  y <- cbind(y, i_index = i_index, t_index = t_index)[delete_index, ]
  X <- cbind(X, i_index = i_index, t_index = t_index)
  # Panel is unbalanced
  estim <- tv_pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 15, verbose = F)
  # Correct nrow
  expect_equal(nrow(estim$alpha_hat), max(t_index))
  # Time periods without support are omitted from the output
  expect_true(any(is.na(estim$alpha_hat[1, , ])))
})
