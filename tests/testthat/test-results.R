test_that("pagfl PLS output", {
  skip_on_cran()
  source(test_path("fixtures", "test_helper.R"))
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  groups_0 <- sim$groups
  alpha_0 <- sim$alpha
  y <- sim$y
  X <- sim$X
  colnames(X) <- c("a", "b")
  data <- as.data.frame(cbind(y = c(y), X))
  # With data
  estim <- pagfl(y ~ a + b, data = data, n_periods = 150, lambda = 5)
  check_pagfl_pls(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_pagfl_output(estim = estim, X = X)
  # With index
  data$i <- as.character(rep(101:120, each = 150))
  data$t <- rep(1:150, 20)
  estim <- pagfl(y ~ a + b, data = data, index = c("i", "t"), lambda = 5)
  check_pagfl_pls(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_pagfl_output(estim = estim, X = X, i_index = data$i, t_index = data$t)
  # with bias correction
  estim <- pagfl(y ~ a + b, data = data, n_periods = 150, lambda = 5, bias_correc = TRUE)
  check_pagfl_output(estim = estim, X = X)
})

test_that("pagfl PGMM output", {
  skip_on_cran()
  source(test_path("fixtures", "test_helper.R"))
  sim <- readRDS(test_path("fixtures", "pagfl_pgmm_sim.rds"))
  groups_0 <- sim$groups
  alpha_0 <- sim$alpha
  y <- sim$y
  X <- sim$X
  Z <- sim$Z
  colnames(X) <- c("a", "b")
  data <- as.data.frame(cbind(y = c(y), X))
  # With data
  estim <- pagfl(y ~ a + b, data = data, n_periods = 150,  method = "PGMM", Z = Z, lambda = 2)
  check_pagfl_pgmm(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_pagfl_output(estim = estim, X = X)
  # With index
  data$i <- as.character(rep(101:120, each = 150))
  data$t <- rep(1:150, 20)
  estim <- pagfl(y ~ a + b, data = data, index = c("i", "t"),  method = "PGMM", Z = Z, lambda = 2)
  check_pagfl_pgmm(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_pagfl_output(estim = estim, X = X, i_index = data$i, t_index = data$t)
})

test_that("tv_pagfl results", {
  skip_on_cran()
  source(test_path("fixtures", "test_helper.R"))
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  groups_0 <- sim$groups
  alpha_0 <- sim$alpha
  y <- sim$y
  X <- sim$X
  colnames(X) <- "a"
  p <- 2
  # With data
  data <- as.data.frame(cbind(y = c(y), X))
  estim <- tv_pagfl(y ~ 1 + a, data = data, n_periods = 100, lambda = 15)
  check_tv_pagfl(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_tv_pagfl_output(estim = estim, X = X)
  # With index
  data$i <- rep(1:10, each = 100)
  data$t <- rep(1:100, 10)
  estim <- tv_pagfl(y ~ 1 + a, data = data, index = c("i", "t"), lambda = 15)
  check_tv_pagfl(estim = estim, groups_0 = groups_0, alpha_0 = alpha_0)
  check_tv_pagfl_output(estim = estim, X = X, i_index = data$i, t_index = data$t)
})


test_that("tv_pagfl Unbalanced panel output", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim.rds"))
  y <- sim$y
  data <- as.data.frame(cbind(y = c(y)))
  data$i <- as.character(rep(1:20, each = 100))
  data$t <- rep(1:100, 20)
  set.seed(1)
  delete_index <- as.logical(rbinom(n = nrow(data), prob = 0.75, size = 1))
  delete_index[2] <- TRUE
  data[!delete_index, "y"] <- NA
  # Ensure at that the first period of all remaining individuals is omitted
  data[data$t %in% c(1, 100), "y"] <- NA
  # Ensure that the second period for all but one group is omitted
  data[data$t == 2 & data$i != 1, ] <- NA
  estim <- tv_pagfl(y ~ 1, data = data, index = c("i", "t"), lambda = 10, verbose = F)
  # Time periods without support at the beginning or end are omitted from the output
  expect_equal(nrow(estim$coefficients$tv), 98)
  # Time periods without support at the beginning or end for some groups are given NA
  expect_equal(sum(is.na(estim$coefficients$tv[1, , ])), 2)
})

test_that("tv_pagfl const_coef", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim_2.rds"))
  y <- sim$y
  X <- sim$X
  colnames(X) <- "a"
  data <- as.data.frame(cbind(y = c(y), X))
  estim <- tv_pagfl(y ~ 1 + a, data = data, n_periods = 100, lambda = 20, const_coef = "a")

  expect_equal(dim(estim$coefficients$const), c(3, 1))
  expect_equal(colnames(estim$coefficients$const), "a")
  expect_equal(ncol(estim$coefficients$tv), 1)
})
