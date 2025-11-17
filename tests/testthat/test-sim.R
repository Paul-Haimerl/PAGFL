test_that("pagfl sim", {
  skip_on_cran()
  N <- 10
  n_periods <- 10
  p <- 2
  K <- 3
  q <- 3
  sim <- sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = K)
  check_pagfl_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  sim_pgmm <- sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, q = q, error_spec = "AR", dynamic = T)
  check_pagfl_sim(sim = sim_pgmm, N = N, n_periods = n_periods, p = p, K = K, q = q)
  alpha_0 <- matrix(rnorm(K * p), nc = p)
  sim_2 <- sim_DGP(
    N = N, n_periods = n_periods, p = p, n_groups = K, group_proportions = c(.6, .2, .2),
    error_spec = "GARCH", alpha_0 = alpha_0
  )
  check_pagfl_sim(sim = sim_2, N = N, n_periods = n_periods, p = p, K = K, alpha_0 = alpha_0)
})

test_that("pagfl sim input", {
  skip_on_cran()
  N <- 10
  n_periods <- 10
  p <- 2
  K <- 3
  q <- 3
  # Wrong K
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = N + 1))
  # Wrong groups
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, group_proportions = c(.6, .4, .1)))
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = 2, group_proportions = c(.6, .2, .2)))
  # Incorrect error spec
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, error_spec = "a"))
  # p = 0
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = 0))
  # wrong alpha
  alpha_0 <- matrix(rnorm(4), nc = 2)
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = 3, alpha_0 = alpha_0))
  alpha_0 <- matrix(c(.1, .5, .8), nc = 1)
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = 3, alpha_0 = alpha_0))
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = 3, alpha_0 = alpha_0, dynamic = T))
  # Wrong q
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, q = p - 1))
  # To much AC
  alpha_0_ar <- matrix(rnorm(6), nc = 2)
  alpha_0_ar[1, 1] <- 1
  expect_error(sim_DGP(N = N, n_periods = n_periods, dynamic = TRUE, alpha_0 = alpha_0_ar))
})

test_that("fuse_time sim", {
  skip_on_cran()
  N <- 10
  n_periods <- 15
  p <- 2
  K <- 3
  # No intercept
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = FALSE, p = p, n_groups = K, group_proportions = c(.6, .2, .2))
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  # Intercept
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = T, p = p, n_groups = K)
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  # Dynamic, no intercept
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = F, p = p, n_groups = K, dynamic = T, sd_error = 2)
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  # Dynamic with intercept
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = TRUE, p = p, n_groups = K, dynamic = T, sd_error = 2)
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  # Dynamic with intercept and extra var
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = TRUE, p = p + 1, n_groups = K, dynamic = T, sd_error = 2)
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p + 1, K = K)
  # With specific coefficients
  locations <- matrix(stats::runif(p * K, .3, .9), ncol = K)
  scales <- matrix(stats::runif(p * K, .01, .09), ncol = K)
  d <- 3
  polynomial_coef <- array(stats::runif(p * d * K, -20, 20), dim = c(p, d, K))
  polynomial_coef <- aperm(apply(polynomial_coef, c(1, 3), function(x) x - mean(x) + 1 / d), c(2, 1, 3))
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = d, locations = locations, scales = scales, polynomial_coef = polynomial_coef)
  check_fuse_time_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
})

test_that("fuse_time sim input", {
  skip_on_cran()
  N <- 10
  n_periods <- 10
  p <- 2
  K <- 3
  # Wrong p
  expect_warning(sim_tv_DGP(N = N, n_periods = n_periods, intercept = TRUE, p = 1, dynamic = T))
  # Wrong coefficients
  locations <- matrix(stats::runif(p * K, .3, .9), ncol = K - 1)
  scales <- matrix(stats::runif(p * K, .01, .09), ncol = K - 1)
  d <- 3
  polynomial_coef <- array(stats::runif(p * d * K, -20, 20), dim = c(p, d, K))
  polynomial_coef <- aperm(apply(polynomial_coef, c(1, 3), function(x) x - mean(x) + 1 / d), c(2, 1, 3))
  expect_error(sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = 2, locations = locations))
  expect_error(sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = 2, polynomial_coef = polynomial_coef))
  expect_error(sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = 2, scales = scales))
  # AR error
  expect_no_error(sim_tv_DGP(error_spec = "AR"))
})

test_that("sim_DGP dyn_panel deprecated argument triggers warning", {
  skip_on_cran()
  expect_warning(sim <-  sim_DGP(dyn_panel = TRUE, N = 6, n_periods = 5, p = 2, n_groups = 2),
                        class = "lifecycle_warning_deprecated"
  )
  expect_type(sim, "list")
  expect_named(sim, c("alpha", "groups", "y", "X", "Z", "data"))
  expect_equal(ncol(sim$alpha), 2)
})
