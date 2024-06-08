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
  # wrong alpha
  alpha_0 <- matrix(rnorm(4), nc = 2)
  expect_error(sim_DGP(N = N, n_periods = n_periods, p = p, n_groups = 3, alpha_0 = alpha_0))
})

test_that("tv_pagfl sim", {
  skip_on_cran()
  N <- 10
  n_periods <- 15
  p <- 2
  K <- 3
  sim <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = T, p = p, n_groups = K, group_proportions = c(.6, .2, .2))
  check_tv_pagfl_sim(sim = sim, N = N, n_periods = n_periods, p = p, K = K)
  sim_2 <- sim_tv_DGP(N = N, n_periods = n_periods, intercept = F, p = p, n_groups = K, dynamic = T, sd_error = 2)
  check_tv_pagfl_sim(sim = sim_2, N = N, n_periods = n_periods, p = p, K = K)
  locations <- matrix(stats::runif(p * K, .3, .9), ncol = K)
  scales <- matrix(stats::runif(p * K, .01, .09), ncol = K)
  d <- 3
  polynomial_coef <- array(stats::runif(p * d * K, -20, 20), dim = c(p, d, K))
  polynomial_coef <- aperm(apply(polynomial_coef, c(1, 3), function(x) x - mean(x) + 1 / d), c(2, 1, 3))
  sim_3 <- sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = d, locations = locations, scales = scales, polynomial_coef = polynomial_coef)
  check_tv_pagfl_sim(sim = sim_3, N = N, n_periods = n_periods, p = p, K = K)
})

test_that("tv_pagfl sim input", {
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
  expect_error(sim_tv_DGP(N = N, n_periods = n_periods, p = p, n_groups = K, d = 2, locations = locations, scales = scales, polynomial_coef = polynomial_coef))
})
