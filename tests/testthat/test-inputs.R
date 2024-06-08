test_that("pagfl inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  colnames(X) <- c("a", "b")
  data <- as.data.frame(cbind(y = c(y), X))
  # Wrong number of time periods
  expect_error(pagfl(y ~ ., data = data, n_periods = 151, lambda = 1, verbose = F))
  # Char matrix for y
  data_star <- as.data.frame(cbind(y = as.character(c(y)), X))
  expect_error(pagfl(y_star ~ ., data = data_star, n_periods = 150, lambda = 1, verbose = F))
  # Wrong index variables
  data$i <- as.character(rep(1:20, each = 150))
  data$t <- rep(1:150, 20)
  expect_error(pagfl(y ~ a + b, data = data, index = c("a", "t"), lambda = 1, verbose = F))
  expect_error(pagfl(y ~ i, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  expect_error(pagfl(y ~ a + b, data = data, index = c("c", "t"), lambda = 1, verbose = F))
  expect_error(pagfl(y ~ a + b, data = data, index = c("i"), lambda = 1, verbose = F))
  # Nonexistent regressor
  expect_error(pagfl(y ~ a + c, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # PGMM but no instruments
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "PGMM", lambda = 1, verbose = F))
  # PLS but instrument
  expect_warning(pagfl(y ~ a + b, data = data, n_periods = 150, method = "PLS", Z = X, lambda = 1, verbose = TRUE))
  # No method
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "A", lambda = 1, verbose = F))
  # Incorrect argument
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), kappa = -1, lambda = 1, verbose = F))
  # No index or n_periods
  expect_error(pagfl(y ~ ., data = data, lambda = 1, verbose = F))
  # No dependent variable
  expect_error(pagfl(~., data = data, lambda = 1, verbose = F))
  # Intercept
  expect_error(pagfl(y ~ as.matrix(rep(1, length(y))), lambda = 1, verbose = F))
})

test_that("Unbalanced panel pagfl", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  data <- as.data.frame(cbind(y = c(y), X))
  data$i <- as.character(rep(1:20, each = 150))
  data$t <- rep(1:150, 20)
  set.seed(1)
  delete_index <- as.logical(rbinom(n = nrow(data), prob = 0.75, size = 1))
  data[delete_index, "y"] <- NA
  expect_no_error(pagfl(y ~ V2 + V3, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  expect_error(pagfl(y[delete_index] ~ X[delete_index, ], n_periods = 150, lambda = 1, verbose = F))
})


test_that("tv_pagfl inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "tv_pagfl_sim.rds"))
  y <- sim$y
  data <- as.data.frame(cbind(y = c(y)))
  data$i <- as.character(rep(1:20, each = 100))
  data$t <- rep(1:100, 20)
  # Wrong number of time periods
  expect_error(tv_pagfl(y ~ 1, data = data, n_periods = 101, lambda = 1, verbose = F))
  # Char matrix for y
  data_star <- data.frame(y = as.character(c(y)))
  expect_error(tv_pagfl(y ~ 1, data = data_star, n_periods = 100, lambda = 1, verbose = F))
  # Wrong index variables
  data$i <- as.character(rep(1:20, each = 100))
  data$t <- rep(1:100, 20)
  data$a <- stats::rnorm(length(y))
  expect_error(tv_pagfl(y ~ 1 + a, data = data, index = c("a", "t"), lambda = 1, verbose = F))
  expect_error(tv_pagfl(y ~ 1, data = data, index = c("c", "t"), lambda = 1, verbose = F))
  # Nonexistent regressor
  expect_error(tv_pagfl(y ~ 1 + b, data = data, index = c("i", "t"), lambda = 1, verbose = F))
})
