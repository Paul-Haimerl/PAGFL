test_that("pagfl inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  X <- sim$X
  colnames(X) <- c("a", "b")
  data <- as.data.frame(cbind(y = c(y), X))
  # Wrong number of time periods
  expect_error(pagfl(y ~ ., data = data, n_periods = 151, lambda = 1, verbose = F))
  # No data
  expect_error(pagfl(y ~ ., n_periods = 151, lambda = 1, verbose = F))
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
  # Index as y
  expect_error(pagfl(i ~ a + b, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # y as part of x
  expect_error(pagfl(y ~ y + b, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # Nonexistent regressor
  expect_error(pagfl(y ~ a + c, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # PGMM but no instruments
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "PGMM", lambda = 1, verbose = F))
  # To few instruments
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "PGMM", Z = as.matrix(X[, -1]), lambda = 1, verbose = F))
  # Plain PGMM
  expect_no_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "PGMM", lambda = 1, Z = X, verbose = F))
  # Plain PGMM and bias correction
  expect_no_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "PGMM", lambda = 1, Z = X, bias_correc = T, verbose = F))
  # No method
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), method = "A", lambda = 1, verbose = F))
  # Incorrect argument
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), kappa = -1, lambda = 1, verbose = F))
  expect_error(pagfl(y ~ a + b, data = data, index = c("i", "t"), rho = -1, lambda = 1, verbose = F))
  # No index or n_periods
  expect_error(pagfl(y ~ ., data = data, lambda = 1, verbose = F))
  # Both index and n_periods
  expect_warning(pagfl(y ~ ., data = data, index = c("i", "t"), n_periods = 150, lambda = 1, verbose = TRUE))
  # No dependent variable
  expect_error(pagfl(~., data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # Large min_group_frac
  expect_warning(pagfl(y ~ ., data = data, index = c("i", "t"), lambda = 1, verbose = TRUE, min_group_frac = .5))
  # 0 ming_group_frac
  expect_no_error(pagfl(y ~ ., data = data, index = c("i", "t"), lambda = 1, verbose = F, min_group_frac = NULL))
  # Negative ming_group_frac
  expect_error(pagfl(y ~ ., data = data, index = c("i", "t"), lambda = 1, verbose = F, min_group_frac = -.5))
  # Force some trivial groups
  sim_smallNk <- readRDS(test_path("fixtures", "pagfl_pls_sim_smallNk.rds"))
  expect_no_error(pagfl(y ~ ., data = sim_smallNk$data, n_periods = 75, lambda = 2, min_group_frac = .3, verbose = F))
  # Intercept
  data_star_2 <- data
  data_star_2$c <- 1
  expect_error(pagfl(y ~ c, data = data_star_2, index = c("i", "t"), lambda = 1, verbose = F))
  # Force only one group
  expect_no_error(pagfl(y ~ ., data = data, index = c("i", "t"), lambda = 1e6, verbose = F))
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
  expect_error(pagfl(y ~ V2 + V3, data = data, n_periods = 150, lambda = 1, verbose = F))
})

test_that("grouped_plm inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "pagfl_pls_sim.rds"))
  y <- sim$y
  data <- as.data.frame(cbind(y = c(y)))
  groups_0 <- sim$groups
  # Wrong number of time periods
  expect_error(grouped_plm(y ~ 1, data = data, groups = groups_0, n_periods = 151, verbose = F))
  # Wrong group vector
  groups_0_star <- c(groups_0, 1)
  expect_error(grouped_plm(y ~ 1, data = data, groups = groups_0_star, n_periods = 150, verbose = F))
})

test_that("fuse_time inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "fuse_time_sim.rds"))
  y <- sim$y
  data <- as.data.frame(cbind(y = c(y)))
  data$i <- as.character(rep(1:20, each = 100))
  data$t <- rep(1:100, 20)
  # Wrong number of time periods
  expect_error(fuse_time(y ~ 1, data = data, n_periods = 101, lambda = 1, verbose = F))
  # Char matrix for y
  data_star <- data.frame(y = as.character(c(y)))
  expect_error(fuse_time(y ~ 1, data = data_star, n_periods = 100, lambda = 1, verbose = F))
  # Wrong index variables
  data$a <- stats::rnorm(length(y))
  expect_error(fuse_time(y ~ 1 + a, data = data, index = c("a", "t"), lambda = 1, verbose = F))
  expect_error(fuse_time(y ~ 1, data = data, index = c("c", "t"), lambda = 1, verbose = F))
  # Nonexistent regressor
  expect_error(fuse_time(y ~ 1 + b, data = data, index = c("i", "t"), lambda = 1, verbose = F))
  # Wrong argument
  expect_error(fuse_time(y ~ 1, data = data, index = c("i", "t"), lambda = 1, verbose = F, d = -1))
  # Const_coef not in data
  expect_error(fuse_time(y ~ 1, data = data, index = c("i", "t"), lambda = 1, verbose = F, const_coef = "a"))
  # Force only one group with const coef
  expect_no_error(fuse_time(y ~ 1 + a, data = data, index = c("i", "t"), lambda = 1e4, verbose = F, const_coef = "a"))
})

test_that("tv_grouped_plm inputs", {
  skip_on_cran()
  sim <- readRDS(test_path("fixtures", "fuse_time_sim.rds"))
  y <- sim$y
  data <- as.data.frame(cbind(y = c(y)))
  groups_0 <- sim$groups
  # Wrong number of time periods
  expect_error(tv_grouped_plm(y ~ 1, data = data, groups = groups_0, n_periods = 101, verbose = F))
  # Wrong group vector
  groups_0_star <- c(groups_0, 1)
  expect_error(tv_grouped_plm(y ~ 1, data = data, groups = groups_0_star, n_periods = 100, verbose = F))
})
