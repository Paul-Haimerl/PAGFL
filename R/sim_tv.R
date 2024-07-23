#' Simulate a Time-varying Panel With a Group Structure in the Slope Coefficients
#'
#' @description Construct a time-varying panel data set subject to a group structure in the slope coefficients with optional \eqn{AR(1)} innovations.
#'
#' @param N the number of cross-sectional units. Default is 50.
#' @param n_periods the number of simulated time periods \eqn{T}. Default is 40.
#' @param intercept logical. If \code{TRUE}, a time-varying intercept is generated.
#' @param p the number of simulated explanatory variables
#' @param n_groups the number of groups \eqn{K}. Default is 3.
#' @param d the polynomial degree used to construct the time-varying coefficients.
#' @param group_proportions a numeric vector of length \code{n_groups} indicating size of each group as a fraction of \eqn{N}. If \code{NULL}, all groups are of size \eqn{N / K}. Default is \code{NULL}.
#' @param error_spec options include
#' \describe{
#' \item{\code{"iid"}}{for \eqn{iid} errors.}
#' \item{\code{"AR"}}{for an \eqn{AR(1)} error process with an autoregressive coefficient of 0.5.}
#' }
#' Default is \code{"iid"}.
#' @param dynamic Logical. If \code{TRUE}, the panel includes one stationary autoregressive lag of \eqn{y_{it}} as a regressor. Default is \code{FALSE}.
#' @param locations a \eqn{p \times K} matrix of location parameters of a logistic distribution function used to construct the time-varying coefficients. If left empty, the location parameters are drawn randomly. Default is \code{NULL}.
#' @param scales a \eqn{p \times K} matrix of scale parameters of a logistic distribution function used to construct the time-varying coefficients. If left empty, the location parameters are drawn randomly. Default is \code{NULL}.
#' @param polynomial_coef a \eqn{p \times d \times K} array of coefficients for a the polynomials used to construct the time-varying coefficients. If left empty, the location parameters are drawn randomly. Default is \code{NULL}.
#' @param sd_error standard deviation of the cross-sectional errors. Default is 1.
#' @param DGP `r lifecycle::badge("deprecated")` the data generating process. Options are
#' \describe{
#'  \item{1}{generates a trend only.}
#'  \item{2}{simulates a trend and an additional exogenous explanatory variable.}
#'  \item{1}{draws a dynamic panel data model with one \eqn{AR} lag.}
#' }
#'
#' @details
#' The scalar dependent variable \eqn{y_{it}} is generated according to the following time-varying grouped panel data model
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{it} x_{it} + u_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{\gamma_i} is an individual fixed effect and \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables.
#' The coefficient vector \eqn{\beta_i = \{\beta_{i1}^\prime, \dots, \beta_{iT}^\prime \}^\prime} is subject to the group pattern
#' \deqn{\beta_i \left( \frac{t}{T} \right) = \sum_{k = 1}^K \alpha_k \left( \frac{t}{T} \right) \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\sup_{v \in [0,1]} \left( \| \alpha_k(v) - \alpha_j(v) \| \right) \neq 0} for any \eqn{k \neq j}, \eqn{k = 1, \dots, K}. The total number of groups \eqn{K} is determined by \code{n_groups}.
#'
#' The predictors are simulated as:
#' \deqn{x_{it,j} = 0.2 \gamma_i + e_{it,j}, \quad \gamma_i,e_{it,j} \sim i.i.d. N(0, 1), \quad j = \{1, \dots, p\},}
#' where \eqn{e_{it,j}} denotes a series of innovations. \eqn{\gamma_i} and \eqn{e_i} are independent of each other.
#'
#' The errors \eqn{u_{it}} feature a \eqn{iid} standard normal distribution.
#'
#' In case \code{locations = NULL}, the location parameters are drawn from \eqn{\sim U[0.3, 0.9]}.
#' In case \code{scales = NULL}, the scale parameters are drawn from \eqn{\sim U[0.01, 0.09]}.
#' In case \code{polynomial_coef = NULL}, the polynomial coefficients are drawn from \eqn{\sim U[-20, 20]} and normalized so that all coefficients of one polynomial sum up to 1.
#' The final coefficient function follows as \eqn{\alpha_k (t/T) = 3 * F(t/T, location, scale) + \sum_{j=1}^d a_j (t/T)^j}, where \eqn{F(\cdot, location, scale)} denotes a cumulative logistic distribution function and \eqn{a_j} reflects a polynomial coefficient.
#'
#' @examples
#' # Simulate a time-varying panel subject to a time trend and a group structure
#' sim <- sim_tv_DGP(N = 20, n_periods = 50, intercept = TRUE, p = 1)
#' y <- sim$y
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{alpha}}{a \eqn{T \times p \times K} array of group-specific time-varying parameters}
#' \item{\code{beta}}{a \eqn{T \times p \times N} array of individual time-varying parameters}
#' \item{\code{groups}}{a vector indicating the group memberships \eqn{(g_1, \dots, g_N)}, where \eqn{g_i = k} if \eqn{i \in} group \eqn{k}.}
#' \item{\code{y}}{a \eqn{NT \times 1} vector of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.}
#' \item{\code{X}}{a \eqn{NT \times p} matrix of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.}
#' \item{\code{data}}{a \eqn{NT \times (p + 1)} data.frame of the outcome and the explanatory variables.}
#' @export
sim_tv_DGP <- function(N = 50, n_periods = 40, intercept = TRUE, p = 1, n_groups = 3, d = 3, dynamic = FALSE, group_proportions = NULL, error_spec = "iid", locations = NULL, scales = NULL, polynomial_coef = NULL, sd_error = 1, DGP = lifecycle::deprecated()) {
  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  error_spec <- match.arg(error_spec, c("iid", "AR"))
  if (lifecycle::is_present(DGP)) {
    lifecycle::deprecate_warn("1.1.0", "sim_tv_DGP(DGP)", "sim_tv_DGP(p)")
    if (DGP == 1) {
      intercept <- TRUE
      p <- 1
    } else if (DGP == 2) {
      intercept <- TRUE
      p <- 2
    } else {
      intercept <- FALSE
      p <- 1
      dynamic <- TRUE
    }
  }

  if (intercept) {
    p_star <- p
    p <- p - 1
  } else {
    p_star <- p
  }
  if (dynamic) {
    p_star <- p_star
    p <- p - 1
  }
  p <- max(p, 0)

  simChecks(
    dyn = TRUE, N = N, n_groups = n_groups, group_proportions = group_proportions, p = p_star, locations = locations,
    scales = scales, polynomial_coef = polynomial_coef, d = d, intercept = intercept, dynamic = dynamic
  )

  #------------------------------#
  #### Generate parameters    ####
  #------------------------------#

  # Draw the functional forms
  K <- n_groups
  if (is.null(locations)) locations <- matrix(stats::runif(p_star * K, .3, .9), ncol = K)
  if (is.null(scales)) scales <- matrix(stats::runif(p_star * K, .01, .09), ncol = K)
  if (is.null(polynomial_coef)) {
    polynomial_coef <- array(stats::runif(p_star * d * K, -20, 20), dim = c(p_star, d, K))
    polynomial_coef <- aperm(apply(polynomial_coef, c(1, 3), function(x) x - mean(x) + 1 / d), c(2, 1, 3))
  }
  alpha_array <- array(NA, dim = c(n_periods, p_star, n_groups))
  for (k in 1:n_groups) {
    trend_mat_star <- trend_fctn(coef_mat = cbind(locations[, k], scales[, k]), n_periods = n_periods)
    alpha_mat <- poly_fctn(coef_mat = polynomial_coef[, , k], n_periods = n_periods)
    alpha_array[, , k] <- alpha_mat + 3 * trend_mat_star
  }

  #------------------------------#
  #### Simulate the groupings ####
  #------------------------------#

  if (is.null(group_proportions)) {
    group_vec <- rep(1:n_groups, length.out = N)
  } else {
    group_vec <- rep(1:n_groups, round(group_proportions * N))
    if (length(group_vec) != N) group_vec <- rep(group_vec, length.out = N)
  }
  groups_raw <- sample(group_vec, N, replace = FALSE)
  groups <- match(groups_raw, unique(groups_raw))
  beta_array <- alpha_array[, , groups]
  if (is.matrix(beta_array)) beta_array <- array(c(beta_array), dim = c(nrow(beta_array), 1, ncol(beta_array)))
  beta_mat <- matrix(aperm(beta_array, c(1, 3, 2)), ncol = dim(beta_array)[2])

  #------------------------------#
  #### Construct X and y      ####
  #------------------------------#

  # Draw the cross-sectional errors
  u <- stats::rnorm(N * n_periods, sd = sd_error)
  if (error_spec == "AR") {
    uList <- split(u, rep(1:((N * n_periods) %/% n_periods), each = n_periods, length.out = N * n_periods))
    u <- simAR(errorList = uList)
  }
  # Draw the fixed-effects
  gamma <- rep(stats::rnorm(N), each = n_periods)
  # Generate the regressors
  X <- matrix(stats::rnorm(N * p * n_periods), ncol = p)
  if (intercept & p > 0) {
    X <- cbind(rep(1, N * n_periods), X)
  } else if (intercept) {
    X <- as.matrix(rep(1, N * n_periods))
  }
  # Construct the observations
  if (!dynamic) {
    y <- rowSums(X * beta_mat) + u + .2 * gamma
  } else {
    y <- rep(0, n_periods * N)
    X_ar <- as.matrix(rep(0, n_periods * N))
    if (p == 0) {
      X <- X_ar
      if (intercept) X <- cbind(rep(1, N * n_periods), X)
    } else {
      X <- cbind(X_ar, X)
    }
    for (i in 1:N) {
      indx <- (i - 1) * n_periods + 1
      y[indx] <- u[indx]
      if (p > 0) {
        y[indx] <- y[indx] + sum(beta_mat[indx, -1] * X[indx, -1])
      }
      for (t in 1:(n_periods - 1)) {
        X[indx + t, 1] <- y[indx + t - 1]
        y[indx + t] <- u[indx + t] + sum(X[indx + t, ] * beta_mat[indx + t, ])
      }
    }
    y <- y + gamma
  }
  data <- data.frame(y = c(y), X)
  return(list(alpha = alpha_array, beta = beta_array, groups = groups, y = y, X = X, data = data))
}

# Logarithmic CDF as a time trend
trend_fctn <- function(coef_mat, n_periods) {
  trends <- apply(coef_mat, 1, function(coefs, n_periods) {
    1 / (1 + exp(-(((1:n_periods) / n_periods) - coefs[1]) / coefs[2]))
  }, n_periods = n_periods)
  return(trends)
}

# Polynomial coefficient functions
poly_fctn <- function(coef_mat, n_periods) {
  if (!is.matrix(coef_mat)) coef_mat <- t(coef_mat)
  d <- ncol(coef_mat)
  beta_mat <- apply(coef_mat, 1, function(coefs, n_periods, d) {
    beta <- rowSums(sapply(1:d, function(x, coefs, n_periods) {
      coefs[x] * ((1:n_periods) / n_periods)^x
    }, coefs, n_periods))
  }, n_periods = n_periods, d = d)
  return(beta_mat)
}
