#' Simulate a Time-varying Panel With a Latent Group Structure
#'
#' @description Construct a time-varying panel data set subject to a latent group structure with three distinct groups, mirroring Su et al. (2019, sec. 6.1).
#'
#' @param N the number of cross-sectional units. Default is 50.
#' @param n_periods the number of simulated time periods \eqn{T}. Default is 40.
#' @param DGP the data generating process. Options are
#' \describe{
#'  \item{1}{generates a trend only, according to DGP 1 in Su et al. (2019, sec. 6.1)}
#'  \item{2}{simulates a trend and an additional exogenous explanatory variable, as described by DGP 2 in Su et al. (2019, sec. 6.1)}
#'  \item{1}{draws a dynamic panel data model with one \eqn{AR} lag, following DGP 3 in Su et al. (2019, sec. 6.1)}
#' }
#'
#' @details
#' The scalar dependent variable \eqn{y_{it}} is driven by the following panel data model:
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{it} x_{it} + u_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect and \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables. The errors \eqn{u_{it}} feature a \eqn{iid} standard normal distribution.
#' The coefficient vector \eqn{\beta_i = \{\beta_{i1}^\prime, \dots, \beta_{iT}^\prime \}^\prime} is subject to the group pattern
#' \deqn{\beta_i = \sum_{k = 1}^3 \alpha_k \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^3 G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k \| \neq \| \alpha_j \|} for any \eqn{k \neq j}.
#'
#' \describe{
#' \item{\code{DGP = 1}}{\eqn{x_{it} = 1} and thus only includes a group-specific trend.}
#' \item{\code{DGP = 2}}{\eqn{x_{it} = (1, x_{it,2})^\prime}, where \eqn{x_{it,2}} is standard normal.}
#' \item{\code{DGP = 3}}{\eqn{x_{it} = y_{it-1}}.}
#' }
#'
#' The parameterization of the functional coefficients and group sizes follow Su et al. (2019, sec. 6.1).
#'
#' @examples
#' # Simulate a time-varying panel subject to a time trend and a latent group structure
#' sim <- sim_tv_DGP(N = 20, n_periods = 50, DGP = 1)
#' y <- sim$y
#' X <- sim$X
#' @references
#' Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. *Journal of Business & Economic Statistics*, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{alpha}}{a \eqn{T \times p \times K} array of group-specific time-varying parameters}
#' \item{\code{beta}}{a \eqn{T \times p \times N} array of individual time-varying parameters}
#' \item{\code{groups}}{a vector indicating the group memberships \eqn{(g_1, \dots, g_N)}, where \eqn{g_i = k} if \eqn{i \in} group \eqn{k}.}
#' \item{\code{y}}{a \eqn{NT \times 1} vector of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.}
#' \item{\code{X}}{a \eqn{NT \times p} matrix of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.}
#' @export
sim_tv_DGP <- function(N = 50, n_periods = 40, DGP = 1) {
  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  simChecks(dyn = TRUE, N = N, n_groups = 3, group_proportions = NULL, p = 1)
  if (!(DGP %in% 1:3)) stop("Select on of the DGPs 1 to 3\n")

  #------------------------------#
  #### Generate parameters    ####
  #------------------------------#

  # Time trends
  locations <- c(.5, .7, .6)
  scales <- c(.1, .05, .05)
  trend_mat <- trend_fctn(coef_mat = cbind(locations, scales), n_periods = n_periods)
  trend_mat[, 2] <- trend_mat[, 2] + poly_fctn(coef_mat = t(c(2, -6, 4)), n_periods = n_periods)
  trend_mat[, 3] <- trend_mat[, 3] + poly_fctn(coef_mat = t(c(4, -8, 4)), n_periods = n_periods)
  trend_mat <- 6 * trend_mat

  # Polynomial coefficients
  alpha_coef <- array(c(
    2, -4, 2, 1, -3, 2, .5, -.5, 0,
    2, -5, 2, 1, -3, 2, .5, -.5, 0
  ), dim = c(3, 3, 2))
  trend_coef <- array(c(
    .6, .7, .4, .1, .04, .07,
    .6, .2, .8, .03, .04, .07
  ), dim = c(3, 2, 2))

  if (DGP == 2) {
    p <- 2
  } else {
    p <- 1
  }

  n_groups <- 3
  alpha_array <- array(NA, dim = c(n_periods, p, n_groups))

  for (k in 1:n_groups) {
    if (DGP == 2) {
      trend_mat_star <- trend_fctn(coef_mat = t(trend_coef[k, , 1]), n_periods = n_periods)
      poly_coef <- t(alpha_coef[, k, 1])
      alpha_mat <- poly_fctn(coef_mat = poly_coef, n_periods = n_periods)
      trend_coef_mat <- alpha_mat + trend_mat_star
      alpha_array[, , k] <- cbind(.5 * trend_mat[, k], 3 * trend_coef_mat)
    } else if (DGP == 3) {
      trend_mat_star <- trend_fctn(coef_mat = t(trend_coef[k, , 2]), n_periods = n_periods)
      poly_coef <- t(alpha_coef[, k, 2])
      alpha_mat <- poly_fctn(coef_mat = poly_coef, n_periods = n_periods)
      trend_coef_mat <- -.5 * alpha_mat + trend_mat_star
      alpha_array[, , k] <- 3 / 2 * trend_coef_mat
    } else {
      alpha_array[, , k] <- trend_mat[, k]
    }
  }

  #------------------------------#
  #### Simulate the groupings ####
  #------------------------------#

  group_vec <- rep(1:n_groups, round(c(.3, .3, .4) * N))
  if (length(group_vec) != N) group_vec <- rep(group_vec, length.out = N)
  groups_raw <- sample(group_vec, N, replace = FALSE)
  groups <- match(groups_raw, unique(groups_raw))
  beta_array <- alpha_array[, , groups]
  if (is.matrix(beta_array)) {
    beta_array <- array(beta_array, dim = c(n_periods, 1, N))
  }
  beta_mat <- matrix(aperm(beta_array, c(1, 3, 2)), ncol = dim(beta_array)[2])

  #------------------------------#
  #### Construct X and y      ####
  #------------------------------#

  # Draw the cross-sectional errors
  u <- stats::rnorm(N * n_periods)
  # Draw the fixed effects
  gamma <- rep(stats::rnorm(N), each = n_periods)
  # Generate the regressors
  if (DGP != 3) {
    gamma <- 0
    if (DGP == 1) {
      X <- as.matrix(rep(1, N * n_periods))
      # Construct the observations
      y <- X * beta_mat + gamma + u
    } else if (DGP == 2) {
      X <- cbind(rep(1, N * n_periods), stats::rnorm(N * n_periods))
      # Construct the observations
      y <- rowSums(X * beta_mat) + gamma + u
    }
  } else {
    y <- c()
    X <- c()
    u <- matrix(u, ncol = N)
    for (i in 1:N) {
      group <- groups[i]
      y_i <- rep(0, n_periods + 1)
      for (t in 2:(n_periods + 1)) {
        y_i[t] <- y_i[t - 1] * alpha_array[t - 1, , group] + u[t - 1, i]
      }
      y_i <- y_i + unique(gamma)[i]
      y <- c(y, y_i[-(n_periods + 1)])
      X <- c(X, y_i[-1])
    }
    y <- as.matrix(y)
    X <- as.matrix(X)
  }

  return(list(alpha = alpha_array, beta = beta_array, groups = groups, y = y, X = X))
}


sim_tv_DGP_rand <- function(N = 50, n_periods = 40, p = 2, n_groups = 3, d = 3, group_proportions = NULL) {
  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  simChecks(dyn = TRUE, N = N, n_groups = n_groups, group_proportions = group_proportions, p = p)

  #------------------------------#
  #### Generate parameters    ####
  #------------------------------#

  # Time trends
  locations <- stats::runif(n_groups, .3, .9)
  scales <- stats::runif(n_groups, .01, .09)
  trend_mat <- trend_fctn(coef_mat = cbind(locations, scales), n_periods = n_periods)
  # Polynomial coefficients
  alpha_array <- array(NA, dim = c(n_periods, p + 1, n_groups))
  for (k in 1:n_groups) {
    locations <- stats::runif(p, .3, .9)
    scales <- stats::runif(p, .01, .09)
    trend_mat_star <- trend_fctn(coef_mat = cbind(locations, scales), n_periods = n_periods)
    poly_coef <- matrix(stats::runif(p * d, -3, 3), ncol = d)
    alpha_mat <- poly_fctn(coef_mat = poly_coef, n_periods = n_periods)
    alpha_array[, , k] <- cbind(trend_mat[, k], alpha_mat + trend_mat_star)
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
  beta_mat <- matrix(aperm(beta_array, c(1, 3, 2)), ncol = dim(beta_array)[2])

  #------------------------------#
  #### Construct X and y      ####
  #------------------------------#

  # Generate the regressors
  X <- cbind(rep(1, N * n_periods), matrix(stats::rnorm(N * p * n_periods), ncol = p))
  # Draw the cross-sectional errors
  u <- stats::rnorm(N * n_periods)
  # Construct the observations
  y <- rowSums(X * beta_mat) + u
  return(list(alpha = alpha_array, beta = beta_array, groups = groups, y = y, X = X))
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
  d <- ncol(coef_mat)
  beta_mat <- apply(coef_mat, 1, function(coefs, n_periods, d) {
    beta <- rowSums(sapply(1:d, function(x, coefs, n_periods) {
      coefs[x] * ((1:n_periods) / n_periods)^x
    }, coefs, n_periods))
  }, n_periods = n_periods, d = d)
  return(beta_mat)
}
