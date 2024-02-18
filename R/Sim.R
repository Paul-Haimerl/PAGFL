#' Simulate a Panel With a Latent Group Structure
#'
#' @description Construct a static or dynamic, exogenous or endogenous panel data set subject to a latent group structure with optional \eqn{AR(1)} or \eqn{GARCH(1,1)} innovations.
#'
#' @param N the number of cross-sectional units. Default is 50.
#' @param n_periods the number of simulated time periods \eqn{T}. Default is 40.
#' @param p the number of explanatory variables. Default is 2.
#' @param n_groups the number of latent groups \eqn{K}. Default is 3.
#' @param group_proportions a numeric vector of length \code{n_groups} indicating the fraction of \eqn{N} each group will contain. If \code{NULL}, all groups are of size \eqn{\frac{N}{K}}. Default is \code{NULL}.
#' @param error_spec the error specification used. Options are
#' \describe{
#' \item{\code{NULL}}{for \eqn{iid} errors.}
#' \item{\code{'AR'}}{for an \eqn{AR(1)} error process with an autoregressive coefficient of 0.5.}
#' \item{\code{'GARCH'}}{for a \eqn{GARCH(1,1)} error process with a 0.05 constant, a 0.05 ARCH and a 0.9 GARCH coefficient.}
#' }
#' Default is \code{NULL}.
#' @param dyn_panel Logical. If \code{TRUE}, the panel includes one stationary autoregressive lag of the dependent variable (see sec. Details for information on the \eqn{AR} coefficient). Default is \code{FALSE}.
#' @param q the number of exogenous instruments when a panel with endogenous regressors is to be simulated. If panel data set with exogenous regressors is supposed to be generated, pass \code{NULL}. Default is \code{NULL}.
#' @param alpha_0 an optional pre-specified \eqn{K \times p} parameter matrix. If \code{NULL}, the coefficients are drawn randomly (see sec. Details). If \code{dyn_panel = TRUE}, the first column represents the stationary \eqn{AR} coefficient. Default is \code{NULL}.
#'
#' @details
#' The scalar dependent variable \eqn{y_{it}} is driven by the following panel data model
#' \deqn{y_{it} = \eta_i + \beta_i^\prime x_{it} + u_{it}, \quad i = \{1, \dots, N\}, \quad t = \{1, \dots, T\}.}
#' \eqn{\eta_i} represents individual fixed effects and \eqn{x_{it} = (x_{it,1}, \dots, x_{it,p})} a \eqn{p \times 1} vector of regressors.
#' The individual slope coefficient vectors \eqn{\beta_i} are subject to a latent group structure \eqn{\beta_i = \sum_{k = 1}^K \alpha_k \bold{1} \{i \in G_k\}}.
#' As a consequence, the group-level coefficients \eqn{\bold{\alpha} = (\alpha^\prime_1, \dots, \alpha^\prime_K)^\prime} follow the partition \eqn{\bold{G}} of \eqn{N} cross-sectional units \eqn{\bold{G} = (G_1, \dots, G_K)} such that \eqn{\cup_{k=1}^K = \{1,\dots,N\}} and \eqn{G_k \cap G_l = \emptyset, \; \alpha_k \neq \alpha_l} for any two groups \eqn{k \neq l} (Mehrabani, 2023, sec. 2.1).
#'
#' If a panel data set with exogenous regressors is generated (set \code{q = NULL}), the \eqn{p} predictors are simulated as:
#' \deqn{x_{it,j} = 0.2 \eta_i + e_{it,j}, \quad \eta_i,e_{it,j} \sim i.i.d. N(0, 1), \quad j = \{1, \dots, p\},}
#' where \eqn{e_{it,j}} denotes a series of innovations. \eqn{\eta_i} and \eqn{e_i} are independent of each other.
#'
#' In case \code{alpha_0 = NULL}, the group-level slope parameters \eqn{\alpha_{k}} are drawn from \eqn{\sim U[-2, 2]}.
#'
#' If a dynamic panel is specified (\code{dyn_panel = TRUE}), the \eqn{AR} coefficients \eqn{\beta^{\text{AR}}_i} are drawn from a uniform distribution with support \eqn{(-1, 1)} and \eqn{x_{it,j} = e_{it,j}}.
#' The individual fixed effects enter the dependent variable via \eqn{(1 - \beta^{\text{AR}}_i) \eta_i} to account for the autoregressive dependency.
#' I refer to Mehrabani (2023, sec 6) for details.
#'
#' When specifying an endogenous panel (set \code{q} to \eqn{q \geq p}), \eqn{e_{it,j}} correlate with the cross-sectional innovations \eqn{u_{it}} by a magnitude of 0.5 to produce endogenous regressors with \eqn{\text{E}(u|X) \neq 0}. However, the endogenous regressors can be accounted for by exploiting the \eqn{q} instruments in \eqn{\bold{Z}}, for which \eqn{\text{E}(u|Z) = 0} holds.
#' The instruments and the first stage coefficients are generated in the same fashion as \eqn{\bold{X}} and \eqn{\bold{\alpha}} when \code{q = NULL}, respectively.
#'
#' The function nests, among other, the DGPs employed in the simulation study of Mehrabani (2023, sec. 6).
#'
#' @examples
#' # Simulate DGP 1 from Mehrabani (2023, sec. 6)
#' alpha_0_DGP1 <- matrix(c(0.4, 1, 1.6, 1.6, 1, 0.4), ncol = 2)
#' DGP1 <- sim_DGP(
#'   N = 50, n_periods = 20, p = 2, n_groups = 3,
#'   group_proportions = c(.4, .3, .3), alpha_0 = alpha_0_DGP1
#' )
#'
#' # Simulate DGP 6 from Mehrabani (2023, sec. 6)
#' alpha_0_DGP6 <- cbind(
#'   c(0.8, 0.6, 0.4, 0.2, -0.2, -0.4, -0.6, -0.8),
#'   c(-4, -3, -2, -1, 1, 2, 3, 4),
#'   c(4, 3, 2, 1, -1, -2, -3, -4)
#' )
#' @references
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{alpha}}{the \eqn{K \times p} matrix of group-specific slope parameters. In case of \code{dyn_panel = TRUE}, the first column holds the \eqn{AR} coefficient.}
#' \item{\code{groups}}{a vector indicating the group memberships.}
#' \item{\code{y}}{a \eqn{NT \times 1} vector of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.}
#' \item{\code{X}}{a \eqn{NT \times p} matrix of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.}
#' \item{\code{Z}}{a \eqn{NT \times q} matrix of instruments , where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1, \dots, z_N)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. In case a panel with exogenous regressors is generated (\code{q = NULL}), \eqn{\bold{Z}} equals \code{NULL}.}
#' @export


sim_DGP <- function(N = 50, n_periods = 40, p = 2, n_groups = 3, group_proportions = NULL, error_spec = NULL, dyn_panel = FALSE, q = NULL, alpha_0 = NULL) {
  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  if (N < n_groups) stop("Number of groups cannot exceed number of observations\n")
  if (p == 0) stop("Include at least one explanatory variable\n")
  if (!is.null(group_proportions)) {
    if (n_groups != length(group_proportions)) stop("Number of groups and group proportions are of different length\n")
    if (sum(group_proportions) != 1) stop("Group proportions must sum to 1\n")
  }
  if (!is.null(error_spec)) {
    if (!(error_spec %in% c("AR", "GARCH"))) stop("The individual error specificaion must be either AR, GARCH or NULL. Use AR in case of an AR(1) serial correlation, GARCH for an GARCH(1,1) innovation and NULL for iid errors.\n")
  }
  if (!is.null(alpha_0)) {
    if (nrow(alpha_0) != n_groups) stop(paste("There are", n_groups, "groups, but only", nrow(alpha_0), "parameter vectors provided"), "\n")
    if (dyn_panel) {
      if (any(abs(alpha_0[, 1]) >= 1)) stop("The AR parameters must be lower than 1 in absolute value\n")
      if (ncol(alpha_0) - 1 != p) stop(paste("There are", ncol(alpha_0), "group-specific parameters provided, but p + 1 =", p + 1, "are required"), "\n")
    } else {
      if (ncol(alpha_0) != p) stop(paste("There are", ncol(alpha_0), "group-specific parameters provided, but p =", p, "are required"), "\n")
    }
  }
  if (!is.null(q)) {
    if (q < p) stop("There must be at least q = p instruments\n")
  }

  #------------------------------#
  #### Generate parameters    ####
  #------------------------------#

  # Separate the AR coefficient from the slope parameters of exogenous explanatory variables
  if (is.null(alpha_0)) {
    alpha <- matrix(stats::runif(p * n_groups, -2, 2), ncol = p)
    if (dyn_panel) {
      alpha_ar <- stats::runif(n_groups, -.99, .99)
      alpha_0 <- cbind(alpha_ar, alpha)
    } else {
      alpha_0 <- alpha
      alpha_ar <- NULL
    }
  } else {
    if (dyn_panel) {
      alpha_ar <- as.matrix(alpha_0[, 1])
      alpha <- as.matrix(alpha_0[, -1])
    } else {
      alpha_ar <- NULL
      alpha <- alpha_0
    }
  }
  colnames(alpha_0) <- paste0("alpha_", 1:ncol(alpha_0))

  #------------------------------#
  #### Simulate the groupings ####
  #------------------------------#

  if (is.null(group_proportions)) {
    group_vec <- rep(1:n_groups, length.out = N)
  } else {
    group_vec <- rep(1:n_groups, round(group_proportions * N))
  }
  groups_raw <- sample(group_vec, N, replace = FALSE)
  groups <- match(groups_raw, unique(groups_raw))
  # Bring the vector into a stacked format (all p explanatory variables for i=1, p for i=2, ...)
  beta_vec <- as.matrix(c(t(alpha[groups, ])))

  #------------------------------#
  #### Construct X and y      ####
  #------------------------------#

  # Draw the cross-sectional errors
  u <- stats::rnorm(N * n_periods)
  # In case of a endogenous panel, correlate the errors with the regressor
  if (!is.null(q)) {
    corr_mat <- diag(p + 1)
    corr_mat[1,][corr_mat[1,] == 0] <- .5
    corr_mat[,1][corr_mat[,1] == 0] <- .5
    eps <- matrix(stats::rnorm(N * n_periods * p), ncol =p)
    U <- cbind(u, eps) %*% chol(corr_mat)
    u <- U[,1]
    eps <- U[,-1]
  } else {
    eps <- NULL
  }
  # Generate the regressors
  XList <- simX(N = N, n_periods = n_periods, p = p, dyn_panel = dyn_panel, eps = eps, q = q)
  X <- XList$X
  Z <- XList$Z
  # Bring into a block-matrix format to make later computations easier
  X_tilde <- buildDiagX(X = X, n_periods = n_periods, N = N, groups = 1:N)
  # Build a more elaborate error process
  if (!is.null(error_spec)) {
    uList <- split(u, rep(1:((N * n_periods) %/% n_periods), each = n_periods, length.out = N * n_periods))
    if (error_spec == "AR") {
      u <- simAR(errorList = uList)
    } else if (error_spec == "GARCH") {
      u <- simGARCH(errorList = uList)
    }
  }
  # Construct the observation
  y <- X_tilde %*% beta_vec + u
  if (dyn_panel) {
    beta_ar_vec <- c(t(alpha_ar[groups]))
    # Adjust the individual FE in case of a dynamic panel
    y <- y + (1 - beta_ar_vec) * rep(stats::rnorm(N), each = n_periods)
    y <- generateARpanel(y = y, N = N, beta = beta_ar_vec)
  }
  return(list(alpha = alpha_0, groups = groups, y = y, X = X, Z = Z))
}

simGARCH <- function(errorList) {
  n_periods <- length(errorList[[1]])
  # Construct a GARCH(1,1) process with coefficients as in Mehrabani (2023, sec. 6)
  errorGARCHList <- lapply(errorList, function(eps) {
    u <- c(eps[1], rep(NA, n_periods - 1))
    h <- rep(0, n_periods)
    for (i in 2:n_periods) {
      h[i] <- .05 + .05 * u[i - 1]^2 + .9 * h[i - 1]
      u[i] <- sqrt(h[i]) * eps[i]
    }
    return(u)
  })
  u <- c(do.call(rbind, errorGARCHList))
  return(u)
}

simAR <- function(errorList) {
  errorARList <- lapply(errorList, function(x) {
    # Create an AR(1) process with with coefficients as in Mehrabani (2023, sec. 6)
    stats::filter(c(0, x), filter = .5, method = "recursive")[-1]
  })
  u <- c(do.call(rbind, errorARList))
  return(u)
}

simX <- function(N, n_periods, p, dyn_panel, eps, q) {
  # Draw individual FE
  eta <- rep(stats::rnorm(N), each = n_periods)
  # Construct the predictor variables
  if (is.null(eps)) {
    # Draw the exogenous regressors
    X <- matrix(stats::rnorm(N * n_periods * p), ncol =p)
    Z <- NULL
  } else {
    # Draw the instruments and loadings
    Z <- matrix(stats::rnorm(N * n_periods * q), ncol =q)
    gamma <- matrix(stats::runif(q * p, -2, 2), ncol =p)
    # Construct the regressor
    X <- Z %*% gamma + eps
  }
  if (!dyn_panel) {
    X <- .2 * eta + X
  }
  return(list(X = X, Z = Z))
}

generateARpanel <- function(y, N, beta) {
  n_periods <- length(y) / N
  # Create a list with one entry per individual
  yList <- split(y, rep(1:((N * n_periods) %/% n_periods), each = n_periods, length.out = N * n_periods))
  # Add the individual AR terms to the panel
  yARMat <- mapply(function(ar, yy) {
    stats::filter(c(0, yy), filter = ar, method = "recursive")[-1]
  }, beta, yList, SIMPLIFY = T)
  return(c(yARMat))
}
