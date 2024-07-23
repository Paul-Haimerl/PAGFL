#' Simulate a Panel With a Group Structure in the Slope Coefficients
#'
#' @description Construct a static or dynamic, exogenous or endogenous panel data set subject to a group structure in the slope coefficients with optional \eqn{AR(1)} or \eqn{GARCH(1,1)} innovations.
#'
#' @param N the number of cross-sectional units. Default is 50.
#' @param n_periods the number of simulated time periods \eqn{T}. Default is 40.
#' @param p the number of explanatory variables. Default is 2.
#' @param n_groups the number of groups \eqn{K}. Default is 3.
#' @param group_proportions a numeric vector of length \code{n_groups} indicating size of each group as a fraction of \eqn{N}. If \code{NULL}, all groups are of size \eqn{N / K}. Default is \code{NULL}.
#' @param error_spec options include
#' \describe{
#' \item{\code{"iid"}}{for \eqn{iid} errors.}
#' \item{\code{"AR"}}{for an \eqn{AR(1)} error process with an autoregressive coefficient of 0.5.}
#' \item{\code{"GARCH"}}{for a \eqn{GARCH(1,1)} error process with a 0.05 constant, a 0.05 ARCH and a 0.9 GARCH coefficient.}
#' }
#' Default is \code{"iid"}.
#' @param dynamic Logical. If \code{TRUE}, the panel includes one stationary autoregressive lag of \eqn{y_{it}} as an explanatory variable (see sec. Details for more information on the \eqn{AR} coefficient). Default is \code{FALSE}.
#' @param dyn_panel `r lifecycle::badge("deprecated")` deprecated and replaced by \code{dynamic}.
#' @param q the number of exogenous instruments when a panel with endogenous regressors is to be simulated. If panel data set with exogenous regressors is supposed to be generated, pass \code{NULL}. Default is \code{NULL}.
#' @param alpha_0 a \eqn{K \times p} matrix of group-specific coefficients. If \code{dynamic = TRUE}, the first column represents the stationary \eqn{AR} coefficient. If \code{NULL}, the coefficients are drawn randomly (see sec. Details). Default is \code{NULL}.
#'
#' @details
#' The scalar dependent variable \eqn{y_{it}} is generated according to the following grouped panel data model
#' \deqn{y_{it} = \gamma_i + \beta_i^\prime x_{it} + u_{it}, \quad i = \{1, \dots, N\}, \quad t = \{1, \dots, T\}.}
#' \eqn{\gamma_i} represents individual fixed effects and \eqn{x_{it}} a \eqn{p \times 1} vector of regressors.
#' The individual slope coefficient vectors \eqn{\beta_i} are subject to a group structure
#' \deqn{\beta_i = \sum_{k = 1}^K \alpha_k \bold{1} \{i \in G_k\},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k - \alpha_j \| \neq 0} for any \eqn{k \neq j}, \eqn{k = 1, \dots, K}. The total number of groups \eqn{K} is determined by \code{n_groups}.
#'
#' If a panel data set with exogenous regressors is generated (set \code{q = NULL}), the explanatory variables are simulated according to
#' \deqn{x_{it,j} = 0.2 \gamma_i + e_{it,j}, \quad \gamma_i,e_{it,j} \sim i.i.d. N(0, 1), \quad j = \{1, \dots, p\},}
#' where \eqn{e_{it,j}} denotes a series of innovations. \eqn{\gamma_i} and \eqn{e_i} are independent of each other.
#'
#' In case \code{alpha_0 = NULL}, the group-level slope parameters \eqn{\alpha_{k}} are drawn from \eqn{\sim U[-2, 2]}.
#'
#' If a dynamic panel is specified (\code{dynamic = TRUE}), the \eqn{AR} coefficients \eqn{\beta^{\text{AR}}_i} are drawn from a uniform distribution with support \eqn{(-1, 1)} and \eqn{x_{it,j} = e_{it,j}}.
#' Moreover, the individual fixed effects enter the dependent variable via \eqn{(1 - \beta^{\text{AR}}_i) \gamma_i} to account for the autoregressive dependency.
#' We refer to Mehrabani (2023, sec 6) for details.
#'
#' When specifying an endogenous panel (set \code{q} to \eqn{q \geq p}), the \eqn{e_{it,j}} correlate with the cross-sectional innovations \eqn{u_{it}} by a magnitude of 0.5 to produce endogenous regressors (\eqn{\text{E}(u|X) \neq 0}). However, the endogenous regressors can be accounted for by exploiting the \eqn{q} instruments in \eqn{\bold{Z}}, for which \eqn{\text{E}(u|Z) = 0} holds.
#' The instruments and the first stage coefficients are generated in the same fashion as \eqn{\bold{X}} and \eqn{\bold{\alpha}} when \code{q = NULL}.
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
#' @references
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{alpha}}{the \eqn{K \times p} matrix of group-specific slope parameters. If \code{dynamic = TRUE}, the first column holds the \eqn{AR} coefficient.}
#' \item{\code{groups}}{a vector indicating the group memberships \eqn{(g_1, \dots, g_N)}, where \eqn{g_i = k} if \eqn{i \in} group \eqn{k}.}
#' \item{\code{y}}{a \eqn{NT \times 1} vector of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.}
#' \item{\code{X}}{a \eqn{NT \times p} matrix of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.}
#' \item{\code{Z}}{a \eqn{NT \times q} matrix of instruments , where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1, \dots, z_N)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. In case a panel with exogenous regressors is generated (\code{q = NULL}), \eqn{\bold{Z}} equals \code{NULL}.}
#' \item{\code{data}}{a \eqn{NT \times (p + 1)} data.frame of the outcome and the explanatory variables.}
#' @export
sim_DGP <- function(N = 50, n_periods = 40, p = 2, n_groups = 3, group_proportions = NULL, error_spec = "iid", dynamic = FALSE, dyn_panel = lifecycle::deprecated(), q = NULL, alpha_0 = NULL) {
  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  if (lifecycle::is_present(dyn_panel)) {
    lifecycle::deprecate_warn("1.1.0", "sim_DGP(dyn_panel)", "sim_DGP(dynamic)")
    dynamic <- dyn_panel
  }

  error_spec <- match.arg(error_spec, c("AR", "GARCH", "iid"))
  simChecks(
    dyn = FALSE, N = N, n_groups = n_groups, group_proportions = group_proportions, error_spec = error_spec,
    alpha_0 = alpha_0, dyn_panel = dynamic, q = q, p = p, dynamic = dynamic
  )

  #------------------------------#
  #### Generate parameters    ####
  #------------------------------#

  if (dynamic) p <- p - 1
  # Separate the AR coefficient from the slope parameters of exogenous explanatory variables
  if (is.null(alpha_0)) {
    alpha <- matrix(stats::runif(p * n_groups, -2, 2), ncol = p)
    if (dynamic) {
      alpha_ar <- stats::runif(n_groups, -.99, .99)
      alpha_0 <- cbind(alpha_ar, alpha)
    } else {
      alpha_0 <- alpha
      alpha_ar <- NULL
    }
  } else {
    if (dynamic) {
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
    if (length(group_vec) != N) group_vec <- rep(group_vec, length.out = N)
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
  # Draw the FE
  gamma <- rep(stats::rnorm(N), each = n_periods)
  # In case of a endogenous panel, correlate the errors with the regressor
  if (!is.null(q)) {
    corr_mat <- diag(p + 1)
    corr_mat[1, ][corr_mat[1, ] == 0] <- .5
    corr_mat[, 1][corr_mat[, 1] == 0] <- .5
    eps <- matrix(stats::rnorm(N * n_periods * p), ncol = p)
    U <- cbind(u, eps) %*% chol(corr_mat)
    u <- U[, 1]
    eps <- U[, -1]
  } else {
    eps <- NULL
  }
  # Generate the regressors
  XList <- simX(N = N, n_periods = n_periods, p = p, dyn_panel = dynamic, eps = eps, q = q)
  X <- XList$X
  Z <- XList$Z
  # Bring into a block-matrix format to make later computations easier
  X_tilde <- buildDiagX_block_dense(X = X, N = N, i_index = rep(1:N, each = n_periods), groups = 1:N)
  # Build a more elaborate error process
  if (error_spec != "iid") {
    uList <- split(u, rep(1:((N * n_periods) %/% n_periods), each = n_periods, length.out = N * n_periods))
    if (error_spec == "AR") {
      u <- simAR(errorList = uList)
    } else if (error_spec == "GARCH") {
      u <- simGARCH(errorList = uList)
    }
  }
  # Construct the observation
  y <- X_tilde %*% beta_vec + u + gamma
  if (dynamic) {
    beta_ar_vec <- c(t(alpha_ar[groups]))
    # Adjust the individual FE in case of a dynamic panel
    y <- y + (1 - beta_ar_vec) * gamma
    ARPanel <- generateARpanel(y = y, N = N, beta = beta_ar_vec)
    y <- ARPanel[, "y"]
    names(y) <- NULL
    X <- cbind(ARPanel[, "X"], X)
    rownames(X) <- NULL
  }
  data <- data.frame(y = c(y), X)
  return(list(alpha = alpha_0, groups = groups, y = y, X = X, Z = Z, data = data))
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
  u <- c(t(do.call(rbind, errorARList)))
  return(u)
}

simX <- function(N, n_periods, p, dyn_panel, eps, q) {
  # Draw individual FE
  eta <- rep(stats::rnorm(N), each = n_periods)
  # Construct the predictor variables
  if (is.null(eps)) {
    # Draw the exogenous regressors
    X <- matrix(stats::rnorm(N * n_periods * p), ncol = p)
    Z <- NULL
  } else {
    # Draw the instruments and loadings
    Z <- matrix(stats::rnorm(N * n_periods * q), ncol = q)
    gamma <- matrix(stats::runif(q * p, -2, 2), ncol = p)
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
  # Create a regressor vector
  xARVec <- unlist(lapply(yList, function(x) c(0, x[-length(x)])))
  out <- cbind(X = xARVec, y = c(yARMat))
  return(out)
}
