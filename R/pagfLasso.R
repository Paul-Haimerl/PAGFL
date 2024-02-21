#' Pairwise Adaptive Group Fused Lasso
#'
#' @description The pairwise adaptive group fused lasso (PAGFL) by Mehrabani (2023) jointly estimates the latent group structure and group-specific slope parameters in a panel data model.
#' It can handle static and dynamic panels, either with or without endogenous regressors.
#'
#' @param y a \eqn{NT \times 1} vector or data.frame of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.
#' @param X a \eqn{NT \times p} matrix or data.frame of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.
#' @param n_periods the number of observed time periods \eqn{T}.
#' @param lambda the tuning parameter governing the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type information criterion selects the best fitting parameter value.
#' @param method the estimation method. Options are
#' \describe{
#' \item{\code{'PLS'}}{for using the penalized least squares (\emph{PLS}) algorithm. We recommend \emph{PLS} in case of (weakly) exogenous regressors (Mehrabani, 2023, sec. 2.2).}
#' \item{\code{'PGMM'}}{for using the penalized Generalized Method of Moments (\emph{PGMM}). \emph{PGMM} is required when instrumenting endogenous regressors (Mehrabani, 2023, sec. 2.3). A matrix \eqn{Z} contains the necessary exogenous instruments.}
#' } Default is \code{'PLS'}.
#' @param Z a \eqn{NT \times q} matrix of exogenous instruments, where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1, \dots, z_N)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. \eqn{Z} is only required when \code{method = 'PGMM'} is selected. When using \code{'PLS'}, either pass \code{NULL} or any matrix \eqn{\bold{Z}} is disregarded. Default is \code{NULL}.
#' @param min_group_frac the minimum group size as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, a hierarchical classifier allocates its members to the remaining groups. Default is 0.05.
#' @param bias_correc logical. If \code{TRUE}, a Split-panel Jackknife bias correction following Dhaene and Jochmans (2015) is applied to the slope parameters. We recommend using this correction when facing a dynamic panel. Default is \code{FALSE}.
#' @param kappa the weight placed on the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is 2000.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is 0.001.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are placed in the same group if the Frobenius norm of their coefficient parameter difference is below this parameter. If left unspecified, the heuristic \eqn{\sqrt{\frac{p}{\sqrt{NT} \log(\log(NT))}}} is used. We recommend the default.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the information criterion that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For \emph{PLS}, the \eqn{\varrho} value is trivial. However, for \emph{PGMM}, small values lead to slow convergence of the algorithm. If left unspecified, the default heuristic \eqn{\varrho = \max(\frac{\sqrt{5NTp}}{\log(NTp)}-7, 1}) is used.
#' @param verbose logical. If \code{TRUE}, a progression bar is printed when iterating over candidate \eqn{\lambda} values and helpful warning messages are shown. Default is \code{TRUE}.
#'
#' @details
#' The \emph{PLS} method minimizes the following criterion:
#' \deqn{\frac{1}{T} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \beta^\prime_i \tilde{x}_{it})^2 + \frac{\lambda}{N} \sum_{1 \leq i} \sum_{i<j \leq N} \dot{w}_{ij} \| \beta_i - \beta_j \|,}
#' where \eqn{\tilde{y}_{it}} is the demeaned dependent variable, \eqn{\tilde{x}_{it}} represents a vector of demeaned weakly exogenous explanatory variables, \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} reflects adaptive penalty weights (see Mehrabani, 2023, eq. 2.6). \eqn{\| \cdot \|} denotes the Frobenius norm.
#' The adaptive weights \eqn{\dot{w}_{ij}} are obtained by a preliminary least squares estimation.
#' The solution \eqn{\hat{\bold{\beta}}} is computed via an iterative alternating direction method of multipliers (\emph{ADMM}) algorithm (see Mehrabani, 2023, sec. 5.1).
#'
#' \emph{PGMM} employs a set of instruments \eqn{\bold{Z}} to control for endogenous regressors. Using \emph{PGMM}, \eqn{\bold{\beta} = (\beta_1^\prime, \dots, \beta_N^\prime)^\prime} is estimated by minimizing:
#' \deqn{\sum^N_{i = 1} \left[ \frac{1}{N} \sum_{t=1}^T z_{it} (\Delta y_{it} - \beta^\prime_i \Delta x_{it}) \right]^\prime W_i \left[\frac{1}{T} \sum_{t=1}^T z_{it}(\Delta y_{it} - \beta^\prime_i \Delta x_{it}) \right] + \frac{\lambda}{N} \sum_{1 \leq i} \sum_{i<j \leq N} \ddot{w}_{ij} \| \beta_i - \beta_j \|.}
#' \eqn{\ddot{w}_{ij}} are obtained by an initial \emph{GMM} estimation. \eqn{\Delta} gives the first differences operator \eqn{\Delta y_{it} = y_{it} - y_{i t-1}}. \eqn{W_i} represents a data-driven \eqn{q \times q} weight matrix. I refer to Mehrabani (2023, eq. 2.10) for more details.
#' \eqn{\bold{\beta}} is again estimated employing an efficient \emph{ADMM} algorithm (Mehrabani, 2023, sec. 5.2).
#'
#' Two individuals are assigned to the same group if \eqn{\| \hat{\beta}_i - \hat{\beta}_j \| \leq \epsilon_{\text{tol}}}, where \eqn{\epsilon_{\text{tol}}} is given by \code{tol_group}.
#'
#' We suggest identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit of 0 and an upper limit that leads to a fully homogenous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.
#'
#' @examples
#' # Simulate a panel with a group structure
#' sim <- sim_DGP(N = 50, n_periods = 80, p = 2, n_groups = 3)
#' y <- sim$y
#' X <- sim$X
#'
#' # Run the PAGFL procedure for a set of candidate tuning parameter values
#' lambda_set <- exp(log(10) * seq(log10(1e-4), log10(10), length.out = 10))
#' estim <- pagfl(y = y, X = X, n_periods = 80, lambda = lambda_set, method = 'PLS')
#' print(estim)
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#'
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{IC}}{the BIC-type information criterion.}
#' \item{\code{lambda}}{the penalization parameter. If multiple \eqn{\lambda} values were passed, the parameter yielding the lowest IC.}
#' \item{\code{alpha_hat}}{a \eqn{K \times p} matrix of the post-Lasso group-specific parameter estimates.}
#' \item{\code{K_hat}}{the estimated total number of groups \eqn{\widehat{K}}.}
#' \item{\code{groups_hat}}{a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k}.}
#' \item{\code{iter}}{the number of executed algorithm iterations.}
#' \item{\code{convergence}}{logical. If \code{TRUE}, convergence was achieved. If \code{FALSE}, \code{max_iter} was reached.}
#' @export
pagfl <- function(y, X, n_periods, lambda, method = 'PLS', Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 2e3, tol_convergence = 1e-3,
                  tol_group = sqrt(p / (sqrt(N * n_periods) * log(log(N * n_periods)))), rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  N <- nrow(y) / n_periods
  p <- ncol(X)
  varrho <- abs(varrho)

  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  checks(N, n_periods, y, X, method, Z, p, min_group_frac, verbose, dyn = FALSE)

  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  # In case of penalized Least Squares, specify an empty instrument matrix Z
  if (method == 'PLS') {
    Z <- matrix()
  } else {
    Z <- as.matrix(Z)
  }

  #------------------------------#
  #### Iterate over lambda    ####
  #------------------------------#

  if (verbose) {
    FUN <- pbapply::pblapply
  } else {
    FUN <- lapply
  }
  lambdaList <- FUN(lambda, function(lam) {
    # Run the algorithm
    estimOutput <- pagfl_algo(
      y = y, X = X, n_periods = n_periods, method = method,
      Z = Z, bias_correc = bias_correc, lambda = lam, kappa = kappa, min_group_frac = min_group_frac,
      max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group, varrho = varrho
    )
    # Compute the Information Criterion
    IC_val <- IC(estimOutput = estimOutput, y = y, X = X, rho = rho, method = method, n_periods = n_periods,
                 N = N)
    return(c(IC = IC_val, lambda = lam, estimOutput))
  })
  # Pick the estimation result with the lowest IC
  IC_vec <- lapply(lambdaList, function(x) x[[1]])
  estim <- lambdaList[[which.min(IC_vec)]]
  colnames(estim$alpha_hat) <- paste0("alpha_", 1:p)
  rownames(estim$alpha_hat) <- paste0("k=", 1:estim$K_hat)
  estim$groups_hat <- c(estim$groups_hat)
  return(estim)
}

#' Grouped Panel Data Oracle Estimator
#'
#' @description Estimation of a panel data model, given a known group structure.
#'
#' @param y a \eqn{NT \times 1} vector or data.frame of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.
#' @param X a \eqn{NT \times p} matrix or data.frame of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.
#' @param n_periods the number of observed time periods \eqn{T}.
#' @param method the estimation method. Options are
#' \describe{
#' \item{\code{'PLS'}}{for using the penalized least squares (\emph{PLS}) algorithm. We recommend \emph{PLS} in case of (weakly) exogenous regressors (Mehrabani, 2023, sec. 2.2).}
#' \item{\code{'PGMM'}}{for using the penalized Generalized Method of Moments (\emph{PGMM}). \emph{PGMM} is required when instrumenting endogenous regressors (Mehrabani, 2023, sec. 2.3). A matrix \eqn{Z} contains the necessary exogenous instruments.}
#' } Default is \code{'PLS'}.
#' @param Z a \eqn{NT \times q} matrix of exogenous instruments, where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1, \dots, z_N)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. \eqn{Z} is only required when \code{method = 'PGMM'} is selected. When using \code{'PLS'}, either pass \code{NULL} or any matrix \eqn{\bold{Z}} is disregarded. Default is \code{NULL}.
#' @param groups_0 the true group vector \eqn{(g_1, \dots, g_N)}, where \eqn{g_i = k} if the unit \eqn{i} is part of group \eqn{k}.
#' @param bias_correc logical. If \code{TRUE}, a Split-panel Jackknife bias correction following Dhaene and Jochmans (2015) is applied to the slope parameters. We recommend using this correction when facing a dynamic panel. Default is \code{FALSE}.
#'
#' @details
#' Given a known group pattern \eqn{\bold{G}^0}, the least squares oracle estimator (pendant to \emph{PLS}) solves
#' \deqn{\hat{\alpha}_k^{\text{oracle}} = \left( \sum_{i \in G_k^0} \sum_{t = 1}^T \tilde{x}_{it} \tilde{x}_{it}^\prime \right)^{-1} \sum_{i \in G_k^0} \sum_{t = 1}^T \tilde{x}_{it} \tilde{y}_{it},}
#' where \eqn{\tilde{y}_{it}} is the demeaned dependent variable and \eqn{\tilde{x}_{it}} represents a vector of demeaned weakly exogenous explanatory variables.
#'
#' The oracle estimates equivalent to \emph{PGMM} are obtained by
#' \deqn{\hat{\alpha}_k^{\text{oracle}} = \left( \tilde{Q}_{k \Delta x}^{(k)^\prime} \tilde{W}_{NT}^{(k)} \tilde{Q}_{k \Delta x}^{(k)} \right)^{-1} \tilde{Q}_{k \Delta x}^{(k)^\prime} \tilde{W}_{NT}^{(k)} \tilde{Q}_{k \Delta y}^{(k)},}
#' where \eqn{\tilde{Q}_{k \Delta x}^{(k)} = T^{-1} \sum_{i \in G_k^0} \sum_{t = 1}^T z_{it} (\Delta x_{it})^\prime}, \eqn{\tilde{Q}_{k \Delta y}^{(k)} = T^{-1} \sum_{i \in G_k^0} \sum_{t = 1}^T z_{it} \Delta y_{it}}, and \eqn{\tilde{W}_{NT}^{(k)}} is a matrix containing weights. \eqn{z_{it}} represents a \eqn{q \times 1} vector of exogenous instruments employed to control for the endogenous \eqn{x_{it}}.
#'
#' @examples
#' # Simulate a panel with a latent group pattern
#' sim <- sim_DGP(N = 50, n_periods = 80, p = 2, n_groups = 3)
#' y <- sim$y
#' X <- sim$X
#' groups_0 <- sim$groups
#'
#' # Estimate the grouped panel data model model, given the true grouping
#' estim <- oracle(y = y, X = X, n_periods = 80, groups_0 = groups_0)
#' print(estim)
#' @author Paul Haimerl
#'
#' @return A \eqn{K \times p} matrix of the group-specific slope coefficients.
#'
#' @export
oracle <- function(y, X, n_periods, groups_0, method = "PLS", Z = NULL, bias_correc = FALSE) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  N <- nrow(y) / n_periods
  p <- ncol(X)

  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  checks(N = N, n_periods = n_periods, y = y, X = X, p = p, min_group_frac = NULL, verbose = FALSE, dyn = FALSE, method = method, Z = Z)

  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  # In case of penalized Least Squares, specify an empty instrument matrix Z
  if (method == 'PLS') {
    Z <- matrix()
  } else {
    Z <- as.matrix(Z)
  }

  # Net out fixed effects
  data <- netFE(y = y, X = X, method = method, n_periods = n_periods, N = N)
  y_tilde <- data[[1]]
  X_tilde <- matrix(data[[2]], ncol = ncol(X))
  if (method == "PGMM")
  {
    Z_tilde = deleteFirstObsMat(Z, n_periods, N, ncol(Z))
    n_periods = n_periods - 1
  } else {
    Z_tilde <- matrix()
  }

  #------------------------------#
  #### Oracle estimation      ####
  #------------------------------#

  alpha_hat <- getAlpha(X_tilde, y_tilde, Z_tilde, method, n_periods, N, p, groups_0)
  if (bias_correc){
    alpha_hat <- spjCorrec(alpha_hat, X, y, Z, n_periods, N, p, groups_0, method)
  }
  colnames(alpha_hat) <- paste0("alpha_", 1:p)
  rownames(alpha_hat) <- paste0("k=", 1:max(groups_0))
  return(alpha_hat)
}
