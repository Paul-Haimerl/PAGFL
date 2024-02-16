#' Apply the Pairwise Adaptive Group Fused Lasso
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
#' @param rho the tuning parameter balancing the fitness and penalty terms in the information criterion that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\sqrt{NT} \log(NT)}{NT}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For \emph{PLS}, the \eqn{\varrho} value is trivial. However, for \emph{PGMM}, small values lead to slow convergence of the algorithm. If left unspecified, the default heuristic \eqn{\varrho = \max(\frac{\sqrt{5NTp}}{\log(NTp)}-7, 1}) is used.
#' @param verbose logical. If \code{TRUE}, a progression bar is printed when iterating over candidate \eqn{\lambda} values and helpful warning messages are shown. Default is \code{TRUE}.
#'
#' @details
#' The \emph{PLS} method minimizes the following criterion:
#' \deqn{\frac{1}{T} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \beta^\prime_i \tilde{x}_{it})^2 + \frac{\lambda}{N} \sum_{1 \leq i} \sum_{i<j \leq N} \dot{w}_{ij} \| \beta_i - \beta_j \|,}
#' where \eqn{\tilde{y}_{it}} is the de-meaned dependent variable, \eqn{\tilde{x}_{it}} represents a vector of de-meaned weakly exogenous explanatory variables, \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} reflects adaptive penalty weights (see Mehrabani, 2023, eq. 2.6). \eqn{\| \cdot \|} denotes the Frobenius norm.
#' The adaptive weights \eqn{\dot{w}_{ij}} are obtained by a preliminary least squares estimation.
#' The solution \eqn{\hat{\beta}} is computed via an iterative alternating direction method of multipliers (\emph{ADMM}) algorithm (see Mehrabani, 2023, sec. 5.1).
#'
#' \emph{PGMM} employs a set of instruments \eqn{Z} to control for endogenous regressors. Using \emph{PGMM}, \eqn{\bold{\beta} = (\beta_1^\prime, \dots, \beta_N^\prime)^\prime} is estimated by minimizing:
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
#' estim <- PAGFL(y = y, X = X, n_periods = 80, lambda = lambda_set, method = 'PLS')
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
#' \item{\code{K_hat}}{the estimated total number of groups.}
#' \item{\code{groups_hat}}{a vector of estimated group memberships.}
#' \item{\code{iter}}{the number of executed algorithm iterations.}
#' \item{\code{convergence}}{logical. If \code{TRUE}, convergence was achieved. If \code{FALSE}, \code{max_iter} was reached.}
#' @export

PAGFL <- function(y, X, n_periods, lambda, method = 'PLS', Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 2e3, tol_convergence = 1e-3,
                  tol_group = sqrt(p / (sqrt(N * n_periods) * log(log(N * n_periods)))), rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  N <- nrow(y) / n_periods
  p <- ncol(X)

  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  if (round(N) != N) stop('n_periods does not match the number of time periods of the dependent variable y\n')
  if (N * n_periods != nrow(X)) stop('The number of time periods of the dependent variable y and the predictor matrix X do not match\n')
  if (!(method %in% c('PLS', 'PGMM'))) stop('The estimation method must be either PLS or PGMM. Use PLS in case of (weakly) exogenous regressors and PGMM for endogenous regressors.\n')
  if (method == 'PGMM') {
    if (is.null(Z)) stop('PGMM requires a matrix of exogenous instruments Z\n')
    if (ncol(Z) < p) stop(paste('Provide at least p =', p, 'exogenous instruments Z'), '\n')
    if (nrow(Z) != nrow(X)) stop('The number of time periods of the instrument matrix Z does do not match the remaining data\n')
  } else {
    if (!is.null(Z) & verbose) warning('The instrument matrix Z is ignored by the PLS estimation algorithm. To instrument endogenous regressors using Z, specify method = PGMM\n')
  }
  if (is.null(min_group_frac)) min_group_frac <- 0
  if (min_group_frac > 1 | min_group_frac < 0) stop('Provide a min_group_frac between 0 and 1\n')
  if (min_group_frac >= .4 & verbose) warning('Large min_group_frac values may lead to all groups falling below the group cardinality threshold, in which case the hierarchical clustering algorithm cannot be employed\n')
  if (any(is.na(X))) stop('The predictor matrix X contains missing values\n')
  if (any(is.na(y))) stop('The dependent variable y contains missing values\n')

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
    estimOutput <- PAGFL_Algo(
      y = y, X = X, n_periods = n_periods, method = method,
      Z = Z, bias_correc = bias_correc, lambda = lam, kappa = kappa, min_group_frac = min_group_frac,
      max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group, varrho = varrho
    )
    estimOutput$groups_hat <- c(estimOutput$groups_hat)
    # Compute the Information Criterion
    IC_val <- IC(estimOutput = estimOutput, y = y, X = X, rho = rho, method = method, n_periods = n_periods,
                 N = N)
    return(c(IC = IC_val, lambda = lam, estimOutput))
  })
  # Pick the estimation result with the lowest IC
  IC_vec <- lapply(lambdaList, function(x) x[[1]])
  return(lambdaList[[which.min(IC_vec)]])
}
