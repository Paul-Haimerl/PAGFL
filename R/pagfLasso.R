#' Pairwise Adaptive Group Fused Lasso
#'
#' @description The pairwise adaptive group fused lasso (PAGFL) by Mehrabani (2023) jointly estimates the latent group structure and group-specific slope parameters in a panel data model.
#' It can handle static and dynamic panels, either with or without endogenous regressors.
#'
#' @param y a \eqn{NT \times 1} vector, matrix, or data.frame of the dependent variable. If no \code{index} variables are explicitly provided, the data must be balanced and ordered in the long format as follows \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} with the scalar \eqn{y_{it}}. If \code{y} is not ordered or not balanced, \code{y} must include two index variables, declaring the cross-sectional unit \eqn{i} and the time period \eqn{t} for each observation.
#' @param X a \eqn{NT \times p} vector, matrix, or data.frame of explanatory variables. If no \code{index} variables are explicitly provided, the data must be balanced and ordered in the long format as follows \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}. If \code{X} is not ordered or not balanced, \code{X} must include two index, declaring the cross-sectional unit \eqn{i} and the time period \eqn{t} for each observation.
#' @param index a character vector holding two strings specifying the variable names that identify the cross-sectional unit and time period for each observation. The first string denotes the individual unit, while the second string represents the time period. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied. However, in such instances no indicator variables must be present in \code{y} or \code{X}.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Else, the data supplied must be a balanced panel data set ordered in the long format.
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
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is 5e3.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is 1e-8.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are placed in the same group if the Frobenius norm of their coefficient parameter difference is below this parameter. Default is 1e-3.
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
#' estim <- pagfl(y = y, X = X, n_periods = 80, lambda = lambda_set, method = "PLS")
#' print(estim)
#'
#' # Lets pass a panel data set with explicit cross-sectional and time indicators
#' i_index <- rep(1:50, each = 80)
#' t_index <- rep(1:80, 50)
#' y <- cbind(y, i_index = i_index, t_index = t_index)
#' X <- cbind(X, i_index = i_index, t_index = t_index)
#' estim <- pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = lambda_set, method = "PLS")
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#'
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @aliases PAGFL
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
pagfl <- function(y, X, index = NULL, n_periods = NULL, lambda, method = "PLS", Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 5e3, tol_convergence = 1e-8,
                  tol_group = 1e-3, rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE) {
  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  prelim_checks(y, X, Z, index, n_periods, method)

  if (!is.null(index)) {
    df <- merge(x = X, y = y, by = index)
    if (method == "PGMM") df <- merge(x = df, y = Z, by = index)
    df <- stats::na.omit(df)
    df <- df[order(df[, index[1]], df[, index[2]]), ]
    if (method == "PGMM") {
      Z <- df[, (ncol(df) - ncol(Z) + 3):ncol(df)]
      Z <- as.matrix(Z)
      df <- df[, -((ncol(df) - ncol(Z) + 1):ncol(df))]
    }
    y <- df[, ncol(df)]
    df <- df[, -ncol(df)]
    X <- df[, !(colnames(df) %in% index)]
    i_index_labs <- df[, index[1]]
    i_index <- as.integer(factor(i_index_labs))
    t_index_labs <- df[, index[2]]
    t_index <- as.integer(factor(t_index_labs))
    n_periods <- length(unique(t_index))
    N <- length(unique(i_index))
  } else {
    N <- NROW(y) / n_periods
    i_index <- i_index_labs <- rep(1:N, each = n_periods)
    t_index <- rep(1:n_periods, N)
  }

  y <- as.matrix(y)
  X <- as.matrix(X)
  p <- ncol(X)
  varrho <- abs(varrho)

  second_checks(N, index, n_periods, y, X, method, Z, p, min_group_frac, verbose, dyn = FALSE)

  # In case of penalized Least Squares, specify an empty instrument matrix Z
  if (method == "PLS") {
    Z <- matrix()
  } else {
    Z <- as.matrix(Z)
  }

  #------------------------------#
  #### Iterate over lambda    ####
  #------------------------------#

  if (verbose & length(lambda) > 1) {
    FUN <- pbapply::pblapply
  } else {
    FUN <- lapply
  }
  lambdaList <- FUN(lambda, function(lam) {
    # Run the algorithm
    estimOutput <- pagfl_algo(
      y = y, X = X, method = method, Z = Z, bias_correc = bias_correc, i_index = i_index,
      t_index = t_index, N = N, lambda = lam, kappa = kappa, min_group_frac = min_group_frac,
      max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group, varrho = varrho
    )
    # Compute the Information Criterion
    IC_val <- IC(estimOutput = estimOutput, y = y, X = X, rho = rho, method = method, N = N, i_index = i_index)
    return(c(IC = IC_val, lambda = lam, estimOutput))
  })
  # Pick the estimation result with the lowest IC
  IC_vec <- lapply(lambdaList, function(x) x[[1]])
  estim <- lambdaList[[which.min(IC_vec)]]
  estim$groups_hat <- c(estim$groups_hat)
  names(estim$groups_hat) <- unique(i_index_labs)
  rownames(estim$alpha_hat) <- paste("Group", 1:estim$K_hat)
  colnames(estim$alpha_hat) <- colnames(X)
  return(estim)
}

#' @export
PAGFL <- function(y, X, n_periods, lambda, method = "PLS", Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 2e3, tol_convergence = 1e-3,
                  tol_group = sqrt(p / (sqrt(N * n_periods) * log(log(N * n_periods)))), rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE) {
  lifecycle::deprecate_warn(when = "1.1.0", what = "PAGFL()", with = "pagfl()", env = asNamespace("PAGFL"))
  N <- NROW(y) / n_periods
  p <- ncol(X)
  estim <- pagfl(
    y = y, X = X, n_periods = n_periods, lambda = lambda, method = method, Z = Z, min_group_frac = min_group_frac, bias_correc = bias_correc,
    kappa = kappa, max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group, rho = rho, varrho = varrho, verbose = verbose
  )
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

  prelim_checks(y = y, X = X, Z = Z, n_periods = n_periods, method = method)

  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  # In case of penalized Least Squares, specify an empty instrument matrix Z
  if (method == "PLS") {
    Z <- matrix()
  } else {
    Z <- as.matrix(Z)
  }

  i_index <- rep(1:N, each = n_periods)

  # Net out fixed effects
  data <- netFE(y = y, X = X, method = method, N = N, i_index = i_index)
  y_tilde <- data[[1]]
  X_tilde <- matrix(data[[2]], ncol = ncol(X))
  if (method == "PGMM") {
    Z_tilde <- deleteObsMat(Z, N, i_index = i_index, first = TRUE)
    n_periods <- n_periods - 1
  } else {
    Z_tilde <- matrix()
  }

  #------------------------------#
  #### Oracle estimation      ####
  #------------------------------#

  alpha_hat <- getAlpha(X_tilde, y_tilde, Z_tilde, method, N, i_index = i_index, p, groups_hat = groups_0)
  if (bias_correc) {
    alpha_hat <- spjCorrec(alpha_hat, X, y, Z, N, i_index = i_index, p = p, groups_0, method)
  }
  rownames(alpha_hat) <- paste("Group", 1:max(groups_0))
  colnames(alpha_hat) <- colnames(X)
  return(alpha_hat)
}
