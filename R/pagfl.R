#' Pairwise Adaptive Group Fused Lasso
#'
#' @description Estimate a panel data model subject to a latent group structure using the pairwise adaptive group fused Lasso (\emph{PAGFL}) by Mehrabani (2023). The \emph{PAGFL} jointly identifies the group structure and group-specific slope parameters.
#' The function supports both static and dynamic panels, with or without endogenous regressors.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, \bold{x}_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables that declare the cross-sectional unit \eqn{i} and the time period \eqn{t} of each observation.
#' @param index a character vector holding two strings. The first string denotes the name of the index variable identifying the cross-sectional unit \eqn{i} and the second string represents the name of the variable declaring the time period \eqn{t}. The data is automatically sorted according to the variables in \code{index}, which may produce errors when the time index is a character variable. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the number of time periods \code{n_periods} is supplied.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{NULL}.
#' @param lambda the tuning parameter determining the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type IC automatically selects the best fitting \eqn{\lambda} value.
#' @param method the estimation method. Options are
#' \describe{
#' \item{\code{"PLS"}}{for using the penalized least squares (\emph{PLS}) algorithm. We recommend \emph{PLS} in case of (weakly) exogenous regressors (Mehrabani, 2023, sec. 2.2).}
#' \item{\code{"PGMM"}}{for using the penalized Generalized Method of Moments (\emph{PGMM}). \emph{PGMM} is required when instrumenting endogenous regressors, in which case a matrix \eqn{\bold{Z}} containing the necessary exogenous instruments must be supplied (Mehrabani, 2023, sec. 2.3).}
#' } Default is \code{"PLS"}.
#' @param Z a \eqn{NT \times q} \code{matrix} or \code{data.frame} of exogenous instruments, where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1^\prime, \dots, z_N^\prime)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. \code{Z} is only required when \code{method = "PGMM"} is selected. When using \code{"PLS"}, the argument can be left empty or it is disregarded. Default is \code{NULL}.
#' @param min_group_frac the minimum group cardinality as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, each of its members is allocated to one of the remaining groups according to the \emph{MSE}. Default is 0.05.
#' @param bias_correc logical. If \code{TRUE}, a Split-panel Jackknife bias correction following Dhaene and Jochmans (2015) is applied to the slope parameters. We recommend using the correction when working with dynamic panels. Default is \code{FALSE}.
#' @param kappa the a non-negative weight used to obtain the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is \eqn{1*10^4}.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is \eqn{1*10^{-8}}.
#' @param tol_group the tolerance limit for within-group differences. Two individuals \eqn{i}, \eqn{j} are assigned to the same group if the Frobenius norm of their coefficient vector difference is below this threshold. Default is \eqn{1*10^{-3}}.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the IC that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For \emph{PLS}, the \eqn{\varrho} value is trivial. However, for \emph{PGMM}, small values lead to slow convergence. If left unspecified, the default heuristic \eqn{\varrho = \max(\frac{\sqrt{5NTp}}{\log(NTp)}-7, 1}) is used.
#' @param verbose logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, certain operations are parallelized across multiple cores. Default is \code{TRUE}.
#' @param ... ellipsis
#'
#' @details
#' Consider the panel data model
#' \deqn{y_{it} = \gamma_i^0 + \bold{\beta}^{0 \prime}_{i} \bold{x}_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i^0} is an individual fixed effect, \eqn{\bold{x}_{it}} is a \eqn{p \times 1} vector of weakly exogenous explanatory variables, and \eqn{\epsilon_{it}} is a zero mean error.
#' The coefficient vector \eqn{\bold{\beta}_i^0} follows the latent group pattern
#' \deqn{\bold{\beta}_i^0 = \sum_{k = 1}^K \bold{\alpha}_k^0 \bold{1} \{i \in G_k^0 \},}
#' with \eqn{\cup_{k = 1}^K G_k^0 = \{1, \dots, N\}}, \eqn{G_k^0 \cap G_j^0 = \emptyset} and \eqn{\| \bold{\alpha}_k^0 - \bold{\alpha}_j^0 \| \neq 0} for any \eqn{k \neq j},  \eqn{k,j = 1, \dots, K}.
#'
#' The \emph{PLS} method jointly estimates the latent group structure and group-specific coefficients by minimizing the criterion
#' \deqn{Q_{NT} (\bold{\beta}, \lambda) = \frac{1}{T} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \bold{\beta}^\prime_i \tilde{\bold{x}}_{it})^2 + \frac{\lambda}{N} \sum_{i = 1}^{N - 1} \sum_{j=i}^N \dot{\omega}_{ij} \| \bold{\beta}_i - \bold{\beta}_j \|}
#' with respect to \eqn{\bold{\beta} = (\bold{\beta}_1^\prime, \dots, \bold{\beta}_N^\prime)^\prime}. \eqn{\tilde{a}_{it} = a_{it} - T^{-1} \sum_{t = 1}^T a_{it}}, \eqn{a = \{y, \bold{x}\}} to concentrate out the individual fixed effects \eqn{\gamma_i^0} (within-transformation). \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{\omega}_{ij}} reflects adaptive penalty weights (see Mehrabani, 2023, eq. 2.6). \eqn{\| \cdot \|} denotes the Frobenius norm.
#' The adaptive weights \eqn{\dot{w}_{ij}} are obtained by a preliminary individual least squares estimation.
#' The criterion function is minimized via an iterative alternating direction method of multipliers (\emph{ADMM}) algorithm (see Mehrabani, 2023, sec. 5.1).
#'
#' \emph{PGMM} employs a set of instruments \code{Z} to control for endogenous regressors. Using \emph{PGMM}, \eqn{\bold{\beta}} is estimated by minimizing
#' \deqn{
#' Q_{NT}(\bold{\beta}, \lambda) = \sum^N_{i = 1} \left[ \frac{1}{N} \sum_{t=1}^T \bold{z}_{it} (\Delta y_{it} - \bold{\beta}^\prime_i \Delta \bold{x}_{it}) \right]^\prime \bold{W}_i \left[\frac{1}{T} \sum_{t=1}^T \bold{z}_{it}(\Delta y_{it} - \bold{\beta}^\prime_i \Delta \bold{x}_{it}) \right]
#' }
#' \deqn{
#' \quad + \frac{\lambda}{N} \sum_{i = 1}^{N - 1} \sum_{j=i+1}^N \ddot{\omega}_{ij} \| \bold{\beta}_i - \bold{\beta}_j \|.
#' }
#' \eqn{\ddot{\omega}_{ij}} are obtained by an initial \emph{GMM} estimation. \eqn{\Delta} gives the first differences operator \eqn{\Delta y_{it} = y_{it} - y_{i t-1}}. \eqn{\bold{W}_i} represents a data-driven \eqn{q \times q} weight matrix. I refer to Mehrabani (2023, eq. 2.10) for more details.
#' Again, the criterion function is minimized using an efficient \emph{ADMM} algorithm (Mehrabani, 2023, sec. 5.2).
#'
#' Two individuals are assigned to the same group if \eqn{\| \hat{\bold{\beta}}_i - \hat{\bold{\beta}}_j \| \leq \epsilon_{\text{tol}}} (and hence \eqn{\hat{\bold{\alpha}}_k = \hat{\bold{\beta}}_i = \hat{\bold{\beta}}_j} for some \eqn{k = 1, \dots, \hat{K}}), where \eqn{\epsilon_{\text{tol}}} is determined by \code{tol_group}. Subsequently, the estimated number of groups \eqn{\hat{K}} and group structure follows by examining the number of distinct elements in \eqn{\hat{\bold{\beta}}}. Given an estimated group structure, it is straightforward to obtain post-Lasso estimates using group-wise least squares or \emph{GMM} (see \code{\link{grouped_plm}}).
#'
#' We recommend identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit close to 0 and an upper limit that leads to a fully homogeneous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.
#'
#' @examples
#' # Simulate a panel with a group structure
#' set.seed(1)
#' sim <- sim_DGP(N = 20, n_periods = 80, p = 2, n_groups = 3)
#' y <- sim$y
#' X <- sim$X
#' df <- cbind(y = c(y), X)
#'
#' # Run the PAGFL procedure
#' estim <- pagfl(y ~ ., data = df, n_periods = 80, lambda = 0.5, method = "PLS")
#' summary(estim)
#'
#' # Lets pass a panel data set with explicit cross-sectional and time indicators
#' i_index <- rep(1:20, each = 80)
#' t_index <- rep(1:80, 20)
#' df <- data.frame(y = c(y), X, i_index = i_index, t_index = t_index)
#' estim <- pagfl(
#'   y ~ .,
#'   data = df, index = c("i_index", "t_index"), lambda = 0.5, method = "PLS"
#' )
#' summary(estim)
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#'
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @aliases PAGFL
#'
#' @return An object of class \code{pagfl} holding
#' \item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as cross-sectional and time indices,}
#' \item{\code{coefficients}}{a \eqn{\hat{K} \times p} matrix of the post-Lasso group-specific parameter estimates,}
#' \item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{\hat{K}} and (ii) a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k},}
#' \item{\code{residuals}}{a vector of residuals of the demeaned model,}
#' \item{\code{fitted}}{a vector of fitted values of the demeaned model,}
#' \item{\code{args}}{a \code{list} of additional arguments,}
#' \item{\code{IC}}{a \code{list} containing (i) the value of the IC, (ii) the employed tuning parameter \eqn{\lambda}, and (iii) the \emph{MSE},}
#' \item{\code{convergence}}{a \code{list} containing (i) a logical variable indicating if convergence was achieved and (ii) the number of executed \emph{ADMM} algorithm iterations,}
#' \item{\code{call}}{the function call.}
#'
#' A \code{pagfl} object has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual}, and \code{coef} S3 methods.
#' @export
pagfl <- function(formula, data, index = NULL, n_periods = NULL, lambda, method = "PLS", Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 1e4, tol_convergence = 1e-8,
                  tol_group = 1e-3, rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE, parallel = TRUE, ...) {
  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  # In case of penalized Least Squares, specify an empty instrument matrix Z
  if (method == "PGMM") {
    if (is.null(Z)) stop("PGMM requires a matrix of exogenous instruments `Z`\n")
  }
  if (method == "PLS") {
    Z <- matrix()
  } else {
    Z <- as.matrix(Z)
  }

  if (is.null(min_group_frac)) min_group_frac <- 0
  formula <- stats::as.formula(formula)
  method <- match.arg(method, c("PLS", "PGMM"))
  prelim_checks(formula, data, Z, index, n_periods, method,
    verbose = verbose, min_group_frac = min_group_frac, max_iter = max_iter, kappa = kappa,
    tol_group = tol_group, tol_convergence = tol_convergence
  )
  # Construct a vector of the dependent variable and a regressor matrix
  # If present, remove the intercept
  data <- as.data.frame(data)
  if (any(all.vars(formula[[3]]) == ".")) {
    data <- stats::na.omit(data)
  } else {
    data <- data[stats::complete.cases(data[, c(index, all.vars(formula))]), ]
  }
  if (!is.null(index)) data <- data[order(data[, index[1]], data[, index[2]]), ]
  if (all(all.vars(formula[[3]]) != ".")) {
    # If present, remove the intercept
    formula <- stats::update(formula, . ~ . - 1)
    X <- stats::model.matrix(formula, data)
    regressor_names <- colnames(X)
  } else {
    # Remove the index variables if present
    if (is.null(index)) {
      data_temp <- data
    } else {
      data_temp <- data[, !(colnames(data) %in% index)]
    }
    X <- stats::model.matrix(formula, data_temp)
    regressor_names <- colnames(X)
    # If present, remove the intercept
    intercept_index <- apply(X, 2, stats::sd) != 0
    X <- X[, intercept_index]
    regressor_names <- regressor_names[intercept_index]
    rm(data_temp, intercept_index)
  }
  X <- as.matrix(X)
  y <- as.matrix(data[[all.vars(formula[[2]])]])
  # Build the final data set
  model_data <- as.data.frame(cbind(c(y), X))
  colnames(model_data)[1] <- all.vars(formula[[2]])

  # Extract or produce index variables
  if (!is.null(index)) {
    i_index_labs <- data[, index[1]]
    i_index <- as.integer(factor(i_index_labs))
    t_index_labs <- data[, index[2]]
    t_index <- as.integer(factor(t_index_labs))
    n_periods <- length(unique(t_index))
    N <- length(unique(i_index))
    model_data <- cbind(model_data, data[, index[1]], data[, index[2]])
    colnames(model_data)[(ncol(model_data) - 1):ncol(model_data)] <- index
  } else {
    N <- NROW(y) / n_periods
    if (round(N) != N) stop("Either `n_periods` is incorrect or an unbalanced panel data set is passed without supplying `index`\n")
    t_index <- t_index_labs <- rep(1:n_periods, N)
    i_index <- i_index_labs <- rep(1:N, each = n_periods)
    model_data$i_index <- i_index
    model_data$t_index <- t_index
    index <- c("i_index", "t_index")
  }
  coef_rownames <- as.character(unique(t_index_labs)[order(unique(t_index))])
  p <- ncol(X)

  second_checks(y = y, X = X, method = method, Z = Z, p = p, min_group_frac = min_group_frac, verbose = verbose, dyn = FALSE, rho = rho, varrho = varrho)

  #------------------------------#
  #### Iterate over lambda    ####
  #------------------------------#

  # Run the algorithm
  lambdaList <- pagfl_routine(
    y = y, X = X, method = method, Z = Z, bias_correc = bias_correc, i_index = i_index,
    t_index = t_index, N = N, lambda_vec = lambda, kappa = kappa, min_group_frac = min_group_frac,
    max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group, varrho = varrho,
    rho = rho, parallel = parallel, verbose = verbose
  )
  # Estimate fixed effects
  fe_vec <- getFE(y = y, i_index = i_index, N = N, method = method)

  #------------------------------#
  #### Prepare the output     ####
  #------------------------------#

  # Pick the estimation result with the lowest IC
  indx <- which.min(lapply(lambdaList, function(x) x$IC$IC))
  out <- lambdaList[[indx]]
  # Add FE to the estimated values
  fitted <- c(out$IC$fitted + fe_vec)
  # Attach names to the output
  out$estimOutput$groups_hat <- c(out$estimOutput$groups_hat)
  names(out$estimOutput$groups_hat) <- unique(i_index_labs)
  rownames(out$estimOutput$alpha_hat) <- paste("Group", 1:out$estimOutput$K_hat)
  colnames(out$estimOutput$alpha_hat) <- regressor_names
  # List with additional arguments
  args <- list(
    formula = formula, labs = list(i = i_index_labs, t = t_index_labs, index = index), method = method, min_group_frac = min_group_frac,
    bias_correc = bias_correc, kappa = kappa, rho = rho, max_iter = max_iter, tol_group = tol_group, varrho = varrho
  )
  # Output
  out <- list(
    model = as.data.frame(data), coefficients = out$estimOutput$alpha_hat, groups = list(n_groups = out$estimOutput$K_hat, groups = out$estimOutput$groups_hat), residuals = c(out$IC$resid), fitted = fitted,
    args = args, IC = list(IC = out$IC$IC, lambda = lambda[indx], msr = out$IC$msr), convergence = list(convergence = out$estimOutput$convergence, iter = out$estimOutput$iter), call = match.call()
  )
  class(out) <- "pagfl"
  return(out)
}

#' @export
PAGFL <- function(y, X, n_periods, lambda, method = "PLS", Z = NULL, min_group_frac = .05, bias_correc = FALSE, kappa = 2, max_iter = 2e3, tol_convergence = 1e-3,
                  tol_group = sqrt(p / (sqrt(N * n_periods) * log(log(N * n_periods)))), rho = .07 * log(N * n_periods) / sqrt(N * n_periods),
                  varrho = max(sqrt(5 * N * n_periods * p) / log(N * n_periods * p) - 7, 1), verbose = TRUE) {
  p <- ncol(X)
  N <- length(c(y)) / n_periods
  lifecycle::deprecate_stop(when = "1.1.0", what = "PAGFL()", with = "pagfl()", env = asNamespace("PAGFL"))
}
