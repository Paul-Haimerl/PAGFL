#' Fused Unobserved group Spline Estimation of TIME varying coefficients
#'
#' @description Estimate a time-varying panel data model subject to a latent group structure using \emph{FUSE-TIME}--Fused Unobserved group Spline Estimation of TIME varying coefficients--by Haimerl et al. (2025). \emph{FUSE-TIME} jointly identifies the latent group structure and group-specific time-varying functional coefficients.
#' The time-varying coefficients are approximated as polynomial B-splines. The function supports both static and dynamic panel data models.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, \bold{x}_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables that declare the cross-sectional unit \eqn{i} and the time period \eqn{t} of each observation.
#' @param index a character vector holding two strings. The first string denotes the name of the index variable identifying the cross-sectional unit \eqn{i} and the second string represents the name of the variable declaring the time period \eqn{t}. The data is automatically sorted according to the variables in \code{index}, which may produce errors when the time index is a character variable. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the number of time periods \code{n_periods} is supplied.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{NULL}.
#' @param lambda the tuning parameter determining the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type IC automatically selects the best fitting \eqn{\lambda} value.
#' @param d the polynomial degree of the B-splines. Default is 3.
#' @param M the number of interior knots of the B-splines. If left unspecified, the default heuristic \eqn{M = \text{floor}((NT)^{\frac{1}{7}} - \log(p))} following Haimerl et al. (2025) is used.
#' @param const_coef a character vector containing the variable names of explanatory variables that enter with time-constant coefficients.
#' @param min_group_frac the minimum group cardinality as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, each of its members is allocated to one of the remaining groups according to the \emph{MSE}. Default is 0.05.
#' @param kappa the a non-negative weight used to obtain the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is \eqn{5*10^4}.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is \eqn{1*10^{-10}}.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are assigned to the same group if the Frobenius norm of their coefficient vector difference is below this threshold. Default is \eqn{1*10^{-3}}.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the IC that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Haimerl et al. (2025) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For the employed penalized sieve estimation \emph{PSE}, the \eqn{\varrho} value is not very influential. We recommend the default 1.
#' @param verbose logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, certain operations are parallelized across multiple cores. Default is \code{TRUE}.
#' @param ... ellipsis
#'
#' @details
#' Consider the grouped time-varying panel data model subject to a latent group structure
#' \deqn{y_{it} = \gamma_i^0 + \bold{\beta}^{0\prime}_{i} (t/T) \bold{x}_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i^0} is an individual fixed effect, \eqn{\bold{x}_{it}} is a \eqn{p \times 1} vector of explanatory variables, and \eqn{\epsilon_{it}} denotes a zero mean error.
#' The \eqn{p}-dimensional coefficient vector \eqn{\bold{\beta}_{i}^0 (t/T)} contains smooth functions of time and follows the latent group pattern
#' \deqn{\bold{\beta}_i^0 \left(\frac{t}{T} \right) = \sum_{k = 1}^K \bold{\alpha}_k^0 \left( \frac{t}{T} \right) \bold{1} \{i \in G_k^0 \},}
#' with \eqn{\cup_{k = 1}^K G_k^0 = \{1, \dots, N\}}, \eqn{G_k^0 \cap G_j^0 = \emptyset} for any \eqn{k \neq j},  \eqn{k,j = 1, \dots, K}.
#'
#' The time-varying coefficient functions are estimated as polynomial B-splines. To this end, let \eqn{\bold{b}(v)} denote a \eqn{M + d +1} vector of polynomial basis functions with the polynomial degree \eqn{d} and \eqn{M} interior knots.
#' Then, \eqn{\bold{\beta}_i^0 (t/T)} is approximated by forming linear combinations of these basis functions \eqn{\bold{\beta}_i^0 (t/T) \approx \bold{\Pi}_i^{0 \prime} \bold{b} (t/T)}, where \eqn{\bold{\Pi}_i^{0}} is a \eqn{(M + d + 1) \times p} matrix of spline control points.
#'
#' To estimate \eqn{\bold{\Pi}_i^{0}}, we project the explanatory variables onto the spline basis system, resulting in the \eqn{(M + d + 1)p \times 1} regressor vector \eqn{\bold{z}_{it} = \bold{x}_{it} \otimes \bold{b}(v)}. Subsequently, the DGP can be reformulated as
#' \deqn{y_{it} = \gamma_i^0 + \bold{\pi}_{i}^{0 \prime} \bold{z}_{it} + u_{it},}
#' where \eqn{\bold{\pi}_i^0 = \text{vec}(\bold{\Pi}_i^0)}, and \eqn{u_{it} = \epsilon_{it} + \eta_{it}} collects the idiosyncratic \eqn{\epsilon_{it}} and the sieve approximation error \eqn{\eta_{it}}.
#'
#' Following Haimerl et al. (2025, sec. 2), \emph{FUSE-TIME} jointly estimates the functional coefficients and the group structure by minimizing the criterion
#' \deqn{F_{NT} (\bold{\pi}, \lambda) = \frac{1}{NT} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \bold{\pi}_{i}^\prime \tilde{\bold{z}}_{it})^2 + \frac{\lambda}{N} \sum_{i = 1}^{N - 1} \sum_{j = i+1}^N \dot{\omega}_{ij} \| \bold{\pi}_i - \bold{\pi}_j \|_2}
#' with respect to \eqn{\bold{\pi} = (\bold{\pi}_i^\prime, \dots, \bold{\pi}_N^\prime)^\prime}. \eqn{\tilde{a}_{it} = a_{it} - T^{-1} \sum^{T}_{t=1} a_{it}}, \eqn{a = \{y, \bold{z}\}} to concentrate out the individual fixed effects \eqn{\gamma_i^0} (within-transformation). \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} denotes adaptive penalty weights which are obtained by a preliminary non-penalized estimation.
#' The criterion function is minimized via an iterative alternating direction method of multipliers (\emph{ADMM}) algorithm (Haimerl et al. 2053, Appendix C).
#'
#' Two individuals are assigned to the same group if \eqn{\| \hat{\bold{\pi}}_i - \hat{\bold{\pi}}_j \|_2 \leq \epsilon_{\text{tol}}} (and hence \eqn{\hat{\bold{\xi}}_k = \hat{\bold{\pi}}_i = \hat{\bold{\pi}}_j} for some \eqn{k = 1, \dots, \hat{K}}), where \eqn{\epsilon_{\text{tol}}} is determined by \code{tol_group}. The time-varying coefficients are then retrieved by taking \eqn{\hat{\bold{\beta}}_i (t/T) = \hat{\bold{\Pi}}_i^\prime \bold{b}(t/T)}, where \eqn{\hat{\bold{\pi}}_i = \text{vec}(\hat{\bold{\Pi}}_i)} (analogously \eqn{\hat{\bold{\alpha}}_k (t/T) = \hat{\bold{\Xi}}_k^\prime \bold{b}(t/T)}, using \eqn{\hat{\bold{\xi}}_k = \text{vec}(\hat{\bold{\Xi}}_k)}).
#'
#' Subsequently, the estimated number of groups \eqn{\hat{K}} and group structure follow by examining the number of distinct elements in \eqn{\hat{\bold{\pi}}}. Given an estimated group structure, it is straightforward to obtain post-Lasso estimates \eqn{\hat{\bold{\alpha}}^p_k (t/T) = \hat{\bold{\Xi}}^{p \prime}_k \bold{b}(t/T)} for each \eqn{k = 1, \dots, \hat{K} } using group-wise least squares (see \code{\link{grouped_tv_plm}}).
#'
#' We recommend choosing a \eqn{\lambda} tuning parameter by passing a logarithmically spaced grid of candidate values with a lower limit close to 0 and an upper limit that leads to a fully homogeneous panel. A BIC-type information criterion then automatically selects the best fitting \eqn{\lambda} value.
#'
#' In case of an unbalanced panel data set, the earliest and latest available observations per group define the start and end-points of the interval on which the group-specific time-varying coefficients are defined.
#'
#' We refer to Haimerl et al. (2025) for more details.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' set.seed(1)
#' sim <- sim_tv_DGP(N = 10, n_periods = 50, intercept = TRUE, p = 1)
#' df <- data.frame(y = c(sim$y))
#'
#' # Run FUSE-TIME
#' estim <- fuse_time(y ~ ., data = df, n_periods = 50, lambda = 10, parallel = FALSE)
#' summary(estim)
#'
#' @references
#' Haimerl, P., Smeekes, S., & Wilms, I. (2025). Estimation of latent group structures in time-varying panel data models. *arXiv preprint arXiv:2503.23165*. \doi{10.48550/arXiv.2503.23165}.
#'
#' @author Paul Haimerl
#'
#' @return An object of class \code{fusetime} holding
#' \item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as cross-sectional and time indices,}
#' \item{\code{coefficients}}{let \eqn{p^{(1)}} denote the number of time-varying coefficients and \eqn{p^{(2)}} the number of time constant parameters. A \code{list} holding (i) a \eqn{T \times p^{(1)} \times \hat{K}} array of the post-Lasso group-specific functional coefficients and (ii) a \eqn{K \times p^{(2)}} matrix of time-constant post-Lasso estimates.}
#' \item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{\hat{K}} and (ii) a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k},}
#' \item{\code{residuals}}{a vector of residuals of the demeaned model,}
#' \item{\code{fitted}}{a vector of fitted values of the demeaned model,}
#' \item{\code{args}}{a \code{list} of additional arguments,}
#' \item{\code{IC}}{a \code{list} containing (i) the value of the IC, (ii) the employed tuning parameter \eqn{\lambda}, and (iii) the \emph{MSE},}
#' \item{\code{convergence}}{a \code{list} containing (i) a logical variable if convergence was achieved and (ii) the number of executed \emph{ADMM} algorithm iterations,}
#' \item{\code{call}}{the function call.}
#'
#' An object of class \code{fusetime} has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual}, and \code{coef} S3 methods.
#' @export
fuse_time <- function(formula, data, index = NULL, n_periods = NULL, lambda, d = 3, M = floor(length(y)^(1 / 7) - log(p)), min_group_frac = .05,
                      const_coef = NULL, kappa = 2, max_iter = 5e4, tol_convergence = 1e-10, tol_group = 1e-3,
                      rho = .04 * log(N * n_periods) / sqrt(N * n_periods), varrho = 1, verbose = TRUE, parallel = TRUE, ...) {
  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  if (is.null(min_group_frac)) min_group_frac <- 0
  formula <- stats::as.formula(formula)
  prelim_checks(formula, data,
    index = index, n_periods = n_periods, const_coef = const_coef, verbose = verbose, min_group_frac = min_group_frac,
    max_iter = max_iter, kappa = kappa, tol_group = tol_group, tol_convergence = tol_convergence
  )
  # Construct a vector of the dependent variable and a regressor matrix
  data <- as.data.frame(data)
  if (any(all.vars(formula[[3]]) == ".")) {
    data <- stats::na.omit(data)
  } else {
    data <- data[stats::complete.cases(data[, c(index, all.vars(formula))]), ]
  }
  if (!is.null(index)) data <- data[order(data[, index[1]], data[, index[2]]), ]
  if (any(all.vars(formula[[3]]) == ".")) {
    if (!is.null(index)) {
      data_temp <- data[, !(colnames(data) %in% index)]
      if (NCOL(data) == 3) {
        data_temp <- as.data.frame(data_temp)
        colnames(data_temp) <- colnames(data)[!(colnames(data) %in% index)]
      }
    } else {
      data_temp <- data
    }
    X <- stats::model.matrix(formula, data_temp)
    rm(data_temp)
    intercept_indicator <- "1"
  } else {
    X <- stats::model.matrix(formula, data)
    intercept_indicator <- ifelse(attr(stats::terms(formula), "intercept") == 1, "1", "-1")
  }
  # Remove second intercepts if present
  sd_vec <- apply(as.matrix(X[, !(colnames(X) %in% "(Intercept)")]), 2, stats::sd)
  sd_vec_ind <- sd_vec != 0
  if (intercept_indicator == "1") sd_vec_ind <- c(TRUE, sd_vec_ind)
  regressor_names <- colnames(X)[sd_vec_ind]
  # Pull the regressors
  X <- as.matrix(X[, sd_vec_ind])
  colnames(X) <- regressor_names
  # Pull the outcome
  y <- as.matrix(data[[all.vars(formula[[2]])]])
  # Build the final data set
  model_data <- as.data.frame(cbind(c(y), X))
  colnames(model_data)[1] <- all.vars(formula[[2]])
  # Check constant coefficients
  if (!all(const_coef %in% regressor_names)) {
    stop(paste(const_coef[!const_coef %in% regressor_names], collapse = ", "), " not part of the formula\n")
  }
  regressor_names <- regressor_names[!(regressor_names %in% const_coef)]

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
  coef_t_index <- unique(t_index_labs)[order(unique(t_index))]
  coef_rownames <- as.character(coef_t_index)
  if (is.character(coef_t_index)) coef_t_index <- as.numeric(factor(coef_t_index))

  # Prepare constant coefficients
  if (!is.null(const_coef)) {
    X_const <- as.matrix(X[, const_coef])
    X <- as.matrix(X[, !(colnames(X) %in% const_coef)])
    p_const <- ncol(X_const)
  } else {
    X_const <- matrix()
    p_const <- 0
  }
  p <- ncol(X)
  M <- floor(max(M, 1))
  d <- floor(d)

  second_checks(
    y = y, X = X, p = p, min_group_frac = min_group_frac,
    verbose = verbose, dyn = TRUE, d = d, M = M, rho = rho, varrho = varrho
  )

  #------------------------------#
  #### Iterate over lambda    ####
  #------------------------------#

  # Run the algorithm
  lambdaList <- fuse_time_routine(
    y = y, X = X, X_const = X_const, d = d, M = M, i_index = i_index, t_index = t_index, N = N, p_const = p_const, lambda_vec = lambda, kappa = kappa,
    min_group_frac = min_group_frac, max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group,
    varrho = varrho, rho = rho, parallel = parallel, verbose = verbose
  )
  # Estimate fixed effects
  fe_vec <- getFE(y = y, i_index = i_index, N = N, method = "PLS")

  # Pick the estimation result with the lowest IC
  indx <- which.min(lapply(lambdaList, function(x) x$IC$IC))
  out <- lambdaList[[indx]]

  #------------------------------#
  #### Output                 ####
  #------------------------------#

  out <- fuse_time_output(
    out = out, fe_vec = fe_vec, p = p, p_const = p_const, n_periods = n_periods, d = d, M = M, coef_rownames = coef_rownames,
    regressor_names = regressor_names, const_coef = const_coef, index = index, i_index_labs = i_index_labs, i_index = i_index,
    t_index = t_index, model_data = model_data
  )
  out$args <- list(
    formula = formula, labs = list(i = i_index_labs, t = t_index_labs, index = index), d = d, M = M, min_group_frac = min_group_frac,
    kappa = kappa, rho = rho, max_iter = max_iter, tol_group = tol_group, varrho = varrho
  )
  out$IC$lambda <- lambda[indx]
  out$call <- match.call()

  class(out) <- "fusetime"
  return(out)
}


#' @export
#' @rdname fuse_time
tv_pagfl <- function(formula, data, index = NULL, n_periods = NULL, lambda, d = 3, M, min_group_frac = .05,
                     const_coef = NULL, kappa = 2, max_iter = 5e4, tol_convergence = 1e-10, tol_group = 1e-3,
                     rho, varrho = 1, verbose = TRUE, parallel = TRUE, ...) {
  lifecycle::deprecate_warn(when = "1.1.4", what = "tv_pagfl()", with = "fuse_time()")
  fuse_time(
    formula = formula, data = data, index = index, n_periods = n_periods, lambda = lambda, d = d, M = M, min_group_frac = min_group_frac,
    const_coef = const_coef, kappa = kappa, max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group,
    rho = rho, varrho = varrho, verbose = verbose, parallel = parallel, ...
  )
}
