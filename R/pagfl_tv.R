#' Time-varying Pairwise Adaptive Group Fused Lasso
#'
#' @description Estimate a time-varying panel data model with a latent group structure using the pairwise adaptive group fused lasso (\emph{time-varying PAGFL}). The \emph{time-varying PAGFL} jointly identifies the latent group structure and group-specific time-varying functional coefficients.
#' The time-varying coefficients are modeled as polynomial B-splines. The function supports both static and dynamic panel data models.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, x_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables that declare the cross-sectional unit \eqn{i} and the time period \eqn{t} of each observation.
#' @param index a character vector holding two strings. The first string denotes the name of the index variable identifying the cross-sectional unit \eqn{i} and the second string represents the name of the variable declaring the time period \eqn{t}. The data is automatically sorted according to the variables in \code{index}, which may produce errors when the time index is a character variable. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{Null}.
#' @param lambda the tuning parameter determining the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type IC automatically selects the best fitting \eqn{\lambda} value.
#' @param d the polynomial degree of the B-splines. Default is 3.
#' @param M the number of interior knots of the B-splines. If left unspecified, the default heuristic \eqn{M = \text{floor}((NT)^{\frac{1}{7}} - \log(p))} is used. Note that \eqn{M} does not include the boundary knots and the entire sequence of knots is of length \eqn{M + d + 1}.
#' @param const_coef a character vector containing the variable names of explanatory variables that enter with time-constant coefficients.
#' @param min_group_frac the minimum group cardinality as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, each of its members is allocated to one of the remaining groups according to the \emph{MSE}. Default is 0.05.
#' @param kappa the a non-negative weight used to obtain the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is \eqn{5*10^4}.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is \eqn{1*10^{-10}}.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are assigned to the same group if the Frobenius norm of their coefficient vector difference is below this threshold. Default is \eqn{1*10^{-3}}.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the IC that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For the employed penalized sieve estimation \emph{PSE}, the \eqn{\varrho} value is trivial. We recommend the default 1.
#' @param verbose logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, certain operations are parallelized across multiple cores. Default is \code{TRUE}.
#' @param ... ellipsis
#'
#' @details
#' Consider the grouped time-varying panel data model
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{i} (t/T) x_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect, \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables, and \eqn{\epsilon_{it}} is a zero mean error.
#' The coefficient vector \eqn{\beta_{i} (t/T)} is subject to the latent group pattern
#' \deqn{\beta_i \left(\frac{t}{T} \right) = \sum_{k = 1}^K \alpha_k \left( \frac{t}{T} \right) \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k - \alpha_j \| \neq 0} for any \eqn{k \neq j}, \eqn{k = 1, \dots, K}.
#'
#' The time-varying coefficient functions are estimated as polynomial B-splines using the penalized sieve-technique. To this end, let \eqn{B(v)} denote a \eqn{M + d +1} vector basis functions, where \eqn{d} denotes the polynomial degree and \eqn{M} the number of interior knots.
#' Then, \eqn{\beta_{i}(t/T)} and \eqn{\alpha_{k}(t/T)} are approximated by forming linear combinations of the basis functions \eqn{\beta_{i} (t/T) \approx \pi_i^\prime B(t/T)} and \eqn{\alpha_{i}(t/T) \approx \xi_k^\prime B(t/T)}, where \eqn{\pi_i} and \eqn{\xi_i} are \eqn{(M + d + 1) \times p} coefficient matrices.
#'
#' The explanatory variables are projected onto the spline basis system, which results in the \eqn{(M + d + 1)p \times 1} vector \eqn{z_{it} = x_{it} \otimes B(v)}. Subsequently, the DGP can be reformulated as
#' \deqn{y_{it} = \gamma_i + z_{it}^\prime \text{vec}(\pi_{i}) + u_{it},}
#' where \eqn{u_{it} = \epsilon_{it} + \eta_{it}} and \eqn{\eta_{it}} reflects a sieve approximation error. We refer to Su et al. (2019, sec. 2) for more details on the sieve technique.
#'
#' Inspired by Su et al. (2019) and Mehrabani (2023), the time-varying PAGFL jointly estimates the functional coefficients and the group structure by minimizing the criterion
#' \deqn{Q_{NT} (\bold{\pi}, \lambda) = \frac{1}{NT} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \tilde{z}_{it}^\prime \text{vec}(\pi_{i}))^2 + \frac{\lambda}{N} \sum_{i = 1}^{N - 1} \sum_{j > i}^N \dot{\omega}_{ij} \| \pi_i - \pi_j \|}
#' with respect to \eqn{\bold{\pi} = (\text{vec}(\pi_i)^\prime, \dots, \text{vec}(\pi_N)^\prime)^\prime}. \eqn{\tilde{a}_{it} = a_{it} - T^{-1} \sum^{T}_{t=1} a_{it}}, \eqn{a = \{y, z\}} to concentrate out the individual fixed effects \eqn{\gamma_i}. \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} denotes adaptive penalty weights which are obtained by a preliminary non-penalized estimation. \eqn{\| \cdot \|} represents the Frobenius norm.
#' The solution criterion function is minimized via the iterative alternating direction method of multipliers (\emph{ADMM}) algorithm proposed by Mehrabani (2023, sec. 5.1).
#'
#' Two individuals are assigned to the same group if \eqn{\| \text{vec} (\hat{\pi}_i - \hat{\pi}_j) \| \leq \epsilon_{\text{tol}}}, where \eqn{\epsilon_{\text{tol}}} is determined by \code{tol_group}. Subsequently, the number of groups follows as the number of distinct elements in \eqn{\hat{\bold{\pi}}}. Given an estimated group structure, it is straightforward to obtain post-Lasso estimates \eqn{\hat{\bold{\xi}}} using group-wise least squares (see \code{\link{grouped_tv_plm}}).
#'
#' We recommend identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit close to 0 and an upper limit that leads to a fully homogeneous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.
#'
#' In case of an unbalanced panel data set, the earliest and latest available observations per group define the start and end-points of the interval on which the group-specific time-varying coefficients are defined.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' set.seed(1)
#' sim <- sim_tv_DGP(N = 10, n_periods = 50, intercept = TRUE, p = 1)
#' df <- data.frame(y = c(sim$y))
#'
#' # Run the time-varying PAGFL
#' estim <- tv_pagfl(y ~ ., data = df, n_periods = 50, lambda = 10, parallel = FALSE)
#' summary(estim)
#'
#' @references
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. *Journal of Business & Economic Statistics*, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
#'
#' @author Paul Haimerl
#'
#' @return An object of class \code{tvpagfl} holding
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
#' An object of class \code{tvpagfl} has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual} and \code{coef} S3 methods.
#' @export
tv_pagfl <- function(formula, data, index = NULL, n_periods = NULL, lambda, d = 3, M = floor(length(y)^(1 / 7) - log(p)), min_group_frac = .05,
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
  data <- stats::na.omit(data)
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
  lambdaList <- tv_pagfl_routine(
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

  out <- tv_pagfl_output(
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

  class(out) <- "tvpagfl"
  return(out)
}
