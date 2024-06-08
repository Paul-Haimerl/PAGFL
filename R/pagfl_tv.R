#' Time-varying Pairwise Adaptive Group Fused Lasso
#'
#' @description The time-varying pairwise adaptive group fused lasso (time-varying \emph{PAGFL}) jointly estimates the latent group structure and group-specific time-varying functional coefficients in a panel data model.
#' The time-varying coefficients are modeled as polynomial B-splines.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, x_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables, declaring the cross-sectional unit \eqn{i} and the time period \eqn{t} for each observation.
#' @param index a character vector holding two strings specifying the variable names that identify the cross-sectional unit and the time period for each observation. The first string denotes the individual unit, while the second string represents the time period. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied. Default is \code{Null}.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{Null}.
#' @param lambda the tuning parameter. \eqn{\lambda} governs the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type IC automatically selects the best fitting parameter value.
#' @param d the polynomial degree of the B-splines. Default is 3.
#' @param M the number of interior knots of the B-splines. If left unspecified, the default heuristic \eqn{M = \text{floor}((NT)^{\frac{1}{7}} - \log(p))} is used. Note that \eqn{M} does not include the boundary knots.
#' @param const_coef a character vector containing the variable names of explanatory variables that are estimated with time-constant coefficients. All of concerning regressors must be named variables in \code{data}.
#' @param min_group_frac the minimum group size as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, a hierarchical classifier allocates its members to the remaining groups. Default is 0.05.
#' @param kappa the a non-negative weight placed on the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is 20,000.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is \eqn{1 * 10^{-10}}.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are assigned to the same group if the Frobenius norm of their coefficient vector difference is below this threshold. Default is 0.001.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the IC that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param varrho the non-negative Lagrangian \emph{ADMM} penalty parameter. For the employed penalized sieve estimation \emph{PSE}, the \eqn{\varrho} value is trivial. We recommend the default 1.
#' @param verbose logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, certain operations are parallelized across multiple cores.
#' @param ... ellipsis
#'
#' @details
#' Consider the grouped time-varying panel data model
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{i} (t/T) x_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect, \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables, and \eqn{\epsilon_{it}} is a zero mean error.
#' The coefficient vector \eqn{\beta_{i} (t/T)} is subject to the group pattern
#' \deqn{\beta_i \left(\frac{t}{T} \right) = \sum_{k = 1}^K \alpha_k \left( \frac{t}{T} \right) \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k \| \neq \| \alpha_j \|} for any \eqn{k \neq M}.

#' \eqn{\beta_i (t/T)}, and \eqn{\alpha_k (t/T)} are estimated as polynomial B-splines using penalized sieve-technique. Let \eqn{\bold{B}(v)} denote a \eqn{M + d +1} vector basis functions, where \eqn{d} denotes the polynomial degree and \eqn{M} the number of interior knots.
#' Then, \eqn{\beta_{i}(t/T)} and \eqn{\alpha_{i}(t/T)} are approximated as \eqn{\beta_{i} (t/T) = \pi_i^\prime \bold{B}(t/T)} and \eqn{\alpha_{i}(t/T) = \xi_i^\prime \bold{B}(t/T)}, respectively. \eqn{\pi_i} and \eqn{\xi_i} are \eqn{(M + d + 1) \times p} coefficient matrices which weigh the individual basis functions.
#' The explanatory variables are projected onto the spline basis system, which results in the \eqn{(M + d + 1)*p \times 1} vector \eqn{z_{it} = x_{it} \otimes \bold{B}(v)}. Subsequently, the DGP can be reformulated as
#' \deqn{y_{it} = \gamma_i + z_{it}^\prime \text{vec}(\pi_{i}) + u_{it},}
#' where \eqn{u_{it} = \epsilon_{it} + \eta_{it}} and \eqn{\eta_{it}} contains the sieve approximation error. I refer to Su et al. (2019, sec. 2) for more details on the sieve technique.
#'
#' Inspired by Su et al. (2019) and Mehrabani (2023), the time-varying PAGFL estimates the functional coefficients and the group structure by minimizing the criterion:
#' \deqn{\frac{1}{NT} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \tilde{z}_{it}^\prime \text{vec}(\pi_{i}))^2 + \frac{\lambda}{N} \sum_{1 \leq i} \sum_{i<j \leq N} \dot{w}_{ij} \| \text{vec}( \pi_i - \pi_j) \|,}
#' where \eqn{\tilde{y}_{it}} is the demeaned dependent variable, and \eqn{\tilde{z}_{it}} is likewise demeaned to concentrate out the individual fixed effects \eqn{\gamma_i}. \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} denotes adaptive penalty weights which are obtained by a preliminary non-penalized estimation. \eqn{\| \cdot \|} represents the Frobenius norm.
#' The solution \eqn{\hat{\bold{\beta}}} is computed via the iterative alternating direction method of multipliers (\emph{ADMM}) algorithm proposed in Mehrabani (2023, sec. 5.1), adapted to accommodate the B-spline coefficient functions.
#'
#' Two individuals are assigned to the same group if \eqn{\| \text{vec} (\hat{\pi}_i - \hat{\pi}_j) \| \leq \epsilon_{\text{tol}}}, where \eqn{\epsilon_{\text{tol}}} is given by \code{tol_group}. Subsequently, the number of groups follows as the number of distinct elements in \eqn{\hat{\bold{\beta}}}. Given an estimated group structure, it is straightforward to obtain post-Lasso estimates using least squares.
#'
#' We suggest identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit of 0 and an upper limit that leads to a fully homogeneous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.
#'
#' In case of an unbalanced panel data set, the earliest and latest available observations out of the entire panel are employed as the start and end-points of the interval on which the time-varying coefficients are defined.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' set.seed(1)
#' sim <- sim_tv_DGP(N = 5, n_periods = 20, intercept = TRUE, p = 1)
#' df <- data.frame(y = c(sim$y))
#'
#' # Run the time-varying PAGFL with only an intercept
#' estim <- tv_pagfl(y ~ 1, data = df, n_periods = 20, lambda = 13, max_iter = 100, parallel = FALSE)
#' summary(estim)
#'
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#'
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. *Journal of Business & Economic Statistics*, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
#'
#' @author Paul Haimerl
#'
#' @return An object of class \code{tvpagfl} holding
#' \item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as individual and time indices,}
#' \item{\code{coefficients}}{a \code{list} holding (i) a \eqn{T \times p^{(1)} \times \hat{K}} array of the post-Lasso group-specific functional coefficients and (ii) a \eqn{K \times p^{(2)}} matrix of time-constant post-Lasso estimates. Let \eqn{p^{(1)}} denote the number of time-varying coefficients and \eqn{p^{(2)}} the number of time constant parameters,}
#' \item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{\hat{K}} and (ii) a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k},}
#' \item{\code{residuals}}{a vector of residuals of the demeaned model,}
#' \item{\code{fitted}}{a vector of fitted values of the demeaned model,}
#' \item{\code{args}}{a \code{list} of additional arguments,}
#' \item{\code{IC}}{a \code{list} containing (i) the value of the IC, (ii) the employed tuning parameter \eqn{\lambda}, and (iii) the mean squared error,}
#' \item{\code{convergence}}{a \code{list} containing (i) a logical variable if convergence was achieved and (ii) the number of executed \emph{ADMM} algorithm iterations,}
#' \item{\code{call}}{the function call.}
#'
#' An object of class \code{tvpagfl} has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual} and \code{coef} S3 methods.
#' @export
tv_pagfl <- function(formula, data, index = NULL, n_periods = NULL, lambda, d = 3, M = floor(length(y)^(1 / 7) - log(p)), min_group_frac = .05,
                     const_coef = NULL, kappa = 2, max_iter = 20e3, tol_convergence = 1e-10, tol_group = 1e-3,
                     rho = .07 * log(N * n_periods) / sqrt(N * n_periods), varrho = 1, verbose = TRUE, parallel = TRUE, ...) {
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
  if (any(all.vars(formula[[3]]) == ".")){
    if (!is.null(index)){
      data_temp <- data[,!(colnames(data) %in% index)]
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
    t_index <- t_index_labs <- rep(1:n_periods, N)
    i_index <- i_index_labs <- rep(1:N, each = n_periods)
    model_data$i_index <- i_index
    model_data$t_index <- t_index
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
    N = N, index = index, n_periods = n_periods, y = y, X = X, p = p, min_group_frac = min_group_frac,
    verbose = verbose, dyn = TRUE, d = d, M = M, rho = rho, varrho = varrho
  )

  #------------------------------#
  #### Iterate over lambda    ####
  #------------------------------#

  # Run the algorithm
  lambdaList <- tv_pagfl_routine(
    y = y, X = X, X_const = X_const, d = d, M = M, i_index = i_index, t_index = t_index, N = N, p_const = p_const, lambda_vec = lambda, kappa = kappa,
    min_group_frac = min_group_frac, max_iter = max_iter, tol_convergence = tol_convergence, tol_group = tol_group,
    varrho = varrho, rho = rho, parallel = parallel
  )
  # Estimate fixed effects
  fe_vec <- getFE(y = y, i_index = i_index, N = N, method = "PLS")

  #------------------------------#
  #### Prepare the output     ####
  #------------------------------#

  # Pick the estimation result with the lowest IC
  indx <- which.min(lapply(lambdaList, function(x) x$IC_list$IC))
  out <- lambdaList[[indx]]
  # Add FE to the estimated values
  fitted <- c(out$IC$fitted + fe_vec)
  # Prep the functional coefficients
  coef_raw <- out$estimOutput$alpha_hat
  # In case of time-constant regressors, separate them from the spline coefficients
  if (p_const > 0) {
    indx_p_c <- (ncol(coef_raw) - p_const + 1):ncol(coef_raw)
    alpha_const <- coef_raw[, indx_p_c]
    alpha_hat <- coef_raw[, -indx_p_c]
    if (out$estimOutput$K_hat == 1) {
      alpha_hat <- t(alpha_hat)
      alpha_const <- t(alpha_const)
    }
    if (p_const == 1) {
      alpha_const <- as.matrix(alpha_const)
    }
  } else {
    alpha_hat <- coef_raw
    alpha_const <- NULL
  }
  B <- bspline_system(1:n_periods, d, seq(1, n_periods, length.out = M + 2), TRUE)

  # Transform the spline coefficients to time-varying functional coefficients
  alpha_mat <- getTVAlpha(xi = alpha_hat, K_hat = out$estimOutput$K_hat, p = p, n_periods = n_periods, B = B)
  # Attach names to the output
  dimnames(alpha_mat) <- list(coef_rownames, regressor_names, paste0("Group ", 1:out$estimOutput$K_hat))
  if (!is.null(const_coef)) {
    dimnames(alpha_const) <- list(paste0("Group ", 1:out$estimOutput$K_hat), const_coef)
  }
  out$estimOutput$groups_hat <- c(out$estimOutput$groups_hat)
  names(out$estimOutput$groups_hat) <- unique(i_index_labs)
  # Omit estimates at the beginning and end of the observational period without support
  if (!is.null(index)) {
    alpha_hat_star <- delete_missing_t(
      i_index = i_index, t_index = t_index, K_hat = out$estimOutput$K_hat, groups_hat = out$estimOutput$groups_hat,
      alpha_hat = alpha_mat
    )
    dimnames(alpha_hat_star) <- dimnames(alpha_mat)
    alpha_mat <- alpha_hat_star
  }
  # De-mean a trend, if present
  if ("(Intercept)" %in% regressor_names) {
    trend_vec <- c(alpha_mat[, which(regressor_names == "(Intercept)"), ])
    trend_vec <- demeanIndVec(trend_vec, N = out$estimOutput$K_hat, i_index = rep(1:out$estimOutput$K_hat, each = n_periods))
    alpha_mat[, which(regressor_names == "(Intercept)"), ] <- array(trend_vec, dim = c(n_periods, 1, out$estimOutput$K_hat))
  }
  # List with additional arguments
  args <- list(
    formula = formula, labs = list(i = i_index_labs, t = t_index_labs), d = d, M = M, min_group_frac = min_group_frac,
    kappa = kappa, max_iter = max_iter, tol_group = tol_group, varrho = varrho
  )

  # Output
  out <- list(
    model = as.data.frame(model_data), coefficients = list(tv = alpha_mat, const = alpha_const), groups = list(n_groups = out$estimOutput$K_hat, groups = out$estimOutput$groups_hat), residuals = c(out$IC$resid), fitted = fitted,
    args = args, IC = list(IC = out$IC$IC, lambda = lambda[indx], msr = out$IC$msr), convergence = list(convergence = out$estimOutput$convergence, iter = out$estimOutput$iter), call = match.call()
  )
  class(out) <- "tvpagfl"
  return(out)
}
