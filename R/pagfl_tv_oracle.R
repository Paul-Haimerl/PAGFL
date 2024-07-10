#' Grouped Time-varying Panel Data Model
#'
#' @description Estimate a time-varying panel data model given an exogenous group structure. Coefficient functions are homogeneous within, but heterogeneous across groups.
#' The time-varying coefficients are modeled as polynomial B-splines.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, x_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables, declaring the cross-sectional unit \eqn{i} and the time period \eqn{t} for each observation.
#' @param groups a numerical or character vector of length \eqn{N} that indicates the group membership.
#' @param index a character vector holding two strings specifying the variable names that identify the cross-sectional unit and the time period for each observation. The first string denotes the individual unit, while the second string represents the time period. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied. Default is \code{Null}.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{Null}.
#' @param d the polynomial degree of the B-splines. Default is 3.
#' @param M the number of interior knots of the B-splines. If left unspecified, the default heuristic \eqn{M = \text{floor}((NT)^{\frac{1}{7}} - \log(p))} is used. Note that \eqn{M} does not include the boundary knots.
#' @param const_coef a character vector containing the variable names of explanatory variables that are estimated with time-constant coefficients. All of concerning regressors must be named variables in \code{data}.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the IC. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
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
#'
#' \eqn{\beta_i (t/T)}, and \eqn{\alpha_k (t/T)} are estimated as polynomial B-splines using penalized sieve-technique. Let \eqn{\bold{B}(v)} denote a \eqn{M + d +1} vector basis functions, where \eqn{d} denotes the polynomial degree and \eqn{M} the number of interior knots.
#' Then, \eqn{\beta_{i}(t/T)} and \eqn{\alpha_{i}(t/T)} are approximated as \eqn{\beta_{i} (t/T) = \pi_i^\prime \bold{B}(t/T)} and \eqn{\alpha_{i}(t/T) = \xi_i^\prime \bold{B}(t/T)}, respectively. \eqn{\pi_i} and \eqn{\xi_i} are \eqn{(M + d + 1) \times p} coefficient matrices which weigh the individual basis functions.
#' The explanatory variables are projected onto the spline basis system, which results in the \eqn{(M + d + 1)*p \times 1} vector \eqn{z_{it} = x_{it} \otimes \bold{B}(v)}. Subsequently, the DGP can be reformulated as
#' \deqn{y_{it} = \gamma_i + z_{it}^\prime \text{vec}(\pi_{i}) + u_{it},}
#' where \eqn{u_{it} = \epsilon_{it} + \eta_{it}} and \eqn{\eta_{it}} contains the sieve approximation error. We refer to Su et al. (2019, sec. 2) for more details on the sieve technique.
#'
#' In case of an unbalanced panel data set, the earliest and latest available observations out of the entire panel are employed as the start and end-points of the interval on which the time-varying coefficients are defined.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' set.seed(1)
#' sim <- sim_tv_DGP(N = 5, n_periods = 20, intercept = TRUE, p = 1)
#' df <- data.frame(y = c(sim$y))
#' groups <- sim$groups
#'
#' # Estimate the time-varying grouped panel data model
#' estim <- tv_pagfl_oracle(y ~ 1, data = df, n_periods = 20, groups = groups)
#' summary(estim)
#'
#' @references
#' Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. *Journal of Business & Economic Statistics*, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
#'
#' @author Paul Haimerl
#'
#' @return An object of class \code{tvpagfl_oracle} holding
#' \item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as individual and time indices,}
#' \item{\code{coefficients}}{a \code{list} holding (i) a \eqn{T \times p^{(1)} \times \hat{K}} array of the post-Lasso group-specific functional coefficients and (ii) a \eqn{K \times p^{(2)}} matrix of time-constant post-Lasso estimates. Let \eqn{p^{(1)}} denote the number of time-varying coefficients and \eqn{p^{(2)}} the number of time constant parameters,}
#' \item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{\hat{K}} and (ii) a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is part of group \eqn{k},}
#' \item{\code{residuals}}{a vector of residuals of the demeaned model,}
#' \item{\code{fitted}}{a vector of fitted values of the demeaned model,}
#' \item{\code{args}}{a \code{list} of additional arguments,}
#' \item{\code{IC}}{a \code{list} containing (i) the value of the IC, and (ii) the mean squared error,}
#' \item{\code{call}}{the function call.}
#'
#' An object of class \code{tvpagfl_oracle} has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual} and \code{coef} S3 methods.
#' @export
tv_pagfl_oracle <- function(formula, data, groups, index = NULL, n_periods = NULL, d = 3, M = floor(length(y)^(1 / 7) - log(p)),
                            const_coef = NULL, rho = .04 * log(N * n_periods) / sqrt(N * n_periods), verbose = TRUE, parallel = TRUE, ...) {
  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  formula <- stats::as.formula(formula)
  prelim_checks(formula, data, index = index, n_periods = n_periods, const_coef = const_coef, verbose = verbose)
  # Construct a vector of the dependent variable and a regressor matrix
  data <- as.data.frame(data)
  data <- stats::na.omit(data)
  if (!is.null(index)) data <- data[order(data[, index[1]], data[, index[2]]), ]
  if (any(all.vars(formula[[3]]) == ".")) {
    if (!is.null(index)) {
      data_temp <- data[, !(colnames(data) %in% index)]
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

  second_checks(y = y, X = X, p = p, verbose = verbose, dyn = TRUE, d = d, M = M, rho = rho)
  if (length(groups) != N) stop("`groups` must have the length N")

  # Get the group structure
  n_groups <- length(unique(groups))
  groups_index <- as.integer(factor(groups))

  #------------------------------#
  #### Estimate               ####
  #------------------------------#

  out <- tv_pagfl_oracle_routine(
    y = y, X = X, X_const = X_const, groups = groups_index, d = d, M = M, i_index = i_index, t_index = t_index, N = N, p_const = p_const, rho = rho, parallel = parallel
  )
  # Estimate fixed effects
  fe_vec <- getFE(y = y, i_index = i_index, N = N, method = "PLS")

  #------------------------------#
  #### Output                 ####
  #------------------------------#

  out <- tv_pagfl_output(
    out = out, fe_vec = fe_vec, p = p, p_const = p_const, n_periods = n_periods, d = d, M = M, coef_rownames = coef_rownames,
    regressor_names = regressor_names, const_coef = const_coef, index = index, i_index_labs = i_index_labs, i_index = i_index,
    t_index = t_index, model_data = model_data
  )
  out <- out[names(out) != "convergence"]
  out$IC <- out$IC[names(out$IC) != "lambda"]
  names(groups) <- unique(i_index_labs)
  out$groups$groups <- groups
  # List with additional arguments
  out$args <- list(
    formula = formula, labs = list(i = i_index_labs, t = t_index_labs), d = d, M = M, rho = rho
  )
  out$call <- match.call()
  class(out) <- "tvpagfl_oracle"
  return(out)
}
