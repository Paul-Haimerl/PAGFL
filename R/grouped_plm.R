#' Grouped Panel Data Model
#'
#' @description Estimate a grouped panel data model given an observed group structure. Slope parameters are homogeneous within groups but heterogeneous across groups.
#' This function supports both static and dynamic panel data models, with or without endogenous regressors.
#'
#' @param formula a formula object describing the model to be estimated.
#' @param data a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, x_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables that declare the cross-sectional unit \eqn{i} and the time period \eqn{t} of each observation.
#' @param groups a numerical or character vector of length \eqn{N} that indicates the group membership of each cross-sectional unit \eqn{i}.
#' @param index a character vector holding two strings. The first string denotes the name of the index variable identifying the cross-sectional unit \eqn{i} and the second string represents the name of the variable declaring the time period \eqn{t}. The data is automatically sorted according to the variables in \code{index}, which may produce errors when the time index is a character variable. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied.
#' @param n_periods the number of observed time periods \eqn{T}. If an \code{index} is passed, this argument can be left empty.
#' @param method the estimation method. Options are
#' \describe{
#' \item{\code{"PLS"}}{for using the penalized least squares (\emph{PLS}) algorithm. We recommend \emph{PLS} in case of (weakly) exogenous regressors (Mehrabani, 2023, sec. 2.2).}
#' \item{\code{"PGMM"}}{for using the penalized Generalized Method of Moments (\emph{PGMM}). \emph{PGMM} is required when instrumenting endogenous regressors, in which case a matrix \eqn{\bold{Z}} containing the necessary exogenous instruments must be supplied (Mehrabani, 2023, sec. 2.3).}
#' } Default is \code{"PLS"}.
#' @param Z a \eqn{NT \times q} \code{matrix} or \code{data.frame} of exogenous instruments, where \eqn{q \geq p}, \eqn{\bold{Z}=(z_1, \dots, z_N)^\prime}, \eqn{z_i = (z_{i1}, \dots, z_{iT})^\prime} and \eqn{z_{it}} is a \eqn{q \times 1} vector. \code{Z} is only required when \code{method = "PGMM"} is selected. When using \code{"PLS"}, the argument can be left empty or it is disregarded. Default is \code{NULL}.
#' @param bias_correc logical. If \code{TRUE}, a Split-panel Jackknife bias correction following Dhaene and Jochmans (2015) is applied to the slope parameters. We recommend using the correction when working with dynamic panels. Default is \code{FALSE}.
#' @param rho a tuning parameter balancing the fitness and penalty terms in the IC. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param verbose logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.
#' @param parallel logical. If \code{TRUE}, certain operations are parallelized across multiple cores. Default is \code{TRUE}.
#' @param ... ellipsis
#'
#' @details
#' Consider the grouped panel data model
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{i} x_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect, \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables, and \eqn{\epsilon_{it}} is a zero mean error.
#' The coefficient vector \eqn{\beta_i} is subject to the observed group pattern
#' \deqn{\beta_i = \sum_{k = 1}^K \alpha_k \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k - \alpha_j \| \neq 0} for any \eqn{k \neq j}, \eqn{k = 1, \dots, K}.
#'
#' Using \emph{PLS}, the group-specific coefficients for group \eqn{k} are obtained via \emph{OLS}
#' \deqn{\hat{\alpha}_k = \left( \sum_{i \in G_k} \sum_{t = 1}^T \tilde{x}_{it} \tilde{x}_{it}^\prime \right)^{-1} \sum_{i \in G_k} \sum_{t = 1}^T \tilde{x}_{it} \tilde{y}_{it},}
#' where \eqn{\tilde{a}_{it} = a_{it} - T^{-1} \sum_{t=1}^T a_{it}}, \eqn{a = \{y, x\}} to concentrate out the individual fixed effects \eqn{\gamma_i} (within-transformation).
#'
#' In case of \emph{PGMM}, the slope coefficients are derived as
#' \deqn{
#' \hat{\alpha}_k = \left( \left[ \sum_{i \in G_k} T^{-1} \sum_{t = 1}^T z_{it} \Delta x_{it} \right]^\prime W_k \left[ \sum_{i \in G_k} T^{-1} \sum_{t = 1}^T z_{it} \Delta x_{it} \right] \right)^{-1}
#' }
#' \deqn{
#' \quad \quad \left[ \sum_{i \in G_k} T^{-1} \sum_{t = 1}^T z_{it} \Delta x_{it} \right]^\prime W_k \left[ \sum_{i \in G_k} T^{-1} \sum_{t = 1}^T z_{it} \Delta y_{it} \right],
#' }
#' where \eqn{W_k} is a \eqn{q \times q} p.d. symmetric weight matrix and \eqn{\Delta} denotes the first difference operator \eqn{\Delta x_{it} = x_{it} - x_{it-1}} (first-difference transformation).
#'
#' @examples
#' # Simulate a panel with a group structure
#' set.seed(1)
#' sim <- sim_DGP(N = 20, n_periods = 80, p = 2, n_groups = 3)
#' y <- sim$y
#' X <- sim$X
#' groups <- sim$groups
#' df <- cbind(y = c(y), X)
#'
#' # Estimate the grouped panel data model
#' estim <- grouped_plm(y ~ ., data = df, groups = groups, n_periods = 80, method = "PLS")
#' summary(estim)
#'
#' # Lets pass a panel data set with explicit cross-sectional and time indicators
#' i_index <- rep(1:20, each = 80)
#' t_index <- rep(1:80, 20)
#' df <- data.frame(y = c(y), X, i_index = i_index, t_index = t_index)
#' estim <- grouped_plm(
#'   y ~ .,
#'   data = df, index = c("i_index", "t_index"), groups = groups, method = "PLS"
#' )
#' summary(estim)
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' @author Paul Haimerl
#'
#' @return An object of class \code{gplm} holding
#' \item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as cross-sectional and time indices,}
#' \item{\code{coefficients}}{a \eqn{K \times p} matrix of the group-specific parameter estimates,}
#' \item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{K} and (ii) a vector of group memberships \eqn{g_1, \dots, g_N)}, where \eqn{g_i = k} if \eqn{i} is assigned to group \eqn{k},}
#' \item{\code{residuals}}{a vector of residuals of the demeaned model,}
#' \item{\code{fitted}}{a vector of fitted values of the demeaned model,}
#' \item{\code{args}}{a \code{list} of additional arguments,}
#' \item{\code{IC}}{a \code{list} containing (i) the value of the IC and (ii) the \emph{MSE},}
#' \item{\code{call}}{the function call.}
#'
#' A \code{gplm} object has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual}, and \code{coef} S3 methods.
#' @export
grouped_plm <- function(formula, data, groups, index = NULL, n_periods = NULL, method = "PLS", Z = NULL, bias_correc = FALSE,
                        rho = .07 * log(N * n_periods) / sqrt(N * n_periods), verbose = TRUE, parallel = TRUE, ...) {
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

  formula <- stats::as.formula(formula)
  method <- match.arg(method, c("PLS", "PGMM"))
  prelim_checks(formula, data, Z, index, n_periods, method, verbose = verbose)
  # Construct a vector of the dependent variable and a regressor matrix
  # If present, remove the intercept
  data <- as.data.frame(data)
  data <- stats::na.omit(data)
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

  if (length(groups) != N) stop("`groups` must have the length N")
  second_checks(y = y, X = X, method = method, Z = Z, p = p, verbose = verbose, dyn = FALSE, rho = rho)

  # Get the group structure
  n_groups <- length(unique(groups))
  groups_index <- as.integer(factor(groups))

  #------------------------------#
  #### Estimate               ####
  #------------------------------#

  # Run the algorithm
  out <- pagfl_oracle_routine(
    y = y, X = X, groups = groups, method = method, Z = Z, bias_correc = bias_correc, i_index = i_index,
    t_index = t_index, N = N, rho = rho, parallel = parallel
  )
  # Estimate fixed effects
  fe_vec <- getFE(y = y, i_index = i_index, N = N, method = method)

  #------------------------------#
  #### Output                 ####
  #------------------------------#

  # Add FE to the estimated values
  fitted <- c(out$IC$fitted + fe_vec)
  # Attach names to the output
  names(groups) <- unique(i_index_labs)
  rownames(out$alpha_hat) <- paste("Group", 1:n_groups)
  colnames(out$alpha_hat) <- regressor_names
  # List with additional arguments
  args <- list(
    formula = formula, labs = list(i = i_index_labs, t = t_index_labs, index = index), method = method, bias_correc = bias_correc, rho = rho
  )
  # Output
  out <- list(
    model = as.data.frame(model_data), coefficients = out$alpha_hat, groups = list(n_groups = n_groups, groups = groups), residuals = c(out$IC$resid), fitted = fitted,
    args = args, IC = list(IC = out$IC$IC, msr = out$IC$msr), call = match.call()
  )
  class(out) <- "gplm"
  return(out)
}
