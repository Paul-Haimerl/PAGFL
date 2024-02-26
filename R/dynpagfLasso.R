#' Time-varying Pairwise Adaptive Group Fused Lasso
#'
#' @description The time-varying pairwise adaptive group fused lasso (PAGFL) jointly estimates the latent group structure and group-specific time-varying functional coefficients in a panel data model.
#' The time-varying coefficients are modeled as polynomial B-splines.
#'
#' @param y a \eqn{NT \times 1} vector or data.frame of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.
#' @param X a \eqn{NT \times p} matrix or data.frame of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.
#' @param n_periods the number of observed time periods \eqn{T}.
#' @param lambda the tuning parameter governing the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type information criterion selects the best fitting parameter value.
#' @param d the degree of the polynomial coefficient spline functions. Default is 3.
#' @param J The number of interior knots of the spline functions. If left unspecified, the default heuristic \eqn{\text{floor}((NT)^{\frac{1}{7}})} is used.
#' @param min_group_frac the minimum group size as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, a hierarchical classifier allocates its members to the remaining groups. Default is 0.05.
#' @param kappa the weight placed on the adaptive penalty weights. Default is 2.
#' @param max_iter the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is 2000.
#' @param tol_convergence the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is 0.001.
#' @param tol_group the tolerance limit for within-group differences. Two individuals are placed in the same group if the Frobenius norm of their coefficient parameter difference is below this parameter. If left unspecified, the heuristic \eqn{\sqrt{\frac{p (J + d + 1)}{\sqrt{NT} \log(\log(NT))}}} is used. We recommend the default.
#' @param rho the tuning parameter balancing the fitness and penalty terms in the information criterion that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.
#' @param verbose logical. If \code{TRUE}, a progression bar is printed when iterating over candidate \eqn{\lambda} values and helpful warning messages are shown. Default is \code{TRUE}.
#'
#' @details
#' Consider the time-varying panel data model
#' \deqn{y_{it} = \gamma_i + \beta^\prime_{it} x_{it} + u_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
#' where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect, and \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables.
#' The coefficient vector \eqn{\beta_i = \{\beta_{i1}^\prime, \dots, \beta_{iT}^\prime \}^\prime} is subject to the group pattern
#' \deqn{\beta_i = \sum_{k = 1}^K \alpha_k \bold{1} \{i \in G_k \},}
#' with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k \| \neq \| \alpha_j \|} for any \eqn{k \neq j}.

#' \eqn{\beta_{it}=\beta_i (\frac{t}{T})}, and \eqn{\alpha_{kt}=\alpha_k (\frac{t}{T})} are estimated as polynomial B-spline functions via penalized sieve estimation (\emph{PSE}). Let \eqn{\bold{B}(v)} denote a system of B-spline basis functions of degree \eqn{d} and \eqn{J} interior knots.
#' Then, \eqn{\beta_{it}} and \eqn{\alpha_{it}} may be approximated as \eqn{\beta_{it}=\pi_i^\prime \bold{B}(v)} and \eqn{\alpha_{it}=\xi_i^\prime \bold{B}(v)}, respectively. \eqn{\pi_i} and \eqn{\xi_i} are \eqn{(J + d + 1) \times p} coefficient matrices which map the fitted functional coefficients on the B-splines.
#' The explanatory variables are projected onto the spline basis system, which results in the \eqn{(J + d + 1)p \times 1} vector \eqn{z_{it} = x_{it} \otimes \bold{B}(v)}. Subsequently, the DGP can be reformulated as
#' \deqn{y_{it} = \gamma_i + z_{it}^\prime \text{vec}(\pi_{i}) + e_{it}.}
#' I refer to Su et al. (2019) for more details on the sieve technique.
#'
#' Inspired by Su et al. (2019) and Mehrabani (2023), the time-varying PAGFL estimates the functional coefficients and the group structure by minimizing the following criterion:
#' \deqn{\frac{1}{NT} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \tilde{z}_{it}^\prime \text{vec}(\pi_{i}))^2 + \frac{\lambda}{N} \sum_{1 \leq i} \sum_{i<j \leq N} \dot{w}_{ij} \| \text{vec}( \pi_i - \pi_j) \|,}
#' where \eqn{\tilde{y}_{it}} is the demeaned dependent variable, and \eqn{\tilde{z}_{it}} is likewise demeaned to concentrate out the individual fixed effects \eqn{\gamma_i}. \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} reflects adaptive penalty weights which are obtained by a preliminary non-penalized estimation. \eqn{\| \cdot \|} denotes the Frobenius norm.
#' The solution \eqn{\hat{\bold{\beta}}} is computed via the iterative alternating direction method of multipliers (\emph{ADMM}) algorithm proposed in Mehrabani (2023, sec. 5.1), adapted to accommodate the B-spline coefficient functions.
#'
#' Two individuals are assigned to the same group if \eqn{\| \text{vec} (\hat{\pi}_i - \hat{\pi}_j) \| \leq \epsilon_{\text{tol}}}, where \eqn{\epsilon_{\text{tol}}} is given by \code{tol_group}.
#'
#' We suggest identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit of 0 and an upper limit that leads to a fully homogenous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' sim <- sim_dyn_DGP(N = 50, n_periods = 50, DGP = 1)
#' y <- sim$y
#' X <- sim$X
#'
#' # Run the time-varying PAGFL
#' estim <- dyn_pagfl(y = y, X = X, n_periods = 50, lambda = 3.5)
#' print(estim)
#' @references
#' Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of fixed-effect models. *The Review of Economic Studies*, 82(3), 991-1030. \doi{10.1093/restud/rdv007}.
#'
#' Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. *Journal of Econometrics*, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.
#'
#' Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. *Journal of Business & Economic Statistics*, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
#'
#' @author Paul Haimerl
#'
#' @return A list holding
#' \item{\code{IC}}{the BIC-type information criterion.}
#' \item{\code{lambda}}{the penalization parameter. If multiple \eqn{\lambda} values were passed, the parameter yielding the lowest IC.}
#' \item{\code{alpha_hat}}{a \eqn{T \times p \times \widehat{K}} array of the post-Lasso group-specific parameter estimates.}
#' \item{\code{K_hat}}{the estimated total number of groups \eqn{\widehat{K}}.}
#' \item{\code{groups_hat}}{a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k}.}
#' \item{\code{iter}}{the number of executed algorithm iterations.}
#' \item{\code{convergence}}{logical. If \code{TRUE}, convergence was achieved. If \code{FALSE}, \code{max_iter} was reached.}
#' @export
dyn_pagfl <- function(y, X, n_periods, lambda, d = 3, J = floor(NROW(y)^(1 / 7)), min_group_frac = .05,
                      kappa = 2, max_iter = 2e3, tol_convergence = 1e-6,
                      tol_group = sqrt(p * (J + d + 1) / (sqrt(N * n_periods) * log(log(N * n_periods)))),
                      rho = .07 * log(N * n_periods) / sqrt(N * n_periods), verbose = TRUE) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  N <- nrow(y) / n_periods
  p <- ncol(X)

  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  checks(
    N = N, n_periods = n_periods, y = y, X = X, p = p, min_group_frac = min_group_frac,
    verbose = verbose, dyn = TRUE, d = d, J = J
  )

  #------------------------------#
  #### Build the spline basis ####
  #------------------------------#

  inner_knots <- seq(1, n_periods, length.out = J)
  B <- splines::bs(1:n_periods, knots = inner_knots, degree = d, intercept = F)
  Z <- buildZ(X = X, B = B, n_periods = n_periods, J = J, d = d, p = p, N = N)

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
    estimOutput <- dyn_pagfl_algo(
      y = y, Z = Z, B = B, d = d, J = J, n_periods = n_periods, lambda = lam, kappa = kappa,
      min_group_frac = min_group_frac, max_iter = max_iter, tol_convergence = tol_convergence,
      tol_group = tol_group, varrho = 1
    )
    # Compute the Information Criterion
    IC_val <- IC(estimOutput = estimOutput, y = y, X = Z, rho = rho, method = "PLS", n_periods = n_periods, N = N)
    # Transform the spline coefficients to time-varying functional coefficients
    estimOutput[[1]] <- getDynAlpha(xi = estimOutput[["alpha_hat"]], K_hat = estimOutput[["K_hat"]], p = p, n_periods = n_periods, B)
    return(c(IC = IC_val, lambda = lam, estimOutput))
  })
  # Pick the estimation result with the lowest IC
  IC_vec <- lapply(lambdaList, function(x) x[[1]])
  estim <- lambdaList[[which.min(IC_vec)]]
  estim$groups_hat <- c(estim$groups_hat)
  if (is.null(colnames(X))){
    column_names <- paste0("alpha_", 1:p)
  } else {
    column_names <- colnames(X)
  }
  dimnames(estim$alpha_hat) <- list(paste0("t=", 1:n_periods), column_names, paste0("k=", 1:estim$K_hat))
  return(estim)
}

#' Grouped Time-varying Panel Data Oracle Estimator
#'
#' @description Estimation of a time-varying panel data model, given a known group structure.
#'
#' @param y a \eqn{NT \times 1} vector or data.frame of the dependent variable, with \eqn{\bold{y}=(y_1, \dots, y_N)^\prime}, \eqn{y_i = (y_{i1}, \dots, y_{iT})^\prime} and the scalar \eqn{y_{it}}.
#' @param X a \eqn{NT \times p} matrix or data.frame of explanatory variables, with \eqn{\bold{X}=(x_1, \dots, x_N)^\prime}, \eqn{x_i = (x_{i1}, \dots, x_{iT})^\prime} and the \eqn{p \times 1} vector \eqn{x_{it}}.
#' @param n_periods the number of observed time periods \eqn{T}.
#' @param groups_0 the true group vector \eqn{(g_1, \dots, g_N)}, where \eqn{g_i = k} if the unit \eqn{i} is part of group \eqn{k}.
#' @param d the degree of the coefficient polynomial spline functions. Default is 3.
#' @param J The number of interior knots of the spline functions. If left unspecified, the default heuristic \eqn{\text{floor}((NT)^{\frac{1}{7}})} is used.
#'
#' @details
#' Given a known group pattern \eqn{\bold{G}^0}, the least squares oracle estimator solves
#' \deqn{\text{vec}(\hat{\xi}_k^{\text{oracle}}) = \left( \sum_{i \in G_k^0} \sum_{t = 1}^T \tilde{z}_{it} \tilde{z}_{it}^\prime \right)^{-1} \sum_{i \in G_k^0} \sum_{t = 1}^T \tilde{z}_{it} \tilde{y}_{it},}
#' where \eqn{\tilde{y}_{it}} is the demeaned dependent variable and \eqn{\tilde{z}_{it}} is a \eqn{(J + d + 1) \times 1} vector that interacts the explanatory variables \eqn{x_{it}} with the spline basis system of degree \eqn{d} and \eqn{J} interior knots \eqn{\bold{B}(v)}, demeaned.
#' The final time-varying coefficients are obtained by \eqn{\hat{\alpha}_{kt}^{\text{oracle}} = \hat{\xi}_k^{\text{oracle}^\prime} \bold{B} \left( \frac{t}{T} \right)}.
#'
#' @examples
#' # Simulate a time-varying panel with a trend and a group pattern
#' sim <- sim_dyn_DGP(N = 20, n_periods = 50, DGP = 1)
#' y <- sim$y
#' X <- sim$X
#' groups_0 <- sim$groups
#'
#' # Estimate the time-varying grouped panel data model model, given the true grouping
#' estim <- dyn_oracle(y = y, X = X, n_periods = 50, groups_0 = groups_0)
#' print(estim)
#' @author Paul Haimerl
#'
#' @return A \eqn{T \times p \times K} array of the group-specific time-varying coefficients.
#'
#' @export
dyn_oracle <- function(y, X, n_periods, groups_0, d = 3, J = floor(NROW(y)^(1 / 7))) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  N <- nrow(y) / n_periods
  p <- ncol(X)

  #------------------------------#
  #### Checks                 ####
  #------------------------------#

  checks(N = N, n_periods = n_periods, y = y, X = X, p = p, min_group_frac = NULL, verbose = FALSE, dyn = TRUE, d = d, J = J)

  #------------------------------#
  #### Build the spline basis ####
  #------------------------------#

  knots <- seq(1, n_periods, length.out = J)
  B <- splines::bs(1:n_periods, knots = knots, degree = d, intercept = FALSE)
  Z <- buildZ(X = X, B = B, n_periods = n_periods, J = J, d = d, p = p, N = N)

  #------------------------------#
  #### Preliminaries          ####
  #------------------------------#

  J_star <- J + d
  p_star <- p * J_star
  # Net out fixed effects
  data <- netFE(y = y, X = Z, method = "PLS", n_periods = n_periods, N = N)
  y_tilde <- data[[1]]
  Z_tilde <- matrix(data[[2]], ncol = ncol(Z))

  #------------------------------#
  #### Oracle estimation      ####
  #------------------------------#

  xi_mat_vec <- getAlpha(Z_tilde, y_tilde, matrix(), "PLS", n_periods, N, p_star, groups_0)
  alpha_hat <- getDynAlpha(xi = xi_mat_vec, K_hat = max(groups_0), p = p, n_periods = n_periods, B)
  dimnames(alpha_hat) <- list(paste0("t=", 1:n_periods), paste0("alpha_", 1:p), paste0("k=", 1:max(groups_0)))
  return(alpha_hat)
}