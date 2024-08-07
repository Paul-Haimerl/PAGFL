% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pagfl_tv.R, R/s3_methods_tv.R
\name{tv_pagfl}
\alias{tv_pagfl}
\alias{summary.tvpagfl}
\alias{formula.tvpagfl}
\alias{df.residual.tvpagfl}
\alias{print.tvpagfl}
\alias{coef.tvpagfl}
\alias{residuals.tvpagfl}
\alias{fitted.tvpagfl}
\title{Time-varying Pairwise Adaptive Group Fused Lasso}
\usage{
tv_pagfl(
  formula,
  data,
  index = NULL,
  n_periods = NULL,
  lambda,
  d = 3,
  M = floor(length(y)^(1/7) - log(p)),
  min_group_frac = 0.05,
  const_coef = NULL,
  kappa = 2,
  max_iter = 50000,
  tol_convergence = 1e-10,
  tol_group = 0.001,
  rho = 0.04 * log(N * n_periods)/sqrt(N * n_periods),
  varrho = 1,
  verbose = TRUE,
  parallel = TRUE,
  ...
)

\method{summary}{tvpagfl}(object, ...)

\method{formula}{tvpagfl}(x, ...)

\method{df.residual}{tvpagfl}(object, ...)

\method{print}{tvpagfl}(x, ...)

\method{coef}{tvpagfl}(object, ...)

\method{residuals}{tvpagfl}(object, ...)

\method{fitted}{tvpagfl}(object, ...)
}
\arguments{
\item{formula}{a formula object describing the model to be estimated.}

\item{data}{a \code{data.frame} or \code{matrix} holding a panel data set. If no \code{index} variables are provided, the panel must be balanced and ordered in the long format \eqn{\bold{Y}=(Y_1^\prime, \dots, Y_N^\prime)^\prime}, \eqn{Y_i = (Y_{i1}, \dots, Y_{iT})^\prime} with \eqn{Y_{it} = (y_{it}, x_{it}^\prime)^\prime}. Conversely, if \code{data} is not ordered or not balanced, \code{data} must include two index variables that declare the cross-sectional unit \eqn{i} and the time period \eqn{t} of each observation.}

\item{index}{a character vector holding two strings. The first string denotes the name of the index variable identifying the cross-sectional unit \eqn{i}, and the second string represents the name of the variable declaring the time period \eqn{t}. In case of a balanced panel data set that is ordered in the long format, \code{index} can be left empty if the the number of time periods \code{n_periods} is supplied.}

\item{n_periods}{the number of observed time periods \eqn{T}. If an \code{index} character vector is passed, this argument can be left empty. Default is \code{Null}.}

\item{lambda}{the tuning parameter determining the strength of the penalty term. Either a single \eqn{\lambda} or a vector of candidate values can be passed. If a vector is supplied, a BIC-type IC automatically selects the best fitting \eqn{\lambda} value.}

\item{d}{the polynomial degree of the B-splines. Default is 3.}

\item{M}{the number of interior knots of the B-splines. If left unspecified, the default heuristic \eqn{M = \text{floor}((NT)^{\frac{1}{7}} - \log(p))} is used. Note that \eqn{M} does not include the boundary knots and the entire sequence of knots is of length \eqn{M + d + 1}.}

\item{min_group_frac}{the minimum group cardinality as a fraction of the total number of individuals \eqn{N}. In case a group falls short of this threshold, each of its members is allocated to one of the remaining groups according to the \emph{MSE}. Default is 0.05.}

\item{const_coef}{a character vector containing the variable names of explanatory variables that enter with time-constant coefficients.}

\item{kappa}{the a non-negative weight used to obtain the adaptive penalty weights. Default is 2.}

\item{max_iter}{the maximum number of iterations for the \emph{ADMM} estimation algorithm. Default is \eqn{5*10^4}.}

\item{tol_convergence}{the tolerance limit for the stopping criterion of the iterative \emph{ADMM} estimation algorithm. Default is \eqn{1*10^{-10}}.}

\item{tol_group}{the tolerance limit for within-group differences. Two individuals are assigned to the same group if the Frobenius norm of their coefficient vector difference is below this threshold. Default is \eqn{1*10^{-3}}.}

\item{rho}{the tuning parameter balancing the fitness and penalty terms in the IC that determines the penalty parameter \eqn{\lambda}. If left unspecified, the heuristic \eqn{\rho = 0.07 \frac{\log(NT)}{\sqrt{NT}}} of Mehrabani (2023, sec. 6) is used. We recommend the default.}

\item{varrho}{the non-negative Lagrangian \emph{ADMM} penalty parameter. For the employed penalized sieve estimation \emph{PSE}, the \eqn{\varrho} value is trivial. We recommend the default 1.}

\item{verbose}{logical. If \code{TRUE}, helpful warning messages are shown. Default is \code{TRUE}.}

\item{parallel}{logical. If \code{TRUE}, certain operations are parallelized across multiple cores. Default is \code{TRUE}.}

\item{...}{ellipsis}

\item{object}{of class \code{tvpagfl}.}

\item{x}{of class \code{tvpagfl}.}
}
\value{
An object of class \code{tvpagfl} holding
\item{\code{model}}{a \code{data.frame} containing the dependent and explanatory variables as well as cross-sectional and time indices,}
\item{\code{coefficients}}{let \eqn{p^{(1)}} denote the number of time-varying coefficients and \eqn{p^{(2)}} the number of time constant parameters. A \code{list} holding (i) a \eqn{T \times p^{(1)} \times \hat{K}} array of the post-Lasso group-specific functional coefficients and (ii) a \eqn{K \times p^{(2)}} matrix of time-constant post-Lasso estimates.}
\item{\code{groups}}{a \code{list} containing (i) the total number of groups \eqn{\hat{K}} and (ii) a vector of estimated group memberships \eqn{(\hat{g}_1, \dots, \hat{g}_N)}, where \eqn{\hat{g}_i = k} if \eqn{i} is assigned to group \eqn{k},}
\item{\code{residuals}}{a vector of residuals of the demeaned model,}
\item{\code{fitted}}{a vector of fitted values of the demeaned model,}
\item{\code{args}}{a \code{list} of additional arguments,}
\item{\code{IC}}{a \code{list} containing (i) the value of the IC, (ii) the employed tuning parameter \eqn{\lambda}, and (iii) the \emph{MSE},}
\item{\code{convergence}}{a \code{list} containing (i) a logical variable if convergence was achieved and (ii) the number of executed \emph{ADMM} algorithm iterations,}
\item{\code{call}}{the function call.}

An object of class \code{tvpagfl} has \code{print}, \code{summary}, \code{fitted}, \code{residuals}, \code{formula}, \code{df.residual} and \code{coef} S3 methods.
}
\description{
Estimate a time-varying panel data model with a latent group structure using the pairwise adaptive group fused lasso (\emph{time-varying PAGFL}). The \emph{time-varying PAGFL} jointly identifies the latent group structure and group-specific time-varying functional coefficients.
The time-varying coefficients are modeled as polynomial B-splines. The function supports both static and dynamic panel data models.
}
\details{
Consider the grouped time-varying panel data model
\deqn{y_{it} = \gamma_i + \beta^\prime_{i} (t/T) x_{it} + \epsilon_{it}, \quad i = 1, \dots, N, \; t = 1, \dots, T,}
where \eqn{y_{it}} is the scalar dependent variable, \eqn{\gamma_i} is an individual fixed effect, \eqn{x_{it}} is a \eqn{p \times 1} vector of explanatory variables, and \eqn{\epsilon_{it}} is a zero mean error.
The coefficient vector \eqn{\beta_{i} (t/T)} is subject to the latent group pattern
\deqn{\beta_i \left(\frac{t}{T} \right) = \sum_{k = 1}^K \alpha_k \left( \frac{t}{T} \right) \bold{1} \{i \in G_k \},}
with \eqn{\cup_{k = 1}^K G_k = \{1, \dots, N\}}, \eqn{G_k \cap G_j = \emptyset} and \eqn{\| \alpha_k - \alpha_j \| \neq 0} for any \eqn{k \neq j}, \eqn{k = 1, \dots, K}.

The time-varying coefficient functions are estimated as polynomial B-splines using the penalized sieve-technique. To this end, let \eqn{B(v)} denote a \eqn{M + d +1} vector basis functions, where \eqn{d} denotes the polynomial degree and \eqn{M} the number of interior knots.
Then, \eqn{\beta_{i}(t/T)} and \eqn{\alpha_{k}(t/T)} are approximated by forming linear combinations of the basis functions \eqn{\beta_{i} (t/T) \approx \pi_i^\prime B(t/T)} and \eqn{\alpha_{i}(t/T) \approx \xi_k^\prime B(t/T)}, where \eqn{\pi_i} and \eqn{\xi_i} are \eqn{(M + d + 1) \times p} coefficient matrices.

The explanatory variables are projected onto the spline basis system, which results in the \eqn{(M + d + 1)p \times 1} vector \eqn{z_{it} = x_{it} \otimes B(v)}. Subsequently, the DGP can be reformulated as
\deqn{y_{it} = \gamma_i + z_{it}^\prime \text{vec}(\pi_{i}) + u_{it},}
where \eqn{u_{it} = \epsilon_{it} + \eta_{it}} and \eqn{\eta_{it}} reflects a sieve approximation error. We refer to Su et al. (2019, sec. 2) for more details on the sieve technique.

Inspired by Su et al. (2019) and Mehrabani (2023), the time-varying PAGFL jointly estimates the functional coefficients and the group structure by minimizing the criterion
\deqn{Q_{NT} (\bold{\pi}, \lambda) = \frac{1}{NT} \sum^N_{i=1} \sum^{T}_{t=1}(\tilde{y}_{it} - \tilde{z}_{it}^\prime \text{vec}(\pi_{i}))^2 + \frac{\lambda}{N} \sum_{i = 1}^{N - 1} \sum_{j > i}^N \dot{\omega}_{ij} \| \pi_i - \pi_j \|}
with respect to \eqn{\bold{\pi} = (\text{vec}(\pi_i)^\prime, \dots, \text{vec}(\pi_N)^\prime)^\prime}. \eqn{\tilde{a}_{it} = a_{it} - T^{-1} \sum^{T}_{t=1} a_{it}}, \eqn{a = \{y, z\}} to concentrate out the individual fixed effects \eqn{\gamma_i}. \eqn{\lambda} is the penalty tuning parameter and \eqn{\dot{w}_{ij}} denotes adaptive penalty weights which are obtained by a preliminary non-penalized estimation. \eqn{\| \cdot \|} represents the Frobenius norm.
The solution criterion function is minimized via the iterative alternating direction method of multipliers (\emph{ADMM}) algorithm proposed by Mehrabani (2023, sec. 5.1).

Two individuals are assigned to the same group if \eqn{\| \text{vec} (\hat{\pi}_i - \hat{\pi}_j) \| \leq \epsilon_{\text{tol}}}, where \eqn{\epsilon_{\text{tol}}} is determined by \code{tol_group}. Subsequently, the number of groups follows as the number of distinct elements in \eqn{\hat{\bold{\pi}}}. Given an estimated group structure, it is straightforward to obtain post-Lasso estimates \eqn{\hat{\bold{\xi}}} using group-wise least squares (see \code{\link{grouped_tv_plm}}).

We recommend identifying a suitable \eqn{\lambda} parameter by passing a logarithmically spaced grid of candidate values with a lower limit close to 0 and an upper limit that leads to a fully homogeneous panel. A BIC-type information criterion then selects the best fitting \eqn{\lambda} value.

In case of an unbalanced panel data set, the earliest and latest available observations per group define the start and end-points of the interval on which the group-specific time-varying coefficients are defined.
}
\examples{
# Simulate a time-varying panel with a trend and a group pattern
set.seed(1)
sim <- sim_tv_DGP(N = 10, n_periods = 50, intercept = TRUE, p = 1)
df <- data.frame(y = c(sim$y))

# Run the time-varying PAGFL
estim <- tv_pagfl(y ~ ., data = df, n_periods = 50, lambda = 10, parallel = FALSE)
summary(estim)

}
\references{
Mehrabani, A. (2023). Estimation and identification of latent group structures in panel data. \emph{Journal of Econometrics}, 235(2), 1464-1482. \doi{10.1016/j.jeconom.2022.12.002}.

Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying panel data models with latent structures. \emph{Journal of Business & Economic Statistics}, 37(2), 334-349. \doi{10.1080/07350015.2017.1340299}.
}
\author{
Paul Haimerl
}
