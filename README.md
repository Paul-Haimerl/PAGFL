
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PAGFL

<!-- badges: start -->

[![CRAN_Version_Badge](http://www.r-pkg.org/badges/version/PAGFL)](https://cran.r-project.org/package=PAGFL)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/PAGFL)](https://cran.r-project.org/package=PAGFL)
[![License_GPLv3_Badge](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/Paul-Haimerl/PAGFL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Paul-Haimerl/PAGFL/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/Paul-Haimerl/PAGFL/graph/badge.svg?token=22WHU5SU63)](https://app.codecov.io/gh/Paul-Haimerl/PAGFL)
<!-- badges: end -->

Unobservable group structures are a common challenge in panel data
analysis. Disregarding group-level heterogeneity can introduce bias.
Conversely, estimating individual coefficients for each cross-sectional
unit is inefficient and may lead to high uncertainty.

This package efficiently addresses the issue of unobservable group
structures by implementing the pairwise adaptive group fused Lasso
(*PAGFL*) by Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002)). *PAGFL* is a
regularizer that identifies latent group structures and estimates
group-specific coefficients in a single step. On top of that, we extend
the PAGFL to time-varying functional coefficients.

The `PAGFL` package makes this powerful procedure easy to use. On top of
that, we extend the `PAGFL` to time-varying functional coefficients.

## Installation

You can install the development version of `PAGFL` (1.1.0) from
[GitHub](https://github.com/) with:

``` r
# install.packages('devtools')
devtools::install_github('Paul-Haimerl/PAGFL')
library(PAGFL)
```

The stable version (1.0.1) is available on CRAN:

    install.packages("PAGFL")

## Data

The `PAGFL` packages includes a function that automatically simulates a
panel data set with a group structure in the slope coefficients:

``` r
# Simulate a simple panel with three distinct groups and two exogenous explanatory variables
set.seed(1)
sim <- sim_DGP(N = 20, n_periods = 150, p = 2, n_groups = 3)
y <- sim$y
X <- sim$X
```

$$y_{it} = \beta_i^\prime x_{it} + \eta_i + u_{it}, \quad i = 1, \dots, N, \quad t = 1, \dots, T,$$
where $y_{it}$ is a scalar dependent variable, $x_{it}$ a $p \times 1$
vector of explanatory variables, and $\eta_i$ reflects a fixed effect.
The slope coefficients are subject to the group structure

$$\beta_{i} = \sum_{k = 1}^K \alpha_k \boldsymbol{1} \{i \in G_k \},$$
with $\cup_{k = 1}^K G_k = \{1, \dots, N \}$, and
$G_k \cap G_j = \emptyset$ as well as $|| \alpha_k \neq \alpha_j ||$ for
any $k \neq j$, $k,j = 1, \dots, K$ (see Mehrabani
[2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 2).

`sim_DGP` also nests, among other, all DGPs employed in the simulation
study of Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 6). I refer
to the documentation of `sim_DGP` or Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 6) for more
details.

## Applying PAGFL

To execute the PAGFL procedure, simply pass the dependent and
independent variables, the number of time periods, and a penalization
parameter $\lambda$.

``` r
estim <- pagfl(y ~ X, n_periods = 150, lambda = 20)
summary(estim)
#> Call:
#> pagfl(formula = y ~ X, n_periods = 150, lambda = 20)
#> 
#> Balanced panel: N = 20, T = 150, obs = 3000
#> 
#> Convergence reached:
#> TRUE (49 iterations)
#> 
#> Information criterion:
#>        IC    lambda 
#>  1.353997 20.000000 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -4.47230 -0.72086 -0.00120  0.76214  4.31838 
#> 
#> 2 groups:
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  1  2  1  1  1  1  2  1  1  2  2  2  1  1  1  1  1  2  1 
#> 
#> Coefficients:
#>                X1       X2
#> Group 1 -0.36838  1.61275
#> Group 2 -0.49489 -1.23534
#> 
#> Residual standard error: 1.15012 on 2978 degrees of freedom
#> Mean squared error 1.31307
#> Multiple R-squared: 0.65845, Adjusted R-squared: 0.65605
```

`pagfl()` returns an object of type `pagfl` which holds

1.  `model`: A `data.frame` containing the dependent and explanatory
    variables as well as individual and time indices (if provided).
2.  `coefficients`: A $K \times p$ matrix of the post-Lasso
    group-specific parameter estimates.
3.  `groups`: A `list` containing (i) the total number of groups
    $\hat{K}$ and (ii) a vector of estimated group memberships
    $(\hat{g}_1, \dots, \hat{g}_N)$, where $\hat{g}_i = k$ if $i$ is
    assigned to group $k$.
4.  `residuals`: A vector of residuals of the demeaned model.
5.  `fitted`: A vector of fitted values of the demeaned model.
6.  `args`: A list of additional arguments.
7.  `IC`: A `list` containing (i) the value of the IC, (ii) the employed
    tuning parameter $\lambda$, and (iii) the mean squared error.
8.  `convergence`: A `list` containing (i) a logical variable if
    convergence was achieved and (ii) the number of executed *ADMM*
    algorithm iterations.
9.  `call`: The function call.

Furthermore, `pagfl` objects can be used in a variety of useful generic
methods like `summary()`, `fitted()`, `resid()`, `df.residual`,
`formula`, and `coef()`.

``` r
estim_fit <- fitted(estim)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Selecting a $\lambda$ value a priori can be tricky. For instance, it
seems like `lambda = 20` is too high since the number of groups $K$ is
underestimated. We suggest iterating over a comprehensive range of
candidate values to trace out the correct model. To specify a suitable
grid, create a logarithmic sequence ranging from 0 to a penalty
parameter that induces an entirely homogeneous model (i.e.,
$\widehat{K} = 1$). The resulting $\lambda$ grid vector can be passed in
place of any specific value, and a BIC IC selects the best-fitting
parameter.

Furthermore, it is also possible to supply a `data.frame` with named
variables and choose a specific formula that selects the variables in
that `data.frame`. If the explanatory variables in `X` are named, these
names also appear in the output.

``` r
colnames(X) <- c("a", "b")
data <- cbind(y = c(y), X)

lambda_set <- exp(log(10) * seq(log10(1e-4), log10(10), length.out = 10))
estim_set <- pagfl(y ~ a + b, data = data, n_periods = 150, lambda = lambda_set)
summary(estim_set)
#> Call:
#> pagfl(formula = y ~ a + b, data = data, n_periods = 150, lambda = lambda_set)
#> 
#> Balanced panel: N = 20, T = 150, obs = 3000
#> 
#> Convergence reached:
#> TRUE (50 iterations)
#> 
#> Information criterion:
#>        IC    lambda 
#> 1.1287693 0.2154435 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.47858 -0.66283 -0.02688  0.72880  3.77812 
#> 
#> 3 groups:
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  1  2  3  1  3  3  2  3  3  2  2  2  1  1  1  3  1  2  3 
#> 
#> Coefficients:
#>                 a        b
#> Group 1 -0.95114  1.61719
#> Group 2 -0.49489 -1.23534
#> Group 3  0.24172  1.61613
#> 
#> Residual standard error: 1.03695 on 2978 degrees of freedom
#> Mean squared error 1.06738
#> Multiple R-squared: 0.72236, Adjusted R-squared: 0.7204
```

When, as above, the specific estimation method is left unspecified,
`pagfl` defaults to penalized Least Squares (*PLS*) `method = 'PLS'`
(Mehrabani, [2023](https://doi.org/10.1016/j.jeconom.2022.12.002),
sec. 2.2). *PLS* is very efficient but requires weakly exogenous
regressors. However, even endogenous predictors can be accounted for by
employing a penalized Generalized Method of Moments (*PGMM*) routine in
combination with exogenous instruments $\boldsymbol{Z}$.

Specify a slightly more elaborate endogenous and dynamic panel data set
and apply *PGMM*. When encountering a dynamic panel data set, we
recommend using a Jackknife bias correction, as proposed by Dhaene and
Jochmans ([2015](https://doi.org/10.1093/restud/rdv007)).

``` r
# Generate a panel where the predictors X correlate with the cross-sectional innovation, 
# but can be instrumented with q = 3 variables in Z. Furthermore, include GARCH(1,1) 
# innovations, an AR lag of the dependent variable, and specific group sizes
sim_endo <- sim_DGP(N = 25, n_periods = 200, p = 2, n_groups = 3, group_proportions = c(0.3, 0.3, 0.4), 
error_spec = 'GARCH', q = 2, dynamic = TRUE)
y_endo <- sim_endo$y
X_endo <- sim_endo$X
Z <- sim_endo$Z

# Note that the method PGMM and the instrument matrix Z needs to be passed
estim_endo <- pagfl(y_endo ~ X_endo, n_periods = 200, lambda = 15, method = 'PGMM', Z = Z, bias_correc = TRUE, max_iter = 20e3)
summary(estim_endo)
#> Call:
#> pagfl(formula = y_endo ~ X_endo, n_periods = 200, lambda = 15, 
#>     method = "PGMM", Z = Z, bias_correc = TRUE, max_iter = 20000)
#> 
#> Balanced panel: N = 25, T = 200, obs = 4975
#> 
#> Convergence reached:
#> FALSE (20000 iterations)
#> 
#> Information criterion:
#>      IC  lambda 
#> 56.2871 15.0000 
#> 
#> Residuals:
#>       Min        1Q    Median        3Q       Max 
#> -97.91933  -2.40370   0.01101   2.41854  77.69418 
#> 
#> 3 groups:
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
#>  1  2  2  1  2  2  2  2  2  1  1  1  1  1  1  1  1  1  1  2  3  1  1  1  1 
#> 
#> Coefficients:
#>                 X1       X2
#> Group 1  -2.01709  1.07102
#> Group 2  -2.08412 -0.93513
#> Group 3 -10.47215  4.24136
#> 
#> Residual standard error: 7.51953 on 4948 degrees of freedom
#> Mean squared error 56.23651
#> Multiple R-squared: -8.37322, Adjusted R-squared: -8.42247
```

Furthermore, `pagfl` lets you select a minimum group size, adjust the
efficiency vs. accuracy trade-off of the iterative estimation algorithm,
and modify a list of further settings. Visit the documentation
`?pagfl()` for more information.

## The Time-varying PAGFL

The development version of the package also includes the functions
`sim_tv_DGP()`and `tv_pagfl()`, which generate and estimate a grouped
panel data models with the time-varying coefficients
$\beta_{it} = \beta_i \left( \frac{t}{T} \right)$. Just like in the
static case, the functional coefficients admit to a group structure
$\beta_{it} = \sum_{k = 1}^K \alpha_k \left( \frac{t}{T} \right) 1 \{i \in G_k \}$.
Following Su et
al. ([2019](https://doi.org/10.1080/07350015.2017.1340299)), the
time-varying coefficients are estimated using polynomial B-spline
functions employing a penalized sieve estimation (*PSE*).

``` r
# Simulate a time-varying panel with a trend and a group pattern
N <- 20
n_periods <- 100
tv_sim <- sim_tv_DGP(N = N, n_periods = n_periods, sd_error = 1, intercept = TRUE, p = 1)
y <- tv_sim$y

tv_estim <- tv_pagfl(y ~ 1, n_periods = n_periods, lambda = 5)
summary(tv_estim)
#> Call:
#> tv_pagfl(formula = y ~ 1, n_periods = n_periods, lambda = 5)
#> 
#> Balanced panel: N = 20, T = 100, obs = 2000
#> 
#> Convergence reached:
#> TRUE (207 iterations)
#> 
#> Information criterion:
#>      IC  lambda 
#> 1.18412 5.00000 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.83243 -0.68348  0.02244  0.65131  2.94155 
#> 
#> 3 groups:
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  1  2  3  3  2  1  3  1  1  2  3  2  1  2  1  3  3  2  3 
#> 
#> Residual standard error: 1.00966 on 1973 degrees of freedom
#> Mean squared error 1.00566
#> Multiple R-squared: 0.73854, Adjusted R-squared: 0.73509
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

`tv_pagfl()` returns an object of class `tvpagfl` which contains

1.  `model`: A `data.frame` containing the dependent and explanatory
    variables as well as individual and time indices (if provided).
2.  `coefficients`: A list holding (i) a
    $T \times p^{(1)} \times \hat{K}$ array of the post-Lasso
    group-specific functional coefficients and (ii) a $K \times p^{(2)}$
    matrix of time-constant parameter estimates (when running a mixed
    time-varying panel data model).
3.  `groups`: A `list` containing (i) the total number of groups
    $\hat{K}$ and (ii) a vector of estimated group memberships
    $(\hat{g}_1, \dots, \hat{g}_N)$, where $\hat{g}_i = k$ if $i$ is
    assigned to group $k$.
4.  `residuals`: A vector of residuals of the demeaned model.
5.  `fitted`: A vector of fitted values of the demeaned model.
6.  `args`: A list of additional arguments.
7.  `IC`: A `list` containing (i) the value of the IC, (ii) the employed
    tuning parameter $\lambda$, and (iii) the mean squared error.
8.  `convergence`: A `list` containing (i) a logical variable if
    convergence was achieved and (ii) the number of executed *ADMM*
    algorithm iterations.
9.  `call`: The function call.

Again, `tvpagfl` objects have generic `summary()`, `fitted()`,
`resid()`, `df.residual`, `formula`, and `coef()` methods.

In empirical applications, it is commonplace to encounter unbalanced
panel data sets. In such instances, time-varying coefficient functions
can be estimated nonetheless. The nonparametric spline functions simply
interpolate missing periods. However, when using unbalanced datasets it
is required to provide explicit indicator variables that declare the
cross-sectional individual and time period each observation belongs to.

Lets delete 30% of observations, add indicator variables, and run
`tv_pagfl()` again.

``` r
# Draw some observations to be omitted
delete_index <- as.logical(rbinom(n = N * n_periods, prob = 0.7, size = 1))
# Construct cross-sectional and time indicator variables
i_index <- rep(1:N, each = n_periods)
t_index <- rep(1:n_periods, N)
data <- cbind(y = c(y), i_index = i_index, t_index = t_index)
# Delete some observations and create a named data.frame
data <- data[delete_index,]
# Apply the time-varying PAGFL to an unbalanced panel
tv_estim_unbalanced <- tv_pagfl(y ~ 1, data = data, index = c("i_index", "t_index"), lambda = 10)
summary(tv_estim_unbalanced)
#> Call:
#> tv_pagfl(formula = y ~ 1, data = data, index = c("i_index", "t_index"), 
#>     lambda = 10)
#> 
#> Unbalanced panel: N = 20, T = 62-78, obs = 1376
#> 
#> Convergence reached:
#> TRUE (133 iterations)
#> 
#> Information criterion:
#>       IC   lambda 
#>  1.16034 10.00000 
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.62126 -0.66947  0.00878  0.63187  3.02198 
#> 
#> 3 groups:
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  1  2  3  3  2  1  3  1  1  2  3  2  1  2  1  3  3  2  3 
#> 
#> Residual standard error: 1.00077 on 1349 degrees of freedom
#> Mean squared error 0.98188
#> Multiple R-squared: 0.74348, Adjusted R-squared: 0.73854
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Furthermore, `tv_pagfl` lets you specify a lot more optionalities than
shown here. For example, it is possible to adjust the polyomial degree
and the number of interior knots in the spline basis system, or estimate
a panel data model with a mix of time-varying and time-constant
coefficients. See `?tv_pagfl()` for details.

## References

- Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of
  fixed-effect models. *The Review of Economic Studies*, 82(3),
  991-1030. DOI:
  [10.1093/restud/rdv007](https://doi.org/10.1093/restud/rdv007)

- Mehrabani, A. (2023). Estimation and identification of latent group
  structures in panel data. *Journal of Econometrics*, 235(2),
  1464-1482. DOI:
  [10.1016/j.jeconom.2022.12.002](https://doi.org/10.1016/j.jeconom.2022.12.002)

- Schumaker, L. (2007). Spline functions: basic theory. *Cambridge
  university press*. DOI:
  [10.1017/CBO9780511618994](https://doi.org/10.1017/CBO9780511618994)

- Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying
  panel data models with latent structures. *Journal of Business &
  Economic Statistics*, 37(2), 334-349. DOI:
  [10.1080/07350015.2017.1340299](https://doi.org/10.1080/07350015.2017.1340299)
