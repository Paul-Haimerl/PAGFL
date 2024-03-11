
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PAGFL

<!-- badges: start -->

[![CRAN_Version_Badge](http://www.r-pkg.org/badges/version/PAGFL)](https://cran.r-project.org/package=PAGFL)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/PAGFL)](https://cran.r-project.org/package=PAGFL)
[![License_GPLv3_Badge](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/Paul-Haimerl/PAGFL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Paul-Haimerl/PAGFL/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

In panel data analysis, unobservable group structures are a common
challenge. Disregarding group-level heterogeneity by assuming an
entirely homogeneous panel can introduce bias. Conversely, estimating
individual coefficients for each cross-sectional unit is inefficient and
may lead to high uncertainty.

Mehrabani ([2023](https://doi.org/10.1016/j.jeconom.2022.12.002))
introduces the pairwise adaptive group fused Lasso (*PAGFL*), a fast
methodology to identify latent group structures and estimate
group-specific coefficients simultaneously.

The `PAGFL` package makes this powerful procedure easy to use.

## Installation

You can install the development version of `PAGFL` (1.1.0) from
[GitHub](https://github.com/) with:

``` r
# install.packages('devtools')
devtools::install_github('Paul-Haimerl/PAGFL')
#> Skipping install of 'PAGFL' from a github remote, the SHA1 (c4514ba5) has not changed since last install.
#>   Use `force = TRUE` to force installation
library(PAGFL)
```

The stable version (1.0.1) is available on CRAN:

    install.packages("PAGFL")

## Data

The `PAGFL` packages includes a function that automatically simulates a
panel with a group structure:

``` r
# Simulate a simple panel with three distinct groups and two exogenous explanatory variables
set.seed(1)
sim <- sim_DGP(N = 50, n_periods = 150, p = 2, n_groups = 3)
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
estim <- pagfl(y = y, X = X, n_periods = 150, lambda = 10)
print(estim)
#> $IC
#> [1] 1.27562
#> 
#> $lambda
#> [1] 10
#> 
#> $alpha_hat
#>               [,1]      [,2]
#> Group 1 -0.3373423  1.615148
#> Group 2 -0.5001927 -1.175167
#> 
#> $K_hat
#> [1] 2
#> 
#> $groups_hat
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#>  1  2  1  1  1  1  1  1  2  2  1  1  1  2  2  2  1  1  2  1  1  1  2  2  1  2 
#> 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
#>  2  1  1  1  1  2  1  1  1  1  1  1  1  2  1  1  1  2  1  2  1  2  2  1 
#> 
#> $iter
#> [1] 124
#> 
#> $convergence
#> [1] TRUE
```

`pagfl` returns a list holding

1.  the value of the IC (see Mehrabani
    [2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 3.4)
2.  the $\lambda$ parameter
3.  the $\widehat{K} \times p$ post-Lasso coefficient matrix
    $\hat{\boldsymbol{\alpha}}^p$
4.  the estimated number of groups $\widehat{K}$
5.  the estimated group structure
6.  the number of executed *ADMM* algorithm iterations
7.  a logical indicator if convergence was achieved

Selecting a $\lambda$ value a priori can be tricky. For instance, it
seems like `lambda = 10` is too high since the number of groups $K$ is
underestimated. As a consequence, we suggest iterating over a
comprehensive range of candidate values. To specify a suitable grid,
create a logarithmic sequence ranging from 0 to a penalty parameter that
induces an entirely homogeneous model (i.e., $\widehat{K} = 1$). The
resulting $\lambda$ grid vector can be passed in place of any specific
value, and a BIC IC selects the best-fitting parameter.

Moreover, if the explanatory variables in `X` are named, those names
also appear in the output.

``` r
colnames(X) <- c("a", "b")
lambda_set <- exp(log(10) * seq(log10(1e-4), log10(10), length.out = 10))
estim_set <- pagfl(y = y, X = X, n_periods = 150, lambda = lambda_set)
print(estim_set)
#> $IC
#> [1] 1.020187
#> 
#> $lambda
#> [1] 0.05994843
#> 
#> $alpha_hat
#>                  a         b
#> Group 1 -0.9874154  1.636026
#> Group 2 -0.5001927 -1.175167
#> Group 3  0.2976462  1.613246
#> 
#> $K_hat
#> [1] 3
#> 
#> $groups_hat
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#>  1  2  1  3  3  3  3  3  2  2  3  3  1  2  2  2  3  3  2  1  3  3  2  2  1  2 
#> 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
#>  2  3  3  1  1  2  1  1  3  3  1  1  1  2  3  1  1  2  1  2  1  2  2  1 
#> 
#> $iter
#> [1] 118
#> 
#> $convergence
#> [1] TRUE
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
sim_endo <- sim_DGP(N = 50, n_periods = 150, p = 2, n_groups = 3, group_proportions = c(0.2, 0.2, 0.6), 
error_spec = 'GARCH', q = 3)
y_endo <- sim_endo$y
X_endo <- sim_endo$X
Z <- sim_endo$Z

# Note that the method PGMM and the instrument matrix Z needs to be supplied
estim_endo <- pagfl(y = y_endo, X = X_endo, n_periods = 150, lambda = 0.05, method = 'PGMM', Z = Z, 
bias_correc = TRUE, max_iter = 8e3)
print(estim_endo)
#> $IC
#> [1] 4.101225
#> 
#> $lambda
#> [1] 0.05
#> 
#> $alpha_hat
#>               [,1]       [,2]
#> Group 1  0.3287358 -1.9741355
#> Group 2 -1.2817828 -1.4823077
#> Group 3  1.3671438 -0.9847914
#> 
#> $K_hat
#> [1] 3
#> 
#> $groups_hat
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#>  1  1  1  2  1  3  3  1  1  2  3  3  1  1  1  1  1  2  1  1  1  1  2  2  2  1 
#> 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 
#>  2  3  1  1  3  2  2  1  3  1  1  3  1  1  1  1  1  1  3  3  3  1  2  3 
#> 
#> $iter
#> [1] 8000
#> 
#> $convergence
#> [1] FALSE
```

Furthermore, `pagfl` lets you select a minimum group size, adjust the
efficiency vs. accuracy trade-off of the iterative estimation algorithm,
and modify a list of further settings. Visit the documentation
`?pagfl()` for more information.

## The Time-varying PAGFL

The development version of the package also includes the functions
`sim_dyn_DGP`and `dyn_pagfl`, which generate and estimate a grouped
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
N <- 50
n_periods <- 50
sim_dyn <- sim_dyn_DGP(N = N, n_periods = n_periods, DGP = 1)
y <- sim_dyn$y
X <- sim_dyn$X

dyn_estim <- dyn_pagfl(y = y, X = X, n_periods = n_periods, lambda = 6)
```

In contrast to the time-constant `pagfl`, `dyn_pagfl` does not return a
post-Lasso coefficient matrix but a $T \times p \times \widehat{K}$
array.

In empirical applications, it is commonplace to encounter unbalanced
panel data sets. In such instances, time-varying coefficient functions
can be estimated nonetheless. The nonparametric spline functions simply
smooth over the time periods with only a limited number or even no
observations. However, it is required to provide explicit indicator
variables that declare the cross-sectional individual and time period
each observation belongs to.

Lets delete a couple of observations, add indicator variables, and run
`dyn_pagfl` again.

``` r
# Draw some observations to be omitted
delete_index_y <- as.logical(rbinom(n = N * n_periods, prob = 0.9, size = 1))
delete_index_X <- as.logical(rbinom(n = N * n_periods, prob = 0.9, size = 1))
# Construct cross-sectional and time indicator variables
i_index <- rep(1:N, each = n_periods)
t_index <- rep(1:n_periods, N)
y <- cbind(y, i_index = i_index, t_index = t_index)
X <- cbind(X, i_index = i_index, t_index = t_index)
# Delte some observations
y <- y[delete_index_y,]
X <- X[delete_index_X,]
# Apply the time-varying PAGFL to an unbalanced panel
dyn_estim_unbalanced <- dyn_pagfl(y = y, X = X, index = c("i_index", "t_index"), lambda = 6)
#> Warning in second_checks(N = N, index = index, n_periods = n_periods, y = y, : The panel data set is unbalanced
```

Furthermore, `dyn_pagfl` lets you specify a lot more optionalities than
shown here. For example, it is possible to adjust the degree of the
spline functions, the number of interior knots in the spline system, or
estimate a panel data model with a mix of time-varying and time-constant
coefficients. See `?dyn_pagfl()` for details.

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
