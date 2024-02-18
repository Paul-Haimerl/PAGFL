
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
introduces the pairwise adaptive group fused Lasso (PAGFL), a fast
methodology to identify latent group structures and estimate
group-specific coefficients simultaneously.

The `PAGFL` package makes this powerful procedure easy to use.

## Installation

You can install the development version of `PAGFL` from
[GitHub](https://github.com/) with:

``` r
# install.packages('devtools')
devtools::install_github('Paul-Haimerl/PAGFL')
#> Downloading GitHub repo Paul-Haimerl/PAGFL@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\phaim\AppData\Local\Temp\Rtmpi2lEj4\remotes5b2828a639ef\Paul-Haimerl-PAGFL-6129f7d/DESCRIPTION' ...     checking for file 'C:\Users\phaim\AppData\Local\Temp\Rtmpi2lEj4\remotes5b2828a639ef\Paul-Haimerl-PAGFL-6129f7d/DESCRIPTION' ...   ✔  checking for file 'C:\Users\phaim\AppData\Local\Temp\Rtmpi2lEj4\remotes5b2828a639ef\Paul-Haimerl-PAGFL-6129f7d/DESCRIPTION' (548ms)
#>       ─  preparing 'PAGFL':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>       ─  checking for empty or unneeded directories
#>       ─  building 'PAGFL_1.0.1.tar.gz'
#>      
#> 
library(PAGFL)
```

The stable version is available on CRAN:

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

$$ y_{it} = \beta_i^\prime x_{it} + \eta_i + u_{it}, \quad i = 1, \dots, N, \quad t = 1, \dots, T,$$
where $y_{it}$ is a scalar dependent variable, $x_{it}$ a $p \times 1$
vector of explanatory variables and $\eta_i$ reflects a fixed effect.
The slope coefficients are subject to the group structure

$$\beta_{it} = \sum_{k = 1}^K \alpha_k \boldsymbol{1}_{\{i \in G_k \}},$$
with $\cup_{k = 1}^K G_k = \{1, \dots, N \}$, and
$G_k \cap G_j = \emptyset$ as well as $|| \alpha_k \neq \alpha_j ||$ for
any $k \neq j$, $k,j = 1, \dots, K$ (see Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 2).

`sim_DGP` also nests, among other, all DGPs employed in the simulation
study of Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 6). I refer
to the documentation of `sim_DGP` or Mehrabani
([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 6) for more
details.

## Applying PAGFL

To execute the PAGFL procedure, simply pass the dependent and
independent variables, the number of time periods and a penalization
parameter $\lambda$.

``` r
estim <- PAGFL(y = y, X = X, n_periods = 150, lambda = 5)
print(estim)
#> $IC
#> [1] 1.020187
#> 
#> $lambda
#> [1] 5
#> 
#> $alpha_hat
#>            [,1]      [,2]
#> [1,] -0.9874154  1.636026
#> [2,] -0.5001927 -1.175167
#> [3,]  0.2976462  1.613246
#> 
#> $K_hat
#> [1] 3
#> 
#> $groups_hat
#>  [1] 1 2 1 3 3 3 3 3 2 2 3 3 1 2 2 2 3 3 2 1 3 3 2 2 1 2 2 3 3 1 1 2 1 1 3 3 1 1
#> [39] 1 2 3 1 1 2 1 2 1 2 2 1
#> 
#> $iter
#> [1] 37
#> 
#> $convergence
#> [1] TRUE
```

`PAGFL` returns a list holding

1.  the value of the IC (see Mehrabani
    ([2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 3.4))
2.  the $\lambda$ parameter
3.  the $\widehat{K} \times p$ post-Lasso coefficient matrix
    $\hat{\boldsymbol{\alpha}}^p$
4.  the estimated number of groups $\widehat{K}$
5.  the estimated group structure
6.  the number of executed *ADMM* algorithm iterations
7.  a logical indicator if convergence was achieved

Selecting a $\lambda$ value a priori can be tricky. Therefore, we
suggest iterating over a comprehensive range of candidate values. To
specify a suitable grid, create a logarithmic sequence ranging from 0 to
a penalty parameter that induces an entirely homogeneous model (i.e.,
$\widehat{K} = 1$). The resulting $\lambda$ grid vector can be passed in
place of any specific value, and a BIC IC selects the best-fitting
parameter.

``` r
lambda_set <- exp(log(10) * seq(log10(1e-4), log10(10), length.out = 10))
estim_set <- PAGFL(y = y, X = X, n_periods = 150, lambda = lambda_set)
print(estim_set)
#> $IC
#> [1] 1.020187
#> 
#> $lambda
#> [1] 0.05994843
#> 
#> $alpha_hat
#>            [,1]      [,2]
#> [1,] -0.9874154  1.636026
#> [2,] -0.5001927 -1.175167
#> [3,]  0.2976462  1.613246
#> 
#> $K_hat
#> [1] 3
#> 
#> $groups_hat
#>  [1] 1 2 1 3 3 3 3 3 2 2 3 3 1 2 2 2 3 3 2 1 3 3 2 2 1 2 2 3 3 1 1 2 1 1 3 3 1 1
#> [39] 1 2 3 1 1 2 1 2 1 2 2 1
#> 
#> $iter
#> [1] 31
#> 
#> $convergence
#> [1] TRUE
```

When, as above, the specific estimation method is left unspecified,
`PAGFL` defaults to penalized Least Squares (*PLS*) (Mehrabani,
[2023](https://doi.org/10.1016/j.jeconom.2022.12.002), sec. 2.2). *PLS*
is very efficient but requires weakly exogenous regressors. However,
even endogenous predictors can be accounted for by employing a penalized
Generalized Method of Moments (*PGMM*) routine in combination with
exogenous instruments $Z$.

Specify a slightly more elaborate endogenous and dynamic panel data set
and apply *PGMM*. When encountering a dynamic panel data set, we
recommend using a Jackknife bias correction, as proposed by Dhaene and
Jochmans ([2015](https://doi.org/10.1093/restud/rdv007)).

``` r
# Generate a panel where the predictors X correlate with the cross-sectional innovation, 
# but can be instrumented with q = 3 variables in Z. Furthermore, include GARCH(1,1) 
# innovations, an AR lag of the dependent variable and specific group sizes
sim_endo <- sim_DGP(N = 50, n_periods = 150, p = 2, n_groups = 3, group_proportions = c(0.2, 0.2, 0.6), 
error_spec = 'GARCH', q = 3)
y_endo <- sim_endo$y
X_endo <- sim_endo$X
Z <- sim_endo$Z

# Note that the method PGMM and the instrument matrix Z needs to be supplied
estim_endo <- PAGFL(y = y_endo, X = X_endo, n_periods = 150, lambda = 0.05, method = 'PGMM', Z = Z, 
bias_correc = TRUE, max_iter = 5e3)
print(estim_endo)
#> $IC
#> [1] 1.923868
#> 
#> $lambda
#> [1] 0.05
#> 
#> $alpha_hat
#>            [,1]       [,2]
#> [1,]  0.3416435 -1.9835855
#> [2,] -1.2817828 -1.4823077
#> [3,]  1.6819715 -0.8715061
#> 
#> $K_hat
#> [1] 3
#> 
#> $groups_hat
#>  [1] 1 1 1 2 1 3 3 1 1 2 1 3 1 1 1 1 1 2 1 1 1 1 2 2 2 1 2 3 1 1 3 2 2 1 3 1 1 1
#> [39] 1 1 1 1 1 1 3 3 3 1 2 3
#> 
#> $iter
#> [1] 5000
#> 
#> $convergence
#> [1] FALSE
```

Furthermore, `PAGFL` lets you select a minimum group size, adjust the
efficiency vs. accuracy trade-off of the iterative estimation algorithm,
and modify a list of further settings. Visit the documentation
`?PAGFL()` for more information.

## References

- Dhaene, G., & Jochmans, K. (2015). Split-panel jackknife estimation of
  fixed-effect models. *The Review of Economic Studies*, 82(3),
  991-1030. DOI:
  [10.1093/restud/rdv007](https://doi.org/10.1093/restud/rdv007)

- Mehrabani, A. (2023). Estimation and identification of latent group
  structures in panel data. *Journal of Econometrics*, 235(2),
  1464-1482. DOI:
  [10.1016/j.jeconom.2022.12.002](https://doi.org/10.1016/j.jeconom.2022.12.002)

- Su, L., Wang, X., & Jin, S. (2019). Sieve estimation of time-varying
  panel data models with latent structures. *Journal of Business &
  Economic Statistics*, 37(2), 334-349. DOI:
  [10.1080/07350015.2017.1340299](https://doi.org/10.1080/07350015.2017.1340299)
