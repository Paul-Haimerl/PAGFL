prelim_checks <- function(formula, data, Z = NULL, index = NULL, n_periods = NULL, method = "PLS", const_coef = NULL, verbose, min_group_frac, max_iter, kappa, tol_group, tol_convergence) {
  const_vec <- c(max_iter, kappa, tol_group, tol_convergence)
  if (any(const_vec <= 0)) stop(paste(c("`max_iter`", "`kappa`", "`tol_group`", "`tol_convergence`")[const_vec <= 0], collapse = ", "), " must be bigger than zero.\n")
  if (any(min_group_frac < 0)) stop("`min_group_frac` must be non-zero.\n")
  if (!(method %in% c("PLS", "PGMM"))) stop("The estimation method must be either `PLS` or `PGMM`. Use `PLS` in case of (weakly) exogenous regressors and `PGMM` for endogenous regressors.\n")
  if (is.null(n_periods) & is.null(index)) stop("Either supply cross-sectional and time index variables or in case of a balanced and ordered panel data set, the number of time periods `n_periods`\n")
  if (!is.null(n_periods) & !is.null(index) & verbose) warning("Both index variables and `n_periods` are supplied. `n_periods` is ignored\n")
  if (length(formula) < 3) stop("Please provide a formula with a dependent and explanatory variables\n")
  if (is.null(n_periods)) {
    if (is.null(data)) stop("Please pass a dataset to `data` when specifying index variables\n")
    if (!all(index %in% colnames(data))) stop(paste(index[!(index %in% colnames(data))], "is not part of `data`"))
    if (length(index) != 2 | !is.character(index)) stop("Please supply a character vector holding two strings as `index`\n")
    if (any(index %in% all.vars(formula[[2]]))) stop("Index variables are used as the dependent variable\n")
    if (any(index %in% all.vars(formula[[3]]))) stop("Index variables are used as the explanatory variables\n")
    if (any(is.na(data)) & verbose) warning(paste(deparse(substitute(data)), "contains missing values"), "\n")
  } else {
    if (any(is.na(data))) stop(paste(deparse(substitute(data)), "contains missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables"), "\n")
  }
  if (!is.null(data)) {
    if (!(all.vars(formula[[2]]) %in% colnames(data))) stop(paste(all.vars(formula[[2]]), "not present in", deparse(substitute(data))), "\n")
  } else {
    if (all.vars(formula[[2]]) %in% colnames(data)[-1]) stop("The dependent variable ", all.vars(formula[[2]]), " is part of ", all.vars(formula[[3]]), "\n")
  }
  if (!is.null(const_coef)) {
    if (any(const_coef %in% all.vars(formula[[2]]))) {
      stop(paste(paste(const_coef[const_coef %in% all.vars(formula[[2]])]), "is passed as the dependent variable and as a constant explanatory variable"), "\n")
    }
  }
  if (method == "PGMM") {
    if (is.null(Z)) stop("PGMM requires a matrix of exogenous instruments `Z`\n")
  }
  if (is.null(min_group_frac)) min_group_frac <- 0
  if (min_group_frac > 1 | min_group_frac < 0) stop("Provide a min_group_frac between 0 and 1\n")
  if (min_group_frac >= .4 & verbose) warning("Large min_group_frac values may lead to all groups falling below the group cardinality threshold, in which case the hierarchical clustering algorithm cannot be employed\n")
}


second_checks <- function(N, index, n_periods, y, X, method = NULL, Z = NULL, p, min_group_frac, verbose, dyn, d = NULL, M = NULL, rho, varrho) {
  const_vec <- c(rho, varrho)
  if (any(const_vec <= 0)) stop(paste(c("`rho`", "`varrho`")[const_vec <= 0], collapse = ", "), " must be greater than zero\n")
  if (is.character(y)) stop("Pass a numerical dependent variable\n")
  if (any(is.na(y))) stop("y contains missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables\n")
  if (any(is.na(y))) stop("The explanatory variables contain missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables\n")
  if (is.null(index)) {
    if (round(N) != N) stop("n_periods does not match the number of time periods present in the dataset or the dataset is unblananced.\n  In case of an unbalanced dataset, specify index variables\n")
  }
  if (p == 0) stop("No explanatory variable is supplied\n")
  if (dyn) {
    checksDyn(d = d, M = M, p = p)
  } else {
    checksStat(X, method, Z, p, verbose)
  }
  # if (!all(all.vars(formula[[3]]) %in% colnames(data))) {
  #   coef_names <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in% colnames(data)]
  #   stop(paste(paste(coef_names, collapse = ", "), "not present in", deparse(substitute(data))), "\n")
  # }
}

# Checks for the static PAGFL
checksStat <- function(X, method, Z, p, verbose) {
  sd_vec <- apply(X, 2, stats::sd)
  if (any(sd_vec == 0)) stop("There must not be an intercept or any time-invariant explanatory variables\n")
  if (method == "PGMM") {
    if (ncol(Z) < p) stop(paste("Provide at least p =", p, "exogenous instruments `Z`"), "\n")
    if (nrow(Z) != nrow(X)) stop("The number of time periods of the instrument matrix `Z` does do not match the remaining data\n")
  } else {
    if (!is.null(Z) & verbose) warning("The instrument matrix `Z` is ignored by the PLS estimation algorithm. To instrument endogenous regressors using `Z`, specify method = `PGMM`\n")
  }
}

# Checks for the dynamic PAGFL
checksDyn <- function(d, M, p) {
  const_vec <- c(d, M)
  if (any(const_vec <= 0)) stop("`M` and `d` must both be positive numerical values\n")
  if (p == 0) stop("At least one coefficient must be time-varying. If all coefficients are time-constant, use pagfl() instead\n")
}

# Checks for simulating data
simChecks <- function(dyn, N, n_groups, group_proportions, error_spec = "iid", alpha_0 = NULL, dyn_panel = FALSE, q = NULL, p, locations = NULL, scales = NULL, polynomial_coef = NULL, d = NULL, dynamic, intercept = FALSE) {
  if (N < n_groups) stop("Number of groups cannot exceed number of observations\n")
  if (p == 0) stop("Include at least one explanatory variable\n")
  if (!is.null(group_proportions)) {
    if (n_groups != length(group_proportions)) stop("Number of groups and group proportions are of different length\n")
    if (sum(group_proportions) != 1) stop("Group proportions must sum to 1\n")
  }
  if (!(error_spec %in% c("iid", "AR", "GARCH"))) stop("The individual error specificaion must be either `iid`, `AR` or `GARCH`. Use `iid` for iid errors, `AR` in case of an AR(1) serial correlation, and `GARCH` for GARCH(1,1) innovations.\n")
  if (!is.null(alpha_0)) {
    if (nrow(alpha_0) != n_groups) stop(paste("There are", n_groups, "groups, but only", nrow(alpha_0), "parameter vectors provided"), "\n")
    if (dyn_panel) {
      if (any(abs(alpha_0[, 1]) >= 1)) stop("The AR parameters must be lower than 1 in absolute value\n")
      if (ncol(alpha_0) - 1 != p) stop(paste("There are", ncol(alpha_0), "group-specific parameters provided, but p + 1 =", p + 1, "are required"), "\n")
    } else {
      if (ncol(alpha_0) != p) stop(paste("There are", ncol(alpha_0), "group-specific parameters provided, but p =", p, "are required"), "\n")
    }
  }
  if (!is.null(q)) {
    if (q < p) stop("There must be at least q = p instruments\n")
  }
  if (dyn){
    if (!is.null(locations)){
      if (!is.matrix(locations) | all(dim(locations) != c(p, n_groups))) stop("`locations` must be a matrix with p rows and one column per group\n")
    }
    if (!is.null(scales)){
      if (!is.matrix(scales) | all(dim(scales) != c(p, n_groups))) stop("`scales` must be a matrix with p rows and one column per group\n")
    }
    if (!is.null(polynomial_coef)){
      if (!is.array(polynomial_coef) | all(dim(polynomial_coef) != c(p, d, n_groups))) stop("`polynomial_coef` must be an array with p rows, d columns ane one plane per group\n")
    }
    if (p < 2 & intercept & dynamic) warning("`p = 1` but both an intercept and a dynamic AR coefficient are specified\n")
  }
}
