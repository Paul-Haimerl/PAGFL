prelim_checks <- function(formula, data, Z = NULL, index = NULL, n_periods = NULL, method = "PLS", const_coef = NULL, verbose, min_group_frac = NULL, max_iter = NULL, kappa = NULL, tol_group = NULL, tol_convergence = NULL) {
  const_vec <- c(max_iter, kappa, tol_group, tol_convergence)
  if (any(const_vec <= 0)) stop(paste(c("`max_iter`", "`kappa`", "`tol_group`", "`tol_convergence`")[const_vec <= 0], collapse = ", "), " must be bigger than zero.\n")
  if (is.null(n_periods) & is.null(index)) stop("Either supply cross-sectional and time index variables or in case of a balanced and ordered panel data set, the number of time periods `n_periods`\n")
  if (!is.null(n_periods) & !is.null(index) & verbose) warning("Both index variables and `n_periods` are supplied. `n_periods` is ignored\n")
  if (is.null(n_periods)) {
    if (!all(index %in% colnames(data))) stop(paste(index[!(index %in% colnames(data))], "is not part of `data`\n"))
    if (length(index) != 2 | !is.character(index)) stop("Please supply a character vector holding two strings as `index`\n")
    if (any(index %in% all.vars(formula[[2]]))) stop("Index variables are used as the dependent variable\n")
    if (any(index %in% all.vars(formula[[3]]))) stop("Index variables are used as the explanatory variables\n")
  } else {
    if (any(is.na(data))) stop(paste(deparse(substitute(data)), "contains missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables"), "\n")
  }
  if (!(all.vars(formula[[2]]) %in% colnames(data))) stop(paste(all.vars(formula[[2]]), "not present in", deparse(substitute(data))), "\n")
  if (all.vars(formula[[2]]) %in% all.vars(formula[[3]])) stop("The dependent variable is also passed as a regressor\n")
  if (is.null(min_group_frac)) min_group_frac <- 0
  if (min_group_frac > 1 | min_group_frac < 0) stop("Provide a min_group_frac between 0 and 1\n")
  if (min_group_frac >= .4 & verbose) warning("Large min_group_frac values may lead to all groups falling below the group cardinality threshold, in which case the hierarchical clustering algorithm cannot be employed\n")
}


second_checks <- function(y, X, method = NULL, Z = NULL, p, min_group_frac = NULL, verbose, dyn, d = NULL, M = NULL, rho, varrho = NULL) {
  const_vec <- c(rho, varrho)
  if (any(const_vec <= 0)) stop(paste(c("`rho`", "`varrho`")[const_vec <= 0], collapse = ", "), " must be greater than zero\n")
  if (is.character(y)) stop("Pass a numerical dependent variable\n")
  if (p == 0) stop("No explanatory variable is supplied\n")
  if (dyn) {
    checksDyn(d = d, p = p)
  } else {
    checksStat(X, method, Z, p, verbose)
  }
}

# Checks for the static PAGFL
checksStat <- function(X, method, Z, p, verbose) {
  sd_vec <- apply(X, 2, stats::sd)
  if (any(sd_vec == 0)) stop("There must not be an intercept or any time-invariant explanatory variables\n")
  if (method == "PGMM") {
    if (ncol(Z) < p) stop(paste("Provide at least p =", p, "exogenous instruments `Z`"), "\n")
    if (nrow(Z) != nrow(X)) stop("The number of time periods of the instrument matrix `Z` does do not match the remaining data\n")
  }
}

# Checks for the dynamic PAGFL
checksDyn <- function(d, p) {
  if (any(d <= 0)) stop("`d` must be a positive numerical value\n")
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
  if (dyn) {
    if (!is.null(locations)) {
      if (!is.matrix(locations) | all(dim(locations) != c(p, n_groups))) stop("`locations` must be a matrix with p rows and one column per group\n")
    }
    if (!is.null(scales)) {
      if (!is.matrix(scales) | all(dim(scales) != c(p, n_groups))) stop("`scales` must be a matrix with p rows and one column per group\n")
    }
    if (!is.null(polynomial_coef)) {
      if (!is.array(polynomial_coef) | any(dim(polynomial_coef) != c(p, d, n_groups))) stop("`polynomial_coef` must be an array with p rows, d columns ane one plane per group\n")
    }
    if (p < 2 & intercept & dynamic) warning("`p = 1` but both an intercept and a dynamic AR coefficient are specified\n")
  }
}
