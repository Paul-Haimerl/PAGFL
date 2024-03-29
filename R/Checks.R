prelim_checks <- function(y, X, Z = NULL, index = NULL, n_periods = NULL, method = "PLS", const_coef = NULL) {
  if (is.null(n_periods) & is.null(index)) stop("Either supply cross-sectional and time index variables or in case of a balanced and ordered panel data set, the number of time periods n_periods\n")
  if (is.null(n_periods)) {
    if (length(index) != 2 | !is.character(index)) stop("Please supply a character vector holding two strings as the index\n")
    if (!all(index %in% colnames(y)) | !all(index %in% colnames(X))) stop("The passed index variables must be both present in y and X\n")
    if (NCOL(as.matrix(y)[, !(colnames(y) %in% index)]) != 1) stop("Please provide a univariate dependent variable y\n")
    if (NCOL(as.matrix(X)[, !(colnames(X) %in% index)]) == 0) stop("Please provide an explanatory variable in addition to the index variables in X\n")
  } else {
    if (NCOL(as.matrix(y)) != 1) stop("Please provide a univariate dependent variable y\n")
    if (NCOL(as.matrix(X)) == 0) stop("Please provide an explanatory variable in X\n")
  }
  if (method == "PGMM") {
    if (is.null(Z)) stop("PGMM requires a matrix of exogenous instruments Z\n")
  }
  if (!is.null(const_coef)) {
    if (!all(const_coef %in% colnames(X))) stop("Explanatory variables that are to be estimated with time-constant coefficients must be named columns in the regressor matrix X")
  }
}


second_checks <- function(N, index, n_periods, y, X, method = NULL, Z = NULL, p, min_group_frac, verbose, dyn, d = NULL, J = NULL) {
  if (is.null(min_group_frac)) min_group_frac <- 0
  if (min_group_frac > 1 | min_group_frac < 0) stop("Provide a min_group_frac between 0 and 1\n")
  if (min_group_frac >= .4 & verbose) warning("Large min_group_frac values may lead to all groups falling below the group cardinality threshold, in which case the hierarchical clustering algorithm cannot be employed\n")
  if (is.null(index)) {
    if (round(N) != N) stop("n_periods does not match the number of time periods of the dependent variable y\n")
    if (N * n_periods != nrow(X)) stop("The number of time periods of the dependent variable y and the predictor matrix X do not match\n")
    if (any(is.na(X))) stop("The predictor matrix X contains missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables\n")
    if (any(is.na(y))) stop("The dependent variable y contains missing values. In order to work with unbalanced panel data sets, supply cross-sectional and temporal indix variables\n")
  } else {
    if (verbose) {
      if (nrow(y) != N * n_periods | nrow(X) != N * n_periods) warning("The panel data set is unbalanced\n")
    }
  }
  if (dyn) {
    checksDyn(d = d, J = J)
  } else {
    checksStat(X, method, Z, p, verbose)
  }
}

# Checks for the static PAGFL
checksStat <- function(X, method, Z, p, verbose) {
  if (!(method %in% c("PLS", "PGMM"))) stop("The estimation method must be either PLS or PGMM. Use PLS in case of (weakly) exogenous regressors and PGMM for endogenous regressors.\n")
  if (method == "PGMM") {
    if (ncol(Z) < p) stop(paste("Provide at least p =", p, "exogenous instruments Z"), "\n")
    if (nrow(Z) != nrow(X)) stop("The number of time periods of the instrument matrix Z does do not match the remaining data\n")
  } else {
    if (!is.null(Z) & verbose) warning("The instrument matrix Z is ignored by the PLS estimation algorithm. To instrument endogenous regressors using Z, specify method = PGMM\n")
  }
}

# Checks for the dynamic PAGFL
checksDyn <- function(d, J) {
  if (as.integer(d) != d) stop("An integer must be passed for d\n")
  if (as.integer(J) != J) stop("An integer must be passed for J\n")
}

# Checks for simulating data
simChecks <- function(dyn, N, n_groups, group_proportions, error_spec = NULL, alpha_0 = NULL, dyn_panel = FALSE, q = NULL, p) {
  if (N < n_groups) stop("Number of groups cannot exceed number of observations\n")
  if (p == 0) stop("Include at least one explanatory variable\n")
  if (!is.null(group_proportions)) {
    if (n_groups != length(group_proportions)) stop("Number of groups and group proportions are of different length\n")
    if (sum(group_proportions) != 1) stop("Group proportions must sum to 1\n")
  }
  if (!is.null(error_spec)) {
    if (!(error_spec %in% c("AR", "GARCH"))) stop("The individual error specificaion must be either AR, GARCH or NULL. Use AR in case of an AR(1) serial correlation, GARCH for an GARCH(1,1) innovation and NULL for iid errors.\n")
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
}
