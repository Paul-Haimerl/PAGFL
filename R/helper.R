tv_pagfl_output <- function(out, fe_vec, p, p_const, n_periods, d, M, coef_rownames, regressor_names, const_coef, index, i_index, i_index_labs, t_index, model_data) {

  # Add FE to the estimated values
  fitted <- c(out$IC$fitted + fe_vec)
  # Prep the functional coefficients
  coef_raw <- out$estimOutput$alpha_hat
  # In case of time-constant regressors, separate them from the spline coefficients
  if (p_const > 0) {
    indx_p_c <- (ncol(coef_raw) - p_const + 1):ncol(coef_raw)
    alpha_const <- coef_raw[, indx_p_c]
    alpha_hat <- coef_raw[, -indx_p_c]
    if (out$estimOutput$K_hat == 1) {
      alpha_hat <- t(alpha_hat)
      alpha_const <- t(alpha_const)
    }
    if (p_const == 1) {
      alpha_const <- as.matrix(alpha_const)
    }
  } else {
    alpha_hat <- coef_raw
    alpha_const <- NULL
  }
  B <- bspline_system(1:n_periods, d, seq(1, n_periods, length.out = M + 2), TRUE)

  # Transform the spline coefficients to time-varying functional coefficients
  alpha_mat <- getTVAlpha(xi = alpha_hat, K_hat = out$estimOutput$K_hat, p = p, n_periods = n_periods, B = B)
  # Attach names to the output
  dimnames(alpha_mat) <- list(coef_rownames, regressor_names, paste0("Group ", 1:out$estimOutput$K_hat))
  if (!is.null(const_coef)) {
    dimnames(alpha_const) <- list(paste0("Group ", 1:out$estimOutput$K_hat), const_coef)
  }
  out$estimOutput$groups_hat <- c(out$estimOutput$groups_hat)
  names(out$estimOutput$groups_hat) <- unique(i_index_labs)
  # Omit estimates at the beginning and end of the observational period without support
  if (!is.null(index)) {
    alpha_hat_star <- delete_missing_t(
      i_index = i_index, t_index = t_index, K_hat = out$estimOutput$K_hat, groups_hat = out$estimOutput$groups_hat,
      alpha_hat = alpha_mat
    )
    dimnames(alpha_hat_star) <- dimnames(alpha_mat)
    alpha_mat <- alpha_hat_star
  }
  # De-mean a trend, if present
  if ("(Intercept)" %in% regressor_names) {
    trend_vec <- c(alpha_mat[, which(regressor_names == "(Intercept)"), ])
    trend_vec <- demeanIndVec(trend_vec, N = out$estimOutput$K_hat, i_index = rep(1:out$estimOutput$K_hat, each = n_periods))
    alpha_mat[, which(regressor_names == "(Intercept)"), ] <- array(trend_vec, dim = c(n_periods, 1, out$estimOutput$K_hat))
  }

  # Output
  out <- list(
    model = as.data.frame(model_data), coefficients = list(tv = alpha_mat, const = alpha_const), groups = list(n_groups = out$estimOutput$K_hat, groups = out$estimOutput$groups_hat), residuals = c(out$IC$resid), fitted = fitted,
    args = list(), IC = list(IC = out$IC$IC, lambda = NULL, msr = out$IC$msr), convergence = list(convergence = out$estimOutput$convergence, iter = out$estimOutput$iter), call = list()
  )
  return(out)
}
