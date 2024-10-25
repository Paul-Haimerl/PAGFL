#' @name grouped_tv_plm
#' @param object of class \code{tv_gplm}.
#' @method summary tv_gplm
#' @export
summary.tv_gplm <- function(object, ...) {
  tmp <- object[c("call", "residuals", "coefficients", "groups", "IC", "args", "model")]
  k_tv <- ncol(tmp$coefficients$tv) * (tmp$args$d + tmp$args$M + 1)
  k_const <- ifelse(is.null(ncol(tmp$coefficients$const)), 0, ncol(tmp$coefficients$const))
  k <- k_tv + k_const
  N <- length(unique(object$args$labs$i))
  i_index <- as.numeric(factor(object$args$labs$i))
  measures_vec <- fitMeasures(N = N, k = k, y = object$model[[1]], i_index = i_index, method = "PLS", msr = tmp$IC$msr)
  out <- c(tmp, r.df = round(measures_vec[1]), r.squared = measures_vec[2], adj.r.squared = measures_vec[3], r.se = measures_vec[4], msr = tmp$IC$msr)
  class(out) <- "summary.tv_gplm"
  return(out)
}

#' @name grouped_tv_plm
#' @param x of class \code{tv_gplm}.
#' @method formula tv_gplm
#' @export
formula.tv_gplm <- function(x, ...) {
  x$args$formula
}

#' @name grouped_tv_plm
#' @param object of class \code{tv_gplm}.
#' @method df.residual tv_gplm
#' @export
df.residual.tv_gplm <- function(object, ...) {
  M <- object$args$M + object$args$d + 1
  df_fe <- length(unique(object$args$labs$i))
  p <- max(ncol(object$coefficients$tv), 0)
  p_const <- max(ncol(object$coefficients$const), 0)
  length(object$args$labs$t) - df_fe - (p * M + p_const) * object$groups$n_groups
}

#' @name grouped_tv_plm
#' @param x of class \code{tv_gplm}.
#' @method print tv_gplm
#' @export
print.tv_gplm <- function(x, ...) {
  cat(paste("Groups:", x$groups$n_groups), "\n")
  cat("\nCall:\n")
  print(x$call)
}

#' @export
print.summary.tv_gplm <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  unique_i <- unique(x$args$labs$i)
  N <- length(unique_i)
  n_periods <- length(unique(x$args$labs$t))
  lab_mat <- cbind(i = x$args$labs$i, t = x$args$labs$t)
  min_max_t <- stats::quantile(sapply(unique_i, function(i) length(lab_mat[lab_mat[, 1] == i, 2])), probs = c(0, 1))
  if (min_max_t[1] == min_max_t[2]) {
    balanced <- "Balanced"
    t_range <- min_max_t[2]
  } else {
    balanced <- "Unbalanced"
    t_range <- paste0(min_max_t, collapse = "-")
  }
  cat(paste0("\n", balanced, " panel: N = ", N, ", T = ", t_range, ", obs = ", length(x$residuals), "\n\n"))
  cat("\nInformation criterion:\n")
  print(round(c(IC = x$IC$IC), 5))
  cat("\nResiduals:\n")
  resid_vec <- x$residuals
  quantile_vec <- round(stats::quantile(resid_vec, probs = c(0, .25, .5, .75, 1)), 5)
  names(quantile_vec) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(quantile_vec)
  if (x$groups$n_groups > 1) {
    cat(paste0("\n", x$groups$n_groups), "groups:\n")
    print(x$groups$groups)
  } else {
    cat("\n1 group\n")
  }
  if (!is.null(x$coefficients$const)) {
    cat("\nConstant coefficients:\n ")
    print(round(x$coefficients$const, 5))
  }
  cat("\nResidual standard error:", round(x$r.se, 5), "on", x$r.df, "degrees of freedom\n")
  cat("Mean squared error:", round(x$IC$msr, 5))
  cat("\nMultiple R-squared:", paste0(round(x$r.squared, 5), ","), "Adjusted R-squared:", round(x$adj.r.squared, 5), "\n")
  # Plot of the functional coefficient
  p <- dim(x$coefficients$tv)[2]
  var_names <- colnames(x$coefficients$tv)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  n_periods <- dim(x$coefficients$tv)[1]
  group_names <- factor(as.numeric(gsub("Group\\s+", "", unlist(dimnames(x$coefficients$tv)[3]))))
  coef_df <- data.frame(
    coef = c(x$coefficients$tv),
    var_name = rep(rep(var_names, each = n_periods), x$groups$n_groups),
    index = rep(1:n_periods, x$groups$n_groups * p),
    Group = rep(group_names, each = p * n_periods)
  )
  legend_position <- ifelse(x$groups$n_groups <= 10, "bottom", "none")
  coef_plot <- gen_coef_plot_tvpagfl(coef_df, legend_position)
  print(coef_plot)
}

#' @name grouped_tv_plm
#' @param object of class \code{tv_gplm}.
#' @method coef tv_gplm
#' @export
coef.tv_gplm <- function(object, ...) {
  tv <- object$coefficients$tv
  if (!is.null(object$coefficients$const)) {
    const <- object$coefficients$const
    tmp <- array(NA, dim = dim(tv) + c(0, ncol(const), 0))
    tmp[, 1:dim(tv)[2], ] <- tv
    const_array <- array(rep(c(t(const)), each = dim(tv)[1]), dim = c(dim(tv)[1], ncol(const), dim(tv)[3]))
    tmp[, (dim(tv)[2] + 1):(dim(tv)[2] + ncol(const)), ] <- const_array
    out <- tmp[, , object$groups$groups]
    dimnames(out) <- list(dimnames(tv)[[1]], c(dimnames(tv)[[2]], colnames(const)), names(object$groups$groups))
  } else {
    out <- tv[, , object$groups$groups]
    if (is.matrix(out)) out <- array(c(out), dim = c(nrow(out), 1, ncol(out)))
    dimnames(out) <- list(dimnames(tv)[[1]], dimnames(tv)[[2]], names(object$groups$groups))
  }
  return(out)
}


#' @name grouped_tv_plm
#' @param object of class \code{tv_gplm}.
#' @method residuals tv_gplm
#' @export
residuals.tv_gplm <- function(object, ...) {
  resid_vec <- object$residuals
  i_index <- object$args$labs$i
  t_index <- object$args$labs$t
  resid_df <- data.frame(
    residuals = resid_vec,
    i_index = i_index,
    t_index = t_index
  )
  colnames(resid_df)[-1] <- object$args$labs$index
  return(resid_df)
}


#' @name grouped_tv_plm
#' @param object of class \code{tv_gplm}.
#' @method fitted tv_gplm
#' @export
fitted.tv_gplm <- function(object, ...) {
  fitted_vec <- object$fitted
  i_index <- object$args$labs$i
  t_index <- object$args$labs$t
  fitted_df <- data.frame(
    fit = fitted_vec,
    i_index = i_index,
    t_index = t_index
  )
  plot_df <- fitted_df
  colnames(fitted_df)[-1] <- object$args$labs$index
  # Plot the fit if feasible
  if (length(unique(i_index)) <= 20) {
    if (!is.numeric(t_index)) {
      suppressWarnings(t_index <- as.numeric(t_index))
      if (all(is.na(t_index))) t_index <- as.integer(factor(object$args$labs$t))
      plot_df$t_index <- t_index
    }
    plot_df$i_index <- as.character(plot_df$i_index)
    plot_df$y <- object$model[[1]]
    plot_df <- plot_df[order(plot_df$i_index), ]
    y_name <- colnames(object$model)[1]
    col_map <- c("red", "black")
    names(col_map) <- c("fit", y_name)
    fit_plot <- gen_fit_plot_tvpagfl(plot_df, y_name, col_map)
    print(fit_plot)
  }
  return(fitted_df)
}
