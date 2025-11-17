#' @name fuse_time
#' @param object of class \code{fusetime}.
#' @method summary fusetime
#' @export
summary.fusetime <- function(object, ...) {
  tmp <- object[c("call", "residuals", "coefficients", "groups", "IC", "convergence", "args", "model")]
  k_tv <- ncol(tmp$coefficients$tv) * (tmp$args$d + tmp$args$M + 1)
  k_const <- ifelse(is.null(ncol(tmp$coefficients$const)), 0, ncol(tmp$coefficients$const))
  k <- k_tv + k_const
  N <- length(unique(object$args$labs$i))
  i_index <- as.numeric(factor(object$args$labs$i))
  measures_vec <- fitMeasures(N = N, k = k, y = object$model[[1]], i_index = i_index, method = "PLS", msr = tmp$IC$msr)
  out <- c(tmp, r.df = round(measures_vec[1]), r.squared = measures_vec[2], adj.r.squared = measures_vec[3], r.se = measures_vec[4], msr = tmp$IC$msr)
  class(out) <- "summary.fusetime"
  return(out)
}

#' @name fuse_time
#' @param x of class \code{fusetime}.
#' @method formula fusetime
#' @export
formula.fusetime <- function(x, ...) {
  x$args$formula
}

#' @name fuse_time
#' @param object of class \code{fusetime}.
#' @method df.residual fusetime
#' @export
df.residual.fusetime <- function(object, ...) {
  M <- object$args$M + object$args$d + 1
  df_fe <- length(unique(object$args$labs$i))
  p <- max(ncol(object$coefficients$tv), 0)
  p_const <- max(ncol(object$coefficients$const), 0)
  length(object$args$labs$t) - df_fe - (p * M + p_const) * object$groups$n_groups
}

#' @name fuse_time
#' @param x of class \code{fusetime}.
#' @method print fusetime
#' @export
print.fusetime <- function(x, ...) {
  cat(paste("Groups:", x$groups$n_groups), "\n")
  cat("\nCall:\n")
  print(x$call)
}

#' @export
print.summary.fusetime <- function(x, ...) {
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
  cat("Convergence reached:\n")
  cat(x$convergence$convergence, paste0("(", x$convergence$iter, " iterations)\n"))
  cat("\nInformation criterion:\n")
  ic_vec <- c(IC = x$IC$IC, lambda = x$IC$lambda)
  print(round(ic_vec, 5))
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
  t_index <- rownames(x$coefficients$tv)
  suppressWarnings(t_index <- as.numeric(t_index))
  if (all(is.na(t_index))) t_index <- rep(1:n_periods, x$groups$n_groups * p)
  coef_df <- data.frame(
    coef = c(x$coefficients$tv),
    var_name = rep(rep(var_names, each = n_periods), x$groups$n_groups),
    index = t_index,
    Group = rep(group_names, each = p * n_periods)
  )
  legend_position <- ifelse(x$groups$n_groups <= 10, "bottom", "none")
  coef_plot <- gen_coef_plot_fusetime(coef_df, legend_position)
  print(coef_plot)
}

gen_coef_plot_fusetime <- function(coef_df, legend_position) {
  coef <- coef_df$coef
  index <- coef_df$index
  Group <- coef_df$Group
  var_name <- coef_df$var_name
  ggplot2::ggplot(coef_df, ggplot2::aes(x = index, y = coef)) +
    ggplot2::geom_line(ggplot2::aes(color = Group), na.rm = TRUE) +
    ggplot2::facet_grid(rows = ggplot2::vars(var_name), scales = "free") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(colour = NA, fill = NA),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1), colour = "black"),
      axis.ticks.length = ggplot2::unit(-6, "pt"),
      legend.position = legend_position,
      strip.background = ggplot2::element_rect(colour = NA, fill = NA),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      legend.margin = ggplot2::margin(t = -25),
      legend.title = ggplot2::element_text(face = "italic", colour = "black", size = ggplot2::rel(1)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1), colour = "black")
    )
}


#' @name fuse_time
#' @param object of class \code{fusetime}.
#' @method coef fusetime
#' @export
coef.fusetime <- function(object, ...) {
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


#' @name fuse_time
#' @param object of class \code{fusetime}.
#' @method residuals fusetime
#' @export
residuals.fusetime <- function(object, ...) {
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


#' @name fuse_time
#' @param object of class \code{fusetime}.
#' @method fitted fusetime
#' @export
fitted.fusetime <- function(object, ...) {
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
    fit_plot <- gen_fit_plot_fusetime(plot_df, y_name, col_map)
    print(fit_plot)
  }
  return(fitted_df)
}

gen_fit_plot_fusetime <- function(plot_df, y_name, col_map) {
  t_index <- plot_df$t_index
  y <- plot_df$y
  fit <- plot_df$fit
  ggplot2::ggplot(plot_df, ggplot2::aes(x = t_index)) +
    ggplot2::geom_line(ggplot2::aes(y = y, color = y_name)) +
    ggplot2::geom_line(ggplot2::aes(y = fit, color = "fit")) +
    ggplot2::facet_wrap(~i_index, scales = "free") +
    ggplot2::xlab("") +
    ggplot2::ylab(y_name) +
    ggplot2::scale_color_manual(values = col_map) +
    ggplot2::labs(colour = "") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(colour = NA, fill = NA),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1), colour = "black"),
      axis.ticks.length = ggplot2::unit(-6, "pt"),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(colour = NA, fill = NA),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      legend.margin = ggplot2::margin(t = -25),
      legend.title = ggplot2::element_text(colour = "black", size = ggplot2::rel(1)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1), colour = "black")
    )
}
