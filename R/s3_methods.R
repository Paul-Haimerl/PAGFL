#' @name pagfl
#' @param x of class \code{pagfl}.
#' @method print pagfl
#' @export
print.pagfl <- function(x, ...) {
  cat(paste("Groups:", x$groups$n_groups), "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
}

#' @name pagfl
#' @param x of class \code{pagfl}.
#' @method formula pagfl
#' @export
formula.pagfl <- function(x, ...) {
  x$args$formula
}

#' @name pagfl
#' @param object of class \code{pagfl}.
#' @method df.residual pagfl
#' @export
df.residual.pagfl <- function(object, ...) {
  length(object$args$labs$t) - length(unique(object$args$labs$i)) - ncol(object$coefficients) * object$groups$n_groups
}

#' @name pagfl
#' @param object of class \code{pagfl}.
#' @method summary pagfl
#' @export
summary.pagfl <- function(object, ...) {
  tmp <- object[c("call", "residuals", "coefficients", "groups", "IC", "convergence", "args", "model")]
  k <- ncol(tmp$coefficients)
  N <- length(unique(object$args$labs$i))
  i_index <- as.numeric(factor(object$args$labs$i))
  measures_vec <- fitMeasures(
    N = N, k = k, y = object$model[[1]], i_index = i_index,
    method = object$args$method, msr = tmp$IC$msr
  )
  out <- c(tmp, r.df = round(measures_vec[1]), r.squared = measures_vec[2], adj.r.squared = measures_vec[3], r.se = measures_vec[4], msr = tmp$IC$msr)
  class(out) <- "summary.pagfl"
  return(out)
}

#' @export
print.summary.pagfl <- function(x, ...) {
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
  cat("\nCoefficients:\n ")
  print(round(x$coefficients, 5))
  cat("\nResidual standard error:", round(x$r.se, 5), "on", x$r.df, "degrees of freedom\n")
  cat("Mean squared error:", round(x$IC$msr, 5))
  cat("\nMultiple R-squared:", paste0(round(x$r.squared, 5), ","), "Adjusted R-squared:", round(x$adj.r.squared, 5), "\n")
}

#' @name pagfl
#' @param object of class \code{pagfl}.
#' @method coef pagfl
#' @export
coef.pagfl <- function(object, ...) {
  coef_mat <- object$coefficients
  groups_hat <- object$groups$groups
  beta_mat <- coef_mat[groups_hat, ]
  row.names(beta_mat) <- names(groups_hat)
  return(beta_mat)
}


#' @name pagfl
#' @param object An object of class \code{pagfl}.
#' @method residuals pagfl
#' @export
residuals.pagfl <- function(object, ...) {
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


#' @name pagfl
#' @param object of class \code{pagfl}.
#' @method fitted pagfl
#' @export
fitted.pagfl <- function(object, ...) {
  fitted_vec <- object$fitted
  i_index <- object$args$labs$i
  t_index <- object$args$labs$t
  if (is.character(t_index)) t_index <- as.numeric(factor(t_index))
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
    fit_plot <- gen_fit_plot_pagfl(plot_df = plot_df, y_name = y_name, col_map = col_map)
    print(fit_plot)
  }
  return(fitted_df)
}

gen_fit_plot_pagfl <- function(plot_df, y_name, col_map) {
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
      panel.background = ggplot2::element_rect(colour = NA, fill = NA),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(colour = "#dcdcdc", fill = "#dcdcdc"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.margin = ggplot2::margin(t = -25),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1))
    )
}
