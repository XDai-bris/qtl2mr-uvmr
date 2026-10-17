qtl2mr_plot_errorbarh <- function(x, lo, hi, y, col = "black", lwd = 2) {
  graphics::segments(lo, y, hi, y, col = col, lwd = lwd)
  tick <- 0.06
  graphics::segments(lo, y - tick, lo, y + tick, col = col, lwd = lwd)
  graphics::segments(hi, y - tick, hi, y + tick, col = col, lwd = lwd)
  graphics::points(x, y, pch = 16, col = col)
}

plot_per_locus_so_forest <- function(per_locus,
                                     estimates = NULL,
                                     global_theta = NULL,
                                     max_loci = 40L,
                                     top_n = NULL,
                                     order_by = c("F", "position", "effect", "se", "weight"),
                                     main = "Per-locus SO estimates",
                                     xlab = "MR effect",
                                     ...) {
  order_by <- match.arg(order_by)
  if (!is.null(top_n)) max_loci <- top_n
  required <- c("marker", "theta", "profile_ci_lower", "profile_ci_upper")
  miss <- setdiff(required, names(per_locus))
  if (length(miss)) stop("per_locus is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  plot_data <- per_locus[is.finite(per_locus$theta), , drop = FALSE]
  if (!nrow(plot_data)) {
    stop("per_locus has no finite theta values to plot", call. = FALSE)
  }
  if (identical(order_by, "F") && "F_X_block" %in% names(plot_data)) {
    plot_data <- plot_data[order(plot_data$F_X_block, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if (identical(order_by, "effect")) {
    plot_data <- plot_data[order(abs(plot_data$theta), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if (identical(order_by, "se") && "se" %in% names(plot_data)) {
    plot_data <- plot_data[order(plot_data$se, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  } else if (identical(order_by, "weight") && "weight" %in% names(plot_data)) {
    plot_data <- plot_data[order(plot_data$weight, decreasing = TRUE), , drop = FALSE]
  } else if (identical(order_by, "position") && "locus_index" %in% names(plot_data)) {
    plot_data <- plot_data[order(plot_data$locus_index), , drop = FALSE]
  }
  if (!is.null(max_loci) && nrow(plot_data) > max_loci) plot_data <- utils::head(plot_data, max_loci)
  plot_data <- plot_data[rev(seq_len(nrow(plot_data))), , drop = FALSE]
  y <- seq_len(nrow(plot_data))
  x_values <- c(
    plot_data$theta,
    plot_data$profile_ci_lower,
    plot_data$profile_ci_upper,
    global_theta,
    if (!is.null(estimates) && "theta" %in% names(estimates)) estimates$theta else numeric()
  )
  xlim <- range(x_values, finite = TRUE)
  if (!all(is.finite(xlim))) stop("no finite MR effect values are available for the forest plot", call. = FALSE)
  if (diff(xlim) <= .Machine$double.eps) {
    span <- max(abs(xlim[[1L]]), 1)
    xlim <- xlim + c(-0.1, 0.1) * span
  }
  pad <- diff(xlim) * 0.08
  xlim <- xlim + c(-pad, pad)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 9, 4, 2))
  graphics::plot(
    plot_data$theta, y,
    xlim = xlim,
    ylim = c(0.5, length(y) + 0.5),
    yaxt = "n",
    xlab = xlab,
    ylab = "",
    pch = 16,
    main = main,
    ...
  )
  graphics::abline(v = 0, lty = 3, col = "grey70")
  if (!is.null(global_theta) && is.finite(global_theta)) graphics::abline(v = global_theta, lty = 2, col = "#00796b", lwd = 2)
  has_interval <- is.finite(plot_data$profile_ci_lower) & is.finite(plot_data$profile_ci_upper)
  if (any(has_interval)) {
    qtl2mr_plot_errorbarh(
      plot_data$theta[has_interval], plot_data$profile_ci_lower[has_interval],
      plot_data$profile_ci_upper[has_interval], y[has_interval], col = "#4c63d2"
    )
  }
  labels <- if ("chr" %in% names(plot_data)) paste0(plot_data$marker, " chr", plot_data$chr) else plot_data$marker
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.72)
  if (!is.null(estimates) && nrow(estimates)) {
    txt <- paste(sprintf("%s: %.3f", estimates$estimator, estimates$theta), collapse = "   ")
    graphics::mtext(txt, side = 3, line = 0.3, cex = 0.75, col = "grey25")
  }
  if (!is.null(max_loci) && nrow(per_locus) > nrow(plot_data)) {
    graphics::mtext(sprintf("Showing %d of %d loci ordered by %s", nrow(plot_data), nrow(per_locus), order_by),
                    side = 1, line = 3.6, cex = 0.72, col = "grey35")
  }
  invisible(plot_data)
}

plot_founder_residual_heatmap <- function(residual_scores,
                                          value = c("t_one_founder_vs_rest", "e", "z_diag"),
                                          row = c("marker", "trait"),
                                          main = "Founder residual diagnostics",
                                          label_threshold = 2,
                                          zlim = NULL,
                                          palette = NULL,
                                          ...) {
  value <- match.arg(value)
  row <- match.arg(row)
  if (!row %in% names(residual_scores)) row <- if ("trait" %in% names(residual_scores)) "trait" else "marker"
  required <- c(row, "founder", value)
  miss <- setdiff(required, names(residual_scores))
  if (length(miss)) stop("residual_scores is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  rows <- unique(as.character(residual_scores[[row]]))
  founders <- unique(as.character(residual_scores$founder))
  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(founders), dimnames = list(rows, founders))
  for (i in seq_len(nrow(residual_scores))) {
    mat[as.character(residual_scores[[row]][[i]]), as.character(residual_scores$founder[[i]])] <- residual_scores[[value]][[i]]
  }
  mat <- mat[rev(seq_len(nrow(mat))), , drop = FALSE]
  if (is.null(zlim)) {
    mx <- max(abs(mat), na.rm = TRUE)
    zlim <- c(-mx, mx)
  }
  if (is.null(palette)) palette <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(6, 8, 4, 5))
  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat),
    col = palette,
    zlim = zlim,
    axes = FALSE,
    xlab = "Founder",
    ylab = "",
    main = main,
    ...
  )
  graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.78)
  graphics::axis(2, at = seq_len(nrow(mat)), labels = rownames(mat), las = 1, cex.axis = 0.72)
  graphics::box()
  if (!is.null(label_threshold) && is.finite(label_threshold)) {
    idx <- which(abs(mat) >= label_threshold, arr.ind = TRUE)
    if (nrow(idx)) {
      graphics::text(idx[, 2], idx[, 1], sprintf("%.1f", mat[idx]), cex = 0.72)
    }
  }
  invisible(mat)
}

plot_founder_ratio_forest <- function(ratio_table,
                                      group = NULL,
                                      xlim = NULL,
                                      main = "Founder diagnostic ratios",
                                      xlab = "Founder-coordinate Wald ratio",
                                      ...) {
  required <- c("coordinate", "theta_diag", "ci_lower", "ci_upper", "weak_denominator")
  miss <- setdiff(required, names(ratio_table))
  if (length(miss)) stop("ratio_table is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  plot_data <- ratio_table
  if (!is.null(group)) {
    group_col <- if ("marker" %in% names(plot_data)) "marker" else if ("trait" %in% names(plot_data)) "trait" else NULL
    if (!is.null(group_col)) plot_data <- plot_data[as.character(plot_data[[group_col]]) == group, , drop = FALSE]
  }
  plot_data <- plot_data[is.finite(plot_data$theta_diag), , drop = FALSE]
  plot_data <- plot_data[rev(seq_len(nrow(plot_data))), , drop = FALSE]
  y <- seq_len(nrow(plot_data))
  if (is.null(xlim)) {
    stable <- plot_data[!plot_data$weak_denominator &
      is.finite(plot_data$ci_lower) & is.finite(plot_data$ci_upper), , drop = FALSE]
    x_for_limit <- if (nrow(stable)) {
      c(stable$ci_lower, stable$ci_upper, stable$theta_diag, plot_data$theta_ref, 0)
    } else {
      c(plot_data$theta_diag, plot_data$theta_ref, 0)
    }
    xlim <- range(x_for_limit, finite = TRUE)
    if (!all(is.finite(xlim))) xlim <- c(-1, 1)
    if (diff(xlim) <= .Machine$double.eps) xlim <- xlim + c(-1, 1)
    pad <- diff(xlim) * 0.15
    xlim <- xlim + c(-pad, pad)
  }
  cols <- ifelse(plot_data$weak_denominator, "grey55", "#00796b")
  pch <- ifelse(plot_data$weak_denominator, 1, 16)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 8, 4, 2))
  x_point <- pmin(pmax(plot_data$theta_diag, xlim[[1L]]), xlim[[2L]])
  graphics::plot(x_point, y, xlim = xlim, ylim = c(0.5, length(y) + 0.5),
                 yaxt = "n", xlab = xlab, ylab = "", pch = pch, col = cols, main = main, ...)
  graphics::abline(v = 0, lty = 3, col = "grey70")
  if ("theta_ref" %in% names(plot_data) && any(is.finite(plot_data$theta_ref))) {
    graphics::abline(v = plot_data$theta_ref[which(is.finite(plot_data$theta_ref))[1L]], lty = 2, col = "grey35")
  }
  lo <- pmax(plot_data$ci_lower, xlim[[1L]])
  hi <- pmin(plot_data$ci_upper, xlim[[2L]])
  ok <- is.finite(lo) & is.finite(hi) & lo <= hi
  graphics::segments(lo[ok], y[ok], hi[ok], y[ok], col = cols[ok], lwd = 2)
  span <- diff(xlim)
  arrow_dx <- span * 0.035
  left_clip <- is.finite(plot_data$ci_lower) & plot_data$ci_lower < xlim[[1L]]
  right_clip <- is.finite(plot_data$ci_upper) & plot_data$ci_upper > xlim[[2L]]
  if (any(left_clip)) {
    graphics::arrows(xlim[[1L]] + arrow_dx, y[left_clip], xlim[[1L]], y[left_clip],
                     length = 0.07, code = 2, col = cols[left_clip], lwd = 2)
  }
  if (any(right_clip)) {
    graphics::arrows(xlim[[2L]] - arrow_dx, y[right_clip], xlim[[2L]], y[right_clip],
                     length = 0.07, code = 2, col = cols[right_clip], lwd = 2)
  }
  graphics::points(x_point, y, pch = pch, col = cols)
  graphics::axis(2, at = y, labels = plot_data$coordinate, las = 1, tick = FALSE)
  invisible(plot_data)
}

