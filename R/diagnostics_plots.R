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

plot_founder_coordinate_wald_ratio_heatmap <- function(ratio_table,
                                                       marker_order = NULL,
                                                       founders = LETTERS[1:8],
                                                       value = c("theta_delta", "theta_diag"),
                                                       zlim = NULL,
                                                       main = "All-loci founder-coordinate Wald ratio deviations",
                                                       weak_col = "grey86",
                                                       palette = NULL,
                                                       label_every = NULL,
                                                       ...) {
  value <- match.arg(value)
  required <- c("marker", "coordinate", "theta_diag", "theta_ref", "weak_denominator")
  miss <- setdiff(required, names(ratio_table))
  if (length(miss)) stop("ratio_table is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- ratio_table[is.finite(ratio_table$theta_diag), , drop = FALSE]
  if (is.null(marker_order)) {
    if ("locus_rank" %in% names(d)) {
      marker_order <- unique(as.character(d$marker[order(d$locus_rank)]))
    } else {
      marker_order <- unique(as.character(d$marker))
    }
  } else {
    marker_order <- as.character(marker_order)
  }
  founders <- as.character(founders)
  d <- d[as.character(d$marker) %in% marker_order & as.character(d$coordinate) %in% founders, , drop = FALSE]
  if (!nrow(d)) stop("No founder-coordinate rows remain after filtering", call. = FALSE)

  mat <- matrix(NA_real_, nrow = length(marker_order), ncol = length(founders),
                dimnames = list(marker_order, founders))
  weak <- matrix(FALSE, nrow = length(marker_order), ncol = length(founders),
                 dimnames = list(marker_order, founders))
  chr <- rep(NA_character_, length(marker_order))
  names(chr) <- marker_order
  for (i in seq_len(nrow(d))) {
    marker <- as.character(d$marker[[i]])
    founder <- as.character(d$coordinate[[i]])
    val <- if (identical(value, "theta_delta")) d$theta_diag[[i]] - d$theta_ref[[i]] else d$theta_diag[[i]]
    mat[marker, founder] <- val
    weak[marker, founder] <- isTRUE(d$weak_denominator[[i]])
    if ("chr" %in% names(d)) chr[[marker]] <- as.character(d$chr[[i]])
  }

  mat_color <- mat
  mat_color[weak] <- NA_real_
  stable_values <- as.numeric(mat_color[is.finite(mat_color)])
  if (is.null(zlim)) {
    if (length(stable_values)) {
      cap <- as.numeric(stats::quantile(abs(stable_values), 0.98, na.rm = TRUE, names = FALSE))
      if (!is.finite(cap) || cap <= 0) cap <- max(abs(stable_values), 1, na.rm = TRUE)
      zlim <- if (identical(value, "theta_delta")) c(-cap, cap) else range(stable_values, finite = TRUE)
    } else {
      zlim <- c(-1, 1)
    }
  }
  if (diff(zlim) <= .Machine$double.eps) zlim <- zlim + c(-1, 1)
  mat_color <- pmax(pmin(mat_color, zlim[[2L]]), zlim[[1L]])
  if (is.null(palette)) {
    palette <- if (identical(value, "theta_delta")) {
      grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
    } else {
      grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(101)
    }
  }

  mat_plot <- mat_color[rev(seq_len(nrow(mat_color))), , drop = FALSE]
  weak_plot <- weak[rownames(mat_plot), , drop = FALSE]
  chr <- chr[rownames(mat_plot)]
  row_labels <- rownames(mat_plot)
  if (any(!is.na(chr))) row_labels <- paste0(row_labels, " chr", chr)
  if (is.null(label_every)) label_every <- max(1L, ceiling(nrow(mat_plot) / 32L))
  y_at <- seq_len(nrow(mat_plot))
  y_labs <- ifelse((y_at - 1L) %% label_every == 0L, row_labels, "")

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(1)
    graphics::par(op)
  }, add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1L), widths = c(5.4, 0.8))
  graphics::par(mar = c(5.5, 7.5, 4, 0.8))
  graphics::image(
    x = seq_len(ncol(mat_plot)),
    y = seq_len(nrow(mat_plot)),
    z = t(mat_plot),
    col = palette,
    zlim = zlim,
    axes = FALSE,
    xlab = "Founder",
    ylab = "",
    main = main,
    ...
  )
  graphics::axis(1, at = seq_len(ncol(mat_plot)), labels = colnames(mat_plot), las = 2)
  graphics::axis(2, at = y_at, labels = y_labs, las = 1, cex.axis = 0.45, tick = FALSE)
  idx <- which(weak_plot, arr.ind = TRUE)
  if (nrow(idx)) {
    graphics::rect(idx[, 2] - 0.5, idx[, 1] - 0.5, idx[, 2] + 0.5, idx[, 1] + 0.5,
                   col = grDevices::adjustcolor(weak_col, alpha.f = 0.9), border = NA)
    graphics::points(idx[, 2], idx[, 1], pch = 4, cex = 0.38, col = "grey35")
  }
  graphics::box()
  graphics::mtext("Grey cells mark weak exposure denominators", side = 1, line = 4.2, cex = 0.72, col = "grey25")

  graphics::par(mar = c(5.5, 0.4, 4, 2.4))
  zseq <- seq(zlim[[1L]], zlim[[2L]], length.out = length(palette))
  graphics::image(x = 1, y = zseq, z = matrix(zseq, nrow = 1L),
                  col = palette, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1, cex.axis = 0.75)
  legend_lab <- if (identical(value, "theta_delta")) "Wald ratio - global theta" else "Founder-coordinate Wald ratio"
  graphics::mtext(legend_lab, side = 4, line = 2.6, cex = 0.78)
  invisible(list(values = mat, weak_denominator = weak, zlim = zlim))
}

plot_founder_residual_null <- function(null_summary,
                                       group = NULL,
                                       main = "Founder residual proportional-null check",
                                       ylab = "Absolute founder residual score",
                                       ...) {
  d <- summarize_founder_residual_null(null_summary)
  required <- c("founder", "observed_abs_t", "null_p95_abs_t")
  miss <- setdiff(required, names(d))
  if (length(miss)) stop("null_summary is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.null(group)) {
    group_col <- if ("marker" %in% names(d)) "marker" else if ("trait" %in% names(d)) "trait" else NULL
    if (!is.null(group_col)) d <- d[as.character(d[[group_col]]) == group, , drop = FALSE]
  }
  d <- d[order(d$observed_abs_t, decreasing = TRUE), , drop = FALSE]
  x <- seq_len(nrow(d))
  ymax <- max(c(d$observed_abs_t, d$null_p95_abs_t, d$null_p99_abs_t %||% NA_real_), na.rm = TRUE)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(6, 5, 4, 2))
  graphics::plot(x, d$observed_abs_t, type = "b", pch = 16, xaxt = "n",
                 ylim = c(0, ymax * 1.08), xlab = "Founder", ylab = ylab, main = main, ...)
  graphics::axis(1, at = x, labels = d$founder, las = 2)
  graphics::points(x, d$null_p95_abs_t, pch = 1, col = "#c46a1b")
  graphics::lines(x, d$null_p95_abs_t, col = "#c46a1b", lty = 2)
  if ("null_p99_abs_t" %in% names(d)) {
    graphics::points(x, d$null_p99_abs_t, pch = 2, col = "#b2182b")
    graphics::lines(x, d$null_p99_abs_t, col = "#b2182b", lty = 3)
  }
  graphics::legend(
    "topright",
    legend = c("Observed", "Null 95%", if ("null_p99_abs_t" %in% names(d)) "Null 99%" else NULL),
    pch = c(16, 1, if ("null_p99_abs_t" %in% names(d)) 2 else NULL),
    lty = c(1, 2, if ("null_p99_abs_t" %in% names(d)) 3 else NULL),
    col = c("black", "#c46a1b", if ("null_p99_abs_t" %in% names(d)) "#b2182b" else NULL),
    bty = "n"
  )
  invisible(d)
}

plot_score_bin_variance <- function(bin_table,
                                    main = "Score-bin residual variance",
                                    xlab = "Score bin",
                                    ylab = "Residual variance",
                                    ...) {
  required <- c("score_bin", "residual_variance", "ci_lower", "ci_upper")
  miss <- setdiff(required, names(bin_table))
  if (length(miss)) stop("bin_table is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  bin_table <- bin_table[order(bin_table$score_bin), , drop = FALSE]
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  ylim <- range(c(bin_table$ci_lower, bin_table$ci_upper, bin_table$residual_variance), finite = TRUE)
  graphics::plot(bin_table$score_bin, bin_table$residual_variance, type = "b", pch = 16,
                 ylim = ylim, xaxt = "n", xlab = xlab, ylab = ylab, main = main, ...)
  graphics::axis(1, at = bin_table$score_bin, labels = bin_table$score_bin_label %||% bin_table$score_bin)
  graphics::segments(bin_table$score_bin, bin_table$ci_lower, bin_table$score_bin, bin_table$ci_upper, lwd = 2, col = "grey35")
  invisible(bin_table)
}

qtl2mr_is_xprs_heteroskedasticity_diagnostic <- function(x) {
  inherits(x, "qtl2mr_xprs_heteroskedasticity_diagnostic")
}

qtl2mr_xprs_plot_inputs <- function(global_tests,
                                    global_bins = NULL,
                                    loo_tests = NULL,
                                    localization_tests = NULL) {
  if (qtl2mr_is_xprs_heteroskedasticity_diagnostic(global_tests)) {
    x <- global_tests
    global_tests <- x$global_tests
    global_bins <- global_bins %||% x$global_bins
    loo_tests <- loo_tests %||% x$loo_tests
    localization_tests <- localization_tests %||% x$localization_tests
  }
  list(
    global_tests = global_tests,
    global_bins = global_bins,
    loo_tests = loo_tests,
    localization_tests = localization_tests
  )
}

qtl2mr_fmt_p <- function(p) {
  if (!is.finite(p)) return("p=NA")
  if (p < 1e-4) sprintf("p=%.2e", p) else sprintf("p=%.4f", p)
}

qtl2mr_xprs_trait_title <- function(global_tests, trait) {
  if (!"trait_role" %in% names(global_tests)) return(paste0(trait, " residual QC"))
  idx <- match(trait, as.character(global_tests$trait))
  role <- if (is.na(idx)) NA_character_ else as.character(global_tests$trait_role[[idx]])
  if (identical(role, "exposure")) return(paste0("Exposure residual QC: ", trait))
  if (identical(role, "outcome")) return(paste0("Outcome residual QC: ", trait))
  paste0(trait, " residual QC")
}

qtl2mr_plot_xprs_global_p_panel <- function(global_tests,
                                            main = "Primary multi-locus local founder exposure score test") {
  required <- c("trait", "joint_p")
  miss <- setdiff(required, names(global_tests))
  if (length(miss)) stop("global_tests is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- global_tests[is.finite(global_tests$joint_p), , drop = FALSE]
  y <- -log10(pmax(d$joint_p, .Machine$double.xmin))
  cols <- ifelse(d$joint_p < 0.01, "#b2182b", "#4c63d2")
  ymax <- max(c(y, -log10(0.01)), finite = TRUE) * 1.18
  bp <- graphics::barplot(y, names.arg = d$trait, col = cols, ylim = c(0, ymax),
                          ylab = "-log10 joint p", main = main)
  graphics::abline(h = -log10(0.01), lty = 2, col = "#c46a1b")
  graphics::text(bp, y, labels = sprintf("p=%.2g", d$joint_p), pos = 3, cex = 0.78)
  invisible(d)
}

qtl2mr_plot_xprs_global_bin_trait_panel <- function(global_bins,
                                                   global_tests,
                                                   trait,
                                                   normalize = TRUE,
                                                   ylim = NULL,
                                                   main = NULL) {
  d <- global_bins[as.character(global_bins$trait) == trait, , drop = FALSE]
  if (!nrow(d)) stop("No global bin rows for trait: ", trait, call. = FALSE)
  d <- d[order(d$score_bin), , drop = FALSE]
  scale <- if (isTRUE(normalize)) max(mean(d$residual_variance, na.rm = TRUE), .Machine$double.eps) else 1
  y <- d$residual_variance / scale
  lo <- d$ci_lower / scale
  hi <- d$ci_upper / scale
  if (is.null(ylim)) {
    ylim <- range(c(y, lo, hi), finite = TRUE)
    pad <- diff(ylim) * 0.08
    if (!is.finite(pad) || pad <= 0) pad <- max(abs(ylim), 1, na.rm = TRUE) * 0.08
    ylim <- ylim + c(-pad, pad)
  }
  gt <- global_tests[match(trait, as.character(global_tests$trait)), , drop = FALSE]
  p_lab <- if (nrow(gt)) qtl2mr_fmt_p(gt$joint_p[[1L]]) else "p=NA"
  q5 <- if (nrow(gt) && "variance_ratio_q5_q1" %in% names(gt)) gt$variance_ratio_q5_q1[[1L]] else NA_real_
  mm <- if (nrow(gt) && "variance_ratio_max_min" %in% names(gt)) gt$variance_ratio_max_min[[1L]] else NA_real_
  subtitle <- sprintf("%s; Q5/Q1=%s; max/min=%s",
                      p_lab,
                      if (is.finite(q5)) sprintf("%.2f", q5) else "NA",
                      if (is.finite(mm)) sprintf("%.2f", mm) else "NA")
  graphics::plot(
    d$score_bin, y,
    type = "n",
    xaxt = "n",
    ylim = ylim,
    xlab = "Summed local founder exposure score quantile",
    ylab = if (isTRUE(normalize)) "Residual variance / mean" else "Residual variance",
    main = main %||% qtl2mr_xprs_trait_title(global_tests, trait)
  )
  if (isTRUE(normalize)) graphics::abline(h = 1, lty = 3, col = "grey65")
  graphics::axis(1, at = d$score_bin, labels = d$score_bin_label %||% paste0("Q", d$score_bin))
  graphics::segments(d$score_bin, lo, d$score_bin, hi, col = "grey35", lwd = 2)
  graphics::lines(d$score_bin, y, type = "b", pch = 16, col = "#00796b", lwd = 2.4)
  graphics::mtext(subtitle, side = 3, line = 0.15, cex = 0.72, col = "grey25")
  invisible(d)
}

qtl2mr_xprs_global_ylim <- function(global_bins, normalize = TRUE) {
  traits <- qtl2mr_score_trait_levels(global_bins)
  yy <- unlist(lapply(traits, function(tr) {
    d <- global_bins[as.character(global_bins$trait) == tr, , drop = FALSE]
    scale <- if (isTRUE(normalize)) max(mean(d$residual_variance, na.rm = TRUE), .Machine$double.eps) else 1
    c(d$residual_variance, d$ci_lower, d$ci_upper) / scale
  }))
  ylim <- range(yy, finite = TRUE)
  pad <- diff(ylim) * 0.08
  if (!is.finite(pad) || pad <= 0) pad <- max(abs(ylim), 1, na.rm = TRUE) * 0.08
  ylim + c(-pad, pad)
}

plot_multilocus_xprs_heteroskedasticity_global <- function(global_tests,
                                                           global_bins = NULL,
                                                           normalize = TRUE,
                                                           traits = NULL,
                                                           main = "Global multi-locus local founder exposure score residual-variance diagnostic",
                                                           ...) {
  inputs <- qtl2mr_xprs_plot_inputs(global_tests, global_bins = global_bins)
  global_tests <- inputs$global_tests
  global_bins <- inputs$global_bins
  if (!is.null(traits)) {
    traits <- as.character(traits)
    global_tests <- global_tests[as.character(global_tests$trait) %in% traits, , drop = FALSE]
    if (!is.null(global_bins)) {
      global_bins <- global_bins[as.character(global_bins$trait) %in% traits, , drop = FALSE]
    }
    if (!nrow(global_tests)) stop("No global xPRS test rows remain after trait filtering", call. = FALSE)
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  if (!is.null(global_bins) && nrow(global_bins)) {
    traits <- qtl2mr_score_trait_levels(global_bins)
    graphics::par(mfrow = c(1L, length(traits)), mar = c(5, 5, 3.5, 1.2), oma = c(0, 0, 2.2, 0))
    ylim <- qtl2mr_xprs_global_ylim(global_bins, normalize = normalize)
    for (tr in traits) {
      qtl2mr_plot_xprs_global_bin_trait_panel(global_bins, global_tests, tr, normalize = normalize, ylim = ylim)
    }
    graphics::mtext(main, outer = TRUE, cex = 1.08, font = 2)
  } else {
    graphics::par(mar = c(5, 5, 4, 2))
    qtl2mr_plot_xprs_global_p_panel(global_tests, main = main)
  }
  invisible(inputs)
}

qtl2mr_plot_xprs_bin_panel <- function(global_bins,
                                       normalize = TRUE,
                                       main = "Multi-locus local founder exposure score-bin variance") {
  if (qtl2mr_is_xprs_heteroskedasticity_diagnostic(global_bins)) {
    global_bins <- global_bins$global_bins
  }
  required <- c("trait", "score_bin", "residual_variance", "ci_lower", "ci_upper")
  miss <- setdiff(required, names(global_bins))
  if (length(miss)) stop("global_bins is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  traits <- qtl2mr_score_trait_levels(global_bins)
  cols <- c("#4c63d2", "#00796b", "#b2182b", "#c46a1b")
  vals <- lapply(seq_along(traits), function(i) {
    d <- global_bins[as.character(global_bins$trait) == traits[[i]], , drop = FALSE]
    d <- d[order(d$score_bin), , drop = FALSE]
    scale <- if (isTRUE(normalize)) max(mean(d$residual_variance, na.rm = TRUE), .Machine$double.eps) else 1
    list(
      trait = traits[[i]],
      x = d$score_bin,
      labels = d$score_bin_label %||% paste0("Q", d$score_bin),
      y = d$residual_variance / scale,
      lo = d$ci_lower / scale,
      hi = d$ci_upper / scale,
      col = cols[[((i - 1L) %% length(cols)) + 1L]]
    )
  })
  ylim <- range(unlist(lapply(vals, function(v) c(v$lo, v$hi, v$y))), finite = TRUE)
  graphics::plot(vals[[1L]]$x, vals[[1L]]$y, type = "n", xaxt = "n", ylim = ylim,
                 xlab = "Local founder exposure score quantile", ylab = if (isTRUE(normalize)) "Variance / mean" else "Residual variance",
                 main = main)
  graphics::axis(1, at = vals[[1L]]$x, labels = vals[[1L]]$labels)
  if (isTRUE(normalize)) graphics::abline(h = 1, lty = 3, col = "grey65")
  for (v in vals) {
    graphics::segments(v$x, v$lo, v$x, v$hi, col = v$col, lwd = 2)
    graphics::lines(v$x, v$y, type = "b", pch = 16, col = v$col, lwd = 2)
  }
  graphics::legend("topleft", legend = traits, col = vapply(vals, `[[`, character(1), "col"),
                   pch = 16, lty = 1, bty = "n")
  invisible(global_bins)
}

plot_multilocus_xprs_bin_variance <- function(global_bins,
                                              normalize = TRUE,
                                              main = "Multi-locus local founder exposure score-bin residual variance",
                                              ...) {
  inputs <- qtl2mr_xprs_plot_inputs(global_bins)
  global_bins <- inputs$global_bins %||% global_bins
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  qtl2mr_plot_xprs_bin_panel(global_bins, normalize = normalize, main = main)
}

qtl2mr_plot_xprs_loo_panel <- function(loo_tests,
                                       global_tests,
                                       main = "Leave-one-locus-out local founder exposure score sensitivity") {
  required <- c("omitted_marker", "trait", "locus_rank", "joint_p")
  miss <- setdiff(required, names(loo_tests))
  if (length(miss)) stop("loo_tests is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  traits <- qtl2mr_score_trait_levels(loo_tests)
  cols <- c("#4c63d2", "#00796b", "#b2182b", "#c46a1b")
  ymax <- max(-log10(pmax(c(loo_tests$joint_p, global_tests$joint_p), .Machine$double.xmin)), finite = TRUE) * 1.1
  graphics::plot(NA, NA, xlim = range(loo_tests$locus_rank, finite = TRUE), ylim = c(0, ymax),
                 xlab = "Omitted exposure locus rank", ylab = "-log10 joint p", main = main)
  for (i in seq_along(traits)) {
    tr <- traits[[i]]
    d <- loo_tests[as.character(loo_tests$trait) == tr, , drop = FALSE]
    d <- d[order(d$locus_rank), , drop = FALSE]
    col <- cols[[((i - 1L) %% length(cols)) + 1L]]
    y <- -log10(pmax(d$joint_p, .Machine$double.xmin))
    graphics::lines(d$locus_rank, y, type = "b", pch = 16, col = col, lwd = 1.6)
    full_p <- global_tests$joint_p[match(tr, as.character(global_tests$trait))]
    if (is.finite(full_p)) graphics::abline(h = -log10(max(full_p, .Machine$double.xmin)), col = col, lty = 2)
    full_y <- -log10(max(full_p, .Machine$double.xmin))
    drop <- full_y - y
    if (any(is.finite(drop))) {
      idx <- utils::head(order(drop, decreasing = TRUE), 3L)
      graphics::text(d$locus_rank[idx], y[idx], labels = d$omitted_marker[idx], pos = 3, cex = 0.62, col = col)
    }
  }
  graphics::legend("topright", legend = traits, col = cols[seq_along(traits)], pch = 16, lty = 1, bty = "n")
  invisible(loo_tests)
}

plot_multilocus_xprs_leave_one_out <- function(loo_tests,
                                               global_tests,
                                               main = "Leave-one-locus-out local founder exposure score heteroskedasticity sensitivity",
                                               ...) {
  if (qtl2mr_is_xprs_heteroskedasticity_diagnostic(loo_tests)) {
    inputs <- qtl2mr_xprs_plot_inputs(loo_tests)
    loo_tests <- inputs$loo_tests
    global_tests <- inputs$global_tests
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 5, 4, 2))
  qtl2mr_plot_xprs_loo_panel(loo_tests, global_tests, main = main)
}

qtl2mr_plot_localization_min_p_panel <- function(tests,
                                                 main = "Per-locus localization") {
  d <- qtl2mr_score_order_tests(tests)
  y <- -log10(pmax(d$joint_p, .Machine$double.xmin))
  cols <- ifelse(d$flag_fdr_0.05 %||% FALSE, "#b2182b",
                 ifelse(d$flag_nominal_0.01 %||% FALSE, "#c46a1b", "#4c63d2"))
  graphics::plot(d$locus_rank, y, pch = ifelse(d$flag_fdr_0.05 %||% FALSE, 17, 1),
                 col = cols, xlab = "Exposure locus rank", ylab = "-log10 joint p", main = main)
  graphics::abline(h = -log10(0.01), lty = 2, col = "#c46a1b")
  invisible(d)
}

plot_multilocus_xprs_heteroskedasticity_dashboard <- function(global_tests,
                                                              global_bins = NULL,
                                                              loo_tests = NULL,
                                                              localization_tests = NULL,
                                                              main = "Multi-locus local founder exposure score residual-variance QC",
                                                              ...) {
  inputs <- qtl2mr_xprs_plot_inputs(
    global_tests,
    global_bins = global_bins,
    loo_tests = loo_tests,
    localization_tests = localization_tests
  )
  global_tests <- inputs$global_tests
  global_bins <- inputs$global_bins
  loo_tests <- inputs$loo_tests
  localization_tests <- inputs$localization_tests
  if (is.null(global_bins) || !nrow(global_bins)) {
    stop("global_bins is required for the xPRS dashboard", call. = FALSE)
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  panels <- 2L + as.integer(!is.null(loo_tests) && nrow(loo_tests)) + as.integer(!is.null(localization_tests) && nrow(localization_tests))
  layout_mat <- switch(
    as.character(panels),
    "2" = matrix(1:2, nrow = 1L),
    "3" = matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE),
    matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE)
  )
  graphics::layout(layout_mat, heights = if (nrow(layout_mat) == 1L) 1 else c(1.08, 1))
  graphics::par(mar = c(4.4, 4.8, 3, 1.2), oma = c(0, 0, 2.4, 0))
  traits <- qtl2mr_score_trait_levels(global_bins)
  ylim <- qtl2mr_xprs_global_ylim(global_bins, normalize = TRUE)
  for (tr in traits[seq_len(min(2L, length(traits)))]) {
    qtl2mr_plot_xprs_global_bin_trait_panel(global_bins, global_tests, tr, normalize = TRUE, ylim = ylim)
  }
  if (!is.null(loo_tests) && nrow(loo_tests)) {
    qtl2mr_plot_xprs_loo_panel(loo_tests, global_tests, main = "Optional: leave-one-locus-out")
  }
  if (!is.null(localization_tests) && nrow(localization_tests)) {
    qtl2mr_plot_localization_min_p_panel(localization_tests, main = "Optional: per-locus localization")
  }
  graphics::mtext(main, outer = TRUE, cex = 1.1, font = 2)
  invisible(inputs)
}

qtl2mr_score_trait_levels <- function(x) {
  traits <- unique(as.character(x$trait))
  role_order <- c("exposure", "outcome")
  if ("trait_role" %in% names(x)) {
    roles <- unique(as.character(x$trait_role))
    ord <- match(roles, role_order)
    traits <- traits[order(match(as.character(x$trait_role[match(traits, x$trait)]), role_order), traits)]
  }
  traits
}

qtl2mr_score_order_tests <- function(tests) {
  tests <- tests[is.finite(tests$joint_p), , drop = FALSE]
  if ("locus_rank" %in% names(tests)) {
    tests <- tests[order(tests$locus_rank, tests$trait), , drop = FALSE]
  } else {
    tests <- tests[order(tests$marker, tests$trait), , drop = FALSE]
    tests$locus_rank <- match(tests$marker, unique(tests$marker))
  }
  tests
}

qtl2mr_plot_score_pvalue_panel <- function(tests,
                                           trait,
                                           main = NULL,
                                           ylab = "-log10 joint p",
                                           show_xlab = TRUE) {
  d <- qtl2mr_score_order_tests(tests[as.character(tests$trait) == trait, , drop = FALSE])
  if (!nrow(d)) stop("No rows for trait: ", trait, call. = FALSE)
  y <- -log10(pmax(d$joint_p, .Machine$double.xmin))
  cols <- ifelse(d$flag_fdr_0.05 %||% FALSE, "#b2182b",
                 ifelse(d$flag_nominal_0.01 %||% FALSE, "#c46a1b", "#4c63d2"))
  pch <- ifelse(d$flag_fdr_0.05 %||% FALSE, 17,
                ifelse(d$flag_nominal_0.01 %||% FALSE, 16, 1))
  ymax <- max(c(y, -log10(0.01)), finite = TRUE) * 1.08
  graphics::plot(
    d$locus_rank, y,
    type = "p",
    pch = pch,
    col = cols,
    xlab = if (show_xlab) "Exposure locus rank" else "",
    ylab = ylab,
    ylim = c(0, ymax),
    main = main %||% paste0(trait, " score-dependent variance scan")
  )
  graphics::abline(h = -log10(0.01), lty = 2, col = "#c46a1b")
  if ("joint_q_bh" %in% names(d) && any(d$flag_fdr_0.05, na.rm = TRUE)) {
    graphics::points(d$locus_rank[d$flag_fdr_0.05], y[d$flag_fdr_0.05], pch = 17, col = "#b2182b")
  }
  graphics::mtext("orange: p < 0.01; red triangles: BH q < 0.05", side = 3, line = 0.15, cex = 0.65, col = "grey35")
  invisible(d)
}

plot_score_heteroskedasticity_pvalue_scan <- function(tests,
                                                      traits = NULL,
                                                      main = "Score-dependent residual variance p-value scan",
                                                      ...) {
  required <- c("marker", "trait", "joint_p")
  miss <- setdiff(required, names(tests))
  if (length(miss)) stop("tests is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  tests <- qtl2mr_score_order_tests(tests)
  if (is.null(traits)) traits <- qtl2mr_score_trait_levels(tests)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(length(traits), 1L), mar = c(4, 5, 3, 1.5), oma = c(0, 0, 2, 0))
  for (tr in traits) {
    qtl2mr_plot_score_pvalue_panel(tests, tr, main = tr, ...)
  }
  graphics::mtext(main, outer = TRUE, cex = 1.05, font = 2)
  invisible(tests)
}

qtl2mr_score_bin_matrix <- function(bins, trait) {
  d <- bins[as.character(bins$trait) == trait, , drop = FALSE]
  if (!nrow(d)) stop("No bin rows for trait: ", trait, call. = FALSE)
  if ("locus_rank" %in% names(d)) {
    loci <- unique(d$marker[order(d$locus_rank)])
  } else {
    loci <- unique(as.character(d$marker))
  }
  bin_ids <- sort(unique(as.integer(d$score_bin)))
  mat <- matrix(NA_real_, nrow = length(loci), ncol = length(bin_ids),
                dimnames = list(loci, paste0("Q", bin_ids)))
  for (i in seq_len(nrow(d))) {
    mat[as.character(d$marker[[i]]), paste0("Q", as.integer(d$score_bin[[i]]))] <- d$residual_variance[[i]]
  }
  row_mean <- rowMeans(mat, na.rm = TRUE)
  mat <- log2(sweep(mat, 1L, pmax(row_mean, .Machine$double.eps), "/"))
  mat
}

qtl2mr_plot_score_heatmap_panel <- function(bins,
                                            trait,
                                            zlim = NULL,
                                            main = NULL,
                                            palette = NULL) {
  mat <- qtl2mr_score_bin_matrix(bins, trait)
  if (is.null(zlim)) {
    mx <- stats::quantile(abs(as.numeric(mat)), 0.98, na.rm = TRUE, names = FALSE)
    if (!is.finite(mx) || mx <= 0) mx <- max(abs(mat), na.rm = TRUE)
    zlim <- c(-mx, mx)
  }
  if (is.null(palette)) palette <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1L, , drop = FALSE]),
    col = palette,
    zlim = zlim,
    axes = FALSE,
    xlab = "Score bin",
    ylab = "Exposure locus rank",
    main = main %||% trait
  )
  graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  labs_at <- pretty(seq_len(nrow(mat)), n = 6)
  labs_at <- labs_at[labs_at >= 1 & labs_at <= nrow(mat)]
  graphics::axis(2, at = nrow(mat) - labs_at + 1L, labels = labs_at, las = 1, cex.axis = 0.8)
  graphics::box()
  invisible(mat)
}

plot_score_bin_variance_heatmap <- function(bins,
                                            traits = NULL,
                                            main = "All-loci score-bin residual variance",
                                            zlim = NULL,
                                            palette = NULL,
                                            ...) {
  required <- c("marker", "trait", "score_bin", "residual_variance")
  miss <- setdiff(required, names(bins))
  if (length(miss)) stop("bins is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  if (is.null(traits)) traits <- qtl2mr_score_trait_levels(bins)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1L, length(traits)), mar = c(4, 5, 3, 2), oma = c(0, 0, 2, 0))
  mats <- lapply(traits, function(tr) qtl2mr_plot_score_heatmap_panel(bins, tr, zlim = zlim, main = tr, palette = palette))
  graphics::mtext(main, outer = TRUE, cex = 1.05, font = 2)
  invisible(mats)
}

qtl2mr_plot_score_ratio_panel <- function(tests,
                                          trait,
                                          main = NULL,
                                          ylab = "Variance ratio max/min",
                                          show_xlab = TRUE) {
  d <- qtl2mr_score_order_tests(tests[as.character(tests$trait) == trait, , drop = FALSE])
  if (!nrow(d)) stop("No rows for trait: ", trait, call. = FALSE)
  y <- d$variance_ratio_max_min
  cols <- ifelse(d$shape_sensitive_qc_warning %||% FALSE, "#b2182b",
                 ifelse(d$flag_nominal_0.01 %||% FALSE, "#c46a1b", "#00796b"))
  ymax <- stats::quantile(y, 0.98, na.rm = TRUE, names = FALSE)
  if (!is.finite(ymax)) ymax <- max(y, na.rm = TRUE)
  ymax <- max(ymax, 1, na.rm = TRUE) * 1.15
  graphics::plot(
    d$locus_rank, pmin(y, ymax),
    type = "p",
    pch = ifelse(d$shape_sensitive_qc_warning %||% FALSE, 17, 16),
    col = cols,
    xlab = if (show_xlab) "Exposure locus rank" else "",
    ylab = ylab,
    ylim = c(1, ymax),
    main = main %||% paste0(trait, " bin variance ratio")
  )
  graphics::abline(h = 1, lty = 3, col = "grey60")
  invisible(d)
}

plot_score_heteroskedasticity_dashboard <- function(tests,
                                                    bins,
                                                    traits = NULL,
                                                    main = "All-loci score-dependent residual variance QC",
                                                    zlim = NULL,
                                                    palette = NULL,
                                                    ...) {
  if (is.null(traits)) traits <- qtl2mr_score_trait_levels(tests)
  if (length(traits) != 2L) {
    stop("Dashboard currently expects exactly two traits", call. = FALSE)
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::layout(matrix(seq_len(6L), nrow = 3L, byrow = TRUE), heights = c(1, 1, 1.55))
  graphics::par(mar = c(3.6, 4.5, 2.6, 1.2), oma = c(0, 0, 2.2, 0))
  qtl2mr_plot_score_pvalue_panel(tests, traits[[1L]], main = paste0(traits[[1L]], ": p-value scan"))
  qtl2mr_plot_score_pvalue_panel(tests, traits[[2L]], main = paste0(traits[[2L]], ": p-value scan"))
  qtl2mr_plot_score_ratio_panel(tests, traits[[1L]], main = paste0(traits[[1L]], ": bin variance ratio"))
  qtl2mr_plot_score_ratio_panel(tests, traits[[2L]], main = paste0(traits[[2L]], ": bin variance ratio"))
  qtl2mr_plot_score_heatmap_panel(bins, traits[[1L]], zlim = zlim, main = paste0(traits[[1L]], ": log2 bin variance ratio"), palette = palette)
  qtl2mr_plot_score_heatmap_panel(bins, traits[[2L]], zlim = zlim, main = paste0(traits[[2L]], ": log2 bin variance ratio"), palette = palette)
  graphics::mtext(main, outer = TRUE, cex = 1.1, font = 2)
  invisible(list(tests = tests, bins = bins))
}

plot_score_bin_variance_grid <- function(tests,
                                         bins,
                                         top_n = 12L,
                                         normalize = TRUE,
                                         main = "Top score-bin residual variance patterns",
                                         ...) {
  required <- c("marker", "trait", "joint_p")
  miss <- setdiff(required, names(tests))
  if (length(miss)) stop("tests is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- tests
  if (!"joint_q_bh" %in% names(d)) d$joint_q_bh <- NA_real_
  if (!"flag_fdr_0.05" %in% names(d)) d$flag_fdr_0.05 <- FALSE
  if (!"flag_nominal_0.01" %in% names(d)) d$flag_nominal_0.01 <- is.finite(d$joint_p) & d$joint_p < 0.01
  flagged <- (d$flag_fdr_0.05 %||% FALSE) | (d$flag_nominal_0.01 %||% FALSE)
  d <- d[flagged & is.finite(d$joint_p), , drop = FALSE]
  if (!nrow(d)) d <- tests[is.finite(tests$joint_p), , drop = FALSE]
  d <- d[order(d$joint_q_bh %||% d$joint_p, d$joint_p), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n) && nrow(d) > top_n) d <- utils::head(d, top_n)
  ncol <- 3L
  nr <- ceiling(nrow(d) / ncol)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(nr, ncol), mar = c(3.2, 4.3, 2.8, 0.8), oma = c(0, 0, 2, 0))
  for (i in seq_len(nrow(d))) {
    bb <- bins[
      as.character(bins$marker) == as.character(d$marker[[i]]) &
        as.character(bins$trait) == as.character(d$trait[[i]]),
      ,
      drop = FALSE
    ]
    bb <- bb[order(bb$score_bin), , drop = FALSE]
    y <- bb$residual_variance
    ylab <- "Variance"
    if (isTRUE(normalize)) {
      y <- y / max(mean(y, na.rm = TRUE), .Machine$double.eps)
      ylab <- "Variance / mean"
    }
    show_ylab <- ((i - 1L) %% ncol) == 0L
    graphics::plot(bb$score_bin, y, type = "b", pch = 16, xaxt = "n",
                   xlab = "Score bin", ylab = if (show_ylab) ylab else "",
                   main = sprintf("%s %s\np=%.3g q=%.3g", d$marker[[i]], d$trait[[i]], d$joint_p[[i]], d$joint_q_bh[[i]]),
                   cex.main = 0.85, cex.lab = 0.82, cex.axis = 0.75,
                   ...)
    graphics::axis(1, at = bb$score_bin, labels = bb$score_bin_label %||% paste0("Q", bb$score_bin))
  }
  if (nrow(d) < nr * ncol) {
    for (i in seq_len(nr * ncol - nrow(d))) graphics::plot.new()
  }
  graphics::mtext(main, outer = TRUE, cex = 1.05, font = 2)
  invisible(d)
}

plot.qtl2mr_xprs_heteroskedasticity_diagnostic <- function(x,
                                                           type = c("global", "dashboard", "loo", "bin", "localization"),
                                                           top_n = 12L,
                                                           ...) {
  type <- match.arg(type)
  switch(
    type,
    global = plot_multilocus_xprs_heteroskedasticity_global(x, ...),
    dashboard = plot_multilocus_xprs_heteroskedasticity_dashboard(x, ...),
    loo = plot_multilocus_xprs_leave_one_out(x, ...),
    bin = plot_multilocus_xprs_bin_variance(x, ...),
    localization = {
      if (is.null(x$localization_tests) || is.null(x$localization_bins)) {
        stop("Diagnostic object has no localization tests/bins", call. = FALSE)
      }
      plot_score_bin_variance_grid(x$localization_tests, x$localization_bins, top_n = top_n, ...)
    }
  )
  invisible(x)
}
