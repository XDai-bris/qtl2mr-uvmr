qtl2mr_fit_from_result <- function(result) {
  if (inherits(result, "qtl2mr_analysis_result")) return(result$pilot$fit)
  if (is.list(result) && !is.null(result$locus_blocks) && !is.null(result$theta)) return(result)
  stop("Expected a qtl2mr_analysis_result or fit_multilocus_founder_mr result", call. = FALSE)
}

qtl2mr_blocks_from_result <- function(result) {
  if (inherits(result, "qtl2mr_analysis_result")) return(result$pilot$locus_blocks)
  if (is.list(result) && !is.null(result$locus_blocks)) return(result$locus_blocks)
  stop("Could not find locus_blocks", call. = FALSE)
}

loo_multilocus_founder_mr <- function(fit,
                                      parallel = FALSE,
                                      allow_cross_cov = TRUE,
                                      ridge = 1e-10,
                                      ...) {
  if (isTRUE(parallel)) {
    stop("parallel LOO is not implemented in v0.1", call. = FALSE)
  }
  fit0 <- qtl2mr_fit_from_result(fit)
  blocks <- fit0$locus_blocks
  if (!is.list(blocks) || !length(blocks)) stop("fit does not contain locus_blocks", call. = FALSE)
  full_theta <- fit0$theta
  if (length(blocks) == 1L) {
    b <- blocks[[1L]]
    out <- data.frame(
      locus_index = 1L,
      marker = b$marker %||% names(blocks)[[1L]],
      chr = b$chr %||% NA_character_,
      n_loci_minus_m = 0L,
      theta_full = full_theta,
      theta_minus_m = NA_real_,
      se_minus_m = NA_real_,
      profile_ci_lower = NA_real_,
      profile_ci_upper = NA_real_,
      Qmin = NA_real_,
      delta_theta = NA_real_,
      abs_delta_theta = NA_real_,
      influence_rank = NA_integer_,
      status = "single_locus_no_refit",
      notes = "LOO MR is undefined for a one-locus fit.",
      stringsAsFactors = FALSE
    )
    class(out) <- c("qtl2mr_loo_mr", "data.frame")
    return(out)
  }
  rows <- lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    marker <- b$marker %||% names(blocks)[[i]] %||% paste0("locus_", i)
    keep <- blocks[-i]
    nm <- names(blocks)
    if (!is.null(nm)) names(keep) <- nm[-i]
    fit_i <- tryCatch(
      fit_multilocus_founder_mr(
        keep,
        allow_cross_cov = allow_cross_cov,
        ridge = ridge,
        compute_independence = FALSE,
        heterogeneity_method = "chisq"
      ),
      error = function(e) e
    )
    if (inherits(fit_i, "error")) {
      return(data.frame(
        locus_index = i,
        marker = marker,
        chr = b$chr %||% NA_character_,
        n_loci_minus_m = length(keep),
        theta_full = full_theta,
        theta_minus_m = NA_real_,
        se_minus_m = NA_real_,
        profile_ci_lower = NA_real_,
        profile_ci_upper = NA_real_,
        Qmin = NA_real_,
        delta_theta = NA_real_,
        abs_delta_theta = NA_real_,
        influence_rank = NA_integer_,
        status = "failed",
        notes = conditionMessage(fit_i),
        stringsAsFactors = FALSE
      ))
    }
    delta <- fit_i$theta - full_theta
    data.frame(
      locus_index = i,
      marker = marker,
      chr = b$chr %||% NA_character_,
      n_loci_minus_m = length(keep),
      theta_full = full_theta,
      theta_minus_m = fit_i$theta,
      se_minus_m = fit_i$se,
      profile_ci_lower = fit_i$profile_ci_lower,
      profile_ci_upper = fit_i$profile_ci_upper,
      Qmin = fit_i$Q,
      delta_theta = delta,
      abs_delta_theta = abs(delta),
      influence_rank = NA_integer_,
      status = "ok",
      notes = "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  ok <- out$status == "ok" & is.finite(out$abs_delta_theta)
  out$influence_rank[ok] <- rank(-out$abs_delta_theta[ok], ties.method = "first")
  out <- out[order(out$influence_rank, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("qtl2mr_loo_mr", "data.frame")
  out
}

plot_loo_multilocus_mr <- function(loo,
                                   top_n = 25L,
                                   show = c("theta", "delta"),
                                   main = "Leave-one-locus-out MR sensitivity",
                                   xlab = NULL,
                                   ...) {
  show <- match.arg(show)
  if (!is.data.frame(loo)) stop("loo must be a data frame returned by loo_multilocus_founder_mr()", call. = FALSE)
  d <- loo[loo$status == "ok" & is.finite(loo$theta_minus_m), , drop = FALSE]
  if (!nrow(d)) stop("No successful LOO rows to plot", call. = FALSE)
  d <- d[order(d$abs_delta_theta, decreasing = TRUE), , drop = FALSE]
  if (!is.null(top_n) && is.finite(top_n)) d <- utils::head(d, as.integer(top_n))
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  y <- seq_len(nrow(d))
  if (identical(show, "theta")) {
    x <- d$theta_minus_m
    lo <- d$profile_ci_lower
    hi <- d$profile_ci_upper
    ref <- d$theta_full[[1L]]
    xlab <- xlab %||% expression(theta[-m])
  } else {
    x <- d$delta_theta
    lo <- d$profile_ci_lower - d$theta_full
    hi <- d$profile_ci_upper - d$theta_full
    ref <- 0
    xlab <- xlab %||% expression(theta[-m] - theta)
  }
  xlim <- range(c(lo, hi, x, ref), finite = TRUE)
  if (!all(is.finite(xlim)) || diff(xlim) <= .Machine$double.eps) xlim <- c(-1, 1)
  pad <- diff(xlim) * 0.1
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 9, 4, 2))
  graphics::plot(x, y, xlim = xlim + c(-pad, pad), ylim = c(0.5, length(y) + 0.5),
                 yaxt = "n", xlab = xlab, ylab = "", pch = 16,
                 col = "#00796b", main = main, ...)
  graphics::abline(v = ref, lty = 2, col = "grey35", lwd = 2)
  qtl2mr_plot_errorbarh(x, lo, hi, y, col = "#4c63d2")
  labels <- paste0(d$marker, if ("chr" %in% names(d)) paste0(" chr", d$chr) else "")
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.72)
  graphics::mtext("Rows ordered by absolute change in global theta", side = 3, line = 0.2, cex = 0.75, col = "grey35")
  invisible(d)
}

refit_qtl2mr_subset <- function(result,
                                subset = NULL,
                                reason = NULL,
                                allow_cross_cov = TRUE,
                                ridge = 1e-10,
                                set_name = "custom") {
  blocks <- qtl2mr_blocks_from_result(result)
  if (is.null(subset)) subset <- names(blocks)
  subset <- as.character(subset)
  if (!length(subset)) stop("subset must contain at least one marker", call. = FALSE)
  missing <- setdiff(subset, names(blocks))
  if (length(missing)) stop("subset markers not found in result: ", paste(missing, collapse = ", "), call. = FALSE)
  keep <- blocks[subset]
  if (!length(keep)) stop("subset leaves no locus blocks to fit", call. = FALSE)
  fit <- fit_multilocus_founder_mr(
    keep,
    allow_cross_cov = allow_cross_cov,
    ridge = ridge,
    compute_independence = TRUE,
    heterogeneity_method = "chisq"
  )
  data.frame(
    set_name = set_name,
    n_loci = length(keep),
    theta = fit$theta,
    se = fit$se,
    profile_ci_lower = fit$profile_ci_lower,
    profile_ci_upper = fit$profile_ci_upper,
    Qmin = fit$Q,
    overlap_delta = fit$theta_overlap_delta %||% NA_real_,
    status = "ok",
    reason = reason %||% "",
    notes = "",
    stringsAsFactors = FALSE
  )
}

qtl2mr_marker_metadata <- function(result) {
  if (inherits(result, "qtl2mr_analysis_result")) {
    selected <- result$pilot$selected_loci %||% data.frame()
    per <- result$pilot$diagnostic_layer$per_locus %||% data.frame()
  } else {
    fit <- qtl2mr_fit_from_result(result)
    selected <- data.frame(
      marker = names(fit$locus_blocks),
      chr = vapply(fit$locus_blocks, function(b) as.character(b$chr %||% NA_character_), character(1)),
      stringsAsFactors = FALSE
    )
    per <- fit$locus_diagnostics %||% data.frame()
  }
  if (!"marker" %in% names(selected)) selected$marker <- unique(per$marker)
  if ("F_X_block" %in% names(per)) {
    selected <- merge(selected, per[, c("marker", "F_X_block"), drop = FALSE], by = "marker", all.x = TRUE, sort = FALSE)
  }
  selected
}

qtl2mr_default_instrument_sets <- function(result,
                                          loo = NULL,
                                          weak_f_threshold = 2,
                                          top_n = c(10L, 25L, 50L),
                                          high_influence_n = 5L) {
  blocks <- qtl2mr_blocks_from_result(result)
  markers <- names(blocks)
  meta <- qtl2mr_marker_metadata(result)
  sets <- list(all_loci = markers)
  if ("F_X_block" %in% names(meta)) {
    ord <- meta$marker[order(meta$F_X_block, decreasing = TRUE, na.last = TRUE)]
    ord <- ord[ord %in% markers]
    for (n in top_n) {
      if (length(ord) >= n) sets[[paste0("top", n, "_by_block_F")]] <- utils::head(ord, n)
    }
    strong <- meta$marker[is.finite(meta$F_X_block) & meta$F_X_block >= weak_f_threshold]
    strong <- strong[strong %in% markers]
    if (length(strong)) sets[[paste0("exclude_F_lt_", weak_f_threshold)]] <- strong
  }
  if ("chr" %in% names(meta)) {
    for (chr in sort(unique(as.character(meta$chr)))) {
      keep <- meta$marker[as.character(meta$chr) != chr]
      keep <- keep[keep %in% markers]
      if (length(keep)) sets[[paste0("exclude_chr", qtl2mr_slug(chr))]] <- keep
    }
  }
  if (!is.null(loo) && is.data.frame(loo) && "marker" %in% names(loo)) {
    influential <- loo$marker[order(loo$abs_delta_theta, decreasing = TRUE, na.last = TRUE)]
    influential <- utils::head(influential[is.finite(loo$abs_delta_theta)], high_influence_n)
    keep <- setdiff(markers, influential)
    if (length(influential) && length(keep)) sets[[paste0("exclude_top", length(influential), "_loo_influence")]] <- keep
  }
  sets
}

fit_instrument_set_sensitivity <- function(result,
                                           sets = NULL,
                                           loo = NULL,
                                           weak_f_threshold = 2,
                                           top_n = c(10L, 25L, 50L),
                                           allow_cross_cov = TRUE,
                                           ridge = 1e-10) {
  if (is.null(sets)) {
    sets <- qtl2mr_default_instrument_sets(result, loo = loo, weak_f_threshold = weak_f_threshold, top_n = top_n)
  }
  if (!is.list(sets) || !length(sets)) stop("sets must be a non-empty named list of marker vectors", call. = FALSE)
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) names(sets) <- paste0("set_", seq_along(sets))
  rows <- lapply(seq_along(sets), function(i) {
    nm <- names(sets)[[i]]
    refit <- tryCatch(
      refit_qtl2mr_subset(
        result = result,
        subset = sets[[i]],
        reason = nm,
        allow_cross_cov = allow_cross_cov,
        ridge = ridge,
        set_name = nm
      ),
      error = function(e) e
    )
    if (inherits(refit, "error")) {
      return(data.frame(
        set_name = nm,
        n_loci = length(sets[[i]]),
        theta = NA_real_,
        se = NA_real_,
        profile_ci_lower = NA_real_,
        profile_ci_upper = NA_real_,
        Qmin = NA_real_,
        overlap_delta = NA_real_,
        status = "failed",
        reason = nm,
        notes = conditionMessage(refit),
        stringsAsFactors = FALSE
      ))
    }
    refit
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("qtl2mr_instrument_sensitivity", "data.frame")
  out
}

plot_instrument_set_sensitivity <- function(sensitivity,
                                            global_theta = NULL,
                                            main = "Instrument-set sensitivity",
                                            xlab = "MR effect",
                                            ...) {
  if (!is.data.frame(sensitivity)) stop("sensitivity must be a data frame", call. = FALSE)
  d <- sensitivity[sensitivity$status == "ok" & is.finite(sensitivity$theta), , drop = FALSE]
  if (!nrow(d)) stop("No successful instrument-set sensitivity rows to plot", call. = FALSE)
  d <- d[order(d$theta), , drop = FALSE]
  y <- seq_len(nrow(d))
  xlim <- range(c(d$profile_ci_lower, d$profile_ci_upper, d$theta, global_theta), finite = TRUE)
  if (!all(is.finite(xlim)) || diff(xlim) <= .Machine$double.eps) xlim <- c(-1, 1)
  pad <- diff(xlim) * 0.12
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(5, 10, 4, 2))
  labels <- paste0(d$set_name, " (M=", d$n_loci, ")")
  margins <- graphics::par("mai")
  margins[[2L]] <- max(margins[[2L]],
    max(graphics::strwidth(labels, units = "inches", cex = 0.76)) + 0.35)
  graphics::par(mai = margins)
  graphics::plot(d$theta, y, xlim = xlim + c(-pad, pad), ylim = c(0.5, length(y) + 0.5),
                 yaxt = "n", xlab = xlab, ylab = "", pch = 16,
                 col = "#00796b", main = main, ...)
  graphics::abline(v = 0, lty = 3, col = "grey70")
  if (!is.null(global_theta) && is.finite(global_theta)) graphics::abline(v = global_theta, lty = 2, col = "grey35", lwd = 2)
  qtl2mr_plot_errorbarh(d$theta, d$profile_ci_lower, d$profile_ci_upper, y, col = "#4c63d2")
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.76)
  invisible(d)
}

summarize_block_first_stage <- function(block_first_stage, threshold = 2) {
  if (!is.data.frame(block_first_stage) || !"F_X_block" %in% names(block_first_stage)) {
    stop("block_first_stage must contain F_X_block", call. = FALSE)
  }
  f <- block_first_stage$F_X_block
  data.frame(
    n_loci = sum(is.finite(f)),
    n_below_threshold = sum(is.finite(f) & f < threshold),
    threshold = threshold,
    min = min(f, na.rm = TRUE),
    q25 = as.numeric(stats::quantile(f, 0.25, na.rm = TRUE)),
    median = stats::median(f, na.rm = TRUE),
    q75 = as.numeric(stats::quantile(f, 0.75, na.rm = TRUE)),
    max = max(f, na.rm = TRUE),
    notes = "Block F_X is the joint Wald statistic of the 7 founder contrasts divided by 7; F_X < 2 is nominal p > 0.05 for the exposure.",
    stringsAsFactors = FALSE
  )
}

plot_block_first_stage_dashboard <- function(block_first_stage,
                                             threshold = 2,
                                             main = "Block first-stage denominator dashboard",
                                             ...) {
  if (!is.data.frame(block_first_stage) || !"F_X_block" %in% names(block_first_stage)) {
    stop("block_first_stage must contain F_X_block", call. = FALSE)
  }
  f <- block_first_stage$F_X_block
  f <- f[is.finite(f)]
  if (!length(f)) stop("No finite F_X_block values", call. = FALSE)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(5, 4, 4, 1))
  graphics::hist(f, breaks = "FD", col = "#80b1aa", border = "white",
                 xlab = "Block F_X", main = main, ...)
  graphics::abline(v = threshold, lty = 2, col = "#b2182b", lwd = 2)
  graphics::plot(sort(f, decreasing = TRUE), pch = 16, col = ifelse(sort(f, decreasing = TRUE) < threshold, "#c46a1b", "#00796b"),
                 xlab = "Rank", ylab = "Block F_X", main = "Ranked F_X")
  graphics::abline(h = threshold, lty = 2, col = "#b2182b", lwd = 2)
  invisible(summarize_block_first_stage(block_first_stage, threshold = threshold))
}
