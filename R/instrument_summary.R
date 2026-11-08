# F report of an instrument table (0.6.1.9004): one row per exposure with the
# candidate counts by drop reason, the block F distribution of the selected
# loci, and how many fitted candidates reach each F break.

qtl2mr_instrument_summary_settings <- c("selected_by", "lod_thresh", "peak_rule", "peakdrop", "min_spacing",
  "rule_version", "F_X_min", "max_loci", "exclude_chr", "min_founder_copies")

qtl2mr_instrument_summary <- function(instruments, F_breaks = c(1.5, 2, 5, 10)) {
  if (!is.data.frame(instruments)) {
    stop("instruments must be the data frame returned by qtl2mr_select_instruments()", call. = FALSE)
  }
  missing <- setdiff(c("exposure", "selected", "drop_reason", "F_X_block"), names(instruments))
  if (length(missing)) {
    stop("instruments is missing column(s) ", paste(missing, collapse = ", "),
      "; use the table returned by qtl2mr_select_instruments()", call. = FALSE)
  }
  if (!is.numeric(F_breaks) || !length(F_breaks) || any(!is.finite(F_breaks)) ||
      anyDuplicated(as.character(F_breaks))) {
    stop("F_breaks must be distinct finite numbers", call. = FALSE)
  }
  selected <- instruments$selected
  if (!is.logical(selected) || anyNA(selected)) stop("instruments$selected must be TRUE or FALSE", call. = FALSE)
  reason <- as.character(instruments$drop_reason)
  reason[is.na(reason)] <- ""
  F_X <- as.numeric(instruments$F_X_block)
  exposure <- as.character(instruments$exposure)
  founder_qc <- any(c("min_founder_copies", "E_min_x") %in% names(instruments)) || any(reason == "founder_qc")
  reasons <- c("non_autosomal", if (founder_qc) "founder_qc", "nonestimable_x", "x_fit_error", "F", "max_loci")
  reasons <- c(reasons, sort(setdiff(unique(reason[nzchar(reason)]), reasons)))
  settings <- intersect(qtl2mr_instrument_summary_settings, names(instruments))
  one <- function(e) {
    rows <- which(exposure == e)
    out <- data.frame(exposure = e, n_candidates = length(rows), stringsAsFactors = FALSE)
    for (r in reasons) out[[paste0("n_dropped_", r)]] <- sum(reason[rows] == r)
    out$n_selected <- sum(selected[rows])
    # Selected loci with a finite block F, summarised as summarize_block_first_stage() does.
    f <- F_X[rows][selected[rows] & is.finite(F_X[rows])]
    s <- if (length(f)) summarize_block_first_stage(data.frame(F_X_block = f)) else NULL
    out$F_min <- if (length(f)) s$min else NA_real_
    out$F_q25 <- if (length(f)) s$q25 else NA_real_
    out$F_median <- if (length(f)) s$median else NA_real_
    out$F_mean <- if (length(f)) mean(f) else NA_real_
    out$F_q75 <- if (length(f)) s$q75 else NA_real_
    out$F_max <- if (length(f)) s$max else NA_real_
    # Every fitted candidate, selected or not.
    fitted <- F_X[rows][is.finite(F_X[rows])]
    for (b in F_breaks) out[[paste0("n_F_ge_", as.character(b))]] <- sum(fitted >= b)
    for (nm in settings) out[[nm]] <- instruments[[nm]][if (length(rows)) rows[[1L]] else NA_integer_]
    out
  }
  exposures <- unique(exposure)
  if (!length(exposures)) {
    template <- one(NA_character_)
    return(template[0L, , drop = FALSE])
  }
  out <- do.call(rbind, lapply(exposures, one))
  rownames(out) <- NULL
  out
}
