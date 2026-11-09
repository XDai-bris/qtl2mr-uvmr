# 0.6.1.9004: qtl2mr_instrument_summary(), the F report of an instrument table.

v064_instrument_table <- function() {
  data.frame(
    exposure = c(rep("x", 7L), rep("w", 3L)),
    marker = paste0("m", 1:10),
    chr = c("X", "1", "1", "2", "2", "3", "4", "X", "1", "2"),
    pos = 1:10, lod = c(9, 8, 7, 6, 5, 4.5, 4.2, 6, 5, 4),
    F_X_block = c(NA, NA, 1.2, 1.5, 2, 4, 12, NA, NA, 0.8),
    n_x = c(NA, NA, rep(100L, 5L), NA, NA, 90L), hsq_x = c(NA, NA, rep(0.3, 5L), NA, NA, 0.2),
    selected = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    drop_reason = c("non_autosomal", "nonestimable_x", "F", "", "", "", "", "non_autosomal", "x_fit_error", "F"),
    locus_rank = c(NA, NA, NA, 1L, 2L, 3L, 4L, NA, NA, NA),
    selected_by = "lod_peakdrop_max", lod_thresh = 4, peak_rule = "peakdrop", peakdrop = 1.5, min_spacing = NA_real_,
    rule_version = "tophit_v1", F_X_min = 1.5, max_loci = NA_integer_, exclude_chr = "",
    stringsAsFactors = FALSE
  )
}

test_that("hand-computed summaries, counts and settings", {
  inst <- v064_instrument_table()
  # A locus dropped by max_loci is fitted: it counts for the F breaks but not the selected-locus summaries.
  inst <- rbind(inst, transform(inst[7L, ], marker = "m11", F_X_block = 3, selected = FALSE, drop_reason = "max_loci",
    locus_rank = NA_integer_))
  s <- qtl2mr_instrument_summary(inst)
  expect_identical(s$exposure, c("x", "w"))
  expect_identical(names(s), c("exposure", "n_candidates", "n_dropped_non_autosomal", "n_dropped_nonestimable_x",
    "n_dropped_x_fit_error", "n_dropped_F", "n_dropped_max_loci", "n_selected", "F_min", "F_q25", "F_median", "F_mean",
    "F_q75", "F_max", "n_F_ge_1.5", "n_F_ge_2", "n_F_ge_5", "n_F_ge_10", "selected_by", "lod_thresh", "peak_rule",
    "peakdrop", "min_spacing", "rule_version", "F_X_min", "max_loci", "exclude_chr"))
  x <- s[1L, ]
  expect_identical(c(x$n_candidates, x$n_dropped_non_autosomal, x$n_dropped_nonestimable_x, x$n_dropped_x_fit_error,
    x$n_dropped_F, x$n_dropped_max_loci, x$n_selected), c(8L, 1L, 1L, 0L, 1L, 1L, 4L))
  # Selected F: 1.5, 2, 4, 12. Type-7 quantiles: q25 = 1.5 + 0.75 * 0.5, q75 = 4 + 0.25 * 8.
  expect_identical(c(x$F_min, x$F_q25, x$F_median, x$F_mean, x$F_q75, x$F_max), c(1.5, 1.875, 3, 4.875, 6, 12))
  # Fitted F: 1.2, 1.5, 2, 4, 12 and 3 (max_loci).
  expect_identical(c(x[["n_F_ge_1.5"]], x$n_F_ge_2, x$n_F_ge_5, x$n_F_ge_10), c(5L, 4L, 1L, 1L))
  expect_identical(x$peak_rule, "peakdrop")
  expect_identical(x$F_X_min, 1.5)
  # No selected locus: NA summaries, zero counts.
  w <- s[2L, ]
  expect_identical(c(w$n_candidates, w$n_dropped_non_autosomal, w$n_dropped_x_fit_error, w$n_dropped_F, w$n_selected),
    c(3L, 1L, 1L, 1L, 0L))
  expect_true(all(is.na(c(w$F_min, w$F_q25, w$F_median, w$F_mean, w$F_q75, w$F_max))))
  expect_identical(c(w[["n_F_ge_1.5"]], w$n_F_ge_2, w$n_F_ge_5, w$n_F_ge_10), rep(0L, 4L))
  # The same quantile definition as the single analysis's first-stage summary.
  ref <- qtl2mr:::summarize_block_first_stage(data.frame(F_X_block = c(1.5, 2, 4, 12)))
  expect_identical(c(x$F_min, x$F_q25, x$F_median, x$F_q75, x$F_max), c(ref$min, ref$q25, ref$median, ref$q75, ref$max))
  # Breaks name the columns; other drop reasons get their own count after the known ones.
  custom <- qtl2mr_instrument_summary(transform(inst, drop_reason = ifelse(marker == "m3", "ld", drop_reason)),
    F_breaks = c(3, 0.5))
  expect_identical(grep("^n_F_ge_|^n_dropped_", names(custom), value = TRUE),
    c("n_dropped_non_autosomal", "n_dropped_nonestimable_x", "n_dropped_x_fit_error", "n_dropped_F",
      "n_dropped_max_loci", "n_dropped_ld", "n_F_ge_3", "n_F_ge_0.5"))
  expect_identical(custom$n_dropped_ld, c(1L, 0L))
  expect_identical(custom$n_F_ge_3, c(3L, 0L))
  expect_identical(custom[["n_F_ge_0.5"]], c(6L, 1L))
  # Empty table: no rows, same columns.
  empty <- qtl2mr_instrument_summary(inst[0L, ])
  expect_identical(nrow(empty), 0L)
  expect_identical(names(empty), names(s))
  for (bad in list(numeric(), c(2, 2), NA_real_, Inf, "2")) {
    expect_error(qtl2mr_instrument_summary(inst, F_breaks = bad), "F_breaks")
  }
  expect_error(qtl2mr_instrument_summary(inst[, setdiff(names(inst), "drop_reason")]), "drop_reason")
  expect_error(qtl2mr_instrument_summary("instruments.csv"), "data frame")
  expect_error(qtl2mr_instrument_summary(transform(inst, selected = NA)), "selected")
})

test_that("on a selection: counts add up, founder QC appears when used, and n_F_ge counts the fitted F", {
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  ex <- qtl2mr_example
  scan <- matrix(ex$marker_map$exposure_lod, ncol = 1L, dimnames = list(ex$marker_map$marker, ex$exposure))
  resources <- suppressWarnings(qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno, ex$covar,
    marker_map = ex$marker_map))
  # The founder QC is on by default since 0.6.1.9006; NULL turns it off.
  all_loci <- suppressWarnings(qtl2mr_select_instruments(resources, scan, lod_thresh = 3, progress = FALSE,
    min_founder_copies = NULL))
  s <- qtl2mr_instrument_summary(all_loci)
  expect_identical(nrow(s), 1L)
  expect_identical(s$n_candidates, nrow(all_loci))
  expect_identical(s$n_selected, sum(all_loci$selected))
  expect_identical(s$n_candidates, s$n_selected + sum(unlist(s[grep("^n_dropped_", names(s))])))
  expect_false("n_dropped_founder_qc" %in% names(s))
  expect_false("min_founder_copies" %in% names(s))
  f <- all_loci$F_X_block[all_loci$selected]
  expect_identical(s$F_median, stats::median(f))
  expect_identical(s$F_mean, mean(f))
  # Selection no longer filters on F (F_X_min was removed in 0.6.1.9006): no
  # locus is dropped for F, and n_F_ge_t counts the selected loci an F >= t
  # subset would keep.
  expect_identical(s$n_dropped_F, 0L)
  expect_true(all(is.na(all_loci$F_X_min) & is.na(all_loci$max_loci)))
  for (threshold in c(1.5, 2, 5, 10)) {
    expect_identical(s[[paste0("n_F_ge_", threshold)]], sum(all_loci$selected & all_loci$F_X_block >= threshold))
  }
  qc <- suppressWarnings(qtl2mr_select_instruments(resources, scan, lod_thresh = 3, progress = FALSE,
    min_founder_copies = 0))
  sq <- qtl2mr_instrument_summary(qc)
  expect_identical(sq$n_dropped_founder_qc, 0L)
  expect_identical(match("n_dropped_founder_qc", names(sq)), match("n_dropped_non_autosomal", names(sq)) + 1L)
  expect_identical(sq$min_founder_copies, 0)
  expect_identical(sq[setdiff(names(sq), c("n_dropped_founder_qc", "min_founder_copies"))], s)
})
