peaks <- function(lod, ...) qtl2mr:::qtl2mr_get_tophits(lod, lod_thresh = 4, ...)

test_that("the peakdrop rule merges plateaus and splits only at a real drop", {
  # A flat plateau around the threshold is one peak, at its maximum.
  expect_identical(peaks(c(1, 4.02, 3.97, 4.04, 3.96, 4.05, 4.01, 3.95, 1)), 6L)
  # Two maxima separated by a drop of 1.6 are two peaks; a drop of 1.4 is one.
  expect_identical(peaks(c(1, 6, 4.4, 6, 1)), c(2L, 4L))
  expect_identical(peaks(c(1, 6, 4.6, 6, 1)), 2L)
  expect_identical(peaks(c(1, 6, 4.6, 6, 1), peakdrop = 1), c(2L, 4L))
  # The second peak must also rise by more than peakdrop above the dip.
  expect_identical(peaks(c(1, 9, 3, 4.4, 1)), 2L)
  expect_identical(peaks(c(1, 9, 3, 4.6, 1)), c(2L, 4L))
  # A higher maximum after a shallow dip moves the peak.
  expect_identical(peaks(c(1, 5, 4.2, 6, 1)), 4L)
  expect_identical(peaks(c(1, 6, 1, 7, 1), peakdrop = Inf), 4L)
})

test_that("peak rules handle edges, missing values, ties and the threshold", {
  for (rule in c("peakdrop", "lod_island")) {
    # Peak regions touching the first or last marker.
    expect_identical(peaks(c(5, 4.5, 1), rule = rule), 1L, info = rule)
    expect_identical(peaks(c(1, 4.5, 5), rule = rule), 3L, info = rule)
    expect_identical(peaks(c(5, 6, 5), rule = rule), 2L, info = rule)
    expect_identical(peaks(5, rule = rule), 1L, info = rule)
    # Missing and non-finite values neither start a peak nor split one.
    expect_identical(peaks(c(1, 5, NA, 5.5, 1), rule = rule), 4L, info = rule)
    expect_identical(peaks(c(1, 6, NaN, 6, 1), rule = rule), 2L, info = rule)
    expect_identical(peaks(c(1, 5, Inf, 4.5, -Inf), rule = rule), 2L, info = rule)
    expect_identical(peaks(c(NA, NA), rule = rule), integer(), info = rule)
    # Nothing above the threshold.
    expect_identical(peaks(c(1, 3.9, 2), rule = rule), integer(), info = rule)
    expect_identical(peaks(numeric(), rule = rule), integer(), info = rule)
    # A tie at the maximum keeps the first marker in map order.
    expect_identical(peaks(c(1, 5, 5, 1), rule = rule), 2L, info = rule)
    # LOD exactly at the threshold counts as above it.
    expect_identical(peaks(c(1, 4, 1), rule = rule), 2L, info = rule)
  }
  # The island rule (the candidates of "spacing") splits a plateau at every dip below the threshold.
  expect_identical(peaks(c(1, 4.02, 3.97, 4.04, 1), rule = "lod_island"), c(4L, 2L))
  for (bad in list(0, -1, NA_real_, "1.5", TRUE, c(1, 2))) {
    expect_error(peaks(c(1, 5), peakdrop = bad), "peakdrop must be")
  }
})

test_that("scan candidates follow the chosen rule within each chromosome", {
  map <- data.frame(marker = c(paste0("m", 1:7), "n1", "n2"), chr = c(rep("1", 7), "2", "2"),
    pos = c(10, 10.2, 10.4, 10.6, 11.5, 20, 30, 5, 6), stringsAsFactors = FALSE)
  # Scan rows deliberately out of map order.
  ord <- c(9, 3, 1, 7, 2, 8, 4, 6, 5)
  scan <- matrix(c(4.5, 3.9, 4.6, 3.8, 4.4, 1, 6, 1, 4.8)[ord], ncol = 1L,
    dimnames = list(map$marker[ord], "x"))
  cand <- function(rule, ...) qtl2mr:::qtl2mr_scan_candidates(scan, map, "x", 4, "X", peak_rule = rule, ...)
  expect_identical(cand("peakdrop")$marker, c("m7", "n2", "m3"))
  # m1 is within 1 Mb of m3; m5 is 1.1 Mb away.
  expect_identical(cand("spacing")$marker, c("m7", "n2", "m3", "m5"))
  expect_identical(cand("spacing", min_spacing = 2)$marker, c("m7", "n2", "m3"))
  expect_identical(cand("peakdrop", peakdrop = 0.5)$marker, c("m7", "n2", "m3", "m1", "m5"))
  expect_identical(cand("peakdrop")$lod, c(6, 4.8, 4.6))
  expect_error(qtl2mr:::qtl2mr_scan_candidates(scan, map[, c("marker", "chr")], "x", 4, "X",
    peak_rule = "spacing"), "needs finite marker positions")
  expect_error(cand("nearest"), "should be one of")
})

test_that("the instrument table records the peak rule", {
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  ex <- qtl2mr_example
  scan <- matrix(c(8, 4.5, 3.9, 4.2), ncol = 1, dimnames = list(ex$marker_map$marker, ex$exposure))
  resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno, ex$covar,
    marker_map = ex$marker_map)
  select <- function(...) qtl2mr_select_instruments(resources, scan, progress = FALSE, ...)
  primary <- select()
  # m4 follows a dip below the threshold but only 0.3 above it, so it joins the m1 peak.
  expect_identical(primary$marker, "m1")
  expect_identical(unique(primary$peak_rule), "peakdrop")
  expect_identical(unique(primary$selected_by), "lod_peakdrop_max")
  expect_identical(unique(primary$peakdrop), 1.5)
  expect_true(all(is.na(primary$min_spacing)))
  expect_identical(unique(primary$rule_version), "tophit_v1")

  spaced <- select(peak_rule = "spacing", min_spacing = 2)
  expect_identical(spaced$marker, c("m1", "m4"))
  expect_identical(unique(spaced$selected_by), "lod_island_max_min_spacing")
  expect_identical(unique(spaced$min_spacing), 2)
  expect_identical(select(peak_rule = "spacing", min_spacing = 3)$marker, "m1")

  expect_identical(names(primary), names(spaced))
  # The legacy rule of 0.6.1 and earlier was removed in 0.6.1.9006.
  expect_error(select(peak_rule = "legacy"), "should be one of")
  expect_identical(names(select(lod_thresh = 50)), names(primary))
  expect_error(select(peakdrop = 0), "peakdrop must be")
  expect_error(select(min_spacing = -1), "min_spacing must be")
  expect_error(select(peak_rule = "nearest"), "should be one of")

  bp <- ex$marker_map
  bp$pos <- bp$pos * 1e6
  bp_resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno, ex$covar, marker_map = bp)
  expect_warning(qtl2mr_select_instruments(bp_resources, scan, peak_rule = "spacing", progress = FALSE),
    "looks like base pairs")
  expect_silent(qtl2mr_select_instruments(bp_resources, scan, progress = FALSE))
})
