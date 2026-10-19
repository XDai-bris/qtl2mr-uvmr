v060_fixture <- function(nonestimable_b2 = FALSE) {
  set.seed(60017)
  n <- 80L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2", "a3")), `2` = make_chr(c("b1", "b2")), X = make_chr("xm1"))
  if (nonestimable_b2) {
    gp[["2"]][, "A", "b2"] <- gp[["2"]][, "A", "b2"] + gp[["2"]][, "H", "b2"]
    gp[["2"]][, "H", "b2"] <- 0
  }
  G <- matrix(rnorm(n * 20L), n)
  K <- tcrossprod(G) / 20 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  x <- 3 * gp[["1"]][, 1L, "a1"] + 2 * gp[["2"]][, 2L, "b2"] + rnorm(n, sd = 0.3)
  w <- 2 * gp[["1"]][, 3L, "a3"] + rnorm(n, sd = 0.3)
  pheno <- data.frame(x = x, w = w, y = 0.6 * x + 0.3 * w + rnorm(n, sd = 0.3), row.names = ids)
  pheno$y[1:6] <- NA_real_
  pheno$x[7:8] <- NA_real_
  covar <- data.frame(sex = rep(0:1, n / 2L), row.names = ids)
  map <- data.frame(marker = c("a1", "a2", "a3", "b1", "b2", "xm1"), chr = c("1", "1", "1", "2", "2", "X"),
    pos = c(1, 2, 3, 1, 2, 1))
  scan <- matrix(c(8, 4.5, 1, 1, 6, 5,
                   1, 1, 7, 4.2, 1, 1), ncol = 2L, dimnames = list(map$marker, c("x", "w")))
  resources <- qtl2mr_prepare_resources(gp, list(`1` = K, `2` = K, X = K), pheno, covar, map)
  list(resources = resources, scan = scan,
    pairs = data.frame(exposure = c("x", "w", "x"), outcome = c("y", "y", "w")))
}

v060_single <- function(resources, instruments, exposure, outcome, ...) {
  input <- qtl2mr_pair_input(resources, instruments, exposure, outcome)
  args <- utils::modifyList(list(make_html = FALSE, make_figures = FALSE, verbose = FALSE, progress = FALSE,
    heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L, run_loo_mr = FALSE,
    run_instrument_sensitivity = FALSE), list(...))
  do.call(run_qtl2mr_analysis, c(list(input, output_dir = tempfile("v060-single-")), args))
}

test_that("instrument selection records every LOD candidate and applies its rules in order", {
  f <- v060_fixture()
  inst <- qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)
  expect_identical(inst$exposure, c("x", "x", "w", "w"))
  expect_identical(inst$marker, c("a1", "b2", "a3", "b1"))
  expect_identical(inst$lod, c(8, 6, 7, 4.2))
  expect_true(all(inst$selected) && all(inst$drop_reason == ""))
  expect_identical(inst$locus_rank, c(1L, 2L, 1L, 2L))
  expect_true(all(is.finite(inst$F_X_block) & inst$F_X_block > 0))
  expect_identical(inst$n_x, c(78L, 78L, 80L, 80L))
  expect_true(all(is.finite(inst$hsq_x)))
  expect_identical(unique(inst$exclude_chr), "X")
  expect_true(all(is.na(inst$F_X_min) & is.na(inst$max_loci)))

  threshold <- sort(inst$F_X_block)[[2L]]
  strict <- qtl2mr_select_instruments(f$resources, f$scan, F_X_min = threshold, progress = FALSE)
  expect_identical(strict$F_X_block, inst$F_X_block)
  expect_identical(strict$selected, inst$F_X_block >= threshold)
  expect_identical(strict$drop_reason, ifelse(inst$F_X_block >= threshold, "", "F"))
  expect_identical(unique(strict$F_X_min), threshold)

  top <- qtl2mr_select_instruments(f$resources, f$scan, max_loci = 1L, progress = FALSE)
  expect_identical(top$drop_reason, c("", "max_loci", "", "max_loci"))
  expect_identical(top$locus_rank, c(1L, NA, 1L, NA))

  with_x <- qtl2mr_select_instruments(f$resources, f$scan, exposures = "x", exclude_chr = NULL, progress = FALSE)
  expect_identical(with_x$marker, c("a1", "b2", "xm1"))
  expect_identical(with_x$drop_reason, c("", "", "non_autosomal"))
  expect_true(is.na(with_x$F_X_block[[3L]]))

  expect_identical(qtl2mr_select_instruments(f$resources, f$scan, n_workers = 2L, progress = FALSE), inst)
  expect_error(qtl2mr_select_instruments(f$resources, f$scan, exposures = "absent"), "missing exposure")
  expect_error(qtl2mr_select_instruments(f$resources, f$scan, F_X_min = -1), "non-negative")
  expect_error(qtl2mr_select_instruments(list(), f$scan), "qtl2mr_prepare_resources")
  empty <- qtl2mr_select_instruments(f$resources, f$scan, lod_thresh = 50, progress = FALSE)
  expect_identical(nrow(empty), 0L)
  expect_identical(names(empty), names(inst))
})

test_that("a locus whose exposure effects cannot be estimated is never selected", {
  f <- v060_fixture(nonestimable_b2 = TRUE)
  inst <- qtl2mr_select_instruments(f$resources, f$scan, exposures = "x", progress = FALSE)
  expect_identical(inst$marker, c("a1", "b2"))
  expect_identical(inst$selected, c(TRUE, FALSE))
  expect_identical(inst$drop_reason, c("", "nonestimable_x"))
  expect_identical(inst$locus_rank, c(1L, NA))
})

test_that("selection F equals the block F of every pair that uses the exposure", {
  f <- v060_fixture()
  inst <- qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)
  batch <- run_qtl2mr_batch(f$resources, f$pairs, inst, progress = FALSE)
  expect_true(all(batch$main_results$status == "ok"))
  joined <- merge(batch$locus_results, inst, by = c("exposure", "marker"), suffixes = c("", ".selected"))
  expect_identical(nrow(joined), 6L)
  expect_identical(joined$F_X_block, joined$F_X_block.selected)
  expect_identical(joined$hsq_x, joined$hsq_x.selected)
  expect_identical(joined$n_x, joined$n_x.selected)
})

test_that("each batch row equals the single analysis of the same pair", {
  f <- v060_fixture()
  inst <- qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)
  batch <- run_qtl2mr_batch(f$resources, f$pairs, inst, progress = FALSE)
  bookkeeping <- c("pair_id", "status", "error_message", "locus_fit_status", "runtime_sec")
  for (i in seq_len(nrow(f$pairs))) {
    single <- qtl2mr_primary_result(v060_single(f$resources, inst, f$pairs$exposure[[i]], f$pairs$outcome[[i]]))
    row <- batch$main_results[i, , drop = FALSE]
    shared <- setdiff(names(row), bookkeeping)
    expect_true(all(shared %in% names(single)))
    for (column in shared) expect_identical(row[[column]], single[[column]], info = paste(i, column))
  }
  # x is fitted on its 78 samples in both of its pairs, y on 74, w on all 80.
  expect_equal(batch$main_results$n_x_median, c(78, 80, 78))
  expect_equal(batch$main_results$n_y_median, c(74, 74, 80))

  simulated <- run_qtl2mr_batch(f$resources, f$pairs[1L, ], inst, heterogeneity_method = "simulation",
    heterogeneity_nsim = 25L, progress = FALSE)
  single <- qtl2mr_primary_result(v060_single(f$resources, inst, "x", "y",
    heterogeneity_method = "simulation", heterogeneity_nsim = 25L))
  expect_identical(simulated$main_results$Q_p_reported, single$Q_p_reported)
  expect_identical(simulated$main_results$Q_p_method, single$Q_p_method)
})

test_that("qtl2mr_pair_input uses the selected instruments of one exposure", {
  f <- v060_fixture()
  inst <- qtl2mr_select_instruments(f$resources, f$scan, max_loci = 1L, progress = FALSE)
  input <- qtl2mr_pair_input(f$resources, inst, "x", "y")
  expect_s3_class(input, "qtl2mr_prepared_input")
  expect_identical(input$selected_loci$marker, "a1")
  expect_identical(c(input$exposure, input$outcome), c("x", "y"))
  expect_null(input$exposure_samples)
  expect_error(qtl2mr_pair_input(f$resources, inst, "y", "x"), "no_instruments")
  expect_error(qtl2mr_pair_input(f$resources, inst, "x", NA_character_), "outcome must be one")
  expect_error(qtl2mr_pair_input(list(), inst, "x", "y"), "qtl2mr_prepare_resources")
})

test_that("the single analysis reports a Steiger check and batch records the flags", {
  f <- v060_fixture()
  r <- f$resources
  pheno <- r$pheno
  # b1 acts on the outcome directly and barely on the exposure.
  pheno$y_pleio <- 0.6 * pheno$x + 3 * r$allele_probs[["2"]][, 4L, "b1"] + rnorm(nrow(pheno), sd = 0.3)
  r <- qtl2mr_prepare_resources(r$allele_probs, r$kinship, pheno, r$covar, r$marker_map)
  inst <- data.frame(exposure = "x", marker = c("a1", "b1"))
  result <- v060_single(r, inst, "x", "y_pleio", run_instrument_sensitivity = TRUE)
  first <- utils::read.csv(file.path(result$output_dir, "data", "block_first_stage.csv"), stringsAsFactors = FALSE)
  expect_identical(first$marker, c("a1", "b1"))
  expect_identical(first$steiger_pass, c(TRUE, FALSE))
  expect_identical(qtl2mr_primary_result(result)$n_loci_steiger_fail, 1L)

  steiger <- result$steiger_sensitivity
  expect_identical(steiger$n_loci_removed, 1L)
  expect_identical(steiger$removed_markers, "b1")
  expect_identical(steiger$status, "ok")
  refit <- refit_qtl2mr_subset(result, "a1")
  expect_identical(steiger$theta_steiger, refit$theta)
  primary <- qtl2mr_primary_result(result)
  expect_identical(steiger$shift_in_primary_ci_widths,
    (refit$theta - primary$theta) / (primary$profile_ci_upper - primary$profile_ci_lower))
  expect_identical(utils::tail(result$instrument_sensitivity$set_name, 1L), "steiger_filtered")
  expect_true(file.exists(file.path(result$output_dir, "data", "steiger_sensitivity.csv")))
  expect_identical("steiger_estimate_outside_ci" %in% qtl2mr_warnings(result)$warning_id,
    isTRUE(steiger$outside_primary_ci))
  expect_match(paste(readLines(file.path(result$output_dir, "summary.md")), collapse = "\n"), "Actionable", fixed = TRUE)

  batch <- run_qtl2mr_batch(r, data.frame(exposure = "x", outcome = "y_pleio"), inst, progress = FALSE)
  expect_identical(batch$main_results$n_loci_steiger_fail, 1L)
  expect_identical(batch$locus_results$steiger_pass, c(TRUE, FALSE))

  clean <- v060_single(f$resources, data.frame(exposure = "x", marker = "a1"), "x", "y",
    run_instrument_sensitivity = TRUE)
  expect_identical(clean$steiger_sensitivity$status, "no_loci_removed")
  expect_identical(clean$steiger_sensitivity$theta_steiger, qtl2mr_primary_result(clean)$theta)
  expect_identical(clean$steiger_sensitivity$shift_in_primary_ci_widths, 0)
})
