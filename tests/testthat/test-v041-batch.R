v041_batch_fixture <- function() {
  set.seed(41017)
  n <- 40L
  ids <- paste0("s", seq_len(n))
  G <- matrix(rexp(n * 8L), nrow = n)
  G <- G / rowSums(G)
  gp <- list(`1` = array(G, c(n, 8L, 1L), list(ids, LETTERS[1:8], "m1")))
  x <- 5 * G[, 1L] + rnorm(n, sd = 0.15)
  w <- 4 * G[, 2L] + rnorm(n, sd = 0.15)
  pheno <- data.frame(x = x, w = w, y = 0.7 * x + 0.2 * w + rnorm(n, sd = 0.2),
    row.names = ids)
  covar <- data.frame(z = rnorm(n), row.names = ids)
  map <- data.frame(marker = "m1", chr = "1", pos = 1)
  resources <- qtl2mr_prepare_resources(gp, NULL, pheno, covar,
    matrix(c(6, 6), nrow = 1L, dimnames = list("m1", c("x", "w"))), map)
  list(resources = resources,
    pairs = data.frame(exposure = c("x", "w"), outcome = "y"),
    options = list(min_samples = 20L, cache_contexts = FALSE, cache_lmm_fits = FALSE))
}

test_that("single and exact batch sample policies differ until masks are matched", {
  f <- v041_batch_fixture()
  r <- f$resources
  r$pheno$y[1:4] <- NA_real_
  input <- as_qtl2mr_prepared_input(r$allele_probs, r$pheno, r$covar,
    marker_map = r$marker_map, selected_loci = "m1", exposure = "x", outcome = "y")
  fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    heterogeneity_method = "none", diagnostic_bootstrap_B = 0L)
  batch <- run_qtl2mr_batch(r, f$pairs[1L, ], fit_options = f$options, progress = FALSE)
  expect_equal(fit$overlap_summary$n_x, 40L)
  expect_equal(fit$overlap_summary$n_y, 36L)
  expect_equal(batch$main_results$n_x_median, 36L)
  expect_equal(batch$main_results$n_y_median, 36L)
  expect_identical(batch$main_results$sample_policy, "pair_common_complete")
  expect_identical(batch$config$sample_policy, "pair_common_complete")
  expect_identical(batch$run_manifest$sample_policy, "pair_common_complete")
  common <- qtl2mr_make_common_complete_samples(r, "x", "y", "1")
  input$exposure_samples <- input$outcome_samples <- common
  matched <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    heterogeneity_method = "none", diagnostic_bootstrap_B = 0L)
  expect_equal(matched$summary$theta, batch$main_results$theta, tolerance = 1e-10)
  expect_equal(matched$summary$se, batch$main_results$se, tolerance = 1e-10)
})

test_that("all batch routes expose covariate exclusions without changing masks", {
  f <- v041_batch_fixture()
  f$resources$pheno$y[1:3] <- NA_real_
  f$resources$covar$z[2:5] <- c(NA_real_, NaN, Inf, -Inf)
  f$resources$covar <- f$resources$covar[-6L, , drop = FALSE]
  settings <- list(list(sample_policy = "pair_common_complete", n_workers = 1L),
    list(sample_policy = "batch_common_complete", n_workers = 1L),
    list(sample_policy = "pair_common_complete", n_workers = 2L),
    list(sample_policy = "batch_common_complete", n_workers = 2L))
  for (setting in settings) {
    path <- tempfile("v041-exclusions-")
    args <- c(list(resources = f$resources, pairs = f$pairs, fit_options = f$options,
      output_dir = path, progress = FALSE), setting)
    batch <- suppressWarnings(do.call(run_qtl2mr_batch, args))
    expect_true(all(batch$main_results$status == "ok"))
    expect_true(all(batch$main_results$n_x_median == 34L))
    expect_true(all(batch$main_results$n_y_median == 34L))
    ledger <- batch$sample_exclusions
    expect_true(all(ledger$n_input_samples == 40L))
    incomplete <- ledger[ledger$reason == "incomplete_covar", ]
    expect_equal(incomplete$n_excluded, c(4L, 4L))
    expect_true(all(incomplete$excluded_ids == "s2;s3;s4;s5"))
    expect_true(all(ledger$n_excluded[ledger$reason == "missing_covar_id"] == 1L))
    expect_true(all(ledger$n_excluded[ledger$reason == "nonfinite_outcome"] == 3L))
    expect_true(all(ledger$counts_overlap))
    warnings <- batch$warnings[batch$warnings$warning_id == "batch_covariate_sample_exclusions", ]
    expect_equal(nrow(warnings), 2L)
    expect_true(all(grepl("5 input samples", warnings$message)))
    expect_true(file.exists(file.path(path, "batch_sample_exclusions.csv")))
    saved <- read.csv(file.path(path, "batch_sample_exclusions.csv"), stringsAsFactors = FALSE)
    expect_equal(saved$n_excluded, ledger$n_excluded)
  }
})

test_that("rounded resource probabilities use audited local copies in cached batches", {
  f <- v041_batch_fixture()
  r <- f$resources
  rounded <- lapply(r$allele_probs, round, digits = 4L)
  expect_error(qtl2mr_prepare_resources(rounded, NULL, r$pheno, r$covar,
    r$qtl_scan, r$marker_map, probability_policy = "strict"), "sum to 1")
  expect_warning(resources <- qtl2mr_prepare_resources(rounded, NULL, r$pheno, r$covar,
    r$qtl_scan, r$marker_map, metadata = list(custom = "preserved")), "rounding")
  original <- resources$allele_probs
  path <- tempfile("v041-rounded-batch-")
  expect_no_warning(batch <- run_qtl2mr_batch(resources, f$pairs, sample_policy = "batch_common_complete",
    fit_options = list(min_samples = 20L), output_dir = path, progress = FALSE))
  expect_true(all(batch$main_results$status == "ok"))
  expect_identical(resources$allele_probs, original)
  expect_identical(resources$metadata$custom, "preserved")
  expect_identical(batch$config$probability_policy, "normalize_rounding")
  expect_true(nrow(batch$probability_normalization) > 0L)
  expect_true(all(batch$probability_normalization$n_normalized > 0L))
  expect_true(file.exists(file.path(path, "batch_probability_normalization.csv")))
  warnings <- batch$warnings[batch$warnings$warning_id == "probability_rounding_normalized", ]
  expect_equal(nrow(warnings), 2L)
  expect_true(all(grepl("batch_probability_normalization.csv", warnings$message, fixed = TRUE)))
  normalized <- rounded
  for (chr in names(normalized)) for (marker in dimnames(normalized[[chr]])[[3L]]) {
    G <- normalized[[chr]][, , marker]
    adjust <- abs(rowSums(G) - 1) > 1e-6 + 16 * .Machine$double.eps
    G[adjust, ] <- G[adjust, , drop = FALSE] / rowSums(G[adjust, , drop = FALSE])
    normalized[[chr]][, , marker] <- G
  }
  explicit <- qtl2mr_prepare_resources(normalized, NULL, r$pheno, r$covar,
    r$qtl_scan, r$marker_map, probability_policy = "strict")
  comparison <- run_qtl2mr_batch(explicit, f$pairs, sample_policy = "batch_common_complete",
    fit_options = list(min_samples = 20L), progress = FALSE)
  expect_equal(batch$main_results$theta, comparison$main_results$theta, tolerance = 1e-10)
  expect_equal(batch$main_results$se, comparison$main_results$se, tolerance = 1e-10)
})
