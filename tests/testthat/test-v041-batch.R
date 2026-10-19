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
  resources <- qtl2mr_prepare_resources(gp, NULL, pheno, covar, map)
  list(resources = resources,
    instruments = data.frame(exposure = c("x", "w"), marker = "m1"),
    pairs = data.frame(exposure = c("x", "w"), outcome = "y"))
}

test_that("single and batch analyses fit the same role-specific samples", {
  f <- v041_batch_fixture()
  r <- f$resources
  r$pheno$y[1:4] <- NA_real_
  input <- as_qtl2mr_prepared_input(r$allele_probs, r$pheno, r$covar,
    marker_map = r$marker_map, selected_loci = "m1", exposure = "x", outcome = "y")
  fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L)
  batch <- run_qtl2mr_batch(r, f$pairs[1L, ], f$instruments, min_samples = 20L, progress = FALSE)
  expect_equal(fit$overlap_summary$n_x, 40L)
  expect_equal(fit$overlap_summary$n_y, 36L)
  expect_equal(batch$main_results$n_x_median, 40L)
  expect_equal(batch$main_results$n_y_median, 36L)
  expect_identical(batch$main_results$theta, fit$summary$theta)
  expect_identical(batch$main_results$se, fit$summary$se)
  expect_identical(batch$main_results$profile_ci_lower, fit$summary$profile_ci_lower)
})

test_that("batch runs report covariate exclusions per trait for every worker count", {
  f <- v041_batch_fixture()
  f$resources$pheno$y[1:3] <- NA_real_
  f$resources$covar$z[2:5] <- c(NA_real_, NaN, Inf, -Inf)
  f$resources$covar <- f$resources$covar[-6L, , drop = FALSE]
  for (n_workers in c(1L, 2L)) {
    path <- tempfile("v041-exclusions-")
    batch <- suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L,
      n_workers = n_workers, output_dir = path, progress = FALSE))
    expect_true(all(batch$main_results$status == "ok"))
    expect_true(all(batch$main_results$n_x_median == 35L))
    expect_true(all(batch$main_results$n_y_median == 34L))
    ledger <- batch$sample_exclusions
    expect_setequal(unique(ledger$trait), c("x", "w", "y"))
    expect_true(all(ledger$n_input_samples == 40L))
    incomplete <- ledger[ledger$reason == "incomplete_covar", ]
    expect_equal(incomplete$n_excluded, c(4L, 4L, 4L))
    expect_true(all(incomplete$excluded_ids == "s2;s3;s4;s5"))
    expect_true(all(ledger$n_excluded[ledger$reason == "missing_covar_id"] == 1L))
    expect_identical(ledger$n_excluded[ledger$reason == "nonfinite_value"], c(0L, 0L, 3L))
    warnings <- batch$warnings[batch$warnings$warning_id == "incomplete_covariates", ]
    expect_equal(nrow(warnings), 2L)
    expect_true(all(grepl("batch_sample_exclusions.csv", warnings$message, fixed = TRUE)))
    saved <- read.csv(file.path(path, "batch_sample_exclusions.csv"), stringsAsFactors = FALSE)
    expect_equal(saved$n_excluded, ledger$n_excluded)
  }
})

test_that("rounded resource probabilities use audited local copies in cached batches", {
  f <- v041_batch_fixture()
  r <- f$resources
  rounded <- lapply(r$allele_probs, round, digits = 4L)
  expect_error(qtl2mr_prepare_resources(rounded, NULL, r$pheno, r$covar,
    r$marker_map, probability_policy = "strict"), "sum to 1")
  expect_warning(resources <- qtl2mr_prepare_resources(rounded, NULL, r$pheno, r$covar,
    r$marker_map, metadata = list(custom = "preserved")), "rounding")
  original <- resources$allele_probs
  path <- tempfile("v041-rounded-batch-")
  expect_no_warning(batch <- run_qtl2mr_batch(resources, f$pairs, f$instruments, min_samples = 20L,
    output_dir = path, progress = FALSE))
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
    r$marker_map, probability_policy = "strict")
  comparison <- run_qtl2mr_batch(explicit, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  expect_equal(batch$main_results$theta, comparison$main_results$theta, tolerance = 1e-10)
  expect_equal(batch$main_results$se, comparison$main_results$se, tolerance = 1e-10)
})
