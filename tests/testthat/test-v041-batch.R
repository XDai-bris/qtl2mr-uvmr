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

v041_batch_hashes <- function(path) {
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE))
  tools::md5sum(files[!dir.exists(files)])
}

test_that("parallel resume rejects changes before touching parent or worker files", {
  f <- v041_batch_fixture()
  for (fast in c(FALSE, TRUE)) for (fail_fast in c(FALSE, TRUE)) {
    path <- tempfile("v041-parallel-integrity-")
    args <- list(resources = f$resources, pairs = f$pairs, fast = fast,
      n_workers = 1L, fit_options = f$options, output_dir = path,
      progress = FALSE, fail_fast = fail_fast)
    first <- do.call(run_qtl2mr_batch_parallel_by_exposure, args)
    expect_true(all(first$main_results$status == "ok"))
    expect_false(is.null(first$resume_fingerprint))
    args$resume <- TRUE
    again <- do.call(run_qtl2mr_batch_parallel_by_exposure, args)
    expect_identical(again$main_results$pair_id, first$main_results$pair_id)
    expect_equal(again$main_results$theta, first$main_results$theta, tolerance = 1e-10)
    expect_equal(nrow(again$main_results), 2L)
    expect_equal(anyDuplicated(again$main_results$pair_id), 0L)

    before <- v041_batch_hashes(path)
    changed <- args
    changed$fit_options$standardize_phenotypes <- FALSE
    expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, changed), "fingerprint changed")
    expect_identical(v041_batch_hashes(path), before)
    changed <- args
    changed$resources$pheno$x[1L] <- changed$resources$pheno$x[1L] + 0.01
    expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, changed), "fingerprint changed")
    expect_identical(v041_batch_hashes(path), before)
    changed$resources$pheno$y[] <- NA_real_
    expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, changed), "fingerprint changed")
    expect_identical(v041_batch_hashes(path), before)
    changed <- args
    changed$sample_retention_min <- 1.01
    expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, changed), "fingerprint changed")
    expect_identical(v041_batch_hashes(path), before)

    targets <- file.path(path, c("batch_main_results.csv", "batch_config.json",
      "by_exposure/x/batch_main_results.csv", "by_exposure/x/batch_resume_state.rds"))
    for (target in targets) {
      original <- readBin(target, "raw", n = file.info(target)$size)
      writeBin(charToRaw("damaged synthetic checkpoint"), target)
      damaged <- v041_batch_hashes(path)
      expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, args), "Cannot resume")
      expect_identical(v041_batch_hashes(path), damaged)
      writeBin(original, target)
    }
  }
})

test_that("all worker fingerprints are checked before the first worker can write", {
  f <- v041_batch_fixture()
  path <- tempfile("v041-worker-preflight-")
  args <- list(resources = f$resources, pairs = f$pairs, fast = FALSE,
    n_workers = 1L, fit_options = f$options, output_dir = path, progress = FALSE)
  do.call(run_qtl2mr_batch_parallel_by_exposure, args)
  child <- file.path(path, "by_exposure", "x", "batch_resume_state.rds")
  state <- readRDS(child)
  state$fingerprint$request <- "different request"
  saveRDS(state, child)
  parent <- file.path(path, "batch_resume_state.rds")
  state <- readRDS(parent)
  state$artifacts[normalizePath(child)] <- unname(tools::md5sum(child))
  saveRDS(state, parent)
  before <- v041_batch_hashes(path)
  args$resume <- TRUE
  expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, args), "fingerprint changed")
  expect_identical(v041_batch_hashes(path), before)
})

test_that("fast-wrapper preflight and retention refusals preserve prior outputs", {
  f <- v041_batch_fixture()
  path <- tempfile("v041-fast-resume-")
  args <- list(resources = f$resources, pairs = f$pairs,
    fit_options = f$options, output_dir = path, progress = FALSE)
  first <- do.call(run_qtl2mr_batch_fast, args)
  args$resume <- TRUE
  again <- do.call(run_qtl2mr_batch_fast, args)
  expect_equal(again$main_results$theta, first$main_results$theta, tolerance = 1e-10)
  before <- v041_batch_hashes(path)
  changed <- args
  changed$resources$pheno$y[] <- NA_real_
  expect_error(do.call(run_qtl2mr_batch_fast, changed), "fingerprint changed")
  expect_identical(v041_batch_hashes(path), before)
  changed <- args
  changed$sample_retention_min <- 1.01
  expect_error(do.call(run_qtl2mr_batch_fast, changed), "Cannot resume")
  expect_identical(v041_batch_hashes(path), before)
})

test_that("Unix parallel resume preserves checkpoints on changed options", {
  skip_if_not(identical(.Platform$OS.type, "unix"))
  f <- v041_batch_fixture()
  path <- tempfile("v041-unix-resume-")
  args <- list(resources = f$resources, pairs = f$pairs, fast = TRUE,
    n_workers = 2L, fit_options = f$options, output_dir = path, progress = FALSE)
  first <- do.call(run_qtl2mr_batch_parallel_by_exposure, args)
  expect_identical(first$run_manifest$parallel_backend, "parallel::mclapply")
  args$resume <- TRUE
  again <- do.call(run_qtl2mr_batch_parallel_by_exposure, args)
  expect_equal(again$main_results$theta, first$main_results$theta, tolerance = 1e-10)
  before <- v041_batch_hashes(path)
  args$fit_options$cross_cov_model <- "zero"
  expect_error(do.call(run_qtl2mr_batch_parallel_by_exposure, args), "fingerprint changed")
  expect_identical(v041_batch_hashes(path), before)
})

test_that("group failures produce one current status per affected pair", {
  f <- v041_batch_fixture()
  f$resources$pheno$z <- rev(f$resources$pheno$y)
  f$pairs <- rbind(f$pairs, data.frame(exposure = "w", outcome = "z"))
  original_runner <- run_qtl2mr_batch
  local_mocked_bindings(run_qtl2mr_batch = function(resources, pairs, ...) {
    if (identical(pairs$exposure[[1L]], "w")) stop("synthetic worker failure")
    original_runner(resources, pairs, ...)
  }, .package = "qtl2mr")
  workers <- if (identical(.Platform$OS.type, "unix")) c(1L, 2L) else 1L
  for (n_workers in workers) {
    batch <- run_qtl2mr_batch_parallel_by_exposure(f$resources, f$pairs,
      fast = FALSE, n_workers = n_workers, fit_options = f$options, progress = FALSE)
    expect_identical(batch$main_results$status, c("ok", "error", "error"))
    expect_equal(batch$run_manifest$n_pairs_requested, 3L)
    expect_equal(batch$run_manifest$n_pairs_ok, 1L)
    expect_equal(batch$run_manifest$n_pairs_error, 2L)
    expect_equal(nrow(batch$errors), 2L)
    expect_setequal(batch$errors$pair_id, c("w__to__y", "w__to__z"))
    expect_true(all(batch$errors$stage == "parallel_exposure_group"))
    expect_error(run_qtl2mr_batch_parallel_by_exposure(f$resources, f$pairs,
      fast = FALSE, n_workers = n_workers, fit_options = f$options,
      progress = FALSE, fail_fast = TRUE), "synthetic worker failure")
  }
})

test_that("resume integrity failures never become ordinary group failures", {
  group <- list(key = "x", pairs = data.frame(exposure = "x", outcome = "y"))
  error <- tryCatch(qtl2mr:::qtl2mr_resume_abort("Cannot resume: synthetic mismatch"), error = identity)
  failure <- qtl2mr:::qtl2mr_parallel_group_failure(group, error)
  expect_error(qtl2mr:::qtl2mr_collect_parallel_group_results(list(failure), list(group),
    fail_fast = FALSE), "Cannot resume")
})

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
  runners <- list(run_qtl2mr_batch, run_qtl2mr_batch_exact_grouped,
    run_qtl2mr_batch_fast, run_qtl2mr_batch_parallel_by_exposure)
  for (runner in runners) {
    path <- tempfile("v041-exclusions-")
    args <- list(resources = f$resources, pairs = f$pairs, fit_options = f$options,
      output_dir = path, progress = FALSE)
    if (identical(runner, run_qtl2mr_batch_parallel_by_exposure)) args$n_workers <- 1L
    batch <- suppressWarnings(do.call(runner, args))
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
  expect_no_warning(batch <- run_qtl2mr_batch_fast(resources, f$pairs,
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
  comparison <- run_qtl2mr_batch_fast(explicit, f$pairs,
    fit_options = list(min_samples = 20L), progress = FALSE)
  expect_equal(batch$main_results$theta, comparison$main_results$theta, tolerance = 1e-10)
  expect_equal(batch$main_results$se, comparison$main_results$se, tolerance = 1e-10)
})
