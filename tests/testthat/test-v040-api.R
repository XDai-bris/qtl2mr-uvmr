v040_fixture <- function() {
  set.seed(6401)
  n <- 64L
  ids <- paste0("s", seq_len(n))
  markers <- paste0("m", 1:3)
  probs <- array(rexp(n * 8L * 3L), c(n, 8L, 3L), list(ids, LETTERS[1:8], markers))
  for (j in 1:3) probs[, , j] <- probs[, , j] / rowSums(probs[, , j])
  x <- 3 * probs[, 1, 1] + 2 * probs[, 2, 3] + rnorm(n, sd = 0.15)
  pheno <- data.frame(x = x, y = 0.7 * x + rnorm(n, sd = 0.2), row.names = ids)
  map <- data.frame(marker = markers, chr = "1", pos = 1:3)
  input <- as_qtl2mr_prepared_input(list(`1` = probs), pheno, marker_map = map,
    selected_loci = c("m1", "m3"), exposure = "x", outcome = "y")
  resources <- qtl2mr_prepare_resources(input$allele_probs, NULL, pheno, NULL,
    qtl_scan = matrix(c(6, 0, 5), ncol = 1, dimnames = list(markers, "x")), marker_map = map)
  list(input = input, resources = resources,
       pairs = data.frame(exposure = "x", outcome = "y"),
       options = list(min_samples = 20L, diagnostic_bootstrap_B = 0L, heterogeneity_method = "none"))
}

test_that("prepared selection is explicit and only autosomal", {
  fixture <- v040_fixture()
  x <- fixture$input
  x$selected_loci <- NULL
  expect_error(qtl2mr:::qtl2mr_prepared_selected_loci(x), "selected_loci is required")
  expect_identical(qtl2mr:::qtl2mr_prepared_selected_loci(x, loci = "m1")$marker, "m1")
  expect_error(qtl2mr:::qtl2mr_prepared_selected_loci(x, loci = c("m1", "m1")), "unique")
  expect_error(qtl2mr:::qtl2mr_prepared_selected_loci(x, loci = -1), "positive integer")
  x$use_all_loci <- TRUE
  expect_equal(nrow(qtl2mr:::qtl2mr_prepared_selected_loci(x)), 3L)
  x <- fixture$input
  names(x$allele_probs) <- "X"
  x$marker_map$chr <- "X"
  expect_error(validate_qtl2mr_input(x), "autosomal")
})

test_that("input probabilities, IDs, and phenotype types fail closed", {
  f <- v040_fixture()
  x <- f$input
  x$pheno$x <- factor(x$pheno$x)
  expect_error(validate_qtl2mr_input(x), "must be numeric")
  x <- f$input
  dimnames(x$allele_probs[[1]])[[1]][2] <- "s1"
  expect_error(validate_qtl2mr_input(x), "unique")
  x <- f$input
  x$allele_probs[[1]] <- x$allele_probs[[1]] * 0.5
  expect_error(validate_qtl2mr_input(x), "sum to 1")
  x <- f$input
  rownames(x$pheno) <- paste0("other", seq_len(nrow(x$pheno)))
  expect_error(validate_qtl2mr_input(x), "no overlapping")
  x <- f$input
  expect_warning(partial <- as_qtl2mr_prepared_input(x$allele_probs, x$pheno[-1, ],
    marker_map = x$marker_map, selected_loci = x$selected_loci, exposure = "x", outcome = "y"), "Sample ID mismatch")
  expect_equal(partial$sample_alignment$n_missing, 1L)
  expect_null(f$resources$covar)
  expect_length(qtl2mr:::qtl2mr_make_common_complete_samples(f$resources, "x", "y", "1"), 64L)
})

test_that("LOD islands stop at chromosome boundaries and native maps work", {
  lod <- matrix(c(1, 2, 5, 6, 7, 5, 1, 1), ncol = 1,
    dimnames = list(c(paste0("a", 1:4), paste0("b", 1:4)), "x"))
  map <- list(`1` = setNames(1:4, paste0("a", 1:4)), `2` = setNames(1:4, paste0("b", 1:4)))
  selected <- qtl2mr_select_loci_from_scan(lod, map, "x", lod_thresh = 4)
  expect_setequal(selected$marker, c("a4", "b1"))
  expect_error(qtl2mr_select_loci_from_scan(lod, list(`1` = map[[1]]), "x"), "missing qtl_scan markers")
})

test_that("missing selection files and misspelled options are not ignored", {
  f <- v040_fixture()
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = list(cross_cov_modle = "zero")), "Unknown option")
  expect_error(suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, sample_rule = list(mdoe = "batch_common_complete"))), "Unknown sample_rule")
  expect_equal(length(qtl2mr:::qtl2mr_make_common_complete_samples(f$resources, "x", "y", "1",
    sample_rule = list(mode = "common_complete", sample_ids = paste0("s", 1:30)))), 30L)
  expect_error(qtl2mr:::run_prepared_qtl2mr(f$input, "x", "y", standadize = FALSE), "Unused arguments")
  f$pairs$selected_loci_file <- file.path(tempdir(), "absent-v040-selected.csv")
  b <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE)
  expect_identical(b$main_results$status, "error")
  expect_match(b$main_results$error_message, "selected_loci_file does not exist")
})

test_that("resume requires exact provenance and leaves one current row per pair", {
  f <- v040_fixture()
  path <- tempfile("v040-batch-")
  b1 <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = path, progress = FALSE)
  expect_identical(b1$main_results$status, "ok")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = path, progress = FALSE), "output already exists")
  b2 <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = path, resume = TRUE, progress = FALSE)
  expect_equal(nrow(b2$main_results), 1L)
  expect_equal(b1$main_results$theta, b2$main_results$theta)
  changed <- f$resources
  changed$pheno$x[1] <- changed$pheno$x[1] + 0.01
  expect_error(run_qtl2mr_batch(changed, f$pairs, fit_options = f$options, output_dir = path, resume = TRUE, progress = FALSE), "fingerprint changed")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = c(f$options, list(cross_cov_model = "zero")), output_dir = path, resume = TRUE, progress = FALSE), "fingerprint changed")
  expect_false(formals(run_qtl2mr_batch)$resume)
  expect_false(formals(run_qtl2mr_batch_fast)$resume)
  legacy <- tempfile("v040-legacy-")
  dir.create(legacy)
  file.copy(file.path(path, "batch_main_results.csv"), file.path(legacy, "batch_main_results.csv"))
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = legacy, resume = TRUE, progress = FALSE), "legacy or incomplete")
  writeLines("modified", file.path(path, "batch_main_results.csv"))
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = path, resume = TRUE, progress = FALSE), "tables changed")
})

test_that("incomplete occupied destinations fail closed without changing files", {
  f <- v040_fixture()
  for (owned in c("batch_selected_loci.csv", "batch_resume_state.rds", "mr_block_cache", "batch_pair_records")) {
    path <- tempfile("v040-incomplete-")
    dir.create(path)
    target <- file.path(path, owned)
    if (owned %in% c("mr_block_cache", "batch_pair_records")) {
      dir.create(target)
      target <- file.path(target, "existing.rds")
    }
    writeLines("preserve this existing content", target)
    hash <- unname(tools::md5sum(target))
    expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options,
      output_dir = path, progress = FALSE), "output already exists")
    expect_error(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options,
      output_dir = path, resume = TRUE, progress = FALSE), "legacy or incomplete|unreadable checkpoint")
    expect_identical(unname(tools::md5sum(target)), hash)
    expect_false(file.exists(file.path(path, "batch_main_results.csv")))
  }
  unrelated <- tempfile("v040-unrelated-")
  dir.create(unrelated)
  note <- file.path(unrelated, "notes.txt")
  writeLines("user notes", note)
  expect_true(qtl2mr:::qtl2mr_check_batch_destination(unrelated))
  result <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = unrelated, progress = FALSE)
  expect_identical(result$main_results$status, "ok")
  expect_identical(readLines(note), "user notes")
})

test_that("cache reuse includes actual phenotype and marker probability values", {
  f <- v040_fixture()
  contexts <- qtl2mr:::qtl2mr_context_cache_new()
  fits <- qtl2mr:::qtl2mr_context_cache_new()
  run <- function(input, cached = TRUE) qtl2mr:::run_prepared_qtl2mr(input, "x", "y",
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none", verbose = FALSE,
    context_cache = if (cached) contexts else NULL, fit_cache = if (cached) fits else NULL)
  original <- run(f$input)
  changed <- f$input
  changed$pheno$y <- rev(changed$pheno$y)
  expect_equal(run(changed)$fit$theta, run(changed, FALSE)$fit$theta, tolerance = 1e-10)
  changed <- f$input
  changed$allele_probs[[1]][, , 1] <- changed$allele_probs[[1]][rev(seq_len(64)), , 1]
  expect_equal(run(changed)$fit$theta, run(changed, FALSE)$fit$theta, tolerance = 1e-10)
})

test_that("mixed file and scan selection requests preserve both pairs", {
  f <- v040_fixture()
  file <- tempfile(fileext = ".csv")
  write.csv(data.frame(marker = "m3", chr = "1"), file, row.names = FALSE)
  pairs <- data.frame(pair_id = c("scan", "file"), exposure = "x", outcome = "y",
                      selected_loci_file = c(NA, file))
  for (policy in c("pair_common_complete", "batch_common_complete")) {
    result <- run_qtl2mr_batch(f$resources, pairs, sample_policy = policy, fit_options = f$options, progress = FALSE)
    expect_equal(result$main_results$status, c("ok", "ok"))
  }
})

test_that("fast common-sample counts exclude nonfinite phenotypes", {
  f <- v040_fixture()
  f$resources$pheno$x[1] <- Inf
  ids <- qtl2mr:::qtl2mr_make_batch_common_complete_samples(f$resources, c("x", "y"), "1")
  expect_length(ids, 63L)
  result <- run_qtl2mr_batch(f$resources, f$pairs, sample_policy = "batch_common_complete",
    fit_options = f$options, progress = FALSE)
  expect_equal(result$sample_qc$n_batch_common, 63L)
  expect_equal(result$main_results$n_x_median, result$main_results$n_y_median)
})

test_that("resuming a mixed batch does not retain duplicate error rows", {
  f <- v040_fixture()
  pairs <- rbind(f$pairs, data.frame(exposure = "x", outcome = "missing"))
  path <- tempfile("v040-mixed-")
  first <- run_qtl2mr_batch(f$resources, pairs, fit_options = f$options, output_dir = path, progress = FALSE)
  resumed <- run_qtl2mr_batch(f$resources, pairs, fit_options = f$options, output_dir = path, resume = TRUE, progress = FALSE)
  expect_equal(first$main_results$status, c("ok", "error"))
  expect_equal(resumed$main_results$status, c("ok", "error"))
  expect_equal(nrow(resumed$errors), 1L)
  expect_equal(anyDuplicated(resumed$main_results$pair_id), 0L)
})

test_that("batch runs preserve valid pairs when another request cannot be planned", {
  f <- v040_fixture()
  pairs <- rbind(f$pairs, data.frame(exposure = "x", outcome = "missing"))
  for (policy in c("pair_common_complete", "batch_common_complete")) for (n_workers in c(1L, 2L)) {
    args <- list(resources = f$resources, pairs = pairs, sample_policy = policy, n_workers = n_workers,
                 fit_options = f$options, progress = FALSE, fail_fast = FALSE)
    result <- do.call(run_qtl2mr_batch, args)
    expect_equal(result$main_results$status, c("ok", "error"))
    expect_equal(length(unique(result$main_results$pair_id)), 2L)
    args$fail_fast <- TRUE
    expect_error(do.call(run_qtl2mr_batch, args), "Outcome trait not found")
  }
})

test_that("scaling metadata uses each locus's fitted samples", {
  f <- v040_fixture()
  fit <- qtl2mr:::run_prepared_qtl2mr(f$input, "x", "y", verbose = FALSE,
    heterogeneity_method = "none", diagnostic_bootstrap_B = 0L)
  expect_identical(fit$summary$effect_units, "SD_y per SD_x")
  expect_equal(fit$summary$exposure_sd, sd(f$input$pheno$x))
  expect_equal(fit$summary$theta_raw_units, fit$summary$theta * sd(f$input$pheno$y) / sd(f$input$pheno$x))
  expect_equal(nrow(fit$phenotype_scaling), 4L)
  expect_true(all(fit$sample_alignment$n_dropped == 0L))
})
