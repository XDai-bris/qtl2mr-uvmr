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
  resources <- qtl2mr_prepare_resources(input$allele_probs, NULL, pheno, NULL, marker_map = map)
  list(input = input, resources = resources,
       instruments = data.frame(exposure = "x", marker = c("m1", "m3")),
       pairs = data.frame(exposure = "x", outcome = "y"),
       options = list(min_samples = 20L, diagnostic_bootstrap_B = 0L, heterogeneity_method = "none"))
}

v040_batch <- function(f, ..., resources = f$resources, pairs = f$pairs, instruments = f$instruments) {
  do.call(run_qtl2mr_batch, c(list(resources, pairs, instruments, progress = FALSE),
    utils::modifyList(f$options, list(...))))
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
  expect_length(qtl2mr:::qtl2mr_complete_samples(f$resources, c("x", "y"), "1"), 64L)
})

test_that("LOD islands stop at chromosome boundaries and native maps work", {
  lod <- matrix(c(1, 2, 5, 6, 7, 5, 1, 1), ncol = 1,
    dimnames = list(c(paste0("a", 1:4), paste0("b", 1:4)), "x"))
  map <- list(`1` = setNames(1:4, paste0("a", 1:4)), `2` = setNames(1:4, paste0("b", 1:4)))
  table <- qtl2mr:::qtl2mr_marker_map_table
  selected <- qtl2mr:::qtl2mr_scan_candidates(lod, table(map), "x", lod_thresh = 4, exclude_chr = "X")
  expect_setequal(selected$marker, c("a4", "b1"))
  expect_error(qtl2mr:::qtl2mr_scan_candidates(lod, table(list(`1` = map[[1]])), "x", 4, "X"), "missing qtl_scan markers")
})

test_that("misspelled options and invalid instrument tables are not ignored", {
  f <- v040_fixture()
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, cross_cov_modle = "zero"), "unused argument")
  expect_error(qtl2mr:::run_prepared_qtl2mr(f$input, "x", "y", standadize = FALSE), "Unused arguments")
  expect_error(v040_batch(f, instruments = file.path(tempdir(), "absent-v040-instruments.csv")),
    "instruments file does not exist")
  expect_error(v040_batch(f, instruments = data.frame(trait = "x", marker = "m1")), "exposure")
  expect_error(v040_batch(f, instruments = data.frame(exposure = "x", marker = c("m1", "m1"))), "more than once")
  expect_error(v040_batch(f, instruments = data.frame(exposure = "x", marker = "m9")), "not in marker_map")
  expect_error(v040_batch(f, instruments = data.frame(exposure = "x", marker = "m1", chr = "2")), "chr does not match")
  expect_error(v040_batch(f, instruments = data.frame(exposure = "x", marker = "m1", selected = NA)), "TRUE or FALSE")
  file <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(exposure = "x", marker = c("m1", "m2", "m3"), selected = c(TRUE, FALSE, TRUE)),
    file, row.names = FALSE)
  from_file <- v040_batch(f, instruments = file)
  expect_identical(from_file$main_results$theta, v040_batch(f)$main_results$theta)
  expect_identical(from_file$locus_results$marker, c("m1", "m3"))
})

test_that("resume requires exact provenance and leaves one current row per pair", {
  f <- v040_fixture()
  path <- tempfile("v040-batch-")
  b1 <- v040_batch(f, output_dir = path)
  expect_identical(b1$main_results$status, "ok")
  expect_error(v040_batch(f, output_dir = path), "output already exists")
  b2 <- v040_batch(f, output_dir = path, resume = TRUE)
  expect_equal(nrow(b2$main_results), 1L)
  expect_equal(b1$main_results$theta, b2$main_results$theta)
  changed <- f$resources
  changed$pheno$x[1] <- changed$pheno$x[1] + 0.01
  expect_error(v040_batch(f, resources = changed, output_dir = path, resume = TRUE), "fingerprint changed")
  expect_error(v040_batch(f, cross_cov_model = "zero", output_dir = path, resume = TRUE), "fingerprint changed")
  expect_false(formals(run_qtl2mr_batch)$resume)
  legacy <- tempfile("v040-legacy-")
  dir.create(legacy)
  file.copy(file.path(path, "batch_main_results.csv"), file.path(legacy, "batch_main_results.csv"))
  expect_error(v040_batch(f, output_dir = legacy, resume = TRUE), "legacy or incomplete")
  writeLines("modified", file.path(path, "batch_main_results.csv"))
  expect_error(v040_batch(f, output_dir = path, resume = TRUE), "tables changed")
})

test_that("incomplete occupied destinations fail closed without changing files", {
  f <- v040_fixture()
  for (owned in c("batch_locus_results.csv", "batch_selected_loci.csv", "batch_resume_state.rds", "mr_block_cache", "batch_pair_records")) {
    path <- tempfile("v040-incomplete-")
    dir.create(path)
    target <- file.path(path, owned)
    if (owned %in% c("mr_block_cache", "batch_pair_records")) {
      dir.create(target)
      target <- file.path(target, "existing.rds")
    }
    writeLines("preserve this existing content", target)
    hash <- unname(tools::md5sum(target))
    expect_error(v040_batch(f, output_dir = path), "output already exists")
    expect_error(v040_batch(f, output_dir = path, resume = TRUE), "legacy or incomplete|unreadable checkpoint")
    expect_identical(unname(tools::md5sum(target)), hash)
    expect_false(file.exists(file.path(path, "batch_main_results.csv")))
  }
  unrelated <- tempfile("v040-unrelated-")
  dir.create(unrelated)
  note <- file.path(unrelated, "notes.txt")
  writeLines("user notes", note)
  expect_true(qtl2mr:::qtl2mr_check_batch_destination(unrelated))
  result <- v040_batch(f, output_dir = unrelated)
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

test_that("each pair uses the instruments of its exposure", {
  f <- v040_fixture()
  f$resources$pheno$x3 <- f$resources$pheno$x
  instruments <- data.frame(exposure = c("x", "x", "x3"), marker = c("m1", "m3", "m3"))
  pairs <- data.frame(pair_id = c("both", "one"), exposure = c("x", "x3"), outcome = "y")
  result <- v040_batch(f, pairs = pairs, instruments = instruments)
  expect_equal(result$main_results$status, c("ok", "ok"))
  expect_equal(result$main_results$n_loci, c(2L, 1L))
  expect_equal(result$locus_results$marker, c("m1", "m3", "m3"))
})

test_that("complete-sample counts exclude nonfinite phenotypes", {
  f <- v040_fixture()
  f$resources$pheno$x[1] <- Inf
  expect_length(qtl2mr:::qtl2mr_complete_samples(f$resources, c("x", "y"), "1"), 63L)
  expect_length(qtl2mr:::qtl2mr_complete_samples(f$resources, "y", "1"), 64L)
  result <- v040_batch(f)
  expect_equal(result$main_results$n_x_median, 63)
  expect_equal(result$main_results$n_y_median, 64)
})

test_that("resuming a mixed batch does not retain duplicate error rows", {
  f <- v040_fixture()
  pairs <- rbind(f$pairs, data.frame(exposure = "x", outcome = "missing"))
  path <- tempfile("v040-mixed-")
  first <- v040_batch(f, pairs = pairs, output_dir = path)
  resumed <- v040_batch(f, pairs = pairs, output_dir = path, resume = TRUE)
  expect_equal(first$main_results$status, c("ok", "error"))
  expect_equal(resumed$main_results$status, c("ok", "error"))
  expect_equal(nrow(resumed$errors), 1L)
  expect_equal(anyDuplicated(resumed$main_results$pair_id), 0L)
})

test_that("batch runs preserve valid pairs when another request cannot be planned", {
  f <- v040_fixture()
  pairs <- rbind(f$pairs, data.frame(exposure = "x", outcome = "missing"))
  for (n_workers in c(1L, 2L)) {
    result <- v040_batch(f, pairs = pairs, n_workers = n_workers)
    expect_equal(result$main_results$status, c("ok", "error"))
    expect_equal(length(unique(result$main_results$pair_id)), 2L)
    expect_error(v040_batch(f, pairs = pairs, n_workers = n_workers, fail_fast = TRUE), "Outcome trait not found")
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
