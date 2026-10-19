v050_batch_fixture <- function() {
  set.seed(50017)
  n <- 48L
  ids <- paste0("s", seq_len(n))
  G1 <- matrix(rexp(n * 8L), nrow = n)
  G1 <- G1 / rowSums(G1)
  G2 <- matrix(rexp(n * 8L), nrow = n)
  G2 <- G2 / rowSums(G2)
  gp <- list(`1` = array(G1, c(n, 8L, 1L), list(ids, LETTERS[1:8], "m1")),
    `2` = array(G2, c(n, 8L, 1L), list(ids, LETTERS[1:8], "m2")))
  x <- 5 * G1[, 1L] + rnorm(n, sd = 0.15)
  w <- 4 * G2[, 2L] + rnorm(n, sd = 0.15)
  pheno <- data.frame(x = x, w = w, y = 0.7 * x + 0.2 * w + rnorm(n, sd = 0.2),
    v = 0.3 * w + rnorm(n, sd = 0.2), row.names = ids)
  pheno$v[1:3] <- NA_real_
  covar <- data.frame(z = rnorm(n), row.names = ids)
  map <- data.frame(marker = c("m1", "m2"), chr = c("1", "2"), pos = c(1, 1))
  scan <- matrix(c(6, 1, 1, 6), nrow = 2L, dimnames = list(c("m1", "m2"), c("x", "w")))
  resources <- qtl2mr_prepare_resources(gp, NULL, pheno, covar, map)
  instruments <- qtl2mr_select_instruments(resources, scan, progress = FALSE)
  list(resources = resources, instruments = instruments,
    pairs = data.frame(exposure = c("x", "x", "w", "w"), outcome = c("y", "v", "y", "v")),
    options = list(min_samples = 20L))
}

v050_hashes <- function(path) {
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE))
  hashes <- tools::md5sum(files[!dir.exists(files)])
  names(hashes) <- sub(paste0("^", path, "/"), "", names(hashes))
  hashes
}

v050_numbers <- function(batch) {
  main <- batch$main_results[order(batch$main_results$pair_id), , drop = FALSE]
  main[, c("pair_id", "status", "theta", "se", "profile_ci_lower", "profile_ci_upper", "Q",
    "n_x_median", "n_y_median"), drop = FALSE]
}

v050_run <- function(f, path, ...) {
  do.call(run_qtl2mr_batch, c(list(f$resources, f$pairs, f$instruments, output_dir = path, progress = FALSE),
    utils::modifyList(f$options, list(...))))
}

test_that("one batch entry point remains and retired arguments are rejected", {
  f <- v050_batch_fixture()
  ns <- getNamespaceExports("qtl2mr")
  for (old in c("run_qtl2mr_batch_fast", "run_qtl2mr_batch_exact_grouped", "run_qtl2mr_batch_parallel_by_exposure",
                "qtl2mr_select_loci_from_scan", "qtl2mr_select_loci_for_traits")) {
    expect_false(old %in% ns)
  }
  expect_true(all(c("qtl2mr_select_instruments", "qtl2mr_pair_input") %in% ns))
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, sample_policy = "batch_common_complete"), "unused argument")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, fit_options = list()), "unused argument")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, force = TRUE), "unused argument")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs), "instruments is required")
  for (column in c("instrument_trait", "lod_thresh", "exclude_chr", "max_loci", "selected_loci_file")) {
    old_pairs <- f$pairs
    old_pairs[[column]] <- "x"
    expect_error(run_qtl2mr_batch(f$resources, old_pairs, f$instruments), "removed in qtl2mr 0.6.0")
  }
  scan <- matrix(1, 2, 1, dimnames = list(c("m1", "m2"), "x"))
  expect_error(qtl2mr_prepare_resources(f$resources$allele_probs, NULL, f$resources$pheno, f$resources$covar, scan),
    "qtl_scan is no longer part of resources")
  batch <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  expect_true(all(batch$main_results$status == "ok"))
  expect_true(all(batch$main_results$Q_p_method == "chisq"))
  expect_true(all(batch$main_results$confidence_set_bounded))
  expect_identical(utils::tail(names(batch$main_results), 2L), c("locus_fit_status", "runtime_sec"))
})

test_that("worker count does not change results", {
  f <- v050_batch_fixture()
  serial <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  parallel <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, n_workers = 2L, progress = FALSE)
  expect_identical(v050_numbers(parallel), v050_numbers(serial))
  if (identical(.Platform$OS.type, "unix")) expect_identical(parallel$run_manifest$parallel_backend, "parallel::mclapply")
})

test_that("resume works in place, after moving the directory and with another worker count", {
  f <- v050_batch_fixture()
  root <- tempfile("v050-resume-")
  path <- file.path(root, "a")
  first <- v050_run(f, path)
  expect_true(dir.exists(file.path(path, "batch_pair_records")))
  expect_length(list.files(file.path(path, "batch_pair_records")), 4L)
  again <- v050_run(f, path, resume = TRUE)
  expect_identical(again$run_manifest$n_pairs_resumed, 4L)
  expect_identical(v050_numbers(again), v050_numbers(first))

  moved <- file.path(root, "moved")
  expect_true(file.rename(path, moved))
  after_move <- v050_run(f, moved, resume = TRUE)
  expect_identical(after_move$run_manifest$n_pairs_resumed, 4L)
  expect_identical(v050_numbers(after_move), v050_numbers(first))

  more_workers <- v050_run(f, moved, resume = TRUE, n_workers = 2L)
  expect_identical(v050_numbers(more_workers), v050_numbers(first))
})

test_that("interrupted runs resume from completed pair records", {
  f <- v050_batch_fixture()
  reference <- v050_numbers(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE))
  interrupt <- function(path, keep) {
    # Emulate a stop after `keep` pairs: tables were never published and the
    # checkpoint carries no table hashes, exactly as during a run.
    state <- readRDS(file.path(path, "batch_resume_state.rds"))
    unlink(c(list.files(path, pattern = "^batch_.*\\.csv$", full.names = TRUE),
      file.path(path, "batch_config.json")))
    records <- sort(list.files(file.path(path, "batch_pair_records"), full.names = TRUE))
    unlink(records[seq_along(records) > keep])
    qtl2mr:::qtl2mr_write_resume_state(path, state$fingerprint)
  }
  for (n_workers in c(1L, 2L)) for (keep in c(0L, 1L, 3L)) {
    path <- tempfile("v050-interrupt-")
    v050_run(f, path, n_workers = n_workers)
    interrupt(path, keep)
    resumed <- v050_run(f, path, resume = TRUE, n_workers = n_workers)
    expect_identical(resumed$run_manifest$n_pairs_resumed, keep)
    expect_identical(v050_numbers(resumed), reference)
    expect_true(file.exists(file.path(path, "batch_main_results.csv")))
    # The interrupted resumed run can itself be interrupted and resumed again.
    interrupt(path, max(keep, 2L))
    expect_identical(v050_numbers(v050_run(f, path, resume = TRUE)), reference)
  }
})

test_that("fail_fast stops at the first invalid pair before writing output", {
  f <- v050_batch_fixture()
  f$pairs <- rbind(f$pairs[1:2, ], data.frame(exposure = "x", outcome = "missing_trait"), f$pairs[3:4, ])
  path <- tempfile("v050-failfast-")
  expect_error(v050_run(f, path, fail_fast = TRUE), "missing_trait")
  expect_false(dir.exists(path))
  resumed <- v050_run(f, path)
  expect_identical(resumed$main_results$status, c("ok", "ok", "error", "ok", "ok"))
})

test_that("changed inputs, options or files are refused without modifying the checkpoint", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-refuse-")
  v050_run(f, path)
  before <- v050_hashes(path)
  changed <- f
  changed$options$standardize_phenotypes <- FALSE
  expect_error(v050_run(changed, path, resume = TRUE), "fingerprint changed")
  changed <- f
  changed$resources$pheno$x[1L] <- changed$resources$pheno$x[1L] + 0.01
  expect_error(v050_run(changed, path, resume = TRUE), "fingerprint changed")
  expect_error(v050_run(f, path, resume = TRUE, heterogeneity_method = "none"), "fingerprint changed")
  changed <- f
  changed$instruments <- data.frame(exposure = c("x", "x", "w"), marker = c("m1", "m2", "m2"))
  expect_error(v050_run(changed, path, resume = TRUE), "fingerprint changed")
  expect_error(v050_run(f, path), "already exists")
  expect_identical(v050_hashes(path), before)

  table <- file.path(path, "batch_main_results.csv")
  writeLines(c(readLines(table), "tampered"), table)
  tampered_table <- v050_hashes(path)
  expect_error(v050_run(f, path, resume = TRUE), "tables changed")
  expect_identical(v050_hashes(path), tampered_table)

  path <- tempfile("v050-refuse-record-")
  v050_run(f, path)
  record_file <- sort(list.files(file.path(path, "batch_pair_records"), full.names = TRUE))[[1L]]
  record <- readRDS(record_file)
  record$payload$main$theta <- record$payload$main$theta + 1
  saveRDS(record, record_file)
  tampered_record <- v050_hashes(path)
  expect_error(v050_run(f, path, resume = TRUE), "changed or does not belong")
  expect_identical(v050_hashes(path), tampered_record)

  writeLines("not an rds", record_file)
  expect_error(v050_run(f, path, resume = TRUE), "unreadable pair record")

  legacy <- tempfile("v050-legacy-")
  dir.create(legacy)
  utils::write.csv(data.frame(pair_id = "x__to__y", status = "ok"), file.path(legacy, "batch_main_results.csv"), row.names = FALSE)
  saveRDS(list(fingerprint = list(schema = "qtl2mr-batch-v040-1")), file.path(legacy, "batch_resume_state.rds"))
  expect_error(v050_run(f, legacy, resume = TRUE), "fingerprint changed")
})

test_that("overwrite replaces only package-owned batch outputs", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-overwrite-")
  v050_run(f, path)
  writeLines("keep me", file.path(path, "notes.txt"))
  dir.create(file.path(path, "by_exposure"))
  changed <- f
  changed$options$standardize_phenotypes <- FALSE
  rerun <- v050_run(changed, path, overwrite = TRUE)
  expect_true(all(rerun$main_results$status == "ok"))
  expect_identical(readLines(file.path(path, "notes.txt")), "keep me")
  expect_false(dir.exists(file.path(path, "by_exposure")))
  expect_identical(v050_numbers(v050_run(changed, path, resume = TRUE)), v050_numbers(rerun))
})

test_that("invalid pairs become error rows", {
  f <- v050_batch_fixture()
  f$pairs <- rbind(f$pairs, data.frame(exposure = c("absent", "x", "y"), outcome = c("y", "x", "x")))
  b <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  expect_identical(b$main_results$status, c(rep("ok", 4L), "error", "error", "error"))
  expect_match(b$main_results$error_message[[5L]], "no_instruments")
  expect_match(b$main_results$error_message[[6L]], "must be different traits")
  expect_match(b$main_results$error_message[[7L]], "no_instruments: instruments has no selected loci for exposure y")
  expect_identical(names(b$main_results), names(qtl2mr:::qtl2mr_empty_batch_table("main_results")))
})

test_that("sample exclusions are reported once per trait of the valid pairs", {
  f <- v050_batch_fixture()
  f$resources$pheno$broken <- NA_real_
  f$resources$pheno$broken[1:20] <- 1
  f$pairs <- rbind(f$pairs, data.frame(exposure = "absent", outcome = "broken"))
  b <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  expect_identical(b$main_results$status, c(rep("ok", 4L), "error"))
  expect_setequal(unique(b$sample_exclusions$trait), c("x", "w", "y", "v"))
  v <- b$sample_exclusions[b$sample_exclusions$trait == "v" & b$sample_exclusions$reason == "nonfinite_value", ]
  expect_identical(v$excluded_ids, "s1;s2;s3")
  expect_equal(b$main_results$n_x_median, c(48, 48, 48, 48, NA))
  expect_equal(b$main_results$n_y_median, c(48, 45, 48, 45, NA))
})

test_that("only package-written tables are checksummed for resume", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-user-files-")
  first <- v050_run(f, path)
  writeLines("a,b", file.path(path, "batch_user_notes.csv"))
  expect_identical(v050_numbers(v050_run(f, path, resume = TRUE)), v050_numbers(first))
  unlink(file.path(path, "batch_user_notes.csv"))
  expect_identical(v050_numbers(v050_run(f, path, resume = TRUE)), v050_numbers(first))
  record <- qtl2mr:::qtl2mr_read_pair_record(path, 1L, first$main_results$pair_id[[1L]],
    readRDS(file.path(path, "batch_resume_state.rds"))$fingerprint)
  expect_identical(record$main$theta, first$main_results$theta[[1L]])
  expect_null(qtl2mr:::qtl2mr_read_pair_record(path, 1L, "another_pair", list()))
})

test_that("the resume fingerprint ignores probabilities at markers the run does not use", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-unused-marker-")
  first <- v050_run(f, path)
  extra <- f
  a <- extra$resources$allele_probs[["1"]]
  b <- array(c(a, a[, , 1L]), c(dim(a)[1:2], 2L), list(dimnames(a)[[1L]], dimnames(a)[[2L]], c("m1", "unused")))
  extra$resources$allele_probs[["1"]] <- b
  instruments <- qtl2mr:::qtl2mr_normalize_instruments(f$instruments, f$resources)
  expect_identical(qtl2mr:::qtl2mr_batch_fingerprint(f$resources, f$pairs, instruments, list(), markers = c("m1", "m2"))$data,
    qtl2mr:::qtl2mr_batch_fingerprint(extra$resources, f$pairs, instruments, list(), markers = c("m1", "m2"))$data)
})

test_that("the batch locus table matches the locus fits and error rows carry NA bounds", {
  f <- v050_batch_fixture()
  f$pairs <- rbind(f$pairs, data.frame(exposure = "absent", outcome = "y"))
  b <- run_qtl2mr_batch(f$resources, f$pairs, f$instruments, min_samples = 20L, progress = FALSE)
  loci <- b$locus_results
  expect_identical(names(loci), names(qtl2mr:::qtl2mr_empty_batch_table("locus_results")))
  expect_true(all(is.finite(loci$F_X_block) & loci$F_X_block > 0))
  expect_true(all(loci$partial_R2_X_block >= 0 & loci$partial_R2_X_block <= 1))
  expect_identical(loci$steiger_pass, loci$partial_R2_X_block >= loci$partial_R2_Y_block)
  expect_true(all(is.na(loci$hsq_x)))
  expect_true(is.na(b$main_results$confidence_set_bounded[b$main_results$status == "error"]))

  single <- qtl2mr:::run_prepared_qtl2mr(qtl2mr_pair_input(f$resources, f$instruments, "x", "y"), "x", "y",
    verbose = FALSE, heterogeneity_method = "none", diagnostic_bootstrap_B = 0L)
  expect_identical(loci$F_X_block[loci$pair_id == "x__to__y"], mr_block_first_stage_f(single$locus_blocks[[1L]]))
  expect_identical(loci$F_X_block[loci$pair_id == "x__to__y"], f$instruments$F_X_block[f$instruments$exposure == "x"])
})

test_that("linked batch outputs are refused before anything is written", {
  skip_on_os("windows")
  f <- v050_batch_fixture()
  root <- tempfile("v050-links-")
  path <- file.path(root, "run")
  v050_run(f, path)
  outside <- file.path(root, "outside.csv")
  table <- file.path(path, "batch_main_results.csv")
  expect_true(file.copy(table, outside))
  expect_true(file.remove(table))
  expect_true(file.symlink(outside, table))
  before <- tools::md5sum(outside)
  expect_error(v050_run(f, path, resume = TRUE), "symbolic links \\(batch_main_results.csv\\)")
  expect_error(v050_run(f, path, overwrite = TRUE), "symbolic links")
  expect_identical(tools::md5sum(outside), before)
  expect_identical(Sys.readlink(table), outside)

  record_dir <- file.path(root, "records")
  path2 <- file.path(root, "run2")
  v050_run(f, path2)
  record <- list.files(file.path(path2, "batch_pair_records"), full.names = TRUE)[[1L]]
  dir.create(record_dir)
  moved <- file.path(record_dir, basename(record))
  expect_true(file.rename(record, moved))
  expect_true(file.symlink(moved, record))
  expect_error(v050_run(f, path2, resume = TRUE), basename(record), fixed = TRUE)
})

test_that("instruments of exposures outside the run do not block a resume", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-unused-exposure-")
  first <- v050_run(f, path)
  unused <- f
  unused$instruments <- data.frame(exposure = c("x", "w", "y"), marker = c("m1", "m2", "m1"))
  again <- v050_run(unused, path, resume = TRUE)
  expect_identical(again$run_manifest$n_pairs_resumed, 4L)
  expect_identical(v050_numbers(again), v050_numbers(first))
})
