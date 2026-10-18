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
  resources <- qtl2mr_prepare_resources(gp, NULL, pheno, covar, scan, map)
  list(resources = resources,
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
  run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, output_dir = path,
    progress = FALSE, ...)
}

test_that("deprecated batch entry points equal run_qtl2mr_batch with mapped options", {
  f <- v050_batch_fixture()
  pair_common <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE)
  batch_common <- run_qtl2mr_batch(f$resources, f$pairs, sample_policy = "batch_common_complete",
    on_nonestimable_locus = "drop", fit_options = f$options, progress = FALSE)
  expect_true(all(pair_common$main_results$status == "ok"))
  expect_false(identical(v050_numbers(pair_common)$theta, v050_numbers(batch_common)$theta))

  expect_warning(grouped <- run_qtl2mr_batch_exact_grouped(f$resources, f$pairs,
    fit_options = f$options, progress = FALSE), "deprecated")
  expect_identical(v050_numbers(grouped), v050_numbers(pair_common))
  expect_warning(fast <- run_qtl2mr_batch_fast(f$resources, f$pairs, fit_options = f$options,
    progress = FALSE), "deprecated")
  expect_identical(v050_numbers(fast), v050_numbers(batch_common))
  for (fast_mode in c(TRUE, FALSE)) {
    expect_warning(par <- run_qtl2mr_batch_parallel_by_exposure(f$resources, f$pairs, fast = fast_mode,
      n_workers = 2L, fit_options = f$options, progress = FALSE), "deprecated")
    expect_identical(v050_numbers(par), v050_numbers(if (fast_mode) batch_common else pair_common))
  }
  expect_warning(legacy <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options,
    sample_rule = list(mode = "batch_common_complete"), progress = FALSE), "deprecated")
  expect_identical(v050_numbers(legacy)$theta, v050_numbers(run_qtl2mr_batch(f$resources, f$pairs,
    sample_policy = "batch_common_complete", fit_options = f$options, progress = FALSE))$theta)
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, sample_polcy = "batch_common_complete"), "Unknown argument")
})

test_that("worker count does not change results", {
  f <- v050_batch_fixture()
  serial <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE)
  parallel <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, n_workers = 2L, progress = FALSE)
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
  reference <- v050_numbers(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE))
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
  expect_error(v050_run(f, path, resume = TRUE, sample_policy = "batch_common_complete"), "fingerprint changed")
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

test_that("force replaces only package-owned batch outputs", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-force-")
  v050_run(f, path)
  writeLines("keep me", file.path(path, "notes.txt"))
  dir.create(file.path(path, "by_exposure"))
  changed <- f
  changed$options$standardize_phenotypes <- FALSE
  rerun <- v050_run(changed, path, force = TRUE)
  expect_true(all(rerun$main_results$status == "ok"))
  expect_identical(readLines(file.path(path, "notes.txt")), "keep me")
  expect_false(dir.exists(file.path(path, "by_exposure")))
  expect_identical(v050_numbers(v050_run(changed, path, resume = TRUE)), v050_numbers(rerun))
})

test_that("batch-common retention gate stops before writing any output", {
  f <- v050_batch_fixture()
  path <- tempfile("v050-gate-")
  expect_error(v050_run(f, path, sample_policy = "batch_common_complete", min_sample_retention = 0.99),
    "min_sample_retention")
  expect_false(dir.exists(path))
  ok <- v050_run(f, path, sample_policy = "batch_common_complete", min_sample_retention = 0.9)
  expect_identical(ok$sample_qc$sample_retention_gate_status, "pass")
  expect_identical(unique(ok$main_results$sample_policy), "batch_common_complete")
  expect_identical(ok$run_manifest$sample_retention_min, 0.9)
})

test_that("invalid pairs become error rows in every sample policy", {
  f <- v050_batch_fixture()
  f$pairs <- rbind(f$pairs, data.frame(exposure = c("absent", "x"), outcome = c("y", "x")))
  for (policy in c("pair_common_complete", "batch_common_complete")) {
    b <- run_qtl2mr_batch(f$resources, f$pairs, sample_policy = policy, fit_options = f$options, progress = FALSE)
    expect_identical(b$main_results$status, c(rep("ok", 4L), "error", "error"))
    expect_match(b$main_results$error_message[[6L]], "must be different traits")
    expect_identical(utils::tail(names(b$main_results), 1L), "sample_policy")
  }
})

test_that("the deprecated sample_rule keeps the 0.4.x behaviour of no retention gate", {
  f <- v050_batch_fixture()
  f$resources$pheno$y[4:12] <- NA_real_
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, sample_policy = "batch_common_complete",
    fit_options = f$options, progress = FALSE), "below min_sample_retention")
  expect_warning(legacy <- run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE,
    sample_rule = list(mode = "batch_common_complete")), "deprecated")
  expect_true(all(legacy$main_results$status == "ok"))
  expect_identical(legacy$run_manifest$sample_retention_min, 0)
  expect_error(suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, sample_policy = "pair_common_complete",
    sample_rule = list(mode = "batch_common_complete"))), "not both")
  # An unnamed value reaches ... only after every formal argument is supplied.
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, list(), "pair_common_complete", "fail", 0.85,
    f$options, 1L, FALSE, FALSE, FALSE, NULL, FALSE, "extra"), "must be named")
})

test_that("sample_ids restrictions are applied in preflight, the gate and the exclusion ledger", {
  f <- v050_batch_fixture()
  keep <- paste0("s", 4:33)  # s1-s3 lack outcome v, so the batch-common set must exclude them
  for (mode in c("common_complete", "batch_common_complete")) {
    b <- suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, fit_options = f$options, progress = FALSE,
      sample_rule = list(mode = mode, sample_ids = keep), min_sample_retention = 0.85))
    expect_true(all(b$main_results$status == "ok"))
    expect_true(all(b$main_results$n_x_median <= 30))
    ledger <- b$sample_exclusions[b$sample_exclusions$reason == "not_in_sample_ids", ]
    expect_true(all(ledger$n_excluded == 18L))
  }
  small <- suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, fit_options = list(min_samples = 40L),
    progress = FALSE, sample_rule = list(mode = "common_complete", sample_ids = keep)))
  expect_true(all(small$errors$stage == "preflight"))
  expect_true(all(grepl("too_few_samples", small$errors$error_message)))
})

test_that("batch-common exclusions only count traits of valid pairs", {
  f <- v050_batch_fixture()
  f$resources$pheno$broken <- NA_real_
  f$resources$pheno$broken[1:20] <- 1
  f$pairs <- rbind(f$pairs, data.frame(exposure = "absent", outcome = "broken"))
  b <- run_qtl2mr_batch(f$resources, f$pairs, sample_policy = "batch_common_complete",
    fit_options = f$options, progress = FALSE)
  expect_identical(b$main_results$status, c(rep("ok", 4L), "error"))
  other <- b$sample_exclusions[b$sample_exclusions$reason == "nonfinite_other_batch_trait" &
    b$sample_exclusions$pair_id == "x__to__y", ]
  expect_identical(other$n_excluded, 3L)
  expect_true(all(b$main_results$n_x_median[1:4] == 45))
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
  expect_identical(qtl2mr:::qtl2mr_batch_fingerprint(f$resources, f$pairs, list(), markers = c("m1", "m2"))$data,
    qtl2mr:::qtl2mr_batch_fingerprint(extra$resources, f$pairs, list(), markers = c("m1", "m2"))$data)
})
