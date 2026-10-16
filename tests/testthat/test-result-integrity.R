qtl2mr_file_bytes <- function(path) readBin(path, what = "raw", n = file.info(path)$size)

qtl2mr_tree_md5 <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  stats::setNames(unname(tools::md5sum(files)),
                  vapply(files, qtl2mr:::qtl2mr_relative_path, character(1), root = path))
}

qtl2mr_test_remove_link_node <- function(path) {
  if (identical(qtl2mr:::qtl2mr_path_node_state(path), "link")) {
    suppressWarnings(fs::link_delete(path))
  }
  invisible(NULL)
}

qtl2mr_prepared_example_run <- function(output_dir, overwrite = FALSE, progress = FALSE) {
  data("qtl2mr_example", package = "qtl2mr")
  prepared <- as_qtl2mr_prepared_input(
    allele_probs = qtl2mr_example$allele_probs,
    pheno = qtl2mr_example$pheno,
    covar = qtl2mr_example$covar,
    kinship = qtl2mr_example$kinship,
    marker_map = qtl2mr_example$marker_map,
    selected_loci = qtl2mr_example$selected_loci,
    exposure = qtl2mr_example$exposure,
    outcome = qtl2mr_example$outcome,
    metadata = qtl2mr_example$metadata
  )
  run_qtl2mr_analysis(
    prepared, adapter = "prepared", loci = 2, output_dir = output_dir,
    overwrite = overwrite, make_figures = FALSE, make_html = TRUE, verbose = FALSE,
    run_qtl2_validation = FALSE, run_heteroskedasticity = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "chisq",
    run_loo_mr = FALSE, run_instrument_sensitivity = FALSE, progress = progress
  )
}

test_that("path-node classifier distinguishes ordinary, absent, and link nodes", {
  root <- tempfile("qtl2mr-path-node-")
  dir.create(root)
  file <- file.path(root, "ordinary.txt")
  writeLines("ordinary", file)
  expect_identical(qtl2mr:::qtl2mr_path_node_state(root), "directory")
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file), "file")
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file.path(root, "absent.txt")), "absent")

  target <- tempfile("qtl2mr-path-node-target-")
  dir.create(target)
  link <- file.path(root, "directory-link")
  if (isTRUE(file.symlink(target, link))) {
    on.exit(qtl2mr_test_remove_link_node(link), add = TRUE)
    expect_identical(qtl2mr:::qtl2mr_path_node_state(link), "link")
    expect_identical(qtl2mr:::qtl2mr_path_node_state(paste0(link, "/")), "link")
  } else {
    skip("directory symlink creation is unavailable on this platform")
  }

  broken <- file.path(root, "broken-link")
  if (isTRUE(file.symlink("missing-target", broken))) {
    on.exit(qtl2mr_test_remove_link_node(broken), add = TRUE)
    expect_identical(qtl2mr:::qtl2mr_path_node_state(broken), "link")
  } else {
    skip("broken symlink creation is unavailable on this platform")
  }
})

test_that("Windows junctions are classified as links whenever they can be created", {
  skip_if_not(identical(.Platform$OS.type, "windows"))
  target <- tempfile("qtl2mr-junction-target-")
  junction <- tempfile("qtl2mr-junction-")
  dir.create(target)
  status <- suppressWarnings(system2(
    "cmd", c("/c", "mklink", "/J", shQuote(junction), shQuote(target)),
    stdout = FALSE, stderr = FALSE
  ))
  if (!identical(as.integer(status), 0L)) skip("junction creation is unavailable on this Windows runner")
  on.exit(qtl2mr_test_remove_link_node(junction), add = TRUE)
  expect_identical(qtl2mr:::qtl2mr_path_node_state(junction), "link")
  expect_identical(qtl2mr:::qtl2mr_path_node_state(paste0(junction, "/")), "link")
})

test_that("unverifiable fs metadata fails closed before writer, transaction, and lite mutations", {
  out_dir <- tempfile("qtl2mr-unverifiable-output-")
  dir.create(out_dir)
  sentinel <- file.path(out_dir, "sentinel.txt")
  writeLines("preserve", sentinel)
  out_key <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  parent_info <- qtl2mr:::qtl2mr_path_node_info
  testthat::local_mocked_bindings(
    qtl2mr_path_node_info = function(path) {
      if (identical(normalizePath(path, winslash = "/", mustWork = FALSE), out_key)) {
        return(data.frame(type = NA_character_))
      }
      parent_info(path)
    },
    .package = "qtl2mr"
  )
  expect_identical(qtl2mr:::qtl2mr_path_node_state(out_dir), "unknown")
  expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                    overwrite = TRUE, make_html = FALSE, make_figures = FALSE),
               "Unable to verify output path safety")
  expect_true(file.exists(sentinel))
})

test_that("NA metadata is absent only after a verified parent listing", {
  root <- tempfile("qtl2mr-na-metadata-")
  dir.create(root)
  candidate <- file.path(root, "candidate")
  original_info <- qtl2mr:::qtl2mr_path_node_info
  testthat::local_mocked_bindings(
    qtl2mr_path_node_info = function(path) {
      if (identical(normalizePath(path, winslash = "/", mustWork = FALSE),
                    normalizePath(candidate, winslash = "/", mustWork = FALSE))) {
        return(data.frame(type = NA_character_))
      }
      original_info(path)
    },
    qtl2mr_path_parent_entries = function(...) stop("injected parent-listing failure", call. = FALSE),
    .package = "qtl2mr"
  )
  expect_identical(qtl2mr:::qtl2mr_path_node_state(candidate), "unknown")
  expect_error(qtl2mr:::qtl2mr_assert_known_path_node(candidate), "Unable to verify output path safety")
})

test_that("unverifiable metadata preserves transaction backup before commit mutation", {
  out_dir <- tempfile("qtl2mr-unverifiable-transaction-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  dir.create(out_dir)
  external <- file.path(out_dir, "external.txt")
  writeLines("external", external)
  old_backup <- transaction$backup_dir
  out_key <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  on.exit({
    if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
    if (dir.exists(old_backup)) unlink(old_backup, recursive = TRUE, force = TRUE)
    if (dir.exists(transaction$stage_dir)) unlink(transaction$stage_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  original_info <- qtl2mr:::qtl2mr_path_node_info
  testthat::local_mocked_bindings(
    qtl2mr_path_node_info = function(path) {
      if (identical(normalizePath(path, winslash = "/", mustWork = FALSE), out_key)) {
        return(data.frame(type = NA_character_))
      }
      original_info(path)
    },
    .package = "qtl2mr"
  )
  expect_identical(qtl2mr:::qtl2mr_path_node_state(out_dir), "unknown")
  expect_error(qtl2mr:::qtl2mr_commit_output_transaction(transaction),
               "Unable to verify output path safety")
  expect_true(file.exists(external))
  expect_true(dir.exists(old_backup))
})

test_that("unverifiable source metadata blocks lite export before publication", {
  source_dir <- tempfile("qtl2mr-unverifiable-lite-source-")
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), source_dir,
                                  make_html = FALSE, make_figures = FALSE)
  source_member <- file.path(source_dir, "result_object.rds")
  source_key <- normalizePath(source_member, winslash = "/", mustWork = FALSE)
  dest_zip <- tempfile("qtl2mr-unverifiable-lite-", fileext = ".zip")
  original_info <- qtl2mr:::qtl2mr_path_node_info
  testthat::local_mocked_bindings(
    qtl2mr_path_node_info = function(path) {
      if (identical(normalizePath(path, winslash = "/", mustWork = FALSE), source_key)) {
        return(data.frame(type = NA_character_))
      }
      original_info(path)
    },
    .package = "qtl2mr"
  )
  expect_error(export_qtl2mr_lite_result(written, dest_zip = dest_zip),
               "Unable to verify result tree safety")
  expect_false(file.exists(dest_zip))
})

test_that("overwrite rejects unknown directories and preserves their contents", {
  out_dir <- tempfile("qtl2mr-unknown-output-")
  dir.create(out_dir)
  sentinel <- file.path(out_dir, "sentinel.txt")
  writeLines("do not replace", sentinel)
  before <- qtl2mr_file_bytes(sentinel)
  expect_error(
    write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir, overwrite = TRUE,
                          make_html = FALSE, make_figures = FALSE),
    "recognized prior qtl2MR"
  )
  expect_identical(qtl2mr_file_bytes(sentinel), before)
})

test_that("nested missing output parents are created only after component checks", {
  root <- tempfile("qtl2mr-nested-output-")
  out_dir <- file.path(root, "a", "b", "result")
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                  make_html = FALSE, make_figures = FALSE)
  expect_s3_class(written, "qtl2mr_analysis_result")
  expect_true(file.exists(file.path(out_dir, "result_object.rds")))
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file.path(root, "a")), "directory")
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file.path(root, "a", "b")), "directory")
})

test_that("nested missing lite ZIP parents are created only after component checks", {
  result_dir <- tempfile("qtl2mr-nested-lite-result-")
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), result_dir,
                                  make_html = TRUE, make_figures = FALSE)
  root <- tempfile("qtl2mr-nested-lite-")
  dest_zip <- file.path(root, "a", "b", "result_lite.zip")
  lite <- export_qtl2mr_lite_result(written, dest_zip = dest_zip, include_rds = FALSE)
  expect_identical(normalizePath(lite), normalizePath(dest_zip))
  expect_true(file.exists(dest_zip))
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file.path(root, "a")), "directory")
  expect_identical(qtl2mr:::qtl2mr_path_node_state(file.path(root, "a", "b")), "directory")
})

test_that("ancestor directory links publish result output at the explicit canonical target", {
  root <- tempfile("qtl2mr-ancestor-link-root-")
  target <- tempfile("qtl2mr-ancestor-link-target-")
  dir.create(root)
  dir.create(file.path(target, "sub"), recursive = TRUE)
  jump <- file.path(root, "jump")
  if (!isTRUE(file.symlink(target, jump))) skip("ancestor directory symlink creation is unavailable on this platform")
  on.exit(qtl2mr_test_remove_link_node(jump), add = TRUE)
  requested <- file.path(jump, "sub", "result")
  canonical <- file.path(normalizePath(target, winslash = "/"), "sub", "result")
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), requested,
                                  make_html = FALSE, make_figures = FALSE)
  actual <- gsub("\\\\", "/", written$output_dir)
  expect_identical(actual, canonical)
  expect_false(identical(actual, gsub("\\\\", "/", qtl2mr:::qtl2mr_trim_node_path(requested))))
  expect_true(file.exists(file.path(canonical, "result_object.rds")))
  expect_false(any(grepl("qtl2mr-stage", list.files(file.path(target, "sub"), all.files = TRUE))))
})

test_that("ancestor directory links publish lite ZIPs at the explicit canonical target", {
  result_dir <- tempfile("qtl2mr-ancestor-link-lite-result-")
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), result_dir,
                                  make_html = TRUE, make_figures = FALSE)
  root <- tempfile("qtl2mr-ancestor-link-lite-root-")
  target <- tempfile("qtl2mr-ancestor-link-lite-target-")
  dir.create(root)
  dir.create(file.path(target, "sub"), recursive = TRUE)
  jump <- file.path(root, "jump")
  if (!isTRUE(file.symlink(target, jump))) skip("ancestor directory symlink creation is unavailable on this platform")
  on.exit(qtl2mr_test_remove_link_node(jump), add = TRUE)
  requested <- file.path(jump, "sub", "result_lite.zip")
  canonical <- file.path(normalizePath(target, winslash = "/"), "sub", "result_lite.zip")
  lite <- export_qtl2mr_lite_result(written, dest_zip = requested, include_rds = FALSE)
  actual <- gsub("\\\\", "/", lite)
  expect_identical(actual, canonical)
  expect_false(identical(actual, gsub("\\\\", "/", qtl2mr:::qtl2mr_trim_node_path(requested))))
  expect_true(file.exists(canonical))
})

test_that("overwrite rejects malformed lookalike results without moving them", {
  out_dir <- tempfile("qtl2mr-lookalike-output-")
  dir.create(file.path(out_dir, "data"), recursive = TRUE)
  writeLines("not an RDS", file.path(out_dir, "result_object.rds"))
  writeLines("not a qtl2MR summary", file.path(out_dir, "summary.md"))
  utils::write.csv(data.frame(
    role = c("model_object", "summary", "data", "data"),
    path_relative = c("result_object.rds", "summary.md", "../escape.csv", "data/warnings.csv"),
    path_absolute = rep(out_dir, 4L), stringsAsFactors = FALSE
  ), file.path(out_dir, "data", "output_manifest.csv"), row.names = FALSE)
  sentinel <- file.path(out_dir, "sentinel.txt")
  writeLines("do not move", sentinel)
  expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                    overwrite = TRUE, make_html = FALSE, make_figures = FALSE),
               "recognized prior qtl2MR")
  expect_true(file.exists(sentinel))
})

test_that("prior-result signature files must physically exist", {
  out_dir <- tempfile("qtl2mr-safe-lookalike-output-")
  dir.create(file.path(out_dir, "data"), recursive = TRUE)
  writeLines("placeholder", file.path(out_dir, "result_object.rds"))
  writeLines("placeholder", file.path(out_dir, "summary.md"))
  utils::write.csv(data.frame(
    role = c("model_object", "summary", "data", "data", "data"),
    path_relative = c("result_object.rds", "summary.md", "data/output_manifest.csv", "data/main_result.csv", "data/warnings.csv"),
    path_absolute = rep(out_dir, 5L), stringsAsFactors = FALSE
  ), file.path(out_dir, "data", "output_manifest.csv"), row.names = FALSE)
  expect_false(qtl2mr:::qtl2mr_is_prior_result_dir(out_dir))
  expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                    overwrite = TRUE, make_html = FALSE, make_figures = FALSE),
               "recognized prior qtl2MR")
})

test_that("directories at required signature paths do not qualify as prior results", {
  out_dir <- tempfile("qtl2mr-directory-signature-output-")
  dir.create(file.path(out_dir, "data", "main_result.csv"), recursive = TRUE)
  dir.create(file.path(out_dir, "data", "warnings.csv"))
  dir.create(file.path(out_dir, "result_object.rds"))
  dir.create(file.path(out_dir, "summary.md"))
  utils::write.csv(data.frame(
    role = c("model_object", "summary", "data", "data", "data"),
    path_relative = c("result_object.rds", "summary.md", "data/output_manifest.csv", "data/main_result.csv", "data/warnings.csv"),
    path_absolute = rep(out_dir, 5L), stringsAsFactors = FALSE
  ), file.path(out_dir, "data", "output_manifest.csv"), row.names = FALSE)
  expect_false(qtl2mr:::qtl2mr_is_prior_result_dir(out_dir))
})

test_that("recognized results replace atomically without stale files", {
  out_dir <- tempfile("qtl2mr-replace-output-")
  old <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                              make_html = TRUE, make_figures = FALSE)
  stale <- file.path(out_dir, "data", "stale_from_previous_run.csv")
  writeLines("stale", stale)
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                  overwrite = TRUE, make_html = TRUE, make_figures = FALSE)
  expect_s3_class(written, "qtl2mr_analysis_result")
  expect_false(file.exists(stale))
  manifest <- utils::read.csv(file.path(out_dir, "data", "output_manifest.csv"), stringsAsFactors = FALSE)
  expect_false(any(grepl("stale_from_previous_run", manifest$path_relative, fixed = TRUE)))
  expect_identical(normalizePath(written$output_dir), normalizePath(out_dir))
  expect_false(any(grepl("qtl2mr-stage|qtl2mr-backup", manifest$path_absolute)))
  expect_true(file.exists(file.path(out_dir, "result_object.rds")))
  expect_false(identical(old$output_dir, written$output_dir) && !file.exists(file.path(out_dir, "summary.md")))
})

test_that("failed writer transaction restores the old result exactly", {
  out_dir <- tempfile("qtl2mr-restore-output-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  sentinel <- file.path(out_dir, "data", "restore-sentinel.txt")
  writeLines("preserve me", sentinel)
  before <- qtl2mr_file_bytes(sentinel)
  before_tree <- qtl2mr_tree_md5(out_dir)
  broken <- make_fake_qtl2mr_analysis_result()
  broken$pilot <- NULL
  expect_error(
    write_qtl2mr_results(broken, out_dir, overwrite = TRUE, make_html = FALSE, make_figures = FALSE)
  )
  expect_true(file.exists(sentinel))
  expect_identical(qtl2mr_file_bytes(sentinel), before)
  expect_identical(qtl2mr_tree_md5(out_dir), before_tree)
  expect_true(file.exists(file.path(out_dir, "result_object.rds")))
  siblings <- list.files(dirname(out_dir), pattern = paste0("^\\.", basename(out_dir), "\\.qtl2mr-(stage|backup)-"))
  expect_length(siblings, 0L)
})

test_that("transaction conflicts preserve an external directory and the recoverable backup", {
  out_dir <- tempfile("qtl2mr-race-output-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  old_sentinel <- file.path(out_dir, "data", "old-sentinel.txt")
  writeLines("old result must survive", old_sentinel)
  old_tree <- qtl2mr_tree_md5(out_dir)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  expect_s3_class(transaction, "qtl2mr_output_transaction")
  expect_false(dir.exists(out_dir))
  dir.create(out_dir)
  external <- file.path(out_dir, "external-sentinel.txt")
  writeLines("external result must survive", external)
  expect_error(qtl2mr:::qtl2mr_commit_output_transaction(transaction), "external destination preserved")
  recovery <- transaction$backup_dir
  expect_identical(transaction$state, "aborted")
  qtl2mr:::qtl2mr_abort_output_transaction(transaction)
  expect_true(file.exists(external))
  expect_true(dir.exists(recovery))
  expect_identical(qtl2mr_tree_md5(recovery), old_tree)
})

test_that("ordinary-file transaction conflict preserves the external node and backup", {
  out_dir <- tempfile("qtl2mr-file-race-output-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  old_tree <- qtl2mr_tree_md5(out_dir)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  writeLines("external file", out_dir)
  conflict <- tryCatch(qtl2mr:::qtl2mr_commit_output_transaction(transaction), error = identity)
  expect_s3_class(conflict, "qtl2mr_output_transaction_conflict")
  expect_match(conditionMessage(conflict), "external destination preserved")
  expect_identical(transaction$state, "aborted")
  expect_true(qtl2mr:::qtl2mr_is_regular_file(out_dir))
  expect_true(dir.exists(transaction$backup_dir))
  expect_identical(qtl2mr_tree_md5(transaction$backup_dir), old_tree)
  expect_false(dir.exists(transaction$stage_dir))
})

test_that("broken-symlink transaction conflict preserves the node and backup", {
  out_dir <- tempfile("qtl2mr-broken-link-race-output-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  if (!isTRUE(file.symlink("missing-external-target", out_dir))) skip("broken symlink creation is unavailable on this platform")
  on.exit({
    qtl2mr_test_remove_link_node(out_dir)
    if (dir.exists(transaction$backup_dir)) unlink(transaction$backup_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  conflict <- tryCatch(qtl2mr:::qtl2mr_commit_output_transaction(transaction), error = identity)
  expect_s3_class(conflict, "qtl2mr_output_transaction_conflict")
  expect_true(qtl2mr:::qtl2mr_is_symlink(out_dir))
  expect_true(dir.exists(transaction$backup_dir))
  expect_identical(transaction$state, "aborted")
})

test_that("published transaction rollback removes only its owned destination", {
  out_dir <- tempfile("qtl2mr-published-rollback-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  old_sentinel <- file.path(out_dir, "data", "old-sentinel.txt")
  writeLines("restore after publish", old_sentinel)
  old_tree <- qtl2mr_tree_md5(out_dir)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  writeLines("staged", file.path(transaction$stage_dir, "staged.txt"))
  qtl2mr:::qtl2mr_commit_output_transaction(transaction, keep_backup = TRUE)
  expect_true(file.exists(file.path(out_dir, ".qtl2mr-transaction-token")))
  qtl2mr:::qtl2mr_abort_output_transaction(transaction)
  expect_identical(qtl2mr_tree_md5(out_dir), old_tree)
  expect_false(dir.exists(transaction$backup_dir))
})

test_that("published new-output rollback removes the token-owned partial publication", {
  out_dir <- tempfile("qtl2mr-new-output-rollback-")
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  writeLines("partial", file.path(transaction$stage_dir, "partial.txt"))
  qtl2mr:::qtl2mr_commit_output_transaction(transaction, keep_backup = TRUE)
  expect_true(qtl2mr:::qtl2mr_transaction_owned(out_dir, transaction))
  expect_null(transaction$backup_dir)
  qtl2mr:::qtl2mr_abort_output_transaction(transaction)
  expect_false(dir.exists(out_dir))
  expect_identical(transaction$state, "aborted")
})

test_that("writer failure after publication rolls a new destination back to absent", {
  out_dir <- tempfile("qtl2mr-writer-post-publish-")
  testthat::local_mocked_bindings(
    qtl2mr_finalize_output_transaction = function(...) stop("injected post-publish failure", call. = FALSE),
    .package = "qtl2mr"
  )
  expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                                    make_html = FALSE, make_figures = FALSE),
               "injected post-publish failure")
  expect_false(dir.exists(out_dir))
})

test_that("token creation failure restores a prior result", {
  out_dir <- tempfile("qtl2mr-token-failure-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  before <- qtl2mr_tree_md5(out_dir)
  testthat::local_mocked_bindings(
    qtl2mr_write_output_transaction_token = function(...) stop("injected token failure", call. = FALSE),
    .package = "qtl2mr"
  )
  expect_error(qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE), "token")
  expect_identical(qtl2mr_tree_md5(out_dir), before)
})

test_that("backup-cleanup failure retains ownership for a recoverable rollback", {
  out_dir <- tempfile("qtl2mr-backup-cleanup-failure-")
  write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), out_dir,
                       make_html = FALSE, make_figures = FALSE)
  old_tree <- qtl2mr_tree_md5(out_dir)
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  writeLines("partial", file.path(transaction$stage_dir, "partial.txt"))
  qtl2mr:::qtl2mr_commit_output_transaction(transaction, keep_backup = TRUE)
  testthat::local_mocked_bindings(
    qtl2mr_remove_output_transaction_backup = function(...) stop("injected backup cleanup failure", call. = FALSE),
    .package = "qtl2mr"
  )
  expect_error(qtl2mr:::qtl2mr_finalize_output_transaction(transaction), "backup cleanup")
  expect_true(qtl2mr:::qtl2mr_transaction_owned(out_dir, transaction))
  qtl2mr:::qtl2mr_abort_output_transaction(transaction)
  expect_identical(qtl2mr_tree_md5(out_dir), old_tree)
})

test_that("token-cleanup failure finalizes without destructive rollback", {
  out_dir <- tempfile("qtl2mr-token-cleanup-failure-")
  transaction <- qtl2mr:::qtl2mr_begin_output_transaction(out_dir, overwrite = TRUE)
  writeLines("complete", file.path(transaction$stage_dir, "complete.txt"))
  qtl2mr:::qtl2mr_commit_output_transaction(transaction, keep_backup = TRUE)
  testthat::local_mocked_bindings(
    qtl2mr_remove_output_transaction_token = function(...) FALSE,
    .package = "qtl2mr"
  )
  expect_warning(qtl2mr:::qtl2mr_finalize_output_transaction(transaction), "token cleanup")
  expect_identical(transaction$state, "finalized")
  expect_true(dir.exists(out_dir))
  expect_true(file.exists(file.path(out_dir, ".qtl2mr-transaction-token")))
})

test_that("empty destinations, repeated overwrite, symlinks, and high-level replacement are safe", {
  empty_dir <- tempfile("qtl2mr-empty-output-")
  dir.create(empty_dir)
  first <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), empty_dir,
                                make_html = FALSE, make_figures = FALSE)
  second <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), empty_dir,
                                 overwrite = TRUE, make_html = FALSE, make_figures = FALSE)
  expect_true(file.exists(file.path(empty_dir, "result_object.rds")))
  expect_identical(normalizePath(first$output_dir), normalizePath(second$output_dir))

  target <- tempfile("qtl2mr-symlink-target-")
  dir.create(target)
  target_sentinel <- file.path(target, "sentinel.txt")
  writeLines("target intact", target_sentinel)
  link <- tempfile("qtl2mr-symlink-output-")
  if (isTRUE(file.symlink(target, link))) {
    on.exit(qtl2mr_test_remove_link_node(link), add = TRUE)
    expect_true(qtl2mr:::qtl2mr_is_symlink(link))
    expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), link,
                                      overwrite = TRUE, make_html = FALSE, make_figures = FALSE), "symlink")
    expect_true(file.exists(target_sentinel))
    expect_error(write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), paste0(link, "/"),
                                      overwrite = TRUE, make_html = FALSE, make_figures = FALSE), "symlink")
    expect_true(file.exists(target_sentinel))
  } else {
    skip("symlink creation is unavailable on this platform")
  }

  high_level_dir <- tempfile("qtl2mr-high-level-output-")
  qtl2mr_prepared_example_run(high_level_dir)
  stale <- file.path(high_level_dir, "data", "stale-high-level.txt")
  writeLines("stale", stale)
  rerun <- qtl2mr_prepared_example_run(high_level_dir, overwrite = TRUE)
  expect_false(file.exists(stale))
  expect_identical(normalizePath(rerun$output_dir), normalizePath(high_level_dir))

  logged_dir <- tempfile("qtl2mr-high-level-progress-")
  logged <- qtl2mr_prepared_example_run(logged_dir, progress = "log")
  progress_log <- file.path(logged$output_dir, "logs", "progress.log")
  expect_true(file.exists(progress_log))
  progress_text <- paste(readLines(progress_log, warn = FALSE), collapse = "\n")
  expect_true(grepl(logged$output_dir, progress_text, fixed = TRUE))
  expect_false(grepl("qtl2mr-stage|qtl2mr-backup", progress_text))
})

test_that("lite ZIP is rebuilt with relative, scrubbed portable artifacts", {
  source_root <- file.path(tempdir(), "qtl2mr-identifying-source-root")
  if (dir.exists(source_root)) unlink(source_root, recursive = TRUE, force = TRUE)
  written <- write_qtl2mr_results(make_fake_qtl2mr_analysis_result(), source_root,
                                  make_html = TRUE, make_figures = FALSE)
  writeLines(c(
    paste("source", normalizePath(source_root)),
    paste0("windows C:", "\\", "Users", "\\Example\\private\\result.csv"),
    paste0("uri file:", "///", "private/tmp/qtl2mr-secret.csv"),
    "web https://example.org/data/resource",
    "root /data/private.txt"
  ), file.path(source_root, "data", "provenance_note.txt"))
  writeLines(paste("log", normalizePath(source_root)), file.path(source_root, "logs", "progress.log"))
  zip_path <- file.path(tempdir(), "qtl2mr-portable-lite.zip")
  if (file.exists(zip_path)) unlink(zip_path)
  lite <- export_qtl2mr_lite_result(written, dest_zip = zip_path, max_file_mb = 10, include_rds = FALSE)
  listing <- utils::unzip(lite, list = TRUE)$Name
  expect_true(all(vapply(listing, qtl2mr:::qtl2mr_is_safe_relative_path, logical(1))))
  expect_false(any(grepl("[.]rds$", listing, ignore.case = TRUE)))
  expect_false(any(grepl("^logs/|progress[.]log$|run[.]log$", listing)))
  expect_false(any(grepl("data/output_manifest[.]csv$", listing)))
  expect_true(any(grepl("LITE_MANIFEST[.]csv$", listing)))
  expect_true(any(grepl("output_manifest[.]html$", listing)))

  unzip_dir <- tempfile("qtl2mr-portable-unzip-")
  dir.create(unzip_dir)
  utils::unzip(lite, exdir = unzip_dir)
  portable_files <- list.files(unzip_dir, recursive = TRUE, full.names = TRUE)
  text_files <- portable_files[tolower(tools::file_ext(portable_files)) %in%
    c("csv", "tsv", "txt", "log", "md", "html", "json", "yaml", "yml", "xml")]
  text <- paste(unlist(lapply(text_files, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl(normalizePath(source_root), text, fixed = TRUE))
  expect_false(qtl2mr:::qtl2mr_has_portable_path_leak(text, normalizePath(source_root), "not-present-lite-root"))
  lite_manifest <- utils::read.csv(file.path(unzip_dir, "LITE_MANIFEST.csv"), stringsAsFactors = FALSE)
  audit <- utils::read.csv(file.path(unzip_dir, "SMOKE_RESULT_SIZE_AUDIT.csv"), stringsAsFactors = FALSE)
  expect_false("path_absolute" %in% names(lite_manifest))
  expect_false("path_absolute" %in% names(audit))
  expect_true(all(vapply(lite_manifest$path_relative, qtl2mr:::qtl2mr_is_safe_relative_path, logical(1))))
  expect_setequal(lite_manifest$path_relative, listing)
  expect_true(grepl("https://example.org/data/resource", text, fixed = TRUE))
  expect_false(grepl("/data/private.txt", text, fixed = TRUE))
  non_self <- setdiff(lite_manifest$path_relative, c("LITE_MANIFEST.csv", "output_manifest.html"))
  expected_size <- round(file.info(file.path(unzip_dir, non_self))$size / 1024, 1)
  observed_size <- lite_manifest$size_kb[match(non_self, lite_manifest$path_relative)]
  expect_true(all(abs(expected_size - observed_size) <= 0.051))
  html <- paste(readLines(file.path(unzip_dir, "output_manifest.html"), warn = FALSE), collapse = "\n")
  expect_true(all(vapply(listing, function(path) grepl(path, html, fixed = TRUE), logical(1))))

  zip_hash <- qtl2mr_file_bytes(lite)
  expect_error(export_qtl2mr_lite_result(written, dest_zip = lite), "already exists")
  expect_identical(qtl2mr_file_bytes(lite), zip_hash)
  expect_error(export_qtl2mr_lite_result(written, dest_zip = file.path(source_root, "inside.zip")), "outside")
  zip_two <- file.path(tempdir(), "qtl2mr-portable-lite-two.zip")
  if (file.exists(zip_two)) unlink(zip_two)
  expect_true(file.exists(export_qtl2mr_lite_result(written, dest_zip = zip_two, include_rds = FALSE)))
})
