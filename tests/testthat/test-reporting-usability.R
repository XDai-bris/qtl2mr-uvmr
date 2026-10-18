test_that("analysis input must be prepared and retired arguments only warn", {
  expect_error(run_qtl2mr_analysis(list(genoprobs = list()), exposure = "x", outcome = "y"),
    "created by as_qtl2mr_prepared_input")
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  ex <- qtl2mr_example
  input <- as_qtl2mr_prepared_input(ex$allele_probs, ex$pheno, ex$covar, ex$kinship,
    ex$marker_map, ex$selected_loci, ex$exposure, ex$outcome)
  expect_warning(
    run_qtl2mr_analysis(input, loci = 1, output_dir = tempfile("qtl2mr-retired-"), adapter = "prepared",
      make_figures = FALSE, make_html = FALSE, verbose = FALSE, progress = FALSE,
      diagnostic_bootstrap_B = 0L, heterogeneity_method = "none",
      run_loo_mr = FALSE, run_instrument_sensitivity = FALSE),
    "adapter\\) is deprecated")
  expect_error(run_qtl2mr_analysis(input, adaptor = "prepared"), "Unknown, unnamed or duplicate")
  # adapter keeps its 0.4.x position: a positional value must not become output_dir.
  out <- tempfile("qtl2mr-positional-")
  expect_warning(result <- run_qtl2mr_analysis(input, "x", "y", "prepared", out, loci = 1,
    make_figures = FALSE, make_html = FALSE, verbose = FALSE, progress = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none",
    run_loo_mr = FALSE, run_instrument_sensitivity = FALSE), "adapter\\) is deprecated")
  expect_identical(normalizePath(result$output_dir), normalizePath(out))
})

test_that("relative output paths use a portable directory boundary", {
  windows_root <- "C:/Temp/Rtmpqtl2/lite"
  expect_equal(
    qtl2mr:::qtl2mr_relative_path(
      paste0(windows_root, "/report.html"),
      windows_root
    ),
    "report.html"
  )
  expect_equal(
    qtl2mr:::qtl2mr_relative_path(
      paste0(windows_root, "/data/output_manifest.csv"),
      windows_root
    ),
    "data/output_manifest.csv"
  )
  expect_false(identical(
    qtl2mr:::qtl2mr_relative_path(
      paste0(windows_root, "-other/report.html"),
      windows_root
    ),
    "report.html"
  ))
})

test_that("standard result exposes accessors, warning taxonomy, and portable manifest schema", {
  result <- make_fake_qtl2mr_analysis_result()
  out_dir <- tempfile("qtl2mr-round2-write-")
  written <- write_qtl2mr_results(result, output_dir = out_dir, make_html = TRUE, make_figures = TRUE)

  expect_s3_class(written, "qtl2mr_analysis_result")
  expect_true(all(c("metadata", "primary", "diagnostics", "tables", "figures", "warnings", "paths", "raw_objects") %in% names(written)))
  expect_equal(qtl2mr_primary_result(written)$theta[[1L]], written$pilot$summary$theta[[1L]])
  expect_true(nrow(qtl2mr_warnings(written)) >= 1L)
  expect_true(nrow(qtl2mr_method_comparison(written)) >= 1L)
  compact_paths <- qtl2mr_output_paths(written)
  full_paths <- qtl2mr_output_paths(written, detail = TRUE)
  expect_true(all(c("item", "path") %in% names(compact_paths)))
  expect_true("report.html" %in% compact_paths$path)
  expect_true("data/output_manifest.csv" %in% compact_paths$path)
  expect_lt(nrow(compact_paths), nrow(full_paths))
  expect_true(all(c("path_relative", "path_absolute") %in% names(full_paths)))

  manifest <- utils::read.csv(file.path(out_dir, "data", "output_manifest.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("role", "path_relative", "path_absolute", "format", "size_kb", "included_in_lite", "notes") %in% names(manifest)))
  expect_true("report.html" %in% manifest$path_relative)
  expect_true("summary.md" %in% manifest$path_relative)
  expect_false(any(grepl("^/", manifest$path_relative)))
  expect_false(any(grepl("^[A-Za-z]:[/\\\\]", manifest$path_relative)))

  warnings <- utils::read.csv(file.path(out_dir, "data", "warnings.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("consequence", "recommended_action", "shown_in_main_report") %in% names(warnings)))
  expect_true(all(warnings$severity %in% c("info", "method_limit", "warning", "error", "technical_qc")))
  if ("raw_space_psd_invalid" %in% warnings$warning_id) {
    raw <- warnings[warnings$warning_id == "raw_space_psd_invalid", , drop = FALSE]
    expect_equal(raw$severity[[1L]], "technical_qc")
    expect_false(isTRUE(raw$shown_in_main_report[[1L]]))
  }

  smry <- readLines(file.path(out_dir, "summary.md"), warn = FALSE)
  expect_false(any(grepl("^Adapter:", smry)))
  expect_true(any(grepl("^Null-model hsq:", smry)))
  expect_true(any(grepl("^Outputs: report.html", smry)))
  main_ids <- smry[grepl("^Actionable warning IDs:", smry)]
  expect_false(any(grepl("raw_space_psd_invalid", main_ids)))

  report <- paste(readLines(file.path(out_dir, "report.html"), warn = FALSE), collapse = "\n")
  img <- regmatches(report, gregexpr('src="[^"]+"', report))[[1L]]
  img <- substring(img, 6L, nchar(img) - 1L)
  expect_gt(length(img), 0L)
  expect_true(all(startsWith(img, "figures/")))
  expect_true(all(file.exists(file.path(out_dir, img))))
  expect_false(grepl("qtl2mr-stage-", report, fixed = TRUE))
  expect_true(all(file.exists(written$figure_manifest$path)))
  expect_identical(written$figures, written$figure_manifest)
  saved <- readRDS(file.path(out_dir, "result_object.rds"))
  expect_identical(saved$figure_manifest, written$figure_manifest)

  relocated <- tempfile("qtl2mr-relocated-report-")
  dir.create(relocated)
  expect_true(file.copy(file.path(out_dir, "report.html"), relocated))
  expect_true(file.copy(file.path(out_dir, "figures"), relocated, recursive = TRUE))
  expect_true(all(file.exists(file.path(relocated, img))))
  expect_false(grepl(normalizePath(out_dir, winslash = "/"), report, fixed = TRUE))
})

test_that("run reports are curated and package coverage is not claimed as run evidence", {
  result <- make_fake_qtl2mr_analysis_result()
  out_dir <- tempfile("qtl2mr-curated-report-")
  written <- write_qtl2mr_results(result, output_dir = out_dir, make_html = TRUE, make_figures = FALSE)

  coverage <- utils::read.csv(file.path(out_dir, "data", "validation_tests.csv"), stringsAsFactors = FALSE)
  expect_true(all(coverage$status == "not_run_in_analysis"))
  expect_false(any(coverage$status == "pass"))

  report <- paste(readLines(file.path(out_dir, "report.html"), warn = FALSE), collapse = "\n")
  expect_match(report, "Analysis overview and primary estimate", fixed = TRUE)
  expect_match(report, "QC and actionable warnings", fixed = TRUE)
  expect_false(grepl("qtl2 validation", report, fixed = TRUE))
  expect_match(report, "data/output_manifest.csv", fixed = TRUE)
  expect_false(grepl("Validation and Test Summary|current local package check", report))
  expect_false(grepl("validation_tests.csv", report, fixed = TRUE))
  expect_false(grepl("<h2>Figures</h2>", report, fixed = TRUE))
  expect_false(grepl("Leave-one-locus-out sensitivity|Instrument-set sensitivity", report))
  expect_false(grepl("<th>source</th>|<th>note</th>|Diagnostic sensitivity estimates", report, ignore.case = TRUE))
  expect_equal(length(gregexpr("Two-step IVW", report, fixed = TRUE)[[1L]][gregexpr("Two-step IVW", report, fixed = TRUE)[[1L]] > 0L]), 1L)

  smry <- readLines(file.path(out_dir, "summary.md"), warn = FALSE)
  expect_true(any(grepl("^Reported estimator:", smry)))
  expect_true(any(grepl("^Actionable warning IDs:", smry)))
  expect_false(any(grepl(normalizePath(out_dir, mustWork = FALSE), smry, fixed = TRUE)))
  expect_lte(length(smry), 21L)
})

test_that("automatic package coverage status is a non-claim", {
  coverage <- qtl2mr:::qtl2mr_validation_tests_table()
  expect_true(all(coverage$status == "not_run_in_analysis"))
  expect_false(any(coverage$status == "pass"))
  legacy_request <- qtl2mr:::qtl2mr_validation_tests_table(check_status = "pass")
  expect_true(all(legacy_request$status == "not_run_in_analysis"))
})

test_that("lite export excludes large RDS files and records decisions", {
  result <- make_fake_qtl2mr_analysis_result()
  out_dir <- tempfile("qtl2mr-lite-src-")
  written <- write_qtl2mr_results(result, output_dir = out_dir, make_html = TRUE, make_figures = TRUE)
  zip_path <- file.path(tempdir(), paste0("qtl2mr-lite-", basename(out_dir), ".zip"))
  lite <- export_qtl2mr_lite_result(written, dest_zip = zip_path, max_file_mb = 0.001, include_rds = FALSE)
  expect_true(file.exists(lite))
  listing <- utils::unzip(lite, list = TRUE)
  expect_true(any(grepl("LITE_MANIFEST.csv$", listing$Name)))
  expect_true(any(grepl("SMOKE_RESULT_SIZE_AUDIT.csv$", listing$Name)))
  expect_false(any(grepl("[.]rds$", listing$Name)))
})

test_that("lite export retains portable report and summary paths", {
  result <- make_fake_qtl2mr_analysis_result()
  out_dir <- tempfile("qtl2mr-lite-portable-src-")
  written <- write_qtl2mr_results(result, output_dir = out_dir, make_html = TRUE, make_figures = FALSE)
  zip_path <- file.path(tempdir(), paste0("qtl2mr-lite-portable-", basename(out_dir), ".zip"))
  lite <- export_qtl2mr_lite_result(written, dest_zip = zip_path, max_file_mb = 10, include_rds = FALSE)
  listing <- utils::unzip(lite, list = TRUE)
  expect_true(any(grepl("report.html$", listing$Name)))
  expect_true(any(grepl("summary.md$", listing$Name)))
  expect_false(any(grepl("^[A-Za-z]:[/\\\\]", listing$Name)))
})
