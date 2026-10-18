test_that("qtl2mr_example is a documented deterministic prepared-input fixture", {
  example_env <- new.env(parent = emptyenv())
  data_index <- utils::data(package = "qtl2mr")
  expect_true("qtl2mr_example" %in% data_index$results[, "Item"])
  utils::data("qtl2mr_example", package = "qtl2mr", envir = example_env)
  example <- example_env$qtl2mr_example

  expect_identical(names(example), c(
    "allele_probs", "pheno", "covar", "kinship", "marker_map",
    "selected_loci", "exposure", "outcome", "true_theta", "metadata"
  ))
  expect_s3_class(example$allele_probs, "calc_genoprob")
  expect_equal(names(example$allele_probs), "1")
  expect_equal(dim(example$allele_probs[["1"]]), c(96L, 8L, 4L))
  expect_equal(dimnames(example$allele_probs[["1"]])[[2L]], LETTERS[1:8])
  expect_equal(dimnames(example$allele_probs[["1"]])[[3L]], paste0("m", 1:4))
  expect_equal(unname(rowSums(example$allele_probs[["1"]][, , "m1"])), rep(1, 96), tolerance = 1e-12)
  expect_equal(rownames(example$pheno), rownames(example$covar))
  expect_identical(rownames(example$pheno), paste0("id", seq_len(96L)))
  expect_equal(rownames(example$pheno), dimnames(example$allele_probs[["1"]])[[1L]])
  expect_null(example$kinship)
  expect_equal(example$selected_loci, c("m1", "m2"))
  expect_equal(example$true_theta, 0.75)
  expect_match(example$metadata$exposure_lod_note, "illustrative", ignore.case = TRUE)
})

test_that("qtl2mr_example supports the public quick workflow and lite export", {
  example_env <- new.env(parent = emptyenv())
  utils::data("qtl2mr_example", package = "qtl2mr", envir = example_env)
  example <- example_env$qtl2mr_example
  prepared <- as_qtl2mr_prepared_input(
    allele_probs = example$allele_probs,
    pheno = example$pheno,
    covar = example$covar,
    kinship = example$kinship,
    marker_map = example$marker_map,
    selected_loci = example$selected_loci,
    exposure = example$exposure,
    outcome = example$outcome,
    metadata = example$metadata
  )
  expect_true(validate_qtl2mr_input(prepared, require_selected_loci = TRUE))

  output_dir <- tempfile("qtl2mr-example-data-")
  result <- run_qtl2mr_analysis(
    prepared,
    loci = 2L,
    output_dir = output_dir,
    overwrite = TRUE,
    make_figures = FALSE,
    make_html = TRUE,
    verbose = FALSE,
    diagnostic_bootstrap_B = 0L,
    heterogeneity_method = "chisq",
    run_loo_mr = FALSE,
    run_instrument_sensitivity = FALSE,
    progress = FALSE
  )
  primary <- qtl2mr_primary_result(result)
  expect_s3_class(result, "qtl2mr_analysis_result")
  expect_equal(primary$n_loci[[1L]], 2L)
  expect_equal(primary$reported_estimator_policy[[1L]], "diagnostic_only")
  expect_true(primary$primary_reliable[[1L]])
  expect_equal(primary$reported_estimator[[1L]], "Star-space qtl2mr")
  expect_true(is.finite(primary$theta_reported[[1L]]))
  expect_true(is.finite(primary$se_reported[[1L]]))

  expected_files <- c("report.html", "summary.md", file.path("data", "main_result.csv"), file.path("data", "output_manifest.csv"))
  expect_true(all(file.exists(file.path(output_dir, expected_files))))
  manifest <- utils::read.csv(file.path(output_dir, "data", "output_manifest.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("report.html", "summary.md") %in% manifest$path_relative))
  expect_false(any(grepl("^/|^[A-Za-z]:[/\\\\]", manifest$path_relative)))

  warning_table <- qtl2mr_warnings(result)
  expect_true("weak_block_first_stage" %in% warning_table$warning_id)
  weak <- warning_table[warning_table$warning_id == "weak_block_first_stage", , drop = FALSE]
  expect_false(any(grepl("Stock-Yogo", weak$message, fixed = TRUE)))

  lite_zip <- tempfile("qtl2mr-example-lite-", fileext = ".zip")
  lite <- export_qtl2mr_lite_result(result, dest_zip = lite_zip, max_file_mb = 10, include_rds = FALSE)
  listing <- utils::unzip(lite, list = TRUE)$Name
  expect_true(file.exists(lite))
  expect_true(any(grepl("(^|/)report[.]html$", listing)))
  expect_true(any(grepl("(^|/)summary[.]md$", listing)))
  expect_true(any(grepl("(^|/)LITE_MANIFEST[.]csv$", listing)))
  expect_false(any(grepl("[.]rds$", listing, ignore.case = TRUE)))
})

test_that("getting-started vignette source is generic and free of local paths", {
  source_vignette <- test_path("..", "..", "vignettes", "getting-started.Rmd")
  vignette <- if (file.exists(source_vignette)) {
    source_vignette
  } else {
    system.file("doc", "getting-started.Rmd", package = "qtl2mr")
  }
  if (!nzchar(vignette) || !file.exists(vignette)) {
    skip("getting-started vignette source is unavailable in this source or installed package")
  }
  text <- paste(readLines(vignette, warn = FALSE), collapse = "\n")
  expect_match(text, "qtl2mr_example", fixed = TRUE)
  local_home_pattern <- paste0("/", "Users", "/")
  token_prefix <- paste0("gh", "p_")
  pat_prefix <- paste0("github", "_pat_")
  expect_false(grepl(paste(local_home_pattern, "[A-Za-z]:[/\\\\]", token_prefix, pat_prefix, sep = "|"), text))
  expect_match(text, "diagnostic_only", fixed = TRUE)
  expect_match(text, "LOCO kinship", fixed = TRUE)
  expect_match(text, "as_qtl2mr_prepared_input", fixed = TRUE)
  expect_match(text, "run_qtl2mr_analysis", fixed = TRUE)
  expect_match(text, "primary <- qtl2mr_primary_result(result)", fixed = TRUE)
  expect_match(text, "warnings <- qtl2mr_warnings(result)", fixed = TRUE)
  expect_match(text, "rownames(actionable_warnings) <- NULL", fixed = TRUE)
  expect_match(text, "knitr::kable", fixed = TRUE)
  expect_match(text, "browseURL(file.path(result_dir, \"report.html\"))", fixed = TRUE)
  expect_false(grepl("(^|[^[:alnum:]_])names\\(qtl2mr_example|head\\(|export_qtl2mr_lite_result|utils::unzip|data.frame\\(\\s*component", text))
  expect_match(text, "prepared_lmm", fixed = TRUE)
  expect_match(text, "not a claim of causal validity", fixed = TRUE)
  expect_lte(length(readLines(vignette, warn = FALSE)), 220L)
})

test_that("package source excludes local paths, credentials and extra datasets", {
  source_root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  package_root <- if (file.exists(file.path(source_root, "DESCRIPTION"))) {
    source_root
  } else {
    system.file(package = "qtl2mr")
  }
  expect_true(nzchar(package_root))
  text_extensions <- c("r", "rd", "rmd", "md", "yaml", "yml", "json", "txt", "html")
  candidates <- list.files(package_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  ignore_file <- file.path(package_root, ".Rbuildignore")
  if (file.exists(ignore_file)) {
    rel <- substring(candidates, nchar(package_root) + 2L)
    ignored <- readLines(ignore_file, warn = FALSE)
    excluded <- vapply(rel, function(path) {
      # Directory exclusions also apply to members underneath that directory.
      parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
      ancestors <- vapply(seq_along(parts), function(i) paste(parts[seq_len(i)], collapse = "/"), character(1))
      any(vapply(ignored, function(rx) any(grepl(rx, ancestors, perl = TRUE, ignore.case = TRUE)), logical(1)))
    }, logical(1))
    candidates <- candidates[!excluded]
  }
  keep <- tolower(tools::file_ext(candidates)) %in% text_extensions |
    basename(candidates) %in% c("DESCRIPTION", "NAMESPACE", "NEWS.md", "README.md")
  text <- paste(vapply(candidates[keep], function(path) {
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }, character(1)), collapse = "\n")
  local_path <- paste0("/", "Users", "/")
  token_prefix <- paste0("gh", "p_")
  pat_prefix <- paste0("github", "_pat_")
  forbidden <- c(local_path, token_prefix, pat_prefix)
  expect_false(any(vapply(forbidden, grepl, logical(1), x = text, ignore.case = TRUE)))
  data_dir <- file.path(package_root, "data")
  if (dir.exists(data_dir)) {
    expect_setequal(list.files(data_dir, pattern = "[.]rda$", ignore.case = TRUE),
      c("qtl2mr_example.rda", "qtl2mr_multilocus_example.rda"))
  }
})

test_that("built installation exposes the getting-started vignette", {
  rendered <- system.file("doc", "getting-started.html", package = "qtl2mr")
  if (!nzchar(rendered)) skip("rendered vignette is created during package build")
  expect_true(file.exists(rendered))
  available <- utils::vignette(package = "qtl2mr")
  expect_true("getting-started" %in% available$results[, "Item"])
})
