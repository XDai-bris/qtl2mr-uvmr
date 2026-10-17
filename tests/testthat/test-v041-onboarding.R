onboarding_vignette_source <- function(name) {
  path <- test_path("..", "..", "vignettes", paste0(name, ".Rmd"))
  if (!file.exists(path)) path <- system.file("doc", paste0(name, ".Rmd"), package = "qtl2mr")
  if (!nzchar(path) || !file.exists(path)) skip("vignette source is unavailable")
  readLines(path, warn = FALSE)
}

test_that("the documented native qtl2 bridge executes without downloaded data", {
  skip_if_not_installed("qtl2")
  source <- onboarding_vignette_source("from-qtl2")
  start <- which(source == "```{r native-qtl2-bridge}")
  expect_length(start, 1L)
  end <- which(seq_along(source) > start & source == "```")[[1L]]
  env <- new.env(parent = globalenv())
  eval(parse(text = source[seq.int(start + 1L, end - 1L)]), envir = env)

  expect_true(qtl2::check_cross2(env$cross_demo))
  expect_s3_class(env$bridge_genoprobs, "calc_genoprob")
  expect_identical(dim(env$bridge_genoprobs[["1"]])[2L], 36L)
  expect_false(attr(env$bridge_genoprobs, "alleleprobs"))
  expect_s3_class(env$bridge_probs, "calc_genoprob")
  expect_true(attr(env$bridge_probs, "alleleprobs"))
  expect_identical(attr(env$bridge_probs, "crosstype"), "do")
  expect_identical(attr(env$bridge_probs, "alleles"), LETTERS[1:8])
  expect_identical(names(env$bridge_probs), c("1", "2"))
  expect_identical(attr(env$bridge_probs, "is_x_chr"), c("1" = FALSE, "2" = FALSE))
  expect_identical(names(env$bridge_kinship), c("1", "2"))
  for (chr in names(env$bridge_probs)) {
    p <- env$bridge_probs[[chr]]
    expect_identical(dim(p), c(64L, 8L, 23L))
    expect_identical(dimnames(p)[[1L]], env$ids)
    expect_identical(dimnames(p)[[3L]], names(env$bridge_map[[chr]]))
    expect_equal(unname(apply(p, c(1L, 3L), sum)), matrix(1, 64L, 23L),
                 tolerance = 1e-12)
    K <- env$bridge_kinship[[chr]]
    expect_identical(dimnames(K), list(env$ids, env$ids))
    expect_equal(K, t(K), tolerance = 1e-12)
    expect_true(all(is.finite(K)))
  }
  # Native single-index subsetting selects individuals, not chromosomes.
  expect_error(suppressWarnings(env$bridge_probs[env$autosomes]))
  expect_s3_class(env$bridge_scan, "scan1")
  expect_identical(colnames(env$bridge_scan), c("x", "y"))
  expect_true(all(is.finite(env$bridge_scan)))
  expect_true(nrow(env$bridge_selected) > 0L)
  expect_true(all(env$bridge_selected$instrument_trait == "x"))
  expect_true(all(env$bridge_selected$marker %in% rownames(env$bridge_scan)))
  expect_true(all(env$bridge_selected$chr %in% c("1", "2")))
  expect_true(validate_qtl2mr_input(env$bridge_prepared, require_selected_loci = TRUE))
  expect_identical(env$bridge_prepared$allele_probs, env$bridge_probs)
})

test_that("onboarding documents sample, units and provenance boundaries", {
  upstream <- paste(onboarding_vignette_source("from-qtl2"), collapse = "\n")
  expect_match(upstream, "allele_probs <- allele_probs[, autosomes]", fixed = TRUE)
  expect_false(grepl("allele_probs <- allele_probs[autosomes]", upstream, fixed = TRUE))
  expect_match(upstream, "role-specific", fixed = TRUE)
  expect_match(upstream, "common complete", fixed = TRUE)
  expect_match(upstream, "identical sample masks", fixed = TRUE)
  getting <- paste(onboarding_vignette_source("getting-started"), collapse = "\n")
  expect_match(getting, "true_theta = 0.75", fixed = TRUE)
  expect_match(getting, "hsq = 0", fixed = TRUE)
  expect_match(getting, "Each draw refits MR", fixed = TRUE)
  multi <- paste(onboarding_vignette_source("multilocus-example"), collapse = "\n")
  expect_match(multi, "not distributed here", fixed = TRUE)
  expect_match(multi, "master seed alone is not a", fixed = TRUE)
  expect_match(multi, "standardized truth = raw truth", fixed = TRUE)
})

test_that("README quick start keeps the documented explicit fast-demo settings", {
  path <- test_path("..", "..", "README.md")
  if (!file.exists(path)) skip("README is unavailable in this installation")
  source <- readLines(path, warn = FALSE)
  start <- which(source == "result <- run_qtl2mr_analysis(")[[1L]]
  end <- which(seq_along(source) > start & source == ")")[[1L]]
  call <- parse(text = source[seq.int(start, end)])[[1L]][[3L]]
  fast <- as.list(call)[-1L]
  settings <- c("adapter", "loci", "overwrite", "make_figures", "make_html",
                "verbose", "run_qtl2_validation", "run_heteroskedasticity",
                "diagnostic_bootstrap_B", "heterogeneity_method", "run_loo_mr",
                "run_instrument_sensitivity", "progress")
  getting <- onboarding_vignette_source("getting-started")
  start <- which(getting == "result <- run_qtl2mr_analysis(")[[1L]]
  end <- which(seq_along(getting) > start & getting == ")")[[1L]]
  call <- parse(text = getting[seq.int(start, end)])[[1L]][[3L]]
  expect_identical(fast[settings], as.list(call)[-1L][settings])
  text <- paste(source, collapse = "\n")
  expect_match(text, "unreleased v0.4.1", fixed = TRUE)
  expect_match(text, "releases/tag/v0.4.0", fixed = TRUE)
  expect_match(text, "pair-common-complete", fixed = TRUE)
})

test_that("the historical extractor fails clearly without its undistributed input", {
  path <- test_path("..", "..", "data-raw", "make_qtl2mr_multilocus_example.R")
  if (!file.exists(path)) skip("data-raw is intentionally excluded from the package")
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  output <- tempfile(fileext = ".rda")
  expect_error(env$extract_qtl2mr_multilocus_example(tempfile(), output),
               "not distributed in this repository", fixed = TRUE)
  expect_false(file.exists(output))
})
