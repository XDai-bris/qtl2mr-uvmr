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
  selected <- env$bridge_instruments[env$bridge_instruments$selected, , drop = FALSE]
  expect_true(nrow(selected) > 0L)
  expect_true(all(selected$exposure == "x"))
  expect_true(all(selected$marker %in% rownames(env$bridge_scan)))
  expect_true(all(selected$chr %in% c("1", "2")))
  expect_true(all(is.finite(selected$F_X_block)))
  expect_identical(env$bridge_prepared$selected_loci$marker, selected$marker)
  expect_true(validate_qtl2mr_input(env$bridge_prepared, require_selected_loci = TRUE))
  expect_identical(env$bridge_prepared$allele_probs, env$bridge_probs)
})

test_that("onboarding documents sample, units and provenance boundaries", {
  upstream <- paste(onboarding_vignette_source("from-qtl2"), collapse = "\n")
  expect_match(upstream, "allele_probs <- allele_probs[, autosomes]", fixed = TRUE)
  expect_false(grepl("allele_probs <- allele_probs[autosomes]", upstream, fixed = TRUE))
  expect_match(upstream, "role-specific", fixed = TRUE)
  expect_match(upstream, "A batch row is identical to the single analysis", fixed = TRUE)
  expect_match(upstream, "identical samples", fixed = TRUE)
  getting <- paste(onboarding_vignette_source("getting-started"), collapse = "\n")
  expect_match(getting, "true_theta = 0.75", fixed = TRUE)
  expect_match(getting, "hsq = 0", fixed = TRUE)
  expect_match(getting, "Each draw refits MR", fixed = TRUE)
  multi <- paste(onboarding_vignette_source("multilocus-example"), collapse = "\n")
  expect_match(multi, "not distributed here", fixed = TRUE)
  expect_match(multi, "master seed alone is not a", fixed = TRUE)
  expect_match(multi, "standardized truth = raw truth", fixed = TRUE)
})

# The R code chunks of a Markdown or R Markdown file.
onboarding_code_chunks <- function(lines) {
  starts <- grep("^```(r|\\{r)", lines)
  lapply(starts, function(s) {
    e <- which(seq_along(lines) > s & lines == "```")[[1L]]
    lines[seq.int(s + 1L, e - 1L)]
  })
}

# Named arguments of calls to the entry functions that the function does not have.
onboarding_unknown_arguments <- function(expr) {
  entries <- c("run_qtl2mr_analysis", "run_qtl2mr_batch", "qtl2mr_select_instruments", "qtl2mr_pair_input",
    "qtl2mr_prepare_resources", "qtl2mr_combine_resources", "qtl2mr_diagnostics")
  if (!is.call(expr)) return(character())
  fn <- expr[[1L]]
  name <- if (is.name(fn)) as.character(fn) else if (is.call(fn) && identical(fn[[1L]], as.name("::"))) {
    as.character(fn[[3L]])
  } else ""
  bad <- character()
  if (name %in% entries) {
    used <- names(as.list(expr)[-1L])
    used <- used[nzchar(used)]
    bad <- paste0(name, "(", setdiff(used, names(formals(get(name, asNamespace("qtl2mr"))))), ")")
    bad <- bad[!endsWith(bad, "()")]
  }
  for (arg in as.list(expr)[-1L]) if (!missing(arg)) bad <- c(bad, onboarding_unknown_arguments(arg))
  bad
}

test_that("README and vignette code pass only arguments the 0.6.1.9006 entries have", {
  # Before 0.6.1.9006 this test required the README and getting-started quick
  # starts to pass the same fast-demo settings, most of which (verbose, the Q
  # and bootstrap settings, run_loo_mr, run_instrument_sensitivity) were removed.
  path <- test_path("..", "..", "README.md")
  if (!file.exists(path)) skip("README is unavailable in this installation")
  source <- readLines(path, warn = FALSE)
  docs <- list(README = source, `getting-started` = onboarding_vignette_source("getting-started"),
    `multilocus-example` = onboarding_vignette_source("multilocus-example"),
    `from-qtl2` = onboarding_vignette_source("from-qtl2"))
  for (doc in names(docs)) {
    chunks <- onboarding_code_chunks(docs[[doc]])
    expect_gt(length(chunks), 0L)
    for (chunk in chunks) {
      exprs <- parse(text = chunk)
      expect_identical(as.character(unlist(lapply(exprs, onboarding_unknown_arguments))), character(), info = doc)
    }
  }
  expect_match(paste(source, collapse = "\n"), "result <- run_qtl2mr_analysis(", fixed = TRUE)
  text <- paste(source, collapse = "\n")
  expect_match(text, "**Current release:** [v0.7.0]", fixed = TRUE)
  expect_match(text, "releases/tag/v0.7.0", fixed = TRUE)
  expect_match(text, "exposure and the outcome on their own complete", fixed = TRUE)
  # Retired names may appear only in the change and upgrade sections.
  current <- sub("## (Changes since|Upgrading from) .*", "", text)
  for (retired in c("run_qtl2mr_batch_fast", "sample_policy", "instrument_rule", "qtl_scan =", "select_loci")) {
    expect_false(grepl(retired, current, fixed = TRUE), info = retired)
  }
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
