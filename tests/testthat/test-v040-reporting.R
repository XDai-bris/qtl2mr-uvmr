test_that("unknown analysis options fail before creating outputs", {
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  d <- qtl2mr_example
  input <- as_qtl2mr_prepared_input(d$allele_probs, d$pheno,
    marker_map = d$marker_map, selected_loci = d$selected_loci,
    exposure = "x", outcome = "y")
  dest <- tempfile("qtl2mr-typo-")
  expect_error(run_qtl2mr_analysis(input, output_dir = dest,
    standardize_phenotpyes = FALSE), "Unknown")
  expect_false(file.exists(dest))
})

test_that("new result exports retain confidence components and interpretation limits", {
  result <- make_fake_qtl2mr_analysis_result()
  out <- tempfile("qtl2mr-v040-report-")
  written <- write_qtl2mr_results(result, output_dir = out,
    make_html = TRUE, make_figures = FALSE)
  expect_true(file.exists(file.path(out, "data", "confidence_set.csv")))
  expect_true(file.exists(file.path(out, "data", "sample_alignment.csv")))
  expect_true(file.exists(file.path(out, "data", "phenotype_scaling.csv")))
  html <- paste(readLines(file.path(out, "report.html")), collapse = "\n")
  expect_match(html, "Numerical success does not establish 95% coverage", fixed = TRUE)
  expect_match(html, "Effect units", fixed = TRUE)
  expect_false(grepl("Reliable reported estimate", html, fixed = TRUE))
})
