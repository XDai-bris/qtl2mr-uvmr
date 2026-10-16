load_multilocus_fixture <- function() {
  e <- new.env(parent = emptyenv())
  utils::data("qtl2mr_multilocus_example", package = "qtl2mr", envir = e)
  e$qtl2mr_multilocus_example
}

test_that("multi-locus fixture contains only the documented synthetic inputs", {
  d <- load_multilocus_fixture()
  expect_named(d, c("allele_probs", "pheno", "covar", "kinship", "marker_map",
    "selected_loci", "exposure", "outcome", "true_theta", "metadata"))
  expect_identical(names(d$pheno), c("x1", "y"))
  expect_identical(dim(d$pheno), c(400L, 2L))
  expect_identical(names(d$covar), c("batch", "age"))
  expect_identical(rownames(d$pheno), sprintf("synthetic_mouse_%03d", 1:400))
  expect_identical(rownames(d$pheno), rownames(d$covar))
  expect_identical(dimnames(d$kinship), list(rownames(d$pheno), rownames(d$pheno)))
  expect_identical(names(d$marker_map), c("marker", "chr", "pos"))
  expect_identical(d$selected_loci, sprintf("synthetic_locus_%03d", 1:40))
  expect_identical(d$selected_loci, d$marker_map$marker)
  expect_identical(d$exposure, "x1")
  expect_identical(d$outcome, "y")
  expect_equal(d$true_theta, -0.25)
  expect_false(d$metadata$standardize_phenotypes)
  expect_match(d$metadata$purpose, "not calibration", fixed = TRUE)
  expect_match(d$metadata$kinship_note, "not LOCO", fixed = TRUE)
  expect_match(d$metadata$selection_note, "no QTL scan", fixed = TRUE)
  expect_match(d$metadata$source_sha256, "^[a-f0-9]{64}$")

  expect_s3_class(d$allele_probs, "calc_genoprob")
  expect_setequal(names(d$allele_probs), as.character(1:19))
  expect_equal(sum(vapply(d$allele_probs, function(a) dim(a)[3L], integer(1))), 40L)
  for (chr in names(d$allele_probs)) {
    a <- d$allele_probs[[chr]]
    expect_identical(dim(a)[1:2], c(400L, 8L))
    expect_identical(dimnames(a)[[1L]], rownames(d$pheno))
    expect_identical(dimnames(a)[[2L]], LETTERS[1:8])
    expect_identical(dimnames(a)[[3L]], d$marker_map$marker[d$marker_map$chr == chr])
    expect_true(all(is.finite(a) & a >= 0 & a <= 1))
    expect_equal(unname(apply(a, c(1L, 3L), sum)),
      matrix(1, 400L, dim(a)[3L]), tolerance = 1e-12)
  }
  family <- outer(1:4, 1:4, function(i, j) 0.3^abs(i - j))
  expect_equal(unname(d$kinship), kronecker(diag(100), family) / 2, tolerance = 1e-14)
  expect_true(all(is.finite(as.matrix(d$pheno))))
  expect_true(all(is.finite(as.matrix(d$covar))))
  text <- paste(capture.output(dput(d$metadata)), collapse = "\n")
  home_prefix <- paste0("/", "Users", "/")
  expect_false(grepl(home_prefix, text, fixed = TRUE))
  expect_false(grepl("[.]rds|[.]rda|[.]git|workshop/", text))
})

test_that("multi-locus example runs the public raw-unit workflow", {
  d <- load_multilocus_fixture()
  input <- as_qtl2mr_prepared_input(
    d$allele_probs, d$pheno, covar = d$covar, kinship = d$kinship,
    marker_map = d$marker_map, selected_loci = d$selected_loci,
    exposure = d$exposure, outcome = d$outcome)
  expect_true(validate_qtl2mr_input(input, require_selected_loci = TRUE))
  result <- run_qtl2mr_analysis(input, output_dir = tempfile("qtl2mr-multilocus-test-"),
    standardize_phenotypes = FALSE, make_figures = FALSE, make_html = TRUE,
    diagnostic_bootstrap_B = 500L, diagnostic_bootstrap_seed = 1L,
    heterogeneity_method = "chisq",
    run_loo_mr = FALSE, run_instrument_sensitivity = FALSE,
    verbose = FALSE, progress = FALSE)
  p <- qtl2mr_primary_result(result)
  expect_equal(p$n_loci, 40L)
  expect_equal(p$n_x_median, 400)
  expect_equal(p$n_y_median, 400)
  expect_equal(p$n_overlap_median, 400)
  expect_identical(p$effect_units, "outcome units per exposure unit")
  expect_identical(p$confidence_set_type, "bounded")
  expect_true(p$numerical_ok)
  expect_identical(p$inference_status, "not_calibrated")
  # Fixed-fixture regression, not a repeated-sampling coverage test.
  expect_equal(p$theta_reported, -0.2750267, tolerance = 1e-3)
  expect_equal(p$reported_ci_lower, -0.5251423, tolerance = 1e-3)
  expect_equal(p$reported_ci_upper, -0.02953669, tolerance = 1e-3)
  expect_equal(unname(coef(lm(y ~ x1 + batch + age, data = cbind(d$pheno, d$covar)))["x1"]),
    0.05, tolerance = 1e-10)
  expect_identical(result$qtl2_validation_status$validation_status, "not_requested")
  methods <- qtl2mr_method_comparison(result)
  expect_identical(methods$method, c("Star-space qtl2MR", "Sigma_xy=0 qtl2MR",
    "Two-step IVW", "Weighted median", "Weighted mode"))
  expect_true(all(is.finite(methods$ci_lower) & is.finite(methods$ci_upper)))
  expect_true(all(methods$ci_lower < methods$ci_upper))
  expect_identical(methods$ci_type, c("profile", "Wald", "Wald",
    "bootstrap-SE Wald", "bootstrap-SE Wald"))
  boot <- result$pilot$diagnostic_layer$estimates
  boot <- boot[boot$estimator %in% c("Weighted median", "Weighted mode"), ]
  expect_true(all(boot$bootstrap_B == 500L & boot$bootstrap_seed == 1L))
  expect_true(all(boot$bootstrap_complete_B == 500L & boot$bootstrap_failed_B == 0L))
  expect_equal(boot$ci_lower, boot$theta - qnorm(.975) * boot$se)
  expect_equal(boot$ci_upper, boot$theta + qnorm(.975) * boot$se)
  expect_equal(methods$estimate[methods$method == "Weighted mode"], -0.2723459669, tolerance = 1e-6)
  exported <- utils::read.csv(file.path(result$output_dir, "data", "method_comparison.csv"))
  expect_equal(exported[, c("method", "ci_lower", "ci_upper", "ci_type")],
    methods[, c("method", "ci_lower", "ci_upper", "ci_type")], ignore_attr = TRUE)
  expect_true(file.exists(file.path(result$output_dir, "report.html")))
  html <- paste(readLines(file.path(result$output_dir, "report.html")), collapse = "\n")
  expect_match(html, "parametric-bootstrap SD", fixed = TRUE)
  expect_match(html, "not profile-likelihood intervals", fixed = TRUE)
  expect_match(html, "CI method", fixed = TRUE)
})

test_that("multi-locus walkthrough is available after vignette building", {
  src <- test_path("..", "..", "vignettes", "multilocus-example.Rmd")
  if (!file.exists(src)) src <- system.file("doc", "multilocus-example.Rmd", package = "qtl2mr")
  if (!nzchar(src) || !file.exists(src)) skip("vignette source not available")
  text <- paste(readLines(src, warn = FALSE), collapse = "\n")
  expect_match(text, "standardize_phenotypes = FALSE", fixed = TRUE)
  expect_match(text, "run_loo_mr = TRUE", fixed = TRUE)
  expect_match(text, "diagnostic_bootstrap_B = 500L", fixed = TRUE)
  expect_match(text, "diagnostic_bootstrap_seed = 1L", fixed = TRUE)
  expect_match(text, "deliberately targeted", fixed = TRUE)
  expect_match(text, "not method", fixed = TRUE)
  expect_lte(length(readLines(src, warn = FALSE)), 180L)
  rendered <- system.file("doc", "multilocus-example.html", package = "qtl2mr")
  if (!nzchar(rendered)) skip("rendered vignette is created during package build")
  expect_true(file.exists(rendered))
  expect_true("multilocus-example" %in% utils::vignette(package = "qtl2mr")$results[, "Item"])
})
