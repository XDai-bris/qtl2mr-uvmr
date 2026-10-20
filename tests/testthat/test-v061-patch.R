v061_fixture <- function() {
  set.seed(61061)
  n <- 90L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2")), `2` = make_chr(c("b1", "b2")))
  G <- matrix(rnorm(n * 15L), n)
  K <- tcrossprod(G) / 15 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  poly <- as.numeric(t(chol(K)) %*% rnorm(n))
  x <- 3 * gp[["1"]][, 1L, "a1"] + 2 * gp[["2"]][, 2L, "b2"] + 0.4 * poly + rnorm(n, sd = 0.3)
  pheno <- data.frame(x = x, y = 0.6 * x + 2 * gp[["1"]][, 4L, "a2"] + 0.4 * poly + rnorm(n, sd = 0.3),
    row.names = ids)
  pheno$y[1:5] <- NA_real_
  list(gp = gp, kinship = list(`1` = K, `2` = K), pheno = pheno,
    covar = data.frame(sex = rep(0:1, n / 2L), row.names = ids),
    map = data.frame(marker = c("a1", "a2", "b1", "b2"), chr = c("1", "1", "2", "2"), pos = c(1, 2, 1, 2)))
}

test_that("raw-unit results do not depend on the units of the traits", {
  f <- v061_fixture()
  fit <- function(sx, sy) {
    pheno <- f$pheno
    pheno$x <- pheno$x * sx
    pheno$y <- pheno$y * sy
    input <- as_qtl2mr_prepared_input(f$gp, pheno, f$covar, f$kinship, f$map, c("a1", "a2", "b2"), "x", "y")
    qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE, standardize_phenotypes = FALSE,
      heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L)
  }
  base <- fit(1, 1)
  expect_true(all(base$hsq$hsq > 0 & base$hsq$hsq < 1))
  expect_true(any(!base$locus_first_stage$steiger_pass))
  for (s in list(c(1e-8, 1), c(1, 1e-8), c(1e-10, 1e-10), c(1e6, 1e-6))) {
    r <- fit(s[[1L]], s[[2L]])
    ratio <- s[[1L]] / s[[2L]]
    # Point estimate and limits agree to the optimiser's precision.
    expect_equal(r$summary$theta * ratio, base$summary$theta, tolerance = 1e-6, info = paste(s, collapse = "/"))
    expect_equal(r$summary$profile_ci_lower * ratio, base$summary$profile_ci_lower, tolerance = 1e-6)
    expect_equal(r$summary$Q, base$summary$Q, tolerance = 1e-7)
    expect_equal(r$locus_first_stage$F_X_block, base$locus_first_stage$F_X_block, tolerance = 1e-7)
    expect_equal(r$locus_first_stage$F_Y_block, base$locus_first_stage$F_Y_block, tolerance = 1e-7)
    expect_identical(r$locus_first_stage$steiger_pass, base$locus_first_stage$steiger_pass)
    expect_equal(r$hsq$hsq, base$hsq$hsq, tolerance = 1e-6)
  }

  # Selection standardizes internally; a raw-unit batch reports the same F.
  pheno <- f$pheno
  pheno$x <- pheno$x * 1e-8
  resources <- qtl2mr_prepare_resources(f$gp, f$kinship, pheno, f$covar, f$map)
  inst <- data.frame(exposure = "x", marker = c("a1", "b2"))
  batch <- run_qtl2mr_batch(resources, data.frame(exposure = "x", outcome = "y"), inst,
    standardize_phenotypes = FALSE, progress = FALSE)
  scan <- matrix(c(6, 1, 1, 5), ncol = 1L, dimnames = list(f$map$marker, "x"))
  selected <- qtl2mr_select_instruments(resources, scan, lod_thresh = 3, progress = FALSE)
  expect_equal(batch$locus_results$F_X_block, selected$F_X_block[match(c("a1", "b2"), selected$marker)],
    tolerance = 1e-8)
})

test_that("allele probabilities must hold the same samples on every chromosome", {
  f <- v061_fixture()
  for (drop_chr in c("1", "2")) {
    gp <- f$gp
    gp[[drop_chr]] <- gp[[drop_chr]][-1L, , , drop = FALSE]
    expect_error(qtl2mr_prepare_resources(gp, NULL, f$pheno, f$covar, f$map),
      "chromosomes must hold the same samples")
    expect_error(as_qtl2mr_prepared_input(gp, f$pheno, f$covar, NULL, f$map, "b2", "x", "y"),
      "chromosomes must hold the same samples")
  }
  reordered <- f$gp
  reordered[["2"]] <- reordered[["2"]][rev(seq_len(dim(reordered[["2"]])[1L])), , , drop = FALSE]
  expect_s3_class(qtl2mr_prepare_resources(reordered, NULL, f$pheno, f$covar, f$map), "qtl2mr_resources")
})

test_that("a failed Steiger refit is logged and warned", {
  f <- v061_fixture()
  input <- as_qtl2mr_prepared_input(f$gp, f$pheno, f$covar, f$kinship, f$map, c("a1", "a2", "b2"), "x", "y")
  local_mocked_bindings(refit_qtl2mr_subset = function(...) stop("refit boom"))
  result <- run_qtl2mr_analysis(input, output_dir = tempfile("v061-steiger-"), make_html = FALSE,
    make_figures = FALSE, verbose = FALSE, progress = FALSE, heterogeneity_method = "chisq",
    diagnostic_bootstrap_B = 0L, run_loo_mr = FALSE)
  expect_identical(result$steiger_sensitivity$status, "failed")
  expect_true("steiger_sensitivity_error" %in% qtl2mr_warnings(result)$warning_id)
  expect_identical(readLines(file.path(result$output_dir, "logs", "steiger_sensitivity_error.log")), "refit boom")
  expect_true(is.finite(qtl2mr_primary_result(result)$theta))
})
