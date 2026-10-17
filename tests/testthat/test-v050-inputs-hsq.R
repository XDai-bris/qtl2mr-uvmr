v050_input <- function(kinship = TRUE) {
  set.seed(50017)
  n <- 60L
  ids <- paste0("s", seq_len(n))
  a <- array(rexp(n * 8L * 2L), c(n, 8L, 2L), dimnames = list(ids, LETTERS[1:8], c("m1", "m2")))
  for (m in seq_len(2L)) a[, , m] <- a[, , m] / rowSums(a[, , m])
  z <- a[, , 1L] %*% make_founder_contrasts(LETTERS[1:8])
  x <- 3 * z[, 1L] + rnorm(n, sd = .2)
  pheno <- data.frame(x = x, y = .7 * x + rnorm(n, sd = .2), row.names = ids)
  G <- matrix(rnorm(n * 30L), n)
  K <- tcrossprod(G) / 60 + diag(.5, n)
  dimnames(K) <- list(ids, ids)
  list(allele_probs = list(`1` = a), pheno = pheno,
    covar = data.frame(age = seq_len(n), row.names = ids),
    kinship = if (kinship) list(`1` = K) else NULL,
    marker_map = data.frame(marker = c("m1", "m2"), chr = "1", pos = 1:2),
    selected_loci = c("m1", "m2"), exposure = "x", outcome = "y")
}

test_that("exposure and outcome must be two different named traits", {
  raw <- v050_input()
  raw$outcome <- "x"
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "must be different traits")
  raw <- v050_input()
  raw$exposure <- c("x", "y")
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "exposure must be one non-empty trait name")
  input <- do.call(as_qtl2mr_prepared_input, v050_input())
  expect_error(run_qtl2mr_analysis(input, exposure = "y", outcome = "y", output_dir = tempfile(), verbose = FALSE, progress = FALSE),
    "must be different traits")
})

test_that("an empty selected-locus set is rejected with a direct message", {
  for (empty in list(character(), data.frame(marker = character(), chr = character()))) {
    raw <- v050_input()
    raw$selected_loci <- empty
    expect_error(do.call(as_qtl2mr_prepared_input, raw), "selected_loci contains no loci")
  }
})

test_that("null-model hsq is reported per role and chromosome and appended to outputs", {
  input <- do.call(as_qtl2mr_prepared_input, v050_input())
  fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "chisq")
  hsq <- fit$hsq
  expect_identical(hsq$role, c("exposure", "outcome"))
  expect_identical(hsq$chr, c("1", "1"))
  expect_identical(hsq$n_loci, c(2L, 2L))
  expect_identical(hsq$markers, c("m1;m2", "m1;m2"))
  expect_true(all(hsq$hsq_source == "reml_null_model"))
  expect_equal(hsq$hsq, c(fit$locus_blocks[[1]]$fit_x$hsq, fit$locus_blocks[[1]]$fit_y$hsq))
  expect_identical(utils::tail(names(fit$summary), 1L), "sample_policy")

  user <- v050_input()
  user$hsq <- 0.3
  fixed <- qtl2mr:::run_prepared_qtl2mr(do.call(as_qtl2mr_prepared_input, user), "x", "y",
    verbose = FALSE, diagnostic_bootstrap_B = 0L, heterogeneity_method = "chisq")
  expect_identical(unique(fixed$hsq$hsq_source), "user_supplied")
  expect_equal(unique(fixed$hsq$hsq), 0.3)

  ols <- qtl2mr:::run_prepared_qtl2mr(do.call(as_qtl2mr_prepared_input, v050_input(kinship = FALSE)),
    "x", "y", verbose = FALSE, diagnostic_bootstrap_B = 0L, heterogeneity_method = "chisq")
  expect_identical(unique(ols$hsq$hsq_source), "no_kinship_ols")
  warnings <- qtl2mr:::qtl2mr_warning_table(ols)
  expect_true("no_kinship_ols" %in% warnings$warning_id)
})

test_that("boundary hsq estimates are flagged and written with the analysis", {
  raw <- v050_input()
  raw$hsq <- list(exposure = 0)
  input <- do.call(as_qtl2mr_prepared_input, raw)
  out <- tempfile("qtl2mr-hsq-")
  result <- run_qtl2mr_analysis(input, output_dir = out, make_figures = FALSE, make_html = TRUE,
    verbose = FALSE, progress = FALSE, heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
    run_loo_mr = FALSE, run_instrument_sensitivity = FALSE)
  written <- utils::read.csv(file.path(out, "data", "hsq.csv"), stringsAsFactors = FALSE)
  expect_identical(written$role, c("exposure", "outcome"))
  expect_identical(written$hsq_source[[1L]], "user_supplied")
  expect_true(any(grepl("^Null-model hsq: ", readLines(file.path(out, "summary.md")))))
  expect_true(any(grepl("Polygenic background", readLines(file.path(out, "report.html")))))

  blocks <- result$pilot$locus_blocks
  blocks[[1]]$fit_y$hsq <- 0
  blocks[[1]]$fit_y$hsq_status <- "boundary_zero"
  flagged <- qtl2mr:::qtl2mr_warning_table(list(summary = result$pilot$summary,
    hsq = qtl2mr:::qtl2mr_hsq_table(blocks, "x", "y")))
  expect_true("hsq_boundary" %in% flagged$warning_id)
  expect_match(flagged$message[flagged$warning_id == "hsq_boundary"], "outcome chr 1")
})

test_that("samples without a covariate row are reported like non-finite covariates", {
  raw <- v050_input()
  raw$covar <- raw$covar[-3L, , drop = FALSE]
  input <- suppressWarnings(do.call(as_qtl2mr_prepared_input, raw))
  expect_warning(fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none"), "Missing covariate rows")
  expect_true(all(fit$sample_alignment$n_missing_covar_id == 1L))
  expect_true("incomplete_covariates" %in% qtl2mr:::qtl2mr_warning_table(fit)$warning_id)
})

test_that("batch warnings describe the shared batch sample policy", {
  input <- do.call(as_qtl2mr_prepared_input, v050_input())
  fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none")
  single <- qtl2mr:::qtl2mr_warning_table(fit)
  batch <- qtl2mr:::qtl2mr_warning_table(fit, sample_policy = "batch_common_complete")
  expect_true("prepared_overlap_engine" %in% single$warning_id)
  expect_false("prepared_overlap_engine" %in% batch$warning_id)
  expect_match(batch$message[batch$warning_id == "batch_shared_samples"], "batch_common_complete")
})
