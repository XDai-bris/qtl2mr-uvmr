make_prepared_smoke_input <- function(n = 96, n_markers = 4, seed = 301) {
  set.seed(seed)
  ids <- paste0("id", seq_len(n))
  founders <- LETTERS[1:8]
  markers <- paste0("m", seq_len(n_markers))
  arr <- array(NA_real_, c(n, 8, n_markers), dimnames = list(ids, founders, markers))
  for (j in seq_len(n_markers)) {
    raw <- matrix(stats::rexp(n * 8), n, 8)
    arr[, , j] <- raw / rowSums(raw)
  }
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  H <- make_founder_contrasts(founders)
  z <- lapply(seq_len(n_markers), function(j) arr[, , j] %*% H)
  x <- as.numeric(1.4 * z[[1]][, "FC1"] + 1.1 * z[[2]][, "FC2"] +
                    0.9 * z[[3]][, "FC3"] + stats::rnorm(n, sd = 0.25))
  y <- as.numeric(0.75 * x + 0.2 * z[[4]][, "FC4"] + stats::rnorm(n, sd = 0.25))
  pheno <- data.frame(x = x, y = y, row.names = ids)
  covar <- data.frame(
    sex = rep(0:1, length.out = n),
    age = as.numeric(scale(seq_len(n))),
    row.names = ids
  )
  marker_map <- data.frame(marker = markers, chr = "1", pos = seq_len(n_markers), stringsAsFactors = FALSE)
  as_qtl2mr_prepared_input(
    allele_probs = gp,
    pheno = pheno,
    covar = covar,
    kinship = NULL,
    marker_map = marker_map,
    selected_loci = markers,
    exposure = "x",
    outcome = "y",
    metadata = list(example = "synthetic")
  )
}

make_prepared_unreliable_so_input <- function(n = 96L, seed = 947L) {
  set.seed(seed)
  ids <- paste0("id", seq_len(n))
  founders <- LETTERS[1:8]
  markers <- paste0("weak", seq_len(3L))
  arr <- array(
    NA_real_,
    c(n, 8L, 3L),
    dimnames = list(ids, founders, markers)
  )
  for (j in seq_len(3L)) {
    raw <- matrix(stats::rexp(n * 8L), n, 8L)
    arr[, , j] <- raw / rowSums(raw)
  }
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  H <- make_founder_contrasts(founders)
  z <- lapply(seq_len(3L), function(j) arr[, , j] %*% H)
  x <- stats::rnorm(n)
  y <- as.numeric(
    2.5 * z[[1L]][, "FC1"] +
      2.2 * z[[2L]][, "FC2"] +
      2.0 * z[[3L]][, "FC3"] +
      stats::rnorm(n, sd = 0.2)
  )
  pheno <- data.frame(x = x, y = y, row.names = ids)
  covar <- data.frame(
    sex = rep(0:1, length.out = n),
    age = as.numeric(scale(seq_len(n))),
    row.names = ids
  )
  marker_map <- data.frame(
    marker = markers,
    chr = "1",
    pos = seq_len(3L),
    stringsAsFactors = FALSE
  )
  as_qtl2mr_prepared_input(
    allele_probs = gp,
    pheno = pheno,
    covar = covar,
    kinship = NULL,
    marker_map = marker_map,
    selected_loci = markers,
    exposure = "x",
    outcome = "y",
    metadata = list(example = "unreliable_so_public_regression")
  )
}

test_that("prepared input constructor validates dimensions and IDs", {
  input <- make_prepared_smoke_input()
  expect_s3_class(input, "qtl2mr_prepared_input")
  expect_true(validate_qtl2mr_input(input))
  bad <- input
  bad$selected_loci <- c("missing_marker")
  expect_error(validate_qtl2mr_input(bad), "not present in allele_probs")
  bad2 <- input
  bad2$sample_ids <- c("id1", "id1")
  expect_error(validate_qtl2mr_input(bad2), "duplicate")
  bad3 <- input
  bad3$outcome <- "missing_trait"
  expect_error(validate_qtl2mr_input(bad3), "Outcome trait not found")
})

test_that("prepared input validates one union kinship source", {
  input <- make_prepared_smoke_input()
  ids <- qtl2mr:::qtl2mr_prepared_sample_ids(input)
  K <- diag(1, length(ids))
  dimnames(K) <- list(ids, ids)
  input$kinship <- list("2" = K)
  expect_error(validate_qtl2mr_input(input), "missing selected chromosome")
  input$kinship <- list(exposure = K, outcome = K)
  expect_error(validate_qtl2mr_input(input), "one kinship source")
})

test_that("prepared overlap engine matches one-sample block fit under full overlap", {
  input <- make_prepared_smoke_input()
  prepared <- qtl2mr:::run_prepared_qtl2mr(
    input,
    exposure = "x",
    outcome = "y",
    loci = 2,
    write_outputs = FALSE,
    verbose = FALSE,
    heterogeneity_method = "chisq"
  )
  selected <- qtl2mr:::qtl2mr_prepared_selected_loci(input, loci = 2)
  x <- qtl2mr:::qtl2mr_prepared_trait(input, "x")
  y <- qtl2mr:::qtl2mr_prepared_trait(input, "y")
  covar <- qtl2mr:::qtl2mr_prepared_covar(input)
  H <- make_founder_contrasts(dimnames(input$allele_probs[[1]])[[2]])
  blocks <- setNames(lapply(selected$marker, function(marker) {
    fit_mr_locus_block(
      input$allele_probs,
      exposure = x,
      outcome = y,
      addcovar = covar,
      contrasts = H,
      marker = marker,
      standardize_phenotypes = TRUE,
      cross_cov_model = "star"
    )
  }), selected$marker)
  joint <- fit_multilocus_founder_mr(blocks, heterogeneity_method = "chisq")
  expect_equal(prepared$fit$theta, joint$theta, tolerance = 1e-8)
  expect_equal(prepared$fit$se, joint$se, tolerance = 1e-8)
  expect_true(all(prepared$overlap_summary$complete_overlap))
})

test_that("prepared overlap engine handles no shared samples as independent samples", {
  input <- make_prepared_smoke_input(n = 96, n_markers = 4)
  input$selected_loci <- "m1"
  input$pheno$x[49:96] <- NA_real_
  input$pheno$y[1:48] <- NA_real_
  prepared <- qtl2mr:::run_prepared_qtl2mr(
    input,
    exposure = "x",
    outcome = "y",
    loci = 1,
    write_outputs = FALSE,
    verbose = FALSE
  )
  block <- prepared$locus_blocks[[1]]
  expect_equal(block$n_overlap, 0L)
  expect_true(all(block$Sigma_xy == 0))
  expect_true(is.na(block$rho_xy))
  expect_true(prepared$overlap_summary$independent_samples[[1L]])
  expect_equal(prepared$overlap_summary$n_x[[1L]], 48L)
  expect_equal(prepared$overlap_summary$n_y[[1L]], 48L)
})

test_that("prepared input runs a tiny high-level smoke", {
  input <- make_prepared_smoke_input()
  outdir <- tempfile("qtl2mr-prepared-result-")
  res <- run_qtl2mr_analysis(
    input,
    loci = 2,
    output_dir = outdir,
    overwrite = TRUE,
    make_figures = FALSE,
    make_html = TRUE,
    verbose = FALSE,
    run_loo_mr = TRUE,
    run_instrument_sensitivity = TRUE
  )
  expect_s3_class(res, "qtl2mr_analysis_result")
  expect_identical(res$pilot$summary$lmm_backend, "ordinary_least_squares")
  expect_identical(res$pilot$backend_provenance$backend, "ordinary_least_squares")
  expect_equal(qtl2mr_primary_result(res)$reported_estimator_policy, "diagnostic_only")
  expect_equal(qtl2mr_primary_result(res)$reported_estimator, "Star-space qtl2mr")
  expect_equal(qtl2mr_primary_result(res)$theta_reported, qtl2mr_primary_result(res)$theta_so)
  expect_true(is.finite(qtl2mr_primary_result(res)$theta))
  expect_true(is.finite(qtl2mr_primary_result(res)$profile_ci_lower))
  expect_true(is.finite(qtl2mr_primary_result(res)$profile_ci_upper))
  expect_equal(qtl2mr_primary_result(res)$n_loci, 2)
  expect_true(file.exists(file.path(outdir, "report.html")))
  expect_false(file.exists(file.path(outdir, "data", "qtl2_validation_summary.csv")))
  expect_true(file.exists(file.path(outdir, "data", "hsq.csv")))
})

test_that("public prepared analysis suppresses unreliable SO unless fallback is explicit", {
  input <- make_prepared_unreliable_so_input()
  common_args <- list(
    input = input,
    overwrite = TRUE,
    make_figures = FALSE,
    verbose = FALSE,
    diagnostic_bootstrap_B = 0L,
    heterogeneity_method = "none",
    run_loo_mr = FALSE,
    run_instrument_sensitivity = FALSE,
    progress = FALSE
  )

  diagnostic_dir <- tempfile("qtl2mr-unreliable-so-diagnostic-only-")
  diagnostic <- do.call(
    run_qtl2mr_analysis,
    c(common_args, list(output_dir = diagnostic_dir, make_html = TRUE))
  )
  reported <- qtl2mr_primary_result(diagnostic)
  expect_false(reported$primary_reliable[[1L]])
  expect_equal(reported$reported_estimator_policy[[1L]], "diagnostic_only")
  expect_equal(reported$reported_estimator[[1L]], "not_reported_unreliable_so")
  expect_equal(reported$reported_estimator_source[[1L]], "none")
  expect_true(all(is.na(unlist(reported[c(
    "theta_reported", "se_reported", "reported_ci_lower", "reported_ci_upper"
  )]))))
  expect_true(is.finite(reported$theta_so[[1L]]))
  expect_true(is.finite(reported$theta_two_step_ivw[[1L]]))
  expect_true("primary_so_unbounded" %in% diagnostic$warnings$warning_id)
  expect_true(file.exists(file.path(diagnostic_dir, "report.html")))

  fallback <- do.call(
    run_qtl2mr_analysis,
    c(
      common_args,
      list(
        output_dir = tempfile("qtl2mr-unreliable-so-auto-fallback-"),
        make_html = FALSE,
        reported_estimator_policy = "auto_fallback"
      )
    )
  )
  fallback_reported <- qtl2mr_primary_result(fallback)
  expect_false(fallback_reported$primary_reliable[[1L]])
  expect_equal(fallback_reported$reported_estimator_policy[[1L]], "auto_fallback")
  expect_equal(fallback_reported$reported_estimator[[1L]], "Two-step IVW")
  expect_equal(fallback_reported$reported_estimator_source[[1L]], "diagnostic_fallback")
  expect_equal(
    fallback_reported$theta_reported[[1L]],
    fallback_reported$theta_two_step_ivw[[1L]]
  )
})

test_that("run progress supports console, silent, and log-only modes", {
  input <- make_prepared_smoke_input()

  out_console <- tempfile("qtl2mr-progress-console-")
  progress_lines <- utils::capture.output(
    run_qtl2mr_analysis(
      input,
      exposure = "x",
      outcome = "y",
      loci = 2,
      output_dir = out_console,
      overwrite = TRUE,
      make_figures = FALSE,
      make_html = FALSE,
      verbose = FALSE,
      progress = TRUE,
      run_loo_mr = FALSE,
      run_instrument_sensitivity = FALSE
    ),
    type = "message"
  )
  expect_true(any(grepl("\\[1/5\\]", progress_lines)))
  expect_true(any(grepl("\\[2/5\\]", progress_lines)))
  expect_true(any(grepl("\\[3/5\\].*leave-one-locus-out", progress_lines)))
  expect_true(any(grepl("\\[4/5\\].*instrument-set", progress_lines)))
  expect_true(any(grepl("\\[5/5\\]", progress_lines)))
  expect_true(file.exists(file.path(out_console, "logs", "progress.log")))

  out_silent <- tempfile("qtl2mr-progress-silent-")
  silent_lines <- utils::capture.output(
    run_qtl2mr_analysis(
      input,
      exposure = "x",
      outcome = "y",
      loci = 2,
      output_dir = out_silent,
      overwrite = TRUE,
      make_figures = FALSE,
      make_html = FALSE,
      verbose = FALSE,
      progress = FALSE,
      run_loo_mr = FALSE,
      run_instrument_sensitivity = FALSE
    ),
    type = "message"
  )
  expect_false(any(grepl("\\[[1-5]/5\\]", silent_lines)))
  expect_false(file.exists(file.path(out_silent, "logs", "progress.log")))

  out_log <- tempfile("qtl2mr-progress-log-")
  log_only_lines <- utils::capture.output(
    run_qtl2mr_analysis(
      input,
      exposure = "x",
      outcome = "y",
      loci = 2,
      output_dir = out_log,
      overwrite = TRUE,
      make_figures = FALSE,
      make_html = FALSE,
      verbose = FALSE,
      progress = "log",
      run_loo_mr = FALSE,
      run_instrument_sensitivity = FALSE
    ),
    type = "message"
  )
  progress_log <- file.path(out_log, "logs", "progress.log")
  expect_false(any(grepl("\\[[1-5]/5\\]", log_only_lines)))
  expect_true(file.exists(progress_log))
  logged <- readLines(progress_log, warn = FALSE)
  expect_true(any(grepl("\\[1/5\\]", logged)))
  expect_true(any(grepl("\\[5/5\\]", logged)))
})
