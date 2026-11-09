# 0.6.1.9006 API clean-up: entry signatures, the fixed fit settings, founder QC on
# by default, and the two new helpers qtl2mr_diagnostics() and
# qtl2mr_combine_resources(). All fixtures are synthetic (traits x, y, w).

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Clean three-chromosome fixture: exact probabilities, no missing covariates, so
# a run signals nothing of its own. x has one instrument on each chromosome.
v066_fixture <- function(n = 80L, seed = 66001L, prefix = "s") {
  set.seed(seed)
  ids <- sprintf("%s%03d", prefix, seq_len(n))
  founders <- LETTERS[1:8]
  make_chr <- function(markers) {
    a <- array(stats::rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, founders, markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    a
  }
  gp <- list(`1` = make_chr(c("m11", "m12")), `2` = make_chr(c("m21", "m22")), `3` = make_chr(c("m31", "m32")))
  G <- matrix(stats::rnorm(n * 20L), n)
  K <- tcrossprod(G) / 20 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  kinship <- list(`1` = K, `2` = K * 0.99, `3` = K * 1.01)
  A <- qtl2mr:::make_founder_contrasts(founders)
  z <- function(chr, m, k) as.numeric(gp[[chr]][, , m] %*% A[, k])
  x <- 1.2 * z("1", "m11", 1) + 1.0 * z("2", "m21", 2) + 0.9 * z("3", "m31", 3) + stats::rnorm(n, sd = 0.4)
  w <- 1.0 * z("1", "m12", 4) + 1.1 * z("2", "m22", 5) + stats::rnorm(n, sd = 0.4)
  y <- 0.5 * x + 0.3 * w + stats::rnorm(n, sd = 0.4)
  pheno <- data.frame(x = x, y = y, w = w, row.names = ids)
  covar <- data.frame(sex = rep(0:1, length.out = n), row.names = ids)
  map <- data.frame(marker = c("m11", "m12", "m21", "m22", "m31", "m32"), chr = rep(c("1", "2", "3"), each = 2L),
    pos = rep(c(10, 20), 3L))
  list(ids = ids, gp = gp, kinship = kinship, pheno = pheno, covar = covar, map = map,
    resources = qtl2mr_prepare_resources(gp, kinship, pheno, covar, map),
    inst = data.frame(exposure = c("x", "x", "x", "w", "w"), marker = c("m11", "m21", "m31", "m12", "m22")))
}

# The founder-QC fixture of test-v063-founder-qc.R (copied: test files do not
# share definitions). Founder H at b2 is carried only by s1-s6, where y is
# missing, so y's own sample has E_min = 0 there and x's has E_min = 6.
v066_fqc_fixture <- function() {
  set.seed(63011)
  n <- 90L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2", "a3")), `2` = make_chr(c("b1", "b2")), `3` = make_chr(c("c1", "c2")),
    X = make_chr("xm1"))
  b2 <- gp[["2"]][, , "b2"]
  b2[, "H"] <- 0
  b2[1:6, "H"] <- 0.5
  b2[1:6, -8L] <- b2[1:6, -8L] / rowSums(b2[1:6, -8L]) * 0.5
  b2 <- b2 / rowSums(b2)
  gp[["2"]][, , "b2"] <- b2
  G <- matrix(rnorm(n * 25L), n)
  K <- tcrossprod(G) / 25 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  kinship <- list(`1` = K, `2` = K * 0.98, `3` = K * 1.01, X = K)
  sex <- rep(0:1, n / 2L)
  x <- 3 * gp[["1"]][, 1L, "a1"] + 2 * gp[["2"]][, 2L, "b2"] + 2 * gp[["2"]][, 8L, "b2"] + rnorm(n, sd = 0.3)
  w <- 2 * gp[["1"]][, 3L, "a3"] + 2 * gp[["3"]][, 4L, "c1"] + rnorm(n, sd = 0.3)
  y <- 0.6 * x + 0.3 * w + rnorm(n, sd = 0.3)
  pheno <- data.frame(x = x, w = w, y = y, row.names = ids)
  pheno$y[1:6] <- NA_real_
  pheno$x[7:8] <- NA_real_
  covar <- data.frame(sex = sex, age = round(runif(n, 50, 90)), row.names = ids)
  covar$age[9L] <- NA_real_
  map <- data.frame(marker = c("a1", "a2", "a3", "b1", "b2", "c1", "c2", "xm1"),
    chr = c("1", "1", "1", "2", "2", "3", "3", "X"), pos = c(1, 2, 3, 1, 2, 1, 2, 1))
  scan <- matrix(c(8, 4.5, 1, 1, 6, 1, 1, 5,
                   1, 1, 7, 4.2, 4.1, 6, 1, 1,
                   5, 1, 4.5, 1, 4.4, 1, 1, 1), ncol = 3L,
    dimnames = list(map$marker, c("x", "w", "y")))
  resources <- suppressWarnings(qtl2mr_prepare_resources(gp, kinship, pheno, covar, map))
  list(gp = gp, kinship = kinship, pheno = pheno, covar = covar, map = map, scan = scan, resources = resources)
}

v066_run <- function(input, ..., dir = tempfile("v066-run-")) {
  run_qtl2mr_analysis(input, output_dir = dir, make_html = FALSE, make_figures = FALSE, progress = FALSE, ...)
}

# The parts of a single-pair result that carry estimates and tables. Left out:
# output_dir, created_at, metadata$created_at, output_manifest and paths (paths,
# sizes and times of this run), pilot$prepared_metadata (records where the input
# came from) and raw_objects (a second copy of pilot). Closures inside the fit
# (objective functions) are dropped: identical() compares their environments.
v066_analysis_view <- function(res) {
  drop_closures <- function(x) rapply(list(x), function(f) NULL, classes = "function", how = "replace")[[1L]]
  pilot <- res$pilot
  pilot$prepared_metadata <- NULL
  list(pilot = drop_closures(pilot), loo_mr = res$loo_mr, instrument_sensitivity = res$instrument_sensitivity,
    steiger_sensitivity = res$steiger_sensitivity, warnings = res$warnings, method_comparison = res$method_comparison,
    first_stage_summary = res$first_stage_summary, tables = res$tables, primary = res$primary,
    diagnostics = drop_closures(res$diagnostics), figure_manifest = res$figure_manifest)
}

# The batch tables apart from runtime_sec, the run manifest (run id, times and
# cache counters), the creation time in the config and the resources metadata.
v066_batch_view <- function(b) {
  b$main_results$runtime_sec <- NULL
  b$config$created_at <- NULL
  b$run_manifest <- NULL
  b$resources_metadata <- NULL
  b
}

v066_md5 <- function(dir) {
  files <- list.files(dir, recursive = TRUE, all.files = TRUE)
  stats::setNames(unname(tools::md5sum(file.path(dir, files))), files)
}

# ---------------------------------------------------------------------------
# Signatures and exports
# ---------------------------------------------------------------------------

test_that("the entry signatures are exactly the 0.6.1.9006 ones, in order and with their defaults", {
  expected <- list(
    run_qtl2mr_analysis = function(input, exposure = NULL, outcome = NULL, output_dir = NULL, overwrite = FALSE,
                                   make_html = TRUE, make_figures = TRUE, progress = TRUE, n_workers = 1L,
                                   subsets = NULL) NULL,
    run_qtl2mr_batch = function(resources, pairs, instruments, min_samples = 50L, min_loci = 1L,
                                min_founder_copies = 5, subsets = NULL, n_workers = 1L, fail_fast = FALSE,
                                resume = FALSE, overwrite = FALSE, output_dir = NULL, progress = TRUE,
                                write_locus_blocks = FALSE) NULL,
    qtl2mr_select_instruments = function(resources, qtl_scan, exposures = NULL, lod_thresh = 4,
                                         peak_rule = c("peakdrop", "spacing"), peakdrop = 1.5, min_spacing = 1,
                                         exclude_chr = "X", min_founder_copies = 5, n_workers = 1L,
                                         progress = TRUE) NULL,
    qtl2mr_pair_input = function(resources, instruments, exposure, outcome, min_founder_copies = 5) NULL,
    qtl2mr_prepare_resources = function(allele_probs, kinship, pheno, covar, marker_map = NULL,
                                        metadata = list()) NULL,
    qtl2mr_combine_resources = function(x, y, kinship = NULL, common_only = FALSE) NULL,
    qtl2mr_diagnostics = function(result, nsim = 1000L, seed = 1L) NULL
  )
  ns <- asNamespace("qtl2mr")
  for (name in names(expected)) {
    actual <- formals(get(name, envir = ns))
    expect_identical(names(actual), names(formals(expected[[name]])), info = name)
    expect_identical(as.list(actual), as.list(formals(expected[[name]])), info = name)
  }
  # The defaults that changed or matter most, spelled out.
  expect_identical(formals(run_qtl2mr_analysis)$progress, TRUE)
  expect_identical(formals(run_qtl2mr_batch)$progress, TRUE)
  expect_identical(formals(qtl2mr_select_instruments)$progress, TRUE)
  for (fn in list(run_qtl2mr_batch, qtl2mr_select_instruments, qtl2mr_pair_input)) {
    expect_identical(formals(fn)$min_founder_copies, 5)
  }
  expect_identical(formals(qtl2mr_diagnostics)$nsim, 1000L)
  expect_identical(formals(qtl2mr_diagnostics)$seed, 1L)
})

test_that("every removed argument is refused as an unused argument", {
  removed_analysis <- list(loci = "m11", verbose = FALSE, on_nonestimable_locus = "drop",
    reported_estimator_policy = "strict_so", heterogeneity_method = "simulation", heterogeneity_nsim = 10L,
    heterogeneity_seed = 2L, diagnostic_bootstrap_B = 10L, diagnostic_bootstrap_seed = 2L,
    standardize_phenotypes = FALSE, cross_cov_model = "raw", run_loo_mr = FALSE, run_instrument_sensitivity = FALSE)
  for (arg in names(removed_analysis)) {
    expect_error(do.call(run_qtl2mr_analysis, c(list(NULL), removed_analysis[arg])), "unused argument", info = arg)
  }
  removed_batch <- list(on_nonestimable_locus = "drop", reported_estimator_policy = "strict_so",
    heterogeneity_method = "simulation", heterogeneity_nsim = 10L, heterogeneity_seed = 2L,
    diagnostic_bootstrap_B = 10L, diagnostic_bootstrap_seed = 2L, standardize_phenotypes = FALSE,
    cross_cov_model = "raw", sample_policy = "pair_common_complete")
  for (arg in names(removed_batch)) {
    expect_error(do.call(run_qtl2mr_batch, c(list(NULL, NULL, NULL), removed_batch[arg])), "unused argument", info = arg)
  }
  expect_error(qtl2mr_select_instruments(NULL, NULL, max_loci = 1L), "unused argument")
  expect_error(qtl2mr_select_instruments(NULL, NULL, F_X_min = 10), "unused argument")
  expect_error(qtl2mr_pair_input(NULL, NULL, "x", "y", sample_policy = "pair_common_complete"), "unused argument")
  expect_error(qtl2mr_prepare_resources(NULL, NULL, NULL, NULL, probability_policy = "strict"), "unused argument")
  # The legacy peak rule is gone; the spacing rule remains.
  f <- v066_fqc_fixture()
  expect_error(qtl2mr_select_instruments(f$resources, f$scan, peak_rule = "legacy", progress = FALSE),
    "should be one of")
  expect_s3_class(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, peak_rule = "spacing",
    progress = FALSE)), "data.frame")
})

test_that("run_qtl2mr_analysis() checks its flags, so a 0.6.1.9005 positional loci is refused", {
  f <- v066_fixture()
  input <- qtl2mr_pair_input(f$resources, f$inst, "x", "y")
  # Checked before anything is written (e54a2c2).
  for (flag in c("overwrite", "make_html", "make_figures")) {
    for (bad in list(NA, 1, "yes", c(TRUE, FALSE))) {
      dir <- tempfile("v066-flag-")
      args <- c(list(input, output_dir = dir, progress = FALSE), stats::setNames(list(bad), flag))
      expect_error(do.call(run_qtl2mr_analysis, args), paste(flag, "must be TRUE or FALSE"),
        info = paste(flag, deparse(bad)))
      expect_false(dir.exists(dir))
    }
  }
  # 0.6.1.9005 took run_qtl2mr_analysis(input, exposure, outcome, output_dir, loci, ...);
  # the fifth position is now overwrite, so an old positional loci is an error.
  for (loci in list(2L, c("m11", "m21"))) {
    dir <- tempfile("v066-positional-")
    expect_error(run_qtl2mr_analysis(input, "x", "y", dir, loci), "overwrite must be TRUE or FALSE")
    expect_false(dir.exists(dir))
  }
})

test_that("as_qtl2mr_prepared_input is internal and the two helpers are exported", {
  exports <- getNamespaceExports("qtl2mr")
  expect_false("as_qtl2mr_prepared_input" %in% exports)
  expect_true(all(c("qtl2mr_combine_resources", "qtl2mr_diagnostics") %in% exports))
  # Still there, unchanged, inside the namespace (qtl2mr_resources_input() uses it).
  expect_true(is.function(qtl2mr:::as_qtl2mr_prepared_input))
  expect_true(all(c("hsq", "exposure_samples", "probability_policy") %in%
    names(formals(qtl2mr:::as_qtl2mr_prepared_input))))
})

# ---------------------------------------------------------------------------
# Fixed settings
# ---------------------------------------------------------------------------

test_that("a batch records the fixed settings in batch_config.json", {
  f <- v066_fixture()
  dir <- tempfile("v066-config-")
  b <- run_qtl2mr_batch(f$resources, data.frame(exposure = c("x", "w"), outcome = "y"), f$inst, min_samples = 20L,
    progress = FALSE, output_dir = dir)
  fixed <- list(on_nonestimable_locus = "fail", reported_estimator_policy = "diagnostic_only",
    heterogeneity_method = "chisq", heterogeneity_nsim = 1000L, heterogeneity_seed = 1L,
    diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L, standardize_phenotypes = TRUE,
    cross_cov_model = "star")
  expect_identical(b$config[names(fixed)], fixed)
  expect_identical(b$config$min_founder_copies, 5)
  expect_null(b$config$sample_policy)
  json <- paste(readLines(file.path(dir, "batch_config.json")), collapse = "")
  for (pattern in c('"on_nonestimable_locus":"fail"', '"reported_estimator_policy":"diagnostic_only"',
                    '"heterogeneity_method":"chisq"', '"heterogeneity_nsim":1000', '"heterogeneity_seed":1',
                    '"diagnostic_bootstrap_B":0', '"diagnostic_bootstrap_seed":1', '"standardize_phenotypes":true',
                    '"cross_cov_model":"star"', '"min_founder_copies":5')) {
    expect_true(grepl(pattern, json, fixed = TRUE), info = pattern)
  }
  expect_false(grepl("sample_policy", json, fixed = TRUE))
  expect_identical(b$main_results$cross_cov_model, c("star", "star"))
  expect_identical(b$main_results$standardize_phenotypes, c(TRUE, TRUE))
  expect_identical(b$main_results$Q_p_method, c("chisq", "chisq"))
})

test_that("a single run uses chi-square Q, no bootstrap, the star model and standardised phenotypes", {
  f <- v066_fixture()
  dir <- tempfile("v066-single-")
  res <- v066_run(qtl2mr_pair_input(f$resources, f$inst, "x", "y"), dir = dir)
  s <- res$pilot$summary
  # Three loci: 0.6.1.9005 defaulted to the simulated Q p-value here (auto).
  expect_identical(s$n_loci, 3L)
  expect_identical(s$Q_p_method, "chisq")
  expect_identical(s$Q_p_reported, s$Q_p_chisq)
  expect_true(is.na(s$Q_p_sim))
  expect_identical(s$diagnostic_bootstrap_B, 0L)
  expect_true(is.na(s$se_weighted_median) && is.na(s$se_weighted_mode))
  expect_true(is.finite(s$theta_weighted_median) && is.finite(s$theta_weighted_mode))
  expect_identical(s$cross_cov_model, "star")
  expect_identical(s$standardize_phenotypes, TRUE)
  expect_identical(s$reported_estimator_policy, "diagnostic_only")
  expect_identical(s$effect_units, "SD_y per SD_x")
  expect_true(is.finite(s$theta_raw_units))
  # Leave-one-locus-out and instrument-set sensitivity always run.
  expect_s3_class(res$loo_mr, "data.frame")
  expect_identical(nrow(res$loo_mr), 3L)
  expect_setequal(res$loo_mr$marker, c("m11", "m21", "m31"))
  expect_s3_class(res$instrument_sensitivity, "data.frame")
  expect_true("all_loci" %in% res$instrument_sensitivity$set_name)
  expect_null(res$loo_error)
  expect_null(res$instrument_sensitivity_error)
  expect_identical(nrow(utils::read.csv(file.path(dir, "data", "loo_mr.csv"))), 3L)
  expect_identical(nrow(utils::read.csv(file.path(dir, "data", "instrument_sensitivity.csv"))),
    nrow(res$instrument_sensitivity))
  # Nonestimable loci fail the pair (on_nonestimable_locus = "fail").
  q <- v066_fqc_fixture()
  input <- qtl2mr_pair_input(q$resources, data.frame(exposure = "x", marker = c("a1", "b2")), "x", "y",
    min_founder_copies = NULL)
  expect_error(suppressWarnings(v066_run(input)), "not estimable at marker b2")
})

test_that("progress = FALSE prints nothing but still writes logs/progress.log", {
  f <- v066_fixture()
  input <- qtl2mr_pair_input(f$resources, f$inst, "x", "y")
  dir <- tempfile("v066-quiet-")
  # HTML report and figures included: no stage prints anything.
  expect_silent(res <- run_qtl2mr_analysis(input, output_dir = dir, progress = FALSE))
  log <- readLines(file.path(dir, "logs", "progress.log"))
  expect_match(log[[1L]], "qtl2mr progress log")
  expect_true(any(grepl("[1/5] Prepare input/output", log, fixed = TRUE)))
  expect_true(any(grepl("[3/5] Run leave-one-locus-out MR: complete", log, fixed = TRUE)))
  expect_true(any(grepl("[5/5] Write outputs: publish staged result", log, fixed = TRUE)))
  # progress = TRUE: the stages as messages and the printed summary, and the same log.
  dir2 <- tempfile("v066-loud-")
  msgs <- character()
  printed <- utils::capture.output(withCallingHandlers(
    run_qtl2mr_analysis(input, output_dir = dir2, make_html = FALSE, make_figures = FALSE, progress = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_true(any(grepl("[3/5] Run leave-one-locus-out MR", msgs, fixed = TRUE)))
  expect_gt(length(printed), 0L)
  # The fitting and sensitivity stages are logged the same way either way (stage 1
  # names the output folder, and stage 5 has no figure line without figures).
  log2 <- readLines(file.path(dir2, "logs", "progress.log"))
  stages <- function(x) {
    x <- sub("^\\[qtl2mr [0-9:]+\\] ", "", x)
    x[grepl("^\\[[2-4]/5\\]", x)]
  }
  expect_gt(length(stages(log)), 5L)
  expect_identical(stages(log2), stages(log))
  expect_error(run_qtl2mr_analysis(input, output_dir = tempfile(), progress = NA), "progress")
})

test_that("the single-pair result equals the internal path with the old explicit chisq / B = 0 arguments", {
  f <- v066_fixture()
  input <- qtl2mr_pair_input(f$resources, f$inst, "x", "y")
  res <- v066_run(input)
  # 0.6.1.9005 run_qtl2mr_analysis(heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L), all else at
  # its defaults, built its settings so and fitted the pair through qtl2mr_fit_prepared().
  old <- qtl2mr:::qtl2mr_fit_settings(on_nonestimable_locus = "fail", reported_estimator_policy = "diagnostic_only",
    heterogeneity_method = "chisq", heterogeneity_nsim = 1000L, heterogeneity_seed = 1L,
    diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L, standardize_phenotypes = TRUE,
    cross_cov_model = "star")
  expect_identical(qtl2mr:::qtl2mr_fixed_fit_settings(), old)
  ref <- qtl2mr:::qtl2mr_fit_prepared(input, "x", "y", old)
  # The same with every argument of run_prepared_qtl2mr() spelled out.
  ref2 <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE, standardize_phenotypes = TRUE,
    cross_cov_model = "star", heterogeneity_method = "chisq", q_sim_n = 1000L, q_sim_seed = 1L,
    diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L, reported_estimator_policy = "diagnostic_only",
    drop_nonestimable_loci = FALSE)
  # write_qtl2mr_results() completes the Q fields of the pilot before it is
  # returned; for a chisq fit this adds nothing to the summary.
  completed <- qtl2mr:::qtl2mr_complete_q_reporting(ref, q_sim_n = NULL)
  expect_identical(res$pilot$summary, ref$summary)
  expect_identical(res$pilot$summary, completed$summary)
  expect_identical(res$pilot$summary, ref2$summary)
  drop_closures <- function(x) rapply(list(x), function(f) NULL, classes = "function", how = "replace")[[1L]]
  expect_identical(drop_closures(res$pilot), drop_closures(completed))
  # And the sensitivity analyses that always run equal those computed from the reference fit.
  expect_identical(res$loo_mr, qtl2mr:::loo_multilocus_founder_mr(completed$fit))
})

# ---------------------------------------------------------------------------
# Founder QC on by default
# ---------------------------------------------------------------------------

test_that("qtl2mr_select_instruments() applies min_founder_copies = 5 by default; NULL turns it off", {
  f <- v066_fqc_fixture()
  def <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  expect_identical(def, suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE,
    min_founder_copies = 5)))
  expect_true(all(c("E_min_x", "min_founder_copies") %in% names(def)))
  expect_identical(unique(def$min_founder_copies), 5)
  # y at b2: no founder H on y's own sample (E_min = 0), dropped by the QC.
  yb2 <- def$exposure == "y" & def$marker == "b2"
  expect_identical(sum(yb2), 1L)
  expect_identical(def$drop_reason[yb2], "founder_qc")
  expect_equal(def$E_min_x[yb2], 0)
  expect_true(is.na(def$locus_rank[yb2]))
  # x at b2 has E_min = 6 (> 5) and is kept.
  xb2 <- def$exposure == "x" & def$marker == "b2"
  expect_equal(def$E_min_x[xb2], 6, tolerance = 1e-12)
  expect_identical(def$drop_reason[xb2], "")
  off <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE, min_founder_copies = NULL))
  expect_false(any(c("E_min_x", "min_founder_copies") %in% names(off)))
  expect_false("founder_qc" %in% off$drop_reason)
  expect_identical(off$marker, def$marker)
  # F_X_min and max_loci are gone; their columns stay, always NA.
  for (tab in list(def, off)) {
    expect_true(all(c("F_X_min", "max_loci") %in% names(tab)))
    expect_true(all(is.na(tab$F_X_min)) && all(is.na(tab$max_loci)))
  }
})

test_that("qtl2mr_pair_input() applies min_founder_copies = 5 by default; NULL turns it off", {
  f <- v066_fqc_fixture()
  inst <- data.frame(exposure = "x", marker = c("a1", "b2"))
  def <- qtl2mr_pair_input(f$resources, inst, "x", "y")
  expect_identical(def, qtl2mr_pair_input(f$resources, inst, "x", "y", min_founder_copies = 5))
  expect_identical(def$selected_loci$marker, "a1")
  qc <- def$metadata$founder_qc_outcome
  expect_identical(qc$min_founder_copies, 5)
  expect_true("E_min_y" %in% names(qc$loci))
  expect_equal(qc$loci$E_min_y[qc$loci$marker == "b2"], 0)
  off <- qtl2mr_pair_input(f$resources, inst, "x", "y", min_founder_copies = NULL)
  expect_identical(off$selected_loci$marker, c("a1", "b2"))
  expect_null(off$metadata$founder_qc_outcome)
  # b2 alone: nothing passes the default QC.
  expect_error(qtl2mr_pair_input(f$resources, data.frame(exposure = "x", marker = "b2"), "x", "y"),
    "too_few_loci: 0 loci remain after the outcome founder QC \\(min_founder_copies = 5\\)")
})

test_that("run_qtl2mr_batch() applies min_founder_copies = 5 by default; NULL turns it off", {
  f <- v066_fqc_fixture()
  inst <- data.frame(exposure = "x", marker = c("a1", "b2"))
  pairs <- data.frame(exposure = "x", outcome = c("y", "w"))
  dir <- tempfile("v066-qc-batch-")
  def <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = dir))
  five <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
    min_founder_copies = 5))
  expect_identical(v066_batch_view(def), v066_batch_view(five))
  expect_identical(def$config$min_founder_copies, 5)
  expect_true(grepl('"min_founder_copies":5', paste(readLines(file.path(dir, "batch_config.json")), collapse = ""),
    fixed = TRUE))
  expect_identical(def$main_results$status, c("ok", "ok"))
  expect_identical(def$main_results$n_loci_dropped_founder_qc, c(1L, 0L))
  xy <- def$locus_results[def$locus_results$pair_id == "x__to__y", ]
  expect_identical(xy$drop_reason, c("", "founder_qc_outcome"))
  expect_equal(xy$E_min_y[[2L]], 0)
  off <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
    min_founder_copies = NULL))
  expect_null(off$config$min_founder_copies)
  expect_false("n_loci_dropped_founder_qc" %in% names(off$main_results))
  expect_false("E_min_y" %in% names(off$locus_results))
  # Without the QC, y's founder effects at b2 are not estimable and the pair fails.
  expect_identical(off$main_results$status, c("error", "ok"))
  expect_match(off$main_results$error_message[[1L]], "not estimable at marker b2")
  # The QC-passing pair is unchanged by the QC.
  expect_identical(off$main_results$theta[[2L]], def$main_results$theta[[2L]])
})

# ---------------------------------------------------------------------------
# qtl2mr_diagnostics()
# ---------------------------------------------------------------------------

test_that("qtl2mr_diagnostics() equals the 0.6.1.9005 simulation and bootstrap path, column for column", {
  f <- v066_fixture()
  input <- qtl2mr_pair_input(f$resources, f$inst, "x", "y")
  dir <- tempfile("v066-diag-")
  res <- v066_run(input, dir = dir)
  expect_gte(res$pilot$summary$n_loci, 2L)
  nsim <- 40L
  seed <- 7L
  dg <- qtl2mr_diagnostics(res, nsim = nsim, seed = seed)
  ref <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE, heterogeneity_method = "simulation",
    q_sim_n = nsim, q_sim_seed = seed, diagnostic_bootstrap_B = nsim, diagnostic_bootstrap_seed = seed)
  expect_s3_class(dg, "data.frame")
  expect_identical(nrow(dg), 1L)
  expect_identical(dg$nsim, nsim)
  expect_identical(dg$seed, seed)
  cols <- setdiff(names(dg), c("nsim", "seed"))
  expect_true(all(cols %in% names(ref$summary)))
  for (col in cols) expect_identical(dg[[col]], ref$summary[[col]], info = col)
  expect_identical(dg$Q_sim_n, nsim)
  expect_true(is.finite(dg$Q_p_sim) && is.finite(dg$se_weighted_median) && is.finite(dg$se_weighted_mode))
  # The point estimates and Q are those of the main result.
  for (col in c("Q", "Q_p_chisq", "theta_weighted_median", "theta_weighted_mode")) {
    expect_identical(dg[[col]], res$pilot$summary[[col]], info = col)
  }
  # Another seed gives another simulation; the same seed reproduces it.
  expect_identical(qtl2mr_diagnostics(res, nsim = nsim, seed = seed), dg)
  other <- qtl2mr_diagnostics(res, nsim = nsim, seed = seed + 1L)
  expect_false(identical(other$se_weighted_median, dg$se_weighted_median))
  # Numeric nsim and seed are taken as whole numbers.
  expect_identical(qtl2mr_diagnostics(res, nsim = 40, seed = 7), dg)
})

test_that("qtl2mr_diagnostics() reads a result directory, writes data/diagnostics.csv and changes nothing else", {
  f <- v066_fixture()
  dir <- tempfile("v066-diag-dir-")
  res <- v066_run(qtl2mr_pair_input(f$resources, f$inst, "x", "y"), dir = dir)
  before <- v066_md5(dir)
  res_before <- res
  dg <- qtl2mr_diagnostics(res, nsim = 30L, seed = 3L)
  # (ii) Object and directory input give the same table (also the RDS path itself).
  expect_identical(qtl2mr_diagnostics(dir, nsim = 30L, seed = 3L), dg)
  expect_identical(qtl2mr_diagnostics(file.path(dir, "result_object.rds"), nsim = 30L, seed = 3L), dg)
  # (iii) The table is written to data/diagnostics.csv.
  path <- file.path(dir, "data", "diagnostics.csv")
  expect_true(file.exists(path))
  written <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_identical(names(written), names(dg))
  num <- names(dg)[vapply(dg, function(v) is.numeric(v) && all(is.finite(v)), logical(1))]
  expect_equal(as.list(written[num]), lapply(as.list(dg[num]), as.numeric), tolerance = 1e-12,
    ignore_attr = TRUE)
  expect_identical(written$ci_type_weighted_median, dg$ci_type_weighted_median)
  # (iv) The result object, main_result.csv and every other file are unchanged.
  expect_identical(res, res_before)
  after <- v066_md5(dir)
  expect_identical(setdiff(names(after), names(before)), "data/diagnostics.csv")
  expect_identical(after[names(before)], before)
  main <- utils::read.csv(file.path(dir, "data", "main_result.csv"), stringsAsFactors = FALSE)
  expect_identical(main$Q_p_method, "chisq")
  expect_identical(main$diagnostic_bootstrap_B, 0L)
  expect_true(is.na(main$se_weighted_median))
  expect_identical(readRDS(file.path(dir, "result_object.rds"))$pilot$summary, res$pilot$summary)
  # No result folder at the recorded path: the table is returned with a warning
  # and nothing is written (e54a2c2).
  moved <- res
  moved$output_dir <- tempfile("v066-absent-")
  expect_warning(absent <- qtl2mr_diagnostics(moved, nsim = 30L, seed = 3L), "no longer holds this run")
  expect_identical(absent, dg)
  expect_false(file.exists(moved$output_dir))
  # A path argument names the folder written to (31454cd), not the folder
  # recorded in the result. A copy, with the original still in place:
  parent <- tempfile("v066-diag-relocated-")
  dir.create(parent)
  expect_true(file.copy(dir, parent, recursive = TRUE))
  copy <- file.path(parent, basename(dir))
  unlink(c(path, file.path(copy, "data", "diagnostics.csv")))
  expect_identical(qtl2mr_diagnostics(copy, nsim = 30L, seed = 3L), dg)
  expect_true(file.exists(file.path(copy, "data", "diagnostics.csv")))
  expect_false(file.exists(path))
  # A moved folder (the recorded one no longer exists), by folder and by RDS path.
  relocated <- file.path(parent, "moved")
  expect_true(file.rename(dir, relocated))
  expect_false(dir.exists(dir))
  expect_identical(qtl2mr_diagnostics(relocated, nsim = 30L, seed = 3L), dg)
  written_moved <- file.path(relocated, "data", "diagnostics.csv")
  expect_true(file.exists(written_moved))
  expect_identical(utils::read.csv(written_moved, stringsAsFactors = FALSE), written)
  unlink(written_moved)
  expect_identical(qtl2mr_diagnostics(file.path(relocated, "result_object.rds"), nsim = 30L, seed = 3L), dg)
  expect_true(file.exists(written_moved))
  expect_false(dir.exists(dir))
  # The object itself, whose recorded folder has moved: a warning, and nothing
  # is written anywhere.
  moved_md5 <- v066_md5(relocated)
  expect_warning(qtl2mr_diagnostics(res, nsim = 30L, seed = 3L), "no longer holds this run")
  expect_false(dir.exists(dir))
  expect_identical(v066_md5(relocated), moved_md5)
  # A data/diagnostics.csv that is a link is refused, and its target is kept.
  if (identical(.Platform$OS.type, "unix")) {
    outside <- tempfile("v066-outside-", fileext = ".csv")
    writeLines("keep", outside)
    unlink(written_moved)
    expect_true(file.symlink(outside, written_moved))
    expect_error(qtl2mr_diagnostics(relocated, nsim = 30L, seed = 3L), "link or not a regular file")
    expect_identical(readLines(outside), "keep")
    expect_identical(Sys.readlink(written_moved), outside)
    expect_identical(list.files(file.path(relocated, "data"), pattern = "[.]tmp$"), character())
  }
})

test_that("qtl2mr_diagnostics() of a result object does not write into a folder another run has overwritten", {
  f <- v066_fixture()
  dir <- tempfile("v066-diag-overwritten-")
  first <- v066_run(qtl2mr_pair_input(f$resources, f$inst, "x", "y"), dir = dir)
  # The same object while its folder still holds it: written there.
  invisible(qtl2mr_diagnostics(first, nsim = 20L))
  expect_true(file.exists(file.path(dir, "data", "diagnostics.csv")))
  # Another analysis (w -> y) overwrites the folder (e54a2c2): the first
  # object's diagnostics are returned with a warning and not written there.
  second <- v066_run(qtl2mr_pair_input(f$resources, f$inst, "w", "y"), dir = dir, overwrite = TRUE)
  expect_false(file.exists(file.path(dir, "data", "diagnostics.csv")))
  expect_warning(dg <- qtl2mr_diagnostics(first, nsim = 20L), "no longer holds this run")
  expect_identical(dg$exposure, "x")
  expect_false(file.exists(file.path(dir, "data", "diagnostics.csv")))
  # The run now in the folder is written there, by object or by folder.
  expect_identical(qtl2mr_diagnostics(second, nsim = 20L)$exposure, "w")
  expect_identical(utils::read.csv(file.path(dir, "data", "diagnostics.csv"))$exposure, "w")
})

test_that("qtl2mr_diagnostics() refuses bad nsim, seed and non-result inputs", {
  f <- v066_fixture()
  dir <- tempfile("v066-diag-bad-")
  res <- v066_run(qtl2mr_pair_input(f$resources, f$inst, "x", "y"), dir = dir)
  for (bad in list(1, 0, -1, 2.5, NA_real_, NA_integer_, Inf, "40", c(10, 20), TRUE, NULL)) {
    expect_error(qtl2mr_diagnostics(res, nsim = bad), "nsim", info = deparse(bad))
  }
  # One draw cannot give a bootstrap SE (e54a2c2): nsim must be at least 2.
  expect_error(qtl2mr_diagnostics(res, nsim = 1L), "nsim must be one whole number >= 2")
  expect_identical(qtl2mr_diagnostics(res, nsim = 2L)$nsim, 2L)
  for (bad in list(1.5, NA_real_, Inf, "1", c(1, 2), TRUE, NULL)) {
    expect_error(qtl2mr_diagnostics(res, seed = bad), "seed", info = deparse(bad))
  }
  not_result <- "must be a run_qtl2mr_analysis\\(\\) result or its output directory"
  expect_error(qtl2mr_diagnostics(list()), not_result)
  expect_error(qtl2mr_diagnostics(res$pilot), not_result)
  expect_error(qtl2mr_diagnostics(data.frame(a = 1)), not_result)
  batch <- run_qtl2mr_batch(f$resources, data.frame(exposure = "x", outcome = "y"), f$inst, min_samples = 20L,
    progress = FALSE)
  expect_error(qtl2mr_diagnostics(batch), not_result)
  expect_error(qtl2mr_diagnostics(tempfile("v066-none-")), "no result_object.rds")
  empty <- tempfile("v066-empty-")
  dir.create(empty)
  expect_error(qtl2mr_diagnostics(empty), "no result_object.rds")
  rds <- tempfile("v066-other-", fileext = ".rds")
  saveRDS(list(a = 1), rds)
  expect_error(qtl2mr_diagnostics(rds), not_result)
  no_blocks <- res
  no_blocks$pilot$locus_blocks <- list()
  expect_error(qtl2mr_diagnostics(no_blocks), "no locus blocks")
  expect_error(qtl2mr_diagnostics(c(dir, dir)), not_result)
  expect_error(qtl2mr_diagnostics(NA_character_), not_result)
})

# ---------------------------------------------------------------------------
# qtl2mr_combine_resources()
# ---------------------------------------------------------------------------

test_that("qtl2mr_combine_resources() refuses inputs that do not describe the same markers, founders and mice", {
  f <- v066_fixture()
  prep <- function(gp = f$gp, kinship = f$kinship, pheno, covar = f$covar, map = f$map) {
    qtl2mr_prepare_resources(gp, kinship, pheno, covar, map)
  }
  rx <- prep(pheno = f$pheno["x"])
  ry <- prep(pheno = f$pheno[c("y", "w")])
  expect_s3_class(qtl2mr_combine_resources(rx, ry), "qtl2mr_resources")
  expect_error(qtl2mr_combine_resources(rx, list()), "x and y must both come from qtl2mr_prepare_resources")
  expect_error(qtl2mr_combine_resources(unclass(rx), ry), "x and y must both come from qtl2mr_prepare_resources")
  # 1. Marker map: a position differs.
  map2 <- f$map
  map2$pos[[3L]] <- 11
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno["y"], map = map2)), "same marker map")
  # 1. Founders: same names in another order.
  gp_f <- lapply(f$gp, function(a) {
    dimnames(a)[[2L]] <- LETTERS[c(2L, 1L, 3:8)]
    a
  })
  expect_error(qtl2mr_combine_resources(rx, prep(gp = gp_f, pheno = f$pheno["y"])), "same founders")
  # 2. A trait in both.
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno[c("x", "y")])), "share trait name\\(s\\): x")
  # 3. Covariate columns.
  covar2 <- f$covar
  names(covar2) <- "sex2"
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno["y"], covar = covar2)), "same covariate columns")
  covar3 <- cbind(f$covar, batch = 1)
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno["y"], covar = covar3)), "same covariate columns")
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno["y"], covar = NULL)),
    "both have covariates or both have none")
  # The same columns in another order are accepted.
  covar4 <- cbind(batch = 1, f$covar)
  rx4 <- prep(pheno = f$pheno["x"], covar = cbind(f$covar, batch = 1))
  expect_s3_class(qtl2mr_combine_resources(rx4, prep(pheno = f$pheno["y"], covar = covar4)), "qtl2mr_resources")
  # 4. A shared sample ID with other probabilities.
  gp_p <- f$gp
  p <- gp_p[["2"]][f$ids[[5L]], , "m22"]
  gp_p[["2"]][f$ids[[5L]], , "m22"] <- rev(p)
  expect_error(qtl2mr_combine_resources(rx, prep(gp = gp_p, pheno = f$pheno["y"])),
    "identical allele_probs \\(chromosome 2\\)")
  # 4. A shared sample ID with other covariates.
  covar5 <- f$covar
  covar5$sex[[7L]] <- 1L - covar5$sex[[7L]]
  expect_error(qtl2mr_combine_resources(rx, prep(pheno = f$pheno["y"], covar = covar5)),
    "same covariate values \\(column sex differs\\)")
  # Numbers are compared by value (e54a2c2): integer and double storage of the
  # same covariates and positions are accepted, and give the same combination.
  covar_dbl <- f$covar
  covar_dbl$sex <- as.double(covar_dbl$sex)
  map_int <- f$map
  map_int$pos <- as.integer(map_int$pos)
  expect_true(is.integer(f$covar$sex) && is.double(f$map$pos))
  by_value <- qtl2mr_combine_resources(rx, prep(pheno = f$pheno[c("y", "w")], covar = covar_dbl, map = map_int))
  expect_identical(by_value$pheno, qtl2mr_combine_resources(rx, ry)$pheno)
  # Each input must hold its samples in one order on every chromosome (e54a2c2).
  gp_o <- f$gp
  gp_o[["2"]] <- gp_o[["2"]][rev(f$ids), , , drop = FALSE]
  ry_o <- prep(gp = gp_o, pheno = f$pheno["y"])
  expect_error(qtl2mr_combine_resources(rx, ry_o),
    "each input must have the same samples, in the same order, on every chromosome .*chromosome 2 differs")
  expect_error(qtl2mr_combine_resources(ry_o, rx), "chromosome 2 differs")
  # Kinship: only one input has one.
  rx_nok <- prep(kinship = NULL, pheno = f$pheno["x"])
  expect_error(qtl2mr_combine_resources(rx_nok, ry), "only one of x and y has a kinship")
  expect_error(qtl2mr_combine_resources(ry, rx_nok), "only one of x and y has a kinship")
  # Kinship: both cover all samples but disagree.
  kin2 <- f$kinship
  kin2[["3"]][1L, 2L] <- kin2[["3"]][2L, 1L] <- kin2[["3"]][1L, 2L] + 1e-9
  expect_error(qtl2mr_combine_resources(rx, prep(kinship = kin2, pheno = f$pheno["y"])),
    "kinships of x and y do not both cover all combined samples with identical values")
  # Kinship: one matrix for all chromosomes against a LOCO list, or other
  # chromosomes, has its own message (e54a2c2).
  expect_error(qtl2mr_combine_resources(rx, prep(kinship = f$kinship[["1"]], pheno = f$pheno["y"])),
    "kinships of x and y have different structures")
  expect_error(qtl2mr_combine_resources(rx, prep(kinship = f$kinship[c("1", "3", "2")], pheno = f$pheno["y"])),
    "kinships of x and y have different structures")
  # common_only must be a flag.
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1)) {
    expect_error(qtl2mr_combine_resources(rx, ry, common_only = bad), "common_only")
  }
})

# x holds samples 1-50 with trait x, y holds samples 31-80 (in reverse order)
# with traits y and w; both carry the kinship of all 80 samples.
v066_split <- function(f) {
  ids <- f$ids
  ix <- ids[1:50]
  iy <- rev(ids[31:80])
  sub <- function(i) lapply(f$gp, function(a) a[i, , , drop = FALSE])
  list(ix = ix, iy = iy,
    x = qtl2mr_prepare_resources(sub(ix), f$kinship, f$pheno[ix, "x", drop = FALSE], f$covar[ix, , drop = FALSE],
      f$map),
    y = qtl2mr_prepare_resources(sub(iy), f$kinship, f$pheno[iy, c("y", "w")], f$covar[iy, , drop = FALSE], f$map))
}

test_that("qtl2mr_combine_resources() fills NA for one-sided samples and records the counts", {
  f <- v066_fixture()
  s <- v066_split(f)
  rc <- qtl2mr_combine_resources(s$x, s$y)
  expect_s3_class(rc, "qtl2mr_resources")
  # Samples: those of x first, then those only in y, in y's order.
  only_y <- setdiff(s$iy, s$ix)
  ids <- c(s$ix, only_y)
  expect_identical(rc$sample_ids, ids)
  expect_identical(only_y, rev(f$ids[51:80]))
  expect_identical(rownames(rc$pheno), ids)
  expect_identical(names(rc$pheno), c("x", "y", "w"))
  expect_identical(rc$pheno[s$ix, "x"], f$pheno[s$ix, "x"])
  expect_true(all(is.na(rc$pheno[only_y, "x"])))
  expect_identical(rc$pheno[s$iy, "y"], f$pheno[s$iy, "y"])
  expect_identical(rc$pheno[s$iy, "w"], f$pheno[s$iy, "w"])
  only_x <- setdiff(s$ix, s$iy)
  expect_true(all(is.na(rc$pheno[only_x, "y"])) && all(is.na(rc$pheno[only_x, "w"])))
  for (chr in names(f$gp)) {
    expect_identical(dimnames(rc$allele_probs[[chr]]), list(ids, LETTERS[1:8], dimnames(f$gp[[chr]])[[3L]]))
    expect_identical(rc$allele_probs[[chr]], f$gp[[chr]][ids, , , drop = FALSE])
  }
  expect_identical(rc$covar, f$covar[ids, , drop = FALSE])
  # Both kinships cover all 80 samples identically: taken from x and y.
  for (chr in names(f$kinship)) expect_identical(rc$kinship[[chr]], f$kinship[[chr]][ids, ids])
  rec <- rc$metadata$combined_resources
  expect_identical(rec$n_samples_x_only, 30L)
  expect_identical(rec$n_samples_y_only, 30L)
  expect_identical(rec$n_samples_overlap, 20L)
  expect_identical(rec$n_samples_combined, 80L)
  expect_false(rec$common_only)
  expect_identical(rec$traits_x, "x")
  expect_identical(rec$traits_y, c("y", "w"))
  expect_identical(rec$kinship_source, "x_and_y")
  expect_true(all(c("x_metadata", "y_metadata") %in% names(rc$metadata)))
  expect_identical(rc$probability_policy, "normalize_rounding")
  # Each trait of a pair is fitted on its own samples.
  res <- v066_run(qtl2mr_pair_input(rc, f$inst, "x", "y"))
  ov <- res$pilot$overlap_summary
  expect_identical(as.integer(ov$n_x), rep(50L, 3L))
  expect_identical(as.integer(ov$n_y), rep(50L, 3L))
  expect_identical(as.integer(ov$n_overlap), rep(20L, 3L))
  for (b in res$pilot$locus_blocks) {
    expect_setequal(b$fit_x$sample_ids, s$ix)
    expect_setequal(b$fit_y$sample_ids, s$iy)
  }
})

test_that("common_only = TRUE keeps only the samples genotyped in both inputs", {
  f <- v066_fixture()
  s <- v066_split(f)
  rc <- qtl2mr_combine_resources(s$x, s$y, common_only = TRUE)
  common <- intersect(s$ix, s$iy)
  expect_identical(common, f$ids[31:50])
  expect_identical(rc$sample_ids, common)
  expect_identical(rownames(rc$pheno), common)
  expect_identical(rownames(rc$covar), common)
  expect_false(anyNA(rc$pheno))
  for (chr in names(f$gp)) {
    expect_identical(rc$allele_probs[[chr]], f$gp[[chr]][common, , , drop = FALSE])
    expect_identical(rc$kinship[[chr]], f$kinship[[chr]][common, common])
  }
  rec <- rc$metadata$combined_resources
  expect_true(rec$common_only)
  expect_identical(rec$n_samples_combined, 20L)
  expect_identical(c(rec$n_samples_x_only, rec$n_samples_y_only, rec$n_samples_overlap), c(30L, 30L, 20L))
  # 20 samples carry about 5 copies of each founder, so the founder QC is off here.
  res <- v066_run(qtl2mr_pair_input(rc, f$inst[f$inst$exposure == "x", ][1:2, ], "x", "y",
    min_founder_copies = NULL))
  expect_identical(as.integer(res$pilot$overlap_summary$n_x), c(20L, 20L))
  expect_identical(as.integer(res$pilot$overlap_summary$n_overlap), c(20L, 20L))
  # No shared samples: refused.
  g <- v066_fixture(n = 40L, seed = 66002L, prefix = "t")
  ry_other <- qtl2mr_prepare_resources(g$gp, g$kinship, g$pheno["y"], g$covar, g$map)
  expect_error(qtl2mr_combine_resources(s$x, ry_other, common_only = TRUE),
    "common_only = TRUE but x and y share no genotyped samples")
})

test_that("kinship comes from x and y only when both cover every sample identically; otherwise it is required", {
  f <- v066_fixture()
  ix <- f$ids[1:50]
  iy <- f$ids[31:80]
  sub <- function(i) lapply(f$gp, function(a) a[i, , , drop = FALSE])
  kin_sub <- function(i) lapply(f$kinship, function(K) K[i, i])
  # Each input carries the kinship of its own samples only.
  rx <- qtl2mr_prepare_resources(sub(ix), kin_sub(ix), f$pheno[ix, "x", drop = FALSE], f$covar[ix, , drop = FALSE],
    f$map)
  ry <- qtl2mr_prepare_resources(sub(iy), kin_sub(iy), f$pheno[iy, "y", drop = FALSE], f$covar[iy, , drop = FALSE],
    f$map)
  expect_error(qtl2mr_combine_resources(rx, ry), "supply kinship computed on the combined probabilities")
  supplied <- kin_sub(f$ids)
  rc <- qtl2mr_combine_resources(rx, ry, kinship = supplied)
  expect_identical(rc$kinship, supplied)
  expect_identical(rc$metadata$combined_resources$kinship_source, "supplied")
  # A supplied kinship is used even when x and y would agree.
  s <- v066_split(f)
  scaled <- lapply(f$kinship, function(K) K * 2)
  expect_identical(qtl2mr_combine_resources(s$x, s$y, kinship = scaled)$kinship, scaled)
  expect_identical(qtl2mr_combine_resources(s$x, s$y)$metadata$combined_resources$kinship_source, "x_and_y")
  # Neither has a kinship: none.
  rx0 <- qtl2mr_prepare_resources(sub(ix), NULL, f$pheno[ix, "x", drop = FALSE], f$covar[ix, , drop = FALSE], f$map)
  ry0 <- qtl2mr_prepare_resources(sub(iy), NULL, f$pheno[iy, "y", drop = FALSE], f$covar[iy, , drop = FALSE], f$map)
  r0 <- qtl2mr_combine_resources(rx0, ry0)
  expect_null(r0$kinship)
  expect_identical(r0$metadata$combined_resources$kinship_source, "none")
})

test_that("splitting resources by trait and combining them gives identical analyses and batches", {
  f <- v066_fixture()
  # qtl2 allele probabilities carry a class and attributes; the combined
  # resources keep those of x (9a347fa).
  attributes(f$gp) <- c(attributes(f$gp), list(crosstype = "do", is_x_chr = c(`1` = FALSE, `2` = FALSE, `3` = FALSE),
    alleles = LETTERS[1:8], alleleprobs = TRUE))
  class(f$gp) <- c("calc_genoprob", "list")
  pheno <- f$pheno
  pheno$x[1:5] <- NA_real_
  pheno$y[76:80] <- NA_real_
  pheno$w[c(10, 20)] <- NA_real_
  whole <- qtl2mr_prepare_resources(f$gp, f$kinship, pheno, f$covar, f$map)
  rx <- qtl2mr_prepare_resources(f$gp, f$kinship, pheno["x"], f$covar, f$map)
  ry <- qtl2mr_prepare_resources(f$gp, f$kinship, pheno[c("y", "w")], f$covar, f$map)
  rc <- qtl2mr_combine_resources(rx, ry)
  expect_s3_class(rc$allele_probs, "calc_genoprob")
  expect_identical(attributes(rc$allele_probs), attributes(f$gp))
  for (nm in c("allele_probs", "kinship", "pheno", "covar", "marker_map", "sample_ids", "probability_policy",
               "sample_alignment")) {
    expect_identical(rc[[nm]], whole[[nm]], info = nm)
  }
  # Single pair: the inputs differ only in metadata, the results not at all.
  qtl2mr:::qtl2mr_cache_reset()
  in_whole <- qtl2mr_pair_input(whole, f$inst, "x", "y")
  in_comb <- qtl2mr_pair_input(rc, f$inst, "x", "y")
  in_whole$metadata <- in_comb$metadata <- NULL
  expect_identical(in_comb, in_whole)
  a <- v066_run(qtl2mr_pair_input(whole, f$inst, "x", "y"))
  qtl2mr:::qtl2mr_cache_reset()
  b <- v066_run(qtl2mr_pair_input(rc, f$inst, "x", "y"))
  expect_identical(v066_analysis_view(b), v066_analysis_view(a))
  expect_identical(b$pilot$summary$theta, a$pilot$summary$theta)
  expect_identical(b$pilot$summary$n_x_median, 75L)
  expect_identical(b$pilot$summary$n_y_median, 75L)
  # Two-pair batch, with the cache cleared before each run.
  pairs <- data.frame(exposure = c("x", "w"), outcome = c("y", "y"))
  qtl2mr:::qtl2mr_cache_reset()
  ba <- run_qtl2mr_batch(whole, pairs, f$inst, min_samples = 20L, progress = FALSE)
  qtl2mr:::qtl2mr_cache_reset()
  bb <- run_qtl2mr_batch(rc, pairs, f$inst, min_samples = 20L, progress = FALSE)
  expect_identical(bb$main_results$status, c("ok", "ok"))
  expect_identical(v066_batch_view(bb), v066_batch_view(ba))
  # The batch row for x -> y carries the estimate of the single analysis.
  expect_identical(bb$main_results$theta[[1L]], b$pilot$summary$theta)
  expect_identical(bb$main_results$Q[[1L]], b$pilot$summary$Q)
})

test_that("splitting by samples equals one resources object with the other samples' phenotypes missing", {
  f <- v066_fixture()
  ids <- f$ids
  sub <- function(i) lapply(f$gp, function(a) a[i, , , drop = FALSE])
  pheno <- f$pheno
  pheno$x[51:80] <- NA_real_
  pheno[1:30, c("y", "w")] <- NA_real_
  whole <- qtl2mr_prepare_resources(f$gp, f$kinship, pheno, f$covar, f$map)
  pairs <- data.frame(exposure = c("x", "w"), outcome = c("y", "y"))
  combine <- function(iy) {
    ix <- ids[1:50]
    qtl2mr_combine_resources(
      qtl2mr_prepare_resources(sub(ix), f$kinship, f$pheno[ix, "x", drop = FALSE], f$covar[ix, , drop = FALSE], f$map),
      qtl2mr_prepare_resources(sub(iy), f$kinship, f$pheno[iy, c("y", "w")], f$covar[iy, , drop = FALSE], f$map))
  }
  qtl2mr:::qtl2mr_cache_reset()
  a <- v066_run(qtl2mr_pair_input(whole, f$inst, "x", "y"))
  qtl2mr:::qtl2mr_cache_reset()
  ba <- run_qtl2mr_batch(whole, pairs, f$inst, min_samples = 20L, progress = FALSE)
  # Samples of y in the same order: the combined data are those of `whole`, and so are the results.
  rc <- combine(ids[31:80])
  for (nm in c("allele_probs", "kinship", "pheno", "covar", "sample_ids", "sample_alignment")) {
    expect_identical(rc[[nm]], whole[[nm]], info = nm)
  }
  qtl2mr:::qtl2mr_cache_reset()
  b <- v066_run(qtl2mr_pair_input(rc, f$inst, "x", "y"))
  expect_identical(v066_analysis_view(b), v066_analysis_view(a))
  qtl2mr:::qtl2mr_cache_reset()
  expect_identical(v066_batch_view(run_qtl2mr_batch(rc, pairs, f$inst, min_samples = 20L, progress = FALSE)),
    v066_batch_view(ba))
  # Samples of y in another order: the combined samples follow it (x's, then y's own), so
  # sums run in another order and the estimates agree to rounding only.
  rr <- combine(rev(ids[31:80]))
  expect_identical(rr$sample_ids, c(ids[1:50], rev(ids[51:80])))
  qtl2mr:::qtl2mr_cache_reset()
  r <- v066_run(qtl2mr_pair_input(rr, f$inst, "x", "y"))
  expect_equal(r$pilot$summary, a$pilot$summary)
})

test_that("two cohorts with no shared samples: the user's kinship, each trait on its own cohort", {
  c1 <- v066_fixture(n = 60L, seed = 66011L, prefix = "a")
  c2 <- v066_fixture(n = 60L, seed = 66012L, prefix = "b")
  rx <- qtl2mr_prepare_resources(c1$gp, c1$kinship, c1$pheno["x"], c1$covar, c1$map)
  ry <- qtl2mr_prepare_resources(c2$gp, c2$kinship, c2$pheno["y"], c2$covar, c2$map)
  expect_error(qtl2mr_combine_resources(rx, ry), "supply kinship computed on the combined probabilities")
  ids <- c(c1$ids, c2$ids)
  # A kinship of all 120 mice (here block-diagonal: the cohorts are unrelated).
  kinship <- stats::setNames(lapply(names(c1$kinship), function(chr) {
    K <- matrix(0, 120L, 120L, dimnames = list(ids, ids))
    K[c1$ids, c1$ids] <- c1$kinship[[chr]]
    K[c2$ids, c2$ids] <- c2$kinship[[chr]]
    K
  }), names(c1$kinship))
  rc <- qtl2mr_combine_resources(rx, ry, kinship = kinship)
  rec <- rc$metadata$combined_resources
  expect_identical(c(rec$n_samples_x_only, rec$n_samples_y_only, rec$n_samples_overlap, rec$n_samples_combined),
    c(60L, 60L, 0L, 120L))
  expect_identical(rec$kinship_source, "supplied")
  expect_identical(rc$sample_ids, ids)
  expect_true(all(is.na(rc$pheno[c2$ids, "x"])) && all(is.na(rc$pheno[c1$ids, "y"])))
  res <- v066_run(qtl2mr_pair_input(rc, c1$inst, "x", "y"))
  s <- res$pilot$summary
  expect_identical(s$n_x_median, 60L)
  expect_identical(s$n_y_median, 60L)
  expect_identical(s$n_overlap_median, 0L)
  ov <- res$pilot$overlap_summary
  expect_identical(as.integer(ov$n_x), rep(60L, 3L))
  expect_identical(as.integer(ov$n_y), rep(60L, 3L))
  expect_identical(as.integer(ov$n_overlap), rep(0L, 3L))
  for (b in res$pilot$locus_blocks) {
    expect_identical(sort(b$fit_x$sample_ids), sort(c1$ids))
    expect_identical(sort(b$fit_y$sample_ids), sort(c2$ids))
  }
  expect_true(is.finite(s$theta))
  # The exposure fit is that of cohort 1 on its own.
  alone <- qtl2mr:::run_prepared_qtl2mr(qtl2mr_pair_input(c1$resources, c1$inst, "x", "w"), "x", "w",
    verbose = FALSE, heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L)
  for (m in names(res$pilot$locus_blocks)) {
    fx <- res$pilot$locus_blocks[[m]]$fit_x
    fa <- alone$locus_blocks[[m]]$fit_x
    expect_identical(fx$sample_ids, fa$sample_ids, info = m)
    expect_identical(fx$hsq, fa$hsq, info = m)
    expect_identical(fx$beta_fc, fa$beta_fc, info = m)
    expect_identical(fx$vcov_fc, fa$vcov_fc, info = m)
  }
})
