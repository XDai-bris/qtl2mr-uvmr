# 0.6.1.9003: opt-in founder-design QC, min_founder_copies. E = min over founders
# of sum_i 2 p_ik on the samples a trait is fitted on, from the probabilities the
# package fits; a locus is dropped when E < min_founder_copies.

fqc_fixture <- function(round_digits = NULL) {
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
  # Founder H at b2 is carried only by s1-s6 (p = 0.5), the samples where y is missing.
  b2 <- gp[["2"]][, , "b2"]
  b2[, "H"] <- 0
  b2[1:6, "H"] <- 0.5
  b2[1:6, -8L] <- b2[1:6, -8L] / rowSums(b2[1:6, -8L]) * 0.5
  b2 <- b2 / rowSums(b2)
  gp[["2"]][, , "b2"] <- b2
  if (!is.null(round_digits)) gp <- lapply(gp, round, round_digits)
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

# Hand calculation: own complete sample (finite trait, complete covariates),
# probabilities normalised as the normalize_rounding policy does.
fqc_hand_E <- function(f, trait, chr, marker) {
  ids <- rownames(f$pheno)[is.finite(f$pheno[[trait]]) & stats::complete.cases(f$covar)]
  P <- f$gp[[chr]][ids, , marker]
  s <- rowSums(P)
  adjust <- abs(s - 1) > 1e-6 + 16 * .Machine$double.eps
  P[adjust, ] <- P[adjust, ] / s[adjust]
  E <- 2 * colSums(P)
  min(E)
}

test_that("E_min is the rarest founder's expected allele count on the fitted samples", {
  for (digits in list(NULL, 5L)) {
    f <- fqc_fixture(round_digits = digits)
    inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, exclude_chr = NULL, progress = FALSE,
      min_founder_copies = 0))
    auto <- inst$chr != "X"
    hand <- mapply(function(e, c, m) fqc_hand_E(f, e, c, m), inst$exposure[auto], inst$chr[auto], inst$marker[auto])
    expect_equal(inst$E_min_x[auto], unname(hand), tolerance = 1e-13)
    expect_true(all(is.na(inst$E_min_x[!auto])))
    # b2: founder H only on s1-s6 (p = 0.5): 6 copies for x, none for y's own sample.
    expect_equal(inst$E_min_x[inst$exposure == "x" & inst$marker == "b2"], fqc_hand_E(f, "x", "2", "b2"))
    if (is.null(digits)) expect_equal(fqc_hand_E(f, "x", "2", "b2"), 6, tolerance = 1e-12)
    # The samples used are exactly those the LMM context is fitted on.
    for (trait in c("x", "y")) {
      ctx <- qtl2mr:::prepare_lmm_context(qtl2mr:::qtl2mr_prepared_trait(f$resources, trait), kinship = f$kinship[["2"]],
        addcovar = as.data.frame(f$covar), sample_ids = f$resources$sample_ids, trait_name = trait, chr = "2")
      expect_identical(qtl2mr:::qtl2mr_fitted_sample_ids(qtl2mr:::qtl2mr_prepared_trait(f$resources, trait),
        f$kinship[["2"]], as.data.frame(f$covar), f$resources$sample_ids, trait), ctx$sample_ids)
    }
  }
})

test_that("off by default: no columns, fields or config entries are added", {
  f <- fqc_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  expect_false(any(c("E_min_x", "min_founder_copies") %in% names(inst)))
  expect_identical(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE,
    min_founder_copies = NULL)), inst)
  input <- qtl2mr_pair_input(f$resources, inst, "x", "w")
  expect_identical(qtl2mr_pair_input(f$resources, inst, "x", "w", min_founder_copies = NULL), input)
  expect_null(input$metadata$founder_qc_outcome)
  b <- suppressWarnings(run_qtl2mr_batch(f$resources, data.frame(exposure = "x", outcome = "w"), inst,
    min_samples = 20L, progress = FALSE))
  expect_false("n_loci_dropped_founder_qc" %in% names(b$main_results))
  expect_false("E_min_y" %in% names(b$locus_results))
  expect_null(b$config$min_founder_copies)
  for (bad in list(-1, NA_real_, c(1, 2), "5", TRUE, Inf)) {
    expect_error(qtl2mr_select_instruments(f$resources, f$scan, min_founder_copies = bad), "min_founder_copies")
    expect_error(run_qtl2mr_batch(f$resources, data.frame(exposure = "x", outcome = "w"), inst,
      min_founder_copies = bad), "min_founder_copies")
    expect_error(qtl2mr_pair_input(f$resources, inst, "x", "w", min_founder_copies = bad), "min_founder_copies")
  }
})

test_that("exposure QC sits after non_autosomal and before fitting, F and max_loci", {
  f <- fqc_fixture()
  E_b2 <- fqc_hand_E(f, "x", "2", "b2")
  off <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, exposures = "x", exclude_chr = NULL,
    progress = FALSE))
  on <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, exposures = "x", exclude_chr = NULL,
    progress = FALSE, min_founder_copies = E_b2 + 0.01))
  expect_identical(on$marker, off$marker)
  expect_identical(on$drop_reason[on$marker == "xm1"], "non_autosomal")
  expect_identical(on$drop_reason[on$marker == "b2"], "founder_qc")
  # A locus dropped by the QC is not fitted.
  expect_true(is.na(on$F_X_block[on$marker == "b2"]) && is.na(on$n_x[on$marker == "b2"]))
  keep <- on$marker != "b2"
  expect_identical(on$F_X_block[keep], off$F_X_block[keep])
  expect_identical(unique(on$min_founder_copies), E_b2 + 0.01)
  expect_identical(match("E_min_x", names(on)), match("hsq_x", names(on)) + 1L)
  # Threshold 0 drops nothing: the table is the QC-off table plus the QC columns.
  zero <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE, min_founder_copies = 0))
  plain <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  expect_identical(zero[names(plain)], plain)
  # max_loci counts loci that pass the QC: for y, b2 has the highest LOD but no founder H on
  # y's own sample, so a1 (the next) is kept with rank 1 and a3 is the one over max_loci.
  scan <- f$scan
  scan["b2", "y"] <- 9
  top <- suppressWarnings(qtl2mr_select_instruments(f$resources, scan, exposures = "y", max_loci = 1L,
    progress = FALSE, min_founder_copies = 5))
  expect_identical(top$marker, c("b2", "a1", "a3"))
  expect_identical(top$drop_reason, c("founder_qc", "", "max_loci"))
  expect_identical(top$locus_rank, c(NA, 1L, NA))
  expect_equal(top$E_min_x[[1L]], 0)
  # F is applied after the QC, to the fitted loci only.
  fmin <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, exposures = "x", progress = FALSE,
    F_X_min = 1e6, min_founder_copies = E_b2 + 0.01))
  expect_identical(fmin$drop_reason, ifelse(fmin$marker == "b2", "founder_qc", "F"))
  expect_identical(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE, n_workers = 2L,
    min_founder_copies = E_b2 + 0.01)), suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan,
    progress = FALSE, min_founder_copies = E_b2 + 0.01)))
})

fqc_single <- function(f, inst, exposure, outcome, ...) {
  suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, exposure, outcome, ...),
    output_dir = tempfile("v063-single-"), make_html = FALSE, make_figures = FALSE, verbose = FALSE,
    heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L, run_loo_mr = FALSE,
    run_instrument_sensitivity = FALSE))
}

fqc_row_equal <- function(batch_row, single_summary) {
  cols <- qtl2mr:::qtl2mr_batch_summary_columns
  a <- as.list(batch_row[1L, cols])
  b <- as.list(single_summary[1L, cols])
  num <- vapply(a, is.numeric, logical(1)) & vapply(b, is.numeric, logical(1))
  a[num] <- lapply(a[num], as.double)
  b[num] <- lapply(b[num], as.double)
  identical(a, b)
}

test_that("the outcome QC drops a locus for that pair only, and the batch row equals the single analysis", {
  f <- fqc_fixture()
  inst <- data.frame(exposure = c("x", "x", "w", "w"), marker = c("a1", "b2", "a3", "c1"))
  pairs <- data.frame(exposure = c("x", "x", "w"), outcome = c("y", "w", "y"))
  for (policy in c("own_complete", "pair_common_complete")) {
    batch <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
      sample_policy = policy, min_founder_copies = 5))
    main <- batch$main_results
    loci <- batch$locus_results
    expect_identical(main$status, rep("ok", 3L))
    expect_identical(main$n_loci_dropped_founder_qc, c(1L, 0L, 0L))
    expect_identical(main$n_loci, c(1L, 2L, 2L))
    xy <- loci[loci$pair_id == "x__to__y", ]
    expect_identical(xy$marker, c("a1", "b2"))
    expect_identical(xy$kept, c(TRUE, FALSE))
    expect_identical(xy$drop_reason, c("", "founder_qc_outcome"))
    expect_equal(xy$E_min_y[[2L]], 0)
    expect_true(is.na(xy$F_X_block[[2L]]))
    # The same locus stays for x -> w, whose outcome sample carries founder H.
    xw <- loci[loci$pair_id == "x__to__w", ]
    expect_identical(xw$kept, c(TRUE, TRUE))
    expect_true(all(xw$E_min_y >= 5))
    expect_identical(batch$config$min_founder_copies, 5)
    for (i in seq_len(nrow(pairs))) {
      single <- fqc_single(f, inst, pairs$exposure[[i]], pairs$outcome[[i]], sample_policy = policy,
        min_founder_copies = 5)
      expect_true(fqc_row_equal(main[i, ], single$pilot$summary))
      qc <- single$pilot$prepared_metadata$founder_qc_outcome
      expect_identical(qc$min_founder_copies, 5)
      expect_identical(qc$loci$E_min_y, loci$E_min_y[loci$pair_id == main$pair_id[[i]]])
    }
  }
  # Too few loci after the QC follows the too-few-loci path, with the QC rows recorded.
  only_b2 <- suppressWarnings(run_qtl2mr_batch(f$resources, data.frame(exposure = "x", outcome = "y"),
    data.frame(exposure = "x", marker = "b2"), min_samples = 20L, progress = FALSE, min_founder_copies = 5))
  expect_identical(only_b2$main_results$status, "error")
  expect_match(only_b2$main_results$error_message, "too_few_loci: 0 loci remain after the outcome founder QC")
  expect_identical(only_b2$main_results$n_loci_dropped_founder_qc, 1L)
  expect_identical(only_b2$locus_results$drop_reason, "founder_qc_outcome")
  expect_error(qtl2mr_pair_input(f$resources, data.frame(exposure = "x", marker = "b2"), "x", "y",
    min_founder_copies = 5), "too_few_loci: 0 loci remain after the outcome founder QC")
  two_loci <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs[1L, ], inst, min_samples = 20L,
    min_loci = 2L, progress = FALSE, min_founder_copies = 5))
  expect_match(two_loci$main_results$error_message, "too_few_loci: 1 loci remain after the outcome founder QC, below min_loci = 2")
})

test_that("the outcome QC runs before fitting, so a locus it drops cannot fail as non-estimable", {
  f <- fqc_fixture()
  inst <- data.frame(exposure = "x", marker = c("a1", "b2"))
  pair <- data.frame(exposure = "x", outcome = "y")
  run <- function(...) suppressWarnings(run_qtl2mr_batch(f$resources, pair, inst, min_samples = 20L, progress = FALSE, ...))
  # Off: founder H is absent from y's own sample at b2, so y's founder effects there are not estimable.
  fail <- run()
  expect_identical(fail$main_results$status, "error")
  expect_match(fail$main_results$error_message, "not estimable at marker b2")
  dropped <- run(on_nonestimable_locus = "drop")
  expect_identical(dropped$locus_results$drop_reason, c("", "nonestimable_outcome"))
  # On: the QC drops b2 for this pair before any fit, under either on_nonestimable_locus.
  for (policy in c("fail", "drop")) {
    qtl2mr:::qtl2mr_cache_reset()
    on <- run(on_nonestimable_locus = policy, min_founder_copies = 5)
    expect_identical(on$main_results$status, "ok")
    expect_identical(on$main_results$n_loci_dropped_nonestimable, 0L)
    expect_identical(on$locus_results$drop_reason, c("", "founder_qc_outcome"))
    # Only x at a1 and b2 and y at a1 were fitted: y at b2 was never even attempted.
    expect_identical(qtl2mr:::qtl2mr_state$computed[["fit"]], 3)
    expect_identical(qtl2mr:::qtl2mr_state$computed[["fit_attempt"]], 3)
    expect_identical(on$main_results$theta, dropped$main_results$theta)
  }
  # A locus that passes the QC (threshold 0: E = 0 is not below 0) is fitted and follows
  # on_nonestimable_locus exactly as with the QC off.
  qtl2mr:::qtl2mr_cache_reset()
  passed_fail <- run(min_founder_copies = 0)
  expect_gt(qtl2mr:::qtl2mr_state$computed[["fit_attempt"]], qtl2mr:::qtl2mr_state$computed[["fit"]])
  expect_identical(passed_fail$main_results$status, "error")
  expect_identical(passed_fail$main_results$error_message, fail$main_results$error_message)
  expect_identical(passed_fail$main_results$n_loci_dropped_founder_qc, 0L)
  passed_drop <- run(on_nonestimable_locus = "drop", min_founder_copies = 0)
  expect_identical(passed_drop$locus_results$drop_reason, c("", "nonestimable_outcome"))
  expect_identical(passed_drop$locus_results$E_min_y, c(passed_drop$locus_results$E_min_y[[1L]], 0))
  expect_identical(passed_drop$main_results$theta, dropped$main_results$theta)
  expect_identical(passed_drop$main_results$n_loci_dropped_nonestimable, 1L)
})

test_that("the QC threshold is part of the batch fingerprint", {
  f <- fqc_fixture()
  inst <- data.frame(exposure = c("x", "x"), marker = c("a1", "b2"))
  pairs <- data.frame(exposure = "x", outcome = c("y", "w"))
  dir <- tempfile("v063-resume-")
  first <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = dir, min_founder_copies = 5))
  again <- suppressWarnings(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = dir, resume = TRUE, min_founder_copies = 5))
  expect_identical(again$run_manifest$n_pairs_resumed, 2L)
  expect_identical(again$locus_results, first$locus_results)
  expect_error(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE, output_dir = dir,
    resume = TRUE, min_founder_copies = 4), "fingerprint changed")
  expect_error(run_qtl2mr_batch(f$resources, pairs, inst, min_samples = 20L, progress = FALSE, output_dir = dir,
    resume = TRUE), "fingerprint changed")
  written <- utils::read.csv(file.path(dir, "batch_locus_results.csv"), stringsAsFactors = FALSE)
  expect_true("E_min_y" %in% names(written))
})
