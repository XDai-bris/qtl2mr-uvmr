# 0.6.1.9002: speed-ups must leave every result unchanged. Each test compares
# with identical(): against the earlier computation, between worker counts,
# between a fresh and a warm session, and between batch and single analyses.

perf_fixture <- function(round_digits = NULL) {
  set.seed(62031)
  n <- 90L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    if (!is.null(round_digits)) a <- round(a, round_digits)
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2", "a3")), `2` = make_chr(c("b1", "b2")), `3` = make_chr(c("c1", "c2")))
  G <- matrix(rnorm(n * 25L), n)
  K <- tcrossprod(G) / 25 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  kinship <- list(`1` = K, `2` = K * 0.98, `3` = K * 1.01)
  sex <- rep(0:1, n / 2L)
  x <- 3 * gp[["1"]][, 1L, "a1"] + 2 * gp[["2"]][, 2L, "b2"] + rnorm(n, sd = 0.3)
  w <- 2 * gp[["1"]][, 3L, "a3"] + 2 * gp[["3"]][, 4L, "c1"] + rnorm(n, sd = 0.3)
  y <- 0.6 * x + 0.3 * w + rnorm(n, sd = 0.3)
  q <- 0.4 * x + rnorm(n, sd = 0.4)
  pheno <- data.frame(x = x, w = w, y = y, q = q, row.names = ids)
  pheno$y[1:6] <- NA_real_
  pheno$x[7:8] <- NA_real_
  # q is measured in one sex only, so sex is constant on its complete samples
  # and every q context warns about a dependent covariate.
  pheno$q[sex == 1L] <- NA_real_
  covar <- data.frame(sex = sex, age = round(runif(n, 50, 90)), row.names = ids)
  map <- data.frame(marker = c("a1", "a2", "a3", "b1", "b2", "c1", "c2"), chr = c("1", "1", "1", "2", "2", "3", "3"),
    pos = c(1, 2, 3, 1, 2, 1, 2))
  scan <- matrix(c(8, 4.5, 1, 1, 6, 1, 1,
                   1, 1, 7, 4.2, 1, 6, 1,
                   5, 1, 4.5, 1, 4.4, 1, 1,
                   4.6, 1, 1, 1, 5.1, 1, 1), ncol = 4L, dimnames = list(map$marker, c("x", "w", "y", "q")))
  resources <- suppressWarnings(qtl2mr_prepare_resources(gp, kinship, pheno, covar, map))
  list(gp = gp, kinship = kinship, pheno = pheno, covar = covar, map = map, scan = scan, resources = resources,
    pairs = data.frame(exposure = c("x", "w", "x", "y", "q", "x"), outcome = c("y", "y", "w", "q", "x", "q")))
}

# qtl2mr_allele_prob_markers() of 0.6.1.9001 (c5a0e52), kept verbatim as the reference.
old_allele_prob_markers <- function(allele_probs, markers_to_validate = NULL,
                                    probability_policy = c("strict", "normalize_rounding")) {
  probability_policy <- match.arg(probability_policy)
  qtl2mr_validate_ids <- qtl2mr:::qtl2mr_validate_ids
  `%||%` <- function(x, y) if (is.null(x)) y else x
  if (!is.list(allele_probs) || !length(allele_probs)) stop("allele_probs must be a non-empty qtl2-style list of 3D arrays", call. = FALSE)
  qtl2mr_validate_ids(names(allele_probs), "allele_probs chromosome names")
  chr_names <- names(allele_probs) %||% as.character(seq_along(allele_probs))
  rows <- list()
  adjustments <- list()
  for (i in seq_along(allele_probs)) {
    arr <- allele_probs[[i]]
    if (length(dim(arr)) != 3L) stop("Each allele_probs element must be a sample x founder x marker array", call. = FALSE)
    dn <- dimnames(arr)
    for (d in seq_len(3L)) qtl2mr_validate_ids(dn[[d]], paste0("allele_probs dimension ", d))
    if (!is.numeric(arr)) stop("allele_probs must be numeric", call. = FALSE)
    first_ids <- dimnames(allele_probs[[1L]])[[1L]]
    if (i > 1L && !setequal(dn[[1L]], first_ids)) stop("allele_probs chromosomes must hold the same samples", call. = FALSE)
    checked <- if (is.null(markers_to_validate)) seq_len(dim(arr)[3L]) else which(dn[[3L]] %in% markers_to_validate)
    for (first in seq.int(1L, max(length(checked), 1L), by = 128L)) {
      if (!length(checked)) break
      chunk <- arr[, , checked[seq.int(first, min(first + 127L, length(checked)))], drop = FALSE]
      if (any(!is.finite(chunk))) stop("allele_probs must contain finite founder allele probabilities", call. = FALSE)
      if (any(chunk < -1e-8 | chunk > 1 + 1e-8)) stop("range", call. = FALSE)
      deviations <- abs(apply(chunk, c(1L, 3L), sum) - 1)
      sum_tolerance <- if (probability_policy == "strict") 1e-6 else dim(arr)[2L] * 5e-5
      if (any(deviations > sum_tolerance + 16 * .Machine$double.eps)) stop("sum", call. = FALSE)
      adjustments[[length(adjustments) + 1L]] <- data.frame(chr = chr_names[[i]],
        marker = dimnames(chunk)[[3L]], n_normalized = colSums(deviations > 1e-6 + 16 * .Machine$double.eps),
        max_abs_sum_deviation = apply(deviations, 2L, max), stringsAsFactors = FALSE)
    }
    if (i > 1L && !identical(dn[[2L]], dimnames(allele_probs[[1L]])[[2L]])) stop("founders", call. = FALSE)
    markers <- dn[[3L]]
    rows[[length(rows) + 1L]] <- data.frame(marker = markers, chr = chr_names[[i]], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  qtl2mr_validate_ids(out$marker, "allele_probs marker names")
  attr(out, "probability_normalization") <- do.call(rbind, adjustments)
  out
}

perf_random_probs <- function(n, markers, digits = NULL, seed = 1L, names_dims = FALSE) {
  set.seed(seed)
  a <- array(rexp(n * 8L * markers), c(n, 8L, markers),
    list(paste0("i", seq_len(n)), LETTERS[1:8], paste0("m", seed, "_", seq_len(markers))))
  for (j in seq_len(markers)) {
    m <- matrix(a[, , j], nrow = n)
    a[, , j] <- m / rowSums(m)
  }
  if (!is.null(digits)) a <- round(a, digits)
  if (names_dims) names(dimnames(a)) <- c("ind", "gen", "mar")
  a
}

# Results without timestamps and output paths, which differ between runs.
perf_strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    x <- x[setdiff(names(x) %||% character(), c("created_at", "output_dir", "paths", "output_manifest", "figure_manifest"))]
    x[] <- lapply(x, perf_strip)
  } else if (is.data.frame(x)) {
    x <- x[, setdiff(names(x), c("runtime_sec", "wall_sec", "mean_wall_sec_per_pair", "run_id", "created_at")), drop = FALSE]
  }
  x
}
`%||%` <- function(x, y) if (is.null(x)) y else x

perf_capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

perf_same <- function(a, b) identical(perf_strip(a), perf_strip(b), ignore.environment = TRUE)

test_that("the probability check gives the 0.6.1.9001 marker table and record", {
  cases <- list(
    list(`1` = perf_random_probs(30, 3), `2` = perf_random_probs(30, 1, seed = 2)),
    list(`1` = perf_random_probs(25, 300, digits = 4, seed = 3)),
    list(`1` = perf_random_probs(12, 129, digits = 5, seed = 4), `5` = perf_random_probs(12, 2, digits = 6, seed = 5)),
    list(`1` = perf_random_probs(1, 4, digits = 4, seed = 6)),
    list(`1` = perf_random_probs(20, 5, digits = 4, seed = 7, names_dims = TRUE))
  )
  cases[[3]][["5"]] <- cases[[3]][["5"]][dimnames(cases[[3]][["1"]])[[1L]], , , drop = FALSE]
  for (probs in cases) {
    for (policy in c("strict", "normalize_rounding")) {
      for (markers in list(NULL, dimnames(probs[[1L]])[[3L]][1L])) {
        old <- tryCatch(old_allele_prob_markers(probs, markers, policy), error = function(e) conditionMessage(e))
        new <- tryCatch(qtl2mr:::qtl2mr_allele_prob_markers(probs, markers, policy), error = function(e) conditionMessage(e))
        if (is.character(old)) expect_true(is.character(new)) else expect_identical(new, old)
      }
    }
  }
  f <- perf_fixture(round_digits = 4)
  w_new <- perf_capture(qtl2mr_prepare_resources(f$gp, f$kinship, f$pheno, f$covar, f$map))
  ref <- list(allele_probs = f$gp, kinship = f$kinship, pheno = as.data.frame(f$pheno), covar = as.data.frame(f$covar),
    marker_map = qtl2mr:::qtl2mr_marker_map_table(f$map), sample_ids = dimnames(f$gp[[1L]])[[1L]],
    metadata = list(), probability_policy = "normalize_rounding")
  ref$marker_map$chr <- as.character(ref$marker_map$chr)
  class(ref) <- c("qtl2mr_resources", "list")
  validate_qtl2mr_input(ref)
  w_ref <- perf_capture(qtl2mr:::qtl2mr_record_probability_policy(ref, warn = TRUE,
    marker_info = old_allele_prob_markers(ref$allele_probs, NULL, "normalize_rounding")))
  ref <- w_ref$value
  ref$sample_alignment <- qtl2mr:::qtl2mr_input_alignment(ref)
  expect_identical(w_new$value, ref)
  expect_identical(w_new$warnings, w_ref$warnings)
  expect_gt(nrow(ref$metadata$probability_normalization$adjustments), 0L)
})

test_that("a kinship that passed validation does not excuse a different one", {
  f <- perf_fixture()
  inst <- data.frame(exposure = "x", marker = "a1")
  expect_silent(validate_qtl2mr_input(qtl2mr_pair_input(f$resources, inst, "x", "y")))
  bad <- f$resources
  bad$kinship[["1"]][1L, 2L] <- bad$kinship[["1"]][1L, 2L] + 1
  expect_error(validate_qtl2mr_input(qtl2mr:::qtl2mr_resources_input(bad, "a1", "x", "y")), "symmetric")
})

test_that("only the requested markers are copied when probabilities are normalized", {
  f <- perf_fixture(round_digits = 4)
  input <- suppressWarnings(qtl2mr_pair_input(f$resources, data.frame(exposure = "x", marker = c("a1", "b2")), "x", "y"))
  copy <- suppressWarnings(qtl2mr:::qtl2mr_prepared_probability_copy(input, c("a1", "b2")))
  adjusted <- unique(copy$metadata$probability_normalization$adjustments$chr)
  expect_true(length(adjusted) > 0L)
  for (chr in adjusted) {
    expect_identical(dimnames(copy$allele_probs[[chr]])[[3L]], intersect(dimnames(f$gp[[chr]])[[3L]], c("a1", "b2")))
  }
  expect_identical(input$allele_probs, f$resources$allele_probs)
})

test_that("Q-only refits search exactly as the full fit", {
  f <- perf_fixture()
  inst <- data.frame(exposure = "x", marker = c("a1", "a2", "b2", "c1"))
  pilot <- qtl2mr:::qtl2mr_fit_prepared(qtl2mr_pair_input(f$resources, inst, "x", "y"), "x", "y",
    qtl2mr:::qtl2mr_fit_settings("fail", "diagnostic_only", "chisq", 0L, 1L, 0L, 1L, TRUE, "star"))
  blocks <- pilot$locus_blocks
  for (k in seq_along(blocks)) {
    b <- blocks[seq_len(k)]
    for (cross in c(TRUE, FALSE)) {
      full <- qtl2mr:::qtl2mr_so_fit(b, cross, 1e-10)
      expect_identical(qtl2mr:::qtl2mr_so_fit(b, cross, 1e-10, q_only = TRUE), list(Q = full$Q))
    }
  }
})

test_that("ordered worker calls signal warnings, messages and the first error as sequential calls do", {
  skip_on_os("windows")
  fun <- function(i) {
    warning("w", i, call. = FALSE)
    if (i == 4L) message("m", i)
    if (i == 5L) stop("stop at ", i, call. = FALSE)
    i^2
  }
  run <- function(workers) {
    old <- qtl2mr:::qtl2mr_set_workers(workers)
    on.exit(qtl2mr:::qtl2mr_set_workers(old), add = TRUE)
    seen <- character()
    value <- tryCatch(withCallingHandlers(qtl2mr:::qtl2mr_ordered_lapply(as.list(1:7), fun),
      warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { seen <<- c(seen, conditionMessage(m)); invokeRestart("muffleMessage") }),
      error = function(e) conditionMessage(e))
    list(value = value, seen = seen)
  }
  expect_identical(run(2L), run(1L))
  ok <- function(workers) {
    old <- qtl2mr:::qtl2mr_set_workers(workers)
    on.exit(qtl2mr:::qtl2mr_set_workers(old), add = TRUE)
    qtl2mr:::qtl2mr_ordered_lapply(stats::setNames(as.list(1:5), letters[1:5]), function(i) i + 0.5)
  }
  expect_identical(ok(3L), ok(1L))
})

test_that("analysis results are identical with one and two workers", {
  skip_on_os("windows")
  f <- perf_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  out <- tempfile("v062-workers-")
  run <- function(exposure, outcome, workers) {
    qtl2mr:::qtl2mr_cache_reset()
    set.seed(11)
    before <- .Random.seed
    r <- perf_capture(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, exposure, outcome),
      output_dir = out, overwrite = TRUE, make_html = FALSE, make_figures = FALSE, verbose = FALSE,
      heterogeneity_nsim = 40L, diagnostic_bootstrap_B = 40L, n_workers = workers))
    r$seed_kept <- identical(.Random.seed, before)
    r
  }
  for (pair in list(c("x", "y"), c("x", "q"))) {
    one <- run(pair[[1L]], pair[[2L]], 1L)
    two <- run(pair[[1L]], pair[[2L]], 2L)
    expect_true(perf_same(one$value, two$value))
    expect_identical(one$warnings, two$warnings)
    expect_true(one$seed_kept && two$seed_kept)
    expect_identical(one$value$pilot$fit$objective(0.37), two$value$pilot$fit$objective(0.37))
    expect_true(is.finite(one$value$pilot$summary$Q_p_sim))
    expect_true(nrow(one$value$loo_mr) > 1L)
  }
  expect_true(any(grepl("Dropping linearly dependent covariate", one$warnings)))
  expect_error(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, "x", "y"), n_workers = 0L), "n_workers")
})

test_that("single analyses do not depend on what the session has cached", {
  f <- perf_fixture()
  run <- function(inst) perf_capture(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, "x", "q"),
    output_dir = tempfile("v062-warm-"), make_html = FALSE, make_figures = FALSE, verbose = FALSE,
    heterogeneity_nsim = 20L, diagnostic_bootstrap_B = 20L, run_instrument_sensitivity = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  cold <- run(inst)
  qtl2mr:::qtl2mr_cache_reset()
  sel <- perf_capture(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  expect_identical(sel$value, inst)
  after_selection <- run(inst)
  again <- run(inst)
  expect_true(perf_same(after_selection$value, cold$value))
  expect_true(perf_same(again$value, cold$value))
  expect_identical(after_selection$warnings, cold$warnings)
  expect_identical(again$warnings, cold$warnings)
  # The second selection in a warm session reuses every fit and warns the same.
  expect_identical(perf_capture(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), sel)
})

perf_batch <- function(f, inst, ...) {
  perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE, ...))
}

test_that("batch results do not depend on the session cache, worker count or wave size", {
  skip_on_os("windows")
  f <- perf_fixture(round_digits = 5)
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  cold <- perf_batch(f, inst)
  warm <- perf_batch(f, inst)
  expect_true(perf_same(warm$value, cold$value))
  expect_identical(warm$warnings, cold$warnings)
  expect_identical(warm$value$run_manifest[c("context_cache_hits", "context_cache_misses", "lmm_fit_cache_hits",
    "lmm_fit_cache_misses")], cold$value$run_manifest[c("context_cache_hits", "context_cache_misses",
    "lmm_fit_cache_hits", "lmm_fit_cache_misses")])
  expect_true(any(grepl("Dropping linearly dependent covariate", cold$warnings)))
  qtl2mr:::qtl2mr_cache_reset()
  invisible(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE)))
  after_selection <- perf_batch(f, inst)
  expect_true(perf_same(after_selection$value, cold$value))
  expect_identical(after_selection$warnings, cold$warnings)
  # Small waves: planning in several steps, with releases in between.
  assign("budget", 2e6, envir = qtl2mr:::qtl2mr_state)
  small <- tryCatch(perf_batch(f, inst), finally = assign("budget", NULL, envir = qtl2mr:::qtl2mr_state))
  expect_true(perf_same(small$value, cold$value))
  expect_identical(small$warnings, cold$warnings)
  # Two workers: groups run in workers (whose warnings are not returned, as before).
  qtl2mr:::qtl2mr_cache_reset()
  forked_old_layout <- perf_batch(f, inst, n_workers = 2L)
  serial <- cold$value
  expect_true(perf_same(forked_old_layout$value[setdiff(names(serial), "run_manifest")],
    serial[setdiff(names(serial), "run_manifest")]))
  # One exposure group with two workers: trait fits are planned in parallel, pairs fitted in order.
  one_group <- f$pairs[f$pairs$exposure == "x", ]
  seq_one <- perf_capture(run_qtl2mr_batch(f$resources, one_group, inst, min_samples = 20L, progress = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  par_one <- perf_capture(run_qtl2mr_batch(f$resources, one_group, inst, min_samples = 20L, progress = FALSE, n_workers = 2L))
  expect_true(perf_same(par_one$value[setdiff(names(par_one$value), "run_manifest")],
    seq_one$value[setdiff(names(seq_one$value), "run_manifest")]))
  keep <- c("parallel_backend", "context_cache_hits", "context_cache_misses", "lmm_fit_cache_hits", "lmm_fit_cache_misses")
  expect_identical(par_one$value$run_manifest[keep], seq_one$value$run_manifest[keep])
  expect_identical(par_one$warnings, seq_one$warnings)
  # Nothing the batch planned stays cached after it.
  reg <- qtl2mr:::qtl2mr_registry()
  traits <- unique(c(one_group$exposure, one_group$outcome))
  left <- Filter(function(id) reg$entries[[id]]$info$trait %in% traits, ls(reg$entries))
  expect_length(left, 0L)
})

test_that("a batch in a fresh R process gives the same results", {
  skip_on_cran()
  skip_on_os("windows")
  f <- perf_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE))
  here <- perf_batch(f, inst)$value
  here_again <- perf_batch(f, inst)$value
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  saveRDS(list(resources = f$resources, pairs = f$pairs, instruments = inst), input)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "suppressPackageStartupMessages(library(qtl2mr))",
    "d <- readRDS(args[[1L]])",
    "b <- suppressWarnings(run_qtl2mr_batch(d$resources, d$pairs, d$instruments, min_samples = 20L, progress = FALSE))",
    "saveRDS(b, args[[2L]])"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), shQuote(input), shQuote(output)),
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)), stdout = FALSE, stderr = FALSE)
  skip_if(status != 0L || !file.exists(output), "child R process unavailable")
  fresh <- readRDS(output)
  expect_true(perf_same(here, fresh))
  expect_true(perf_same(here_again, fresh))
})

test_that("batch rows equal single analyses under both sample policies", {
  f <- perf_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE))
  cols <- qtl2mr:::qtl2mr_batch_summary_columns
  results <- list()
  for (policy in c("own_complete", "pair_common_complete", "own_complete")) {
    batch <- perf_batch(f, inst, sample_policy = policy)$value
    expect_identical(batch$main_results$status, rep("ok", nrow(f$pairs)))
    for (i in seq_len(nrow(f$pairs))) {
      single <- suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, f$pairs$exposure[[i]],
        f$pairs$outcome[[i]], sample_policy = policy), output_dir = tempfile("v062-policy-"), make_html = FALSE,
        make_figures = FALSE, verbose = FALSE, heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L,
        run_loo_mr = FALSE, run_instrument_sensitivity = FALSE))
      row <- batch$main_results[i, cols]
      ref <- single$pilot$summary[1L, cols]
      rownames(row) <- rownames(ref) <- NULL
      expect_identical(row, ref)
    }
    results[[length(results) + 1L]] <- batch
  }
  own <- results[[1L]]$main_results
  common <- results[[2L]]$main_results
  expect_true(perf_same(results[[3L]], results[[1L]]))
  q_pairs <- f$pairs$exposure == "q" | f$pairs$outcome == "q"
  expect_identical(common$n_x_median, common$n_y_median)
  expect_true(all(common$n_x_median <= pmin(own$n_x_median, own$n_y_median)))
  expect_true(any(own$n_x_median[q_pairs] != common$n_x_median[q_pairs]))
  expect_false(identical(own$theta, common$theta))
  expect_identical(results[[2L]]$config$sample_policy, "pair_common_complete")
  expect_null(results[[1L]]$config$sample_policy)
  input <- qtl2mr_pair_input(f$resources, inst, "x", "q", sample_policy = "pair_common_complete")
  expect_identical(input$exposure_samples, qtl2mr:::qtl2mr_complete_samples(f$resources, c("x", "q"), c("1", "2")))
  expect_identical(input$outcome_samples, input$exposure_samples)
  expect_null(qtl2mr_pair_input(f$resources, inst, "x", "q")$exposure_samples)
  expect_error(qtl2mr_pair_input(f$resources, inst, "x", "q", sample_policy = "batch_common_complete"), "should be one of")
  # q has 45 complete samples of its own, fewer together with any partner.
  too_few <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 45L, progress = FALSE,
    sample_policy = "pair_common_complete"))$value
  expect_true(all(grepl("too_few_samples: .* samples complete for both traits",
    too_few$main_results$error_message[q_pairs])))
  expect_identical(too_few$main_results$status[!q_pairs], rep("ok", sum(!q_pairs)))
  enough <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 45L, progress = FALSE))$value
  expect_identical(enough$main_results$status, rep("ok", nrow(f$pairs)))
})

test_that("per-locus blocks are written on request and refit to the batch estimates", {
  f <- perf_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE))
  plain_dir <- tempfile("v062-plain-")
  plain <- perf_batch(f, inst, output_dir = plain_dir)$value
  expect_false(dir.exists(file.path(plain_dir, "batch_locus_blocks")))
  expect_false(file.exists(file.path(plain_dir, "batch_locus_blocks_manifest.csv")))
  expect_null(plain$locus_blocks_manifest)
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, write_locus_blocks = TRUE), "requires output_dir")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, write_locus_blocks = NA), "TRUE or FALSE")

  dir <- tempfile("v062-blocks-")
  for (policy in c("own_complete", "pair_common_complete")) {
    batch <- perf_batch(f, inst, output_dir = dir, write_locus_blocks = TRUE, sample_policy = policy,
      overwrite = TRUE)$value
    if (identical(policy, "own_complete")) {
      expect_true(perf_same(batch[setdiff(names(batch), c("config", "locus_blocks_manifest"))],
        plain[setdiff(names(plain), "config")]))
    }
    manifest <- utils::read.csv(file.path(dir, "batch_locus_blocks_manifest.csv"), stringsAsFactors = FALSE)
    expect_identical(manifest$pair_id, batch$main_results$pair_id)
    expect_identical(batch$locus_blocks_manifest$file, manifest$file)
    main <- batch$main_results
    loci <- batch$locus_results
    for (i in seq_len(nrow(manifest))) {
      blocks <- readRDS(file.path(dir, manifest$file[[i]]))
      expect_identical(attr(blocks, "qtl2mr_locus_blocks")$sample_policy, policy)
      expect_true(all(c("marker", "chr", "pos", "beta_x", "beta_y", "Sigma_x", "Sigma_y", "Sigma_xy", "F_X_block",
        "F_Y_block", "partial_R2_X_block", "partial_R2_Y_block", "steiger_pass", "df_x", "df_y", "n_x", "n_y",
        "n_overlap", "hsq_x", "hsq_y") %in% names(blocks[[1L]])))
      refit <- fit_multilocus_founder_mr(blocks, allow_cross_cov = TRUE, ridge = 1e-10, heterogeneity_method = "chisq")
      j <- match(manifest$pair_id[[i]], main$pair_id)
      expect_identical(refit$theta, main$theta[[j]])
      expect_identical(refit$Q, main$Q[[j]])
      expect_identical(refit$Q_p_chisq, main$Q_p_chisq[[j]])
      expect_identical(refit$independence$theta, main$theta_independence[[j]])
      l <- loci[loci$pair_id == manifest$pair_id[[i]], ]
      expect_identical(unname(vapply(blocks, `[[`, numeric(1), "F_X_block")), l$F_X_block)
      expect_identical(unname(vapply(blocks, `[[`, logical(1), "steiger_pass")), l$steiger_pass)
      expect_identical(unname(vapply(blocks, `[[`, integer(1), "n_y")), l$n_y)
      expect_identical(as.numeric(blocks[[1L]]$fit_x$df), blocks[[1L]]$df_x)
    }
  }
  # Resume completes the same run and rewrites the block files from the pair records.
  partial <- f$pairs[1:3, ]
  rdir <- tempfile("v062-resume-")
  first <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, write_locus_blocks = TRUE))$value
  unlink(file.path(rdir, "batch_locus_blocks", "pair_000002.rds"))
  unlink(file.path(rdir, "batch_pair_records", "pair_000004.rds"))
  resumed <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, write_locus_blocks = TRUE, resume = TRUE))$value
  expect_identical(resumed$run_manifest$n_pairs_resumed, nrow(f$pairs) - 1L)
  expect_true(perf_same(resumed[setdiff(names(resumed), "run_manifest")], first[setdiff(names(first), "run_manifest")]))
  expect_true(file.exists(file.path(rdir, "batch_locus_blocks", "pair_000002.rds")))
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, resume = TRUE), "fingerprint changed")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE, output_dir = rdir),
    "already exists")
  cleared <- perf_capture(run_qtl2mr_batch(f$resources, partial, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, overwrite = TRUE))$value
  expect_false(dir.exists(file.path(rdir, "batch_locus_blocks")))
  expect_false(file.exists(file.path(rdir, "batch_locus_blocks_manifest.csv")))
  expect_identical(cleared$main_results$theta, first$main_results$theta[1:3])
})

test_that("the session cache stays within its budget and keeps pinned entries", {
  qtl2mr:::qtl2mr_cache_reset()
  on.exit(qtl2mr:::qtl2mr_cache_reset(), add = TRUE)
  assign("budget", 5e5, envir = qtl2mr:::qtl2mr_state)
  value <- function(i) rep(as.numeric(i), 1e4)
  qtl2mr:::qtl2mr_registry_put("fit", "k1", value(1), info = list(trait = "a", chr = "1", marker = "m1"))
  qtl2mr:::qtl2mr_registry_pin("fit:k1")
  for (i in 2:12) qtl2mr:::qtl2mr_registry_put("fit", paste0("k", i), value(i), info = list(trait = "a", chr = "1", marker = paste0("m", i)))
  reg <- qtl2mr:::qtl2mr_registry()
  expect_lte(reg$bytes, 5e5)
  expect_false(is.null(qtl2mr:::qtl2mr_registry_get("fit", "k1")))
  expect_false(is.null(qtl2mr:::qtl2mr_registry_get("fit", "k12")))
  expect_null(qtl2mr:::qtl2mr_registry_get("fit", "k2"))
  expect_equal(reg$bytes, sum(vapply(ls(reg$entries), function(id) reg$entries[[id]]$bytes, numeric(1))))
})
