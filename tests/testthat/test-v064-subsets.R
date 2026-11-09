# 0.6.1.9004: instrument subsets. One run fits each pair on its instruments
# (the maximal set) and refits the multilocus stage on every requested subset
# of the fitted loci, with the run's settings and without refitting the LMMs.

# Three chromosomes with LOCO kinship; exposure x has loci of graded strength,
# y has a direct effect at b2 (which fails Steiger for x -> y), and w is a
# second exposure. y and w have missing values, so own and pair-common samples
# differ.
v064_fixture <- function() {
  set.seed(64001)
  n <- 120L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2")), `2` = make_chr(c("b1", "b2")), `3` = make_chr(c("c1", "c2")))
  G <- matrix(rnorm(n * 25L), n)
  K <- tcrossprod(G) / 25 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  kinship <- list(`1` = K, `2` = K * 0.98, `3` = K * 1.01)
  x <- 5 * gp[["1"]][, 1L, "a1"] + 2.5 * gp[["2"]][, 2L, "b1"] + 1.2 * gp[["3"]][, 3L, "c1"] +
    0.8 * gp[["2"]][, 5L, "b2"] + rnorm(n, sd = 0.4)
  w <- 4 * gp[["3"]][, 6L, "c2"] + 2 * gp[["1"]][, 7L, "a2"] + rnorm(n, sd = 0.4)
  y <- 0.5 * x + 0.3 * w + 4 * gp[["2"]][, 4L, "b2"] + rnorm(n, sd = 0.4)
  pheno <- data.frame(x = x, w = w, y = y, row.names = ids)
  pheno$y[1:4] <- NA_real_
  pheno$w[5:7] <- NA_real_
  covar <- data.frame(sex = rep(0:1, n / 2L), row.names = ids)
  map <- data.frame(marker = c("a1", "a2", "b1", "b2", "c1", "c2"), chr = c("1", "1", "2", "2", "3", "3"),
    pos = c(1, 30, 1, 30, 1, 30))
  resources <- suppressWarnings(qtl2mr_prepare_resources(gp, kinship, pheno, covar, map))
  instruments <- data.frame(exposure = c(rep("x", 6L), "w", "w", "w"),
    marker = c("a1", "b1", "c1", "b2", "a2", "c2", "c2", "a2", "b1"), stringsAsFactors = FALSE)
  list(resources = resources, instruments = instruments,
    pairs = data.frame(exposure = c("x", "x", "w"), outcome = c("y", "w", "y"), stringsAsFactors = FALSE))
}

# F5 (one locus for x), F2 (four), F2_SF (Steiger within F2), zero (equals the
# maximal set), huge (no locus), SF (Steiger alone; equals F2_SF here).
v064_subsets <- list(F5 = list(F_X_min = 5), F2 = list(F_X_min = 2), F2_SF = list(F_X_min = 2, steiger = TRUE),
  zero = list(F_X_min = 0), huge = list(F_X_min = 1e6), SF = list(steiger = TRUE))

v064_batch <- function(f, ..., pairs = f$pairs, instruments = f$instruments) {
  suppressWarnings(run_qtl2mr_batch(f$resources, pairs, instruments, min_samples = 20L, progress = FALSE, ...))
}

v064_drop_runtime <- function(d) d[, setdiff(names(d), "runtime_sec"), drop = FALSE]

# Rows are compared with identical() column by column. Integer and double
# storage of the same value are not distinguished: rbind() of rows promotes an
# integer median (odd number of loci) to double when another row has an even
# number, as for main_results since 0.6.1.9001.
v064_same_values <- function(a, b, cols) {
  a <- as.list(a[1L, cols, drop = FALSE])
  b <- as.list(b[1L, cols, drop = FALSE])
  num <- vapply(a, is.numeric, logical(1)) & vapply(b, is.numeric, logical(1))
  a[num] <- lapply(a[num], as.double)
  b[num] <- lapply(b[num], as.double)
  ok <- vapply(cols, function(k) identical(unname(a[[k]]), unname(b[[k]])), logical(1))
  if (!all(ok)) message("differing columns: ", paste(cols[!ok], collapse = ", "))
  all(ok)
}

# Columns of a subset row that come from the fit: everything but the subset's
# own name, rule text and copy marker (status and error_message are "ok" and ""
# in both; n_loci is the subset's locus count and the separate run's n_loci).
v064_fit_columns <- function(table) setdiff(names(table), c("pair_id", "subset", "rules", "same_loci_as"))

test_that("off by default: no table, file, column or config entry; invalid specifications are refused", {
  f <- v064_fixture()
  dir_off <- tempfile("v064-off-")
  off <- v064_batch(f, output_dir = dir_off)
  expect_identical(v064_drop_runtime(v064_batch(f, subsets = NULL)$main_results), v064_drop_runtime(off$main_results))
  expect_null(off$subset_results)
  expect_false("subset_results" %in% names(off))
  expect_false(any(startsWith(names(off$locus_results), "in_subset_")))
  expect_null(off$config$subsets)
  expect_false(file.exists(file.path(dir_off, "batch_subset_results.csv")))
  expect_false(grepl("subsets", paste(readLines(file.path(dir_off, "batch_config.json")), collapse = "")))
  single <- suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, f$instruments, "x", "y"),
    output_dir = tempfile("v064-single-off-"), make_html = FALSE, make_figures = FALSE, progress = FALSE,
    subsets = NULL))
  expect_false("subset_results" %in% names(single))
  expect_false(file.exists(file.path(single$output_dir, "data", "subset_results.csv")))

  bad <- list(
    "F2", list(), list(list(F_X_min = 1)), stats::setNames(list(list(F_X_min = 1)), ""),
    list(`a b` = list(F_X_min = 1)), list(`F-2` = list(F_X_min = 1)), list(all = list(F_X_min = 1)),
    list(a = list(F_X_min = 1), a = list(F_X_min = 2)), list(a = list()), list(a = c(F_X_min = 2)),
    list(a = list(foo = 1)), list(a = list(ld = 0.2)), list(a = list(2)),
    list(a = list(F_X_min = 1, F_X_min = 2)), list(a = list(F_X_min = -1)), list(a = list(F_X_min = "2")),
    list(a = list(F_X_min = NA_real_)), list(a = list(F_X_min = c(1, 2))), list(a = list(F_X_min = Inf)),
    list(a = list(F_X_min = TRUE)), list(a = list(steiger = NA)), list(a = list(steiger = "yes")),
    list(a = list(steiger = 1)), list(a = list(steiger = c(TRUE, TRUE))), data.frame(F_X_min = 1)
  )
  for (spec in bad) {
    expect_error(qtl2mr:::qtl2mr_check_subsets(spec), "subset")
    # Refused before anything is fitted or written.
    dir <- tempfile("v064-bad-")
    expect_error(run_qtl2mr_batch(f$resources, f$pairs, f$instruments, output_dir = dir, progress = FALSE,
      subsets = spec), "subset")
    expect_false(dir.exists(dir))
    dir <- tempfile("v064-bad-single-")
    expect_error(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, f$instruments, "x", "y"), output_dir = dir,
      progress = FALSE, subsets = spec), "subset")
    expect_false(dir.exists(dir))
  }
  # Accepted forms are normalised: rules in the order F_X_min, steiger; whole numbers as double.
  expect_identical(qtl2mr:::qtl2mr_check_subsets(list(S.1_b = list(steiger = TRUE, F_X_min = 2L))),
    list(S.1_b = list(F_X_min = 2, steiger = TRUE)))
  expect_identical(qtl2mr:::qtl2mr_subset_rules_text(list(steiger = FALSE)), "none")
})

test_that("membership follows the pair's own F_X_block and steiger_pass; subsets nest; too few loci", {
  f <- v064_fixture()
  b <- v064_batch(f, subsets = v064_subsets)
  loci <- b$locus_results
  sub <- b$subset_results
  expect_identical(nrow(sub), nrow(f$pairs) * length(v064_subsets))
  expect_identical(sub$pair_id, rep(b$main_results$pair_id, each = length(v064_subsets)))
  expect_identical(sub$subset, rep(names(v064_subsets), times = nrow(f$pairs)))
  expect_identical(utils::tail(names(loci), length(v064_subsets)), paste0("in_subset_", names(v064_subsets)))
  expect_identical(names(sub), qtl2mr:::qtl2mr_batch_subset_columns())
  for (nm in names(v064_subsets)) {
    r <- v064_subsets[[nm]]
    hand <- (is.null(r$F_X_min) | (is.finite(loci$F_X_block) & loci$F_X_block >= (r$F_X_min %||% 0))) &
      (!isTRUE(r$steiger) | loci$steiger_pass %in% TRUE)
    expect_identical(loci[[paste0("in_subset_", nm)]], hand)
    counts <- tapply(hand, factor(loci$pair_id, levels = b$main_results$pair_id), sum)
    expect_identical(sub$n_loci[sub$subset == nm], as.integer(counts))
  }
  # Nesting: F5 within F2; F2_SF within F2; Steiger drops b2 (F_Y >> F_X) for x -> y.
  expect_true(all(!loci$in_subset_F5 | loci$in_subset_F2))
  expect_true(all(!loci$in_subset_F2_SF | loci$in_subset_F2))
  expect_false(loci$in_subset_SF[loci$pair_id == "x__to__y" & loci$marker == "b2"])
  expect_identical(sub$n_loci[sub$pair_id == "x__to__y"], c(1L, 4L, 3L, 6L, 0L, 3L))
  expect_identical(sub$rules[1:6], c("F_X_min>=5", "F_X_min>=2", "F_X_min>=2;steiger", "F_X_min>=0",
    "F_X_min>=1e+06", "steiger"))
  # Copies: zero has every fitted locus (all); SF equals F2_SF for the x pairs; for w -> y F2 and F2_SF equal F5.
  expect_identical(sub$same_loci_as[sub$subset == "zero"], rep("all", 3L))
  expect_identical(sub$same_loci_as[sub$pair_id == "w__to__y"], c(NA, "F5", "F5", "all", NA, "F5"))
  expect_identical(sub$same_loci_as[sub$pair_id == "x__to__y"], c(NA, NA, NA, "all", NA, "F2_SF"))
  # The maximal-set copy carries the main row's estimates.
  zero <- sub[sub$subset == "zero", ]
  for (i in seq_len(nrow(b$main_results))) {
    expect_true(v064_same_values(zero[i, ], b$main_results[i, ], v064_fit_columns(sub)))
  }
  # Too few loci: NA estimates, a message, no errors row.
  huge <- sub[sub$subset == "huge", ]
  expect_identical(huge$status, rep("too_few_loci", 3L))
  expect_match(huge$error_message, "too_few_loci: 0 loci pass the subset rules, below min_loci = 1")
  expect_true(all(is.na(huge$theta) & is.na(huge$se) & is.na(huge$Q)))
  expect_identical(nrow(b$errors), 0L)
  expect_true(all(sub$status[sub$subset != "huge"] == "ok"))
  # min_loci applies to subsets too.
  two <- v064_batch(f, min_loci = 2L, subsets = v064_subsets["F5"])
  expect_identical(two$subset_results$status, c("too_few_loci", "too_few_loci", "ok"))
  expect_match(two$subset_results$error_message[[1L]], "1 loci pass the subset rules, below min_loci = 2")
  # Config records the normalised specification.
  expect_identical(b$config$subsets, v064_subsets)
})

test_that("Steiger NA is dropped by the steiger rule and kept by F rules", {
  first <- data.frame(marker = c("m1", "m2", "m3", "m4"), F_X_block = c(3, NA, 2.5, Inf),
    steiger_pass = c(TRUE, NA, NA, FALSE))
  members <- qtl2mr:::qtl2mr_subset_members
  expect_identical(members(first, list(steiger = TRUE)), c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(members(first, list(F_X_min = 2)), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(members(first, list(F_X_min = 2, steiger = TRUE)), c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(members(first, list(steiger = FALSE)), rep(TRUE, 4L))
  expect_identical(members(first, list(F_X_min = 0)), c(TRUE, FALSE, TRUE, FALSE))
  # On a fitted pair: a locus whose steiger_pass is NA leaves the Steiger subset
  # but stays in an F subset, and the refit equals a refit on the remaining loci.
  f <- v064_fixture()
  settings <- qtl2mr:::qtl2mr_fit_settings("fail", "diagnostic_only", "chisq", 0L, 1L, 0L, 1L, TRUE, "star")
  pilot <- suppressWarnings(qtl2mr:::qtl2mr_fit_prepared(qtl2mr_pair_input(f$resources, f$instruments, "x", "w"),
    "x", "w", settings))
  expect_true(pilot$locus_first_stage$steiger_pass[pilot$locus_first_stage$marker == "a1"])
  pilot$locus_first_stage$steiger_pass[pilot$locus_first_stage$marker == "a1"] <- NA
  out <- qtl2mr:::qtl2mr_pair_subsets(pilot, list(SF = list(steiger = TRUE), F2 = list(F_X_min = 2)), settings)
  expect_false(out$membership$SF[out$markers == "a1"])
  expect_true(out$membership$F2[out$markers == "a1"])
  keep <- out$membership$SF
  refit <- qtl2mr:::qtl2mr_subset_refit(pilot, keep, settings)$summary
  expect_identical(out$rows$theta[[1L]], refit$theta)
  expect_identical(out$rows$n_loci[[1L]], sum(keep))
})

test_that("pairs without a fit get pair_not_fitted rows; a failing subset leaves the main row intact", {
  f <- v064_fixture()
  pairs <- rbind(f$pairs, data.frame(exposure = "y", outcome = "x"))
  b <- v064_batch(f, pairs = pairs, subsets = v064_subsets[c("F2", "SF")])
  expect_identical(b$main_results$status, c("ok", "ok", "ok", "error"))
  sub <- b$subset_results
  expect_identical(nrow(sub), 8L)
  expect_identical(sub$status[sub$pair_id == "y__to__x"], rep("pair_not_fitted", 2L))
  expect_identical(sub$rules[sub$pair_id == "y__to__x"], c("F_X_min>=2", "steiger"))
  expect_true(all(is.na(sub$n_loci[sub$pair_id == "y__to__x"])))
  # A pair that fails while fitting.
  real_fit <- qtl2mr:::qtl2mr_fit_prepared
  local_mocked_bindings(qtl2mr_fit_prepared = function(input, exposure, outcome, ...) {
    if (identical(outcome, "w")) stop("injected pair failure")
    real_fit(input, exposure, outcome, ...)
  })
  failed <- v064_batch(f, subsets = v064_subsets[c("F2", "SF")])
  expect_identical(failed$main_results$status, c("ok", "error", "ok"))
  expect_identical(failed$subset_results$status, c("ok", "ok", "pair_not_fitted", "pair_not_fitted", "ok", "ok"))
  expect_identical(failed$errors$stage, "pair")
  expect_true(all(is.na(failed$locus_results$in_subset_F2[failed$locus_results$pair_id == "x__to__w"])))
})

test_that("a subset refit error is isolated, copied subsets share it, and other rows are unaffected", {
  f <- v064_fixture()
  subsets <- c(v064_subsets, list(F2_again = list(F_X_min = 2)))
  reference <- v064_batch(f, subsets = subsets)
  real_refit <- qtl2mr:::qtl2mr_subset_refit
  local_mocked_bindings(qtl2mr_subset_refit = function(pilot, keep, ...) {
    if (sum(keep) == 4L) stop("injected subset failure")
    real_refit(pilot, keep, ...)
  })
  b <- v064_batch(f, subsets = subsets)
  expect_identical(v064_drop_runtime(b$main_results), v064_drop_runtime(reference$main_results))
  expect_identical(b$locus_results, reference$locus_results)
  failing <- b$subset_results$n_loci == 4L & b$subset_results$subset %in% c("F2", "F2_again")
  expect_identical(sum(failing), 4L)
  expect_identical(b$subset_results$status[failing], rep("error", 4L))
  expect_identical(b$subset_results$error_message[failing], rep("injected subset failure", 4L))
  expect_true(all(is.na(b$subset_results$theta[failing])))
  expect_identical(b$subset_results$same_loci_as[failing & b$subset_results$subset == "F2_again"], c("F2", "F2"))
  expect_identical(b$subset_results[!failing, ], reference$subset_results[!failing, ])
  expect_identical(b$errors$stage, rep("subset", 4L))
  expect_identical(b$errors$pair_id, c("x__to__y", "x__to__y", "x__to__w", "x__to__w"))
  expect_identical(b$errors$error_message, rep(c("subset F2: injected subset failure",
    "subset F2_again: injected subset failure"), 2L))
  expect_identical(b$run_manifest$n_pairs_error, 0L)
})

test_that("a copied subset row equals an explicit refit", {
  f <- v064_fixture()
  full <- v064_batch(f, subsets = v064_subsets)
  # SF alone is refitted (it copies F2_SF in the full specification).
  alone <- v064_batch(f, subsets = v064_subsets["SF"])
  cols <- v064_fit_columns(full$subset_results)
  for (id in full$main_results$pair_id) {
    a <- full$subset_results[full$subset_results$pair_id == id & full$subset_results$subset == "SF", ]
    b <- alone$subset_results[alone$subset_results$pair_id == id, ]
    expect_true(is.na(b$same_loci_as) || identical(b$same_loci_as, "all"))
    expect_true(v064_same_values(a, b, cols))
  }
  # The maximal-set copy equals refitting every fitted block.
  settings <- qtl2mr:::qtl2mr_fit_settings("fail", "diagnostic_only", "chisq", 1000L, 1L, 0L, 1L, TRUE, "star")
  pilot <- suppressWarnings(qtl2mr:::qtl2mr_fit_prepared(qtl2mr_pair_input(f$resources, f$instruments, "x", "y"),
    "x", "y", settings))
  refit <- qtl2mr:::qtl2mr_subset_refit(pilot, rep(TRUE, length(pilot$locus_blocks)), settings)$summary
  keep <- setdiff(names(pilot$summary), qtl2mr:::qtl2mr_subset_dropped_summary)
  expect_true(v064_same_values(refit, pilot$summary, keep))
  zero <- full$subset_results[full$subset_results$pair_id == "x__to__y" & full$subset_results$subset == "zero", ]
  expect_true(v064_same_values(zero, refit, intersect(names(zero), keep)))
})

test_that("resume restores subset rows and refuses a changed specification; workers do not change rows", {
  f <- v064_fixture()
  dir <- tempfile("v064-resume-")
  first <- v064_batch(f, output_dir = dir, subsets = v064_subsets)
  expect_true(file.exists(file.path(dir, "batch_subset_results.csv")))
  written <- utils::read.csv(file.path(dir, "batch_subset_results.csv"), stringsAsFactors = FALSE)
  expect_identical(names(written), names(first$subset_results))
  expect_identical(written$status, first$subset_results$status)
  expect_true(all(paste0("in_subset_", names(v064_subsets)) %in%
    names(utils::read.csv(file.path(dir, "batch_locus_results.csv")))))
  expect_match(paste(readLines(file.path(dir, "batch_config.json")), collapse = ""),
    "\"subsets\":{\"F5\":{\"F_X_min\":5},", fixed = TRUE)
  again <- v064_batch(f, output_dir = dir, resume = TRUE, subsets = v064_subsets)
  expect_identical(again$run_manifest$n_pairs_resumed, 3L)
  expect_identical(again$subset_results, first$subset_results)
  expect_identical(again$locus_results, first$locus_results)
  # An interrupted run: one pair record removed, then resumed.
  unlink(file.path(dir, "batch_pair_records", "pair_000002.rds"))
  partial <- v064_batch(f, output_dir = dir, resume = TRUE, subsets = v064_subsets)
  expect_identical(partial$run_manifest$n_pairs_resumed, 2L)
  expect_identical(partial$subset_results, first$subset_results)
  expect_identical(partial$errors, first$errors)
  for (changed in list(v064_subsets[-1L], rev(v064_subsets), NULL,
                       utils::modifyList(v064_subsets, list(F2 = list(F_X_min = 2.5))))) {
    expect_error(v064_batch(f, output_dir = dir, resume = TRUE, subsets = changed), "fingerprint changed")
  }
  dir_off <- tempfile("v064-resume-off-")
  v064_batch(f, output_dir = dir_off)
  expect_error(v064_batch(f, output_dir = dir_off, resume = TRUE, subsets = v064_subsets), "fingerprint changed")
  # overwrite = TRUE removes a subset table left by an earlier run.
  v064_batch(f, output_dir = dir, overwrite = TRUE)
  expect_false(file.exists(file.path(dir, "batch_subset_results.csv")))

  # With the simulated Q p-value and the bootstrap (no longer batch options since
  # 0.6.1.9006), so that the refits draw random numbers.
  serial <- qtl2mr_with_settings(v064_batch(f, subsets = v064_subsets), heterogeneity_method = "simulation",
    heterogeneity_nsim = 20L, diagnostic_bootstrap_B = 10L)
  forked <- qtl2mr_with_settings(v064_batch(f, subsets = v064_subsets, n_workers = 2L),
    heterogeneity_method = "simulation", heterogeneity_nsim = 20L, diagnostic_bootstrap_B = 10L)
  # The mock took effect: simulated Q p-values and a 10-draw bootstrap.
  expect_identical(unique(serial$main_results$Q_p_method), "simulation")
  expect_true(all(is.finite(serial$main_results$Q_p_sim)))
  boot <- serial$method_comparison$bootstrap_B[serial$method_comparison$method == "Weighted median"]
  expect_true(length(boot) > 0L && all(boot == 10L))
  expect_identical(forked$subset_results, serial$subset_results)
  expect_identical(forked$locus_results, serial$locus_results)
  expect_identical(forked$errors, serial$errors)
  expect_identical(v064_drop_runtime(forked$main_results), v064_drop_runtime(serial$main_results))
  if (identical(.Platform$OS.type, "unix")) expect_identical(forked$run_manifest$parallel_backend, "parallel::mclapply")
})

test_that("subset refits leave the caller's random number stream as without subsets", {
  f <- v064_fixture()
  # The refits draw random numbers only with the simulated Q p-value and the
  # bootstrap, which the entries no longer take (0.6.1.9006): the batch with the
  # 0.6.1.9005 simulation settings below, the single analysis with its 0.6.1.9005
  # defaults ("auto", which simulates for several loci) and the counts below.
  draw_after <- function(...) {
    set.seed(99L)
    b <- qtl2mr_with_settings(v064_batch(f, ...), heterogeneity_method = "simulation",
      heterogeneity_nsim = 10L, diagnostic_bootstrap_B = 5L)
    expect_identical(unique(b$main_results$Q_p_method), "simulation")  # the mock took effect
    stats::runif(3L)
  }
  expect_identical(draw_after(subsets = v064_subsets), draw_after())
  single_after <- function(...) {
    set.seed(7L)
    r <- qtl2mr_with_settings(suppressWarnings(run_qtl2mr_analysis(
      qtl2mr_pair_input(f$resources, f$instruments, "x", "y"),
      output_dir = tempfile("v064-rng-"), make_html = FALSE, make_figures = FALSE, progress = FALSE, ...)),
      heterogeneity_method = "auto", heterogeneity_nsim = 10L, diagnostic_bootstrap_B = 5L)
    expect_identical(r$pilot$summary$Q_p_method, "simulation")  # the mock took effect
    stats::runif(3L)
  }
  expect_identical(single_after(subsets = v064_subsets[c("F2", "SF")]), single_after())
})

test_that("a failure anywhere in the subset step gives error rows and leaves the main outputs unchanged", {
  f <- v064_fixture()
  subsets <- v064_subsets[c("F2", "SF")]
  reference <- v064_batch(f, subsets = subsets)
  off <- v064_batch(f)
  main_cols <- names(off$locus_results)
  real_pair <- qtl2mr:::qtl2mr_pair_subsets
  # Row assembly fails for one pair: its rows lack a column the batch table needs.
  b <- with_mocked_bindings(v064_batch(f, subsets = subsets), qtl2mr_pair_subsets = function(pilot, ...) {
    out <- real_pair(pilot, ...)
    if (identical(pilot$summary$outcome[[1L]], "w")) out$rows$theta <- NULL
    out
  })
  expect_identical(v064_drop_runtime(b$main_results), v064_drop_runtime(reference$main_results))
  expect_identical(v064_drop_runtime(b$main_results), v064_drop_runtime(off$main_results))
  expect_identical(b$locus_results[main_cols], off$locus_results)
  bad <- b$subset_results$pair_id == "x__to__w"
  expect_identical(b$subset_results$status[bad], c("error", "error"))
  expect_match(b$subset_results$error_message[bad], "undefined columns")
  expect_true(all(is.na(b$subset_results$n_loci[bad]) & is.na(b$subset_results$theta[bad])))
  expect_identical(b$subset_results$rules[bad], c("F_X_min>=2", "steiger"))
  expect_identical(b$subset_results[!bad, ], reference$subset_results[!bad, ])
  expect_true(all(is.na(unlist(b$locus_results[b$locus_results$pair_id == "x__to__w", c("in_subset_F2", "in_subset_SF")]))))
  expect_identical(b$locus_results[b$locus_results$pair_id != "x__to__w", ],
    reference$locus_results[reference$locus_results$pair_id != "x__to__w", ])
  expect_identical(b$errors$stage, c("subset", "subset"))
  expect_identical(b$errors$pair_id, c("x__to__w", "x__to__w"))
  expect_match(b$errors$error_message, "^subset (F2|SF): .*undefined columns")
  expect_identical(b$run_manifest$n_pairs_error, 0L)
  # The whole step fails for every pair; pair records are still written.
  dir <- tempfile("v064-stepfail-")
  all_bad <- with_mocked_bindings(v064_batch(f, subsets = subsets, output_dir = dir),
    qtl2mr_pair_subsets = function(...) stop("injected subset step failure"))
  expect_identical(v064_drop_runtime(all_bad$main_results), v064_drop_runtime(off$main_results))
  expect_identical(all_bad$locus_results[main_cols], off$locus_results)
  expect_identical(all_bad$subset_results$status, rep("error", 6L))
  expect_identical(all_bad$subset_results$error_message, rep("injected subset step failure", 6L))
  expect_identical(names(all_bad$subset_results), qtl2mr:::qtl2mr_batch_subset_columns())
  expect_identical(nrow(all_bad$errors), 6L)
  expect_length(list.files(file.path(dir, "batch_pair_records")), 3L)
  expect_identical(utils::read.csv(file.path(dir, "batch_subset_results.csv"), stringsAsFactors = FALSE)$status,
    rep("error", 6L))
  # Single analysis: error rows with the message, main outputs as without subsets.
  run_single <- function(...) suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, f$instruments, "x", "y"),
    output_dir = tempfile("v064-stepfail-single-"), make_html = FALSE, make_figures = FALSE, progress = FALSE, ...))
  single_off <- run_single()
  single_ok <- run_single(subsets = subsets)
  single_bad <- with_mocked_bindings(run_single(subsets = subsets),
    qtl2mr_pair_subsets = function(...) stop("injected subset step failure"))
  sub <- single_bad$subset_results
  expect_identical(names(sub), names(single_ok$subset_results))
  expect_identical(sub$subset, names(subsets))
  expect_identical(sub$status, c("error", "error"))
  expect_identical(sub$error_message, rep("injected subset step failure", 2L))
  expect_true(all(is.na(sub$theta) & is.na(sub$n_loci)))
  expect_true(file.exists(file.path(single_bad$output_dir, "data", "subset_results.csv")))
  expect_identical(single_bad$pilot$summary, single_off$pilot$summary)
  for (file in c("main_result.csv", "instrument_sensitivity.csv", "steiger_sensitivity.csv", "loo_mr.csv",
                 "method_comparison.csv", "block_first_stage.csv", "warnings.csv", "q_calibration.csv")) {
    expect_identical(unname(tools::md5sum(file.path(single_bad$output_dir, "data", file))),
      unname(tools::md5sum(file.path(single_off$output_dir, "data", file))), label = file)
  }
})

test_that("O1 locus blocks hold the maximal set whether or not subsets are requested", {
  f <- v064_fixture()
  a <- tempfile("v064-blocks-a-")
  b <- tempfile("v064-blocks-b-")
  v064_batch(f, output_dir = a, write_locus_blocks = TRUE)
  v064_batch(f, output_dir = b, write_locus_blocks = TRUE, subsets = v064_subsets)
  files <- list.files(file.path(a, "batch_locus_blocks"))
  expect_identical(files, list.files(file.path(b, "batch_locus_blocks")))
  for (file in files) {
    expect_identical(readRDS(file.path(a, "batch_locus_blocks", file)), readRDS(file.path(b, "batch_locus_blocks", file)))
  }
})

# Key equivalence: each subset row equals the main row of a separate run whose
# instruments are only that subset's loci (for that pair's exposure), with the
# same settings. Compared columns: every subset_results column except subset,
# rules and same_loci_as (which a main row does not have); main-only columns
# (n_loci_selected, n_loci_dropped_nonestimable, locus_fit_status, runtime_sec)
# are not in subset_results.
v064_batch_equivalence <- function(f, subsets, ...) {
  b <- v064_batch(f, subsets = subsets, ...)
  cols <- v064_fit_columns(b$subset_results)
  checked <- 0L
  for (i in seq_len(nrow(b$subset_results))) {
    row <- b$subset_results[i, ]
    if (!identical(row$status, "ok")) next
    loci <- b$locus_results[b$locus_results$pair_id == row$pair_id, ]
    markers <- loci$marker[loci[[paste0("in_subset_", row$subset)]] %in% TRUE]
    pair <- f$pairs[match(row$pair_id, b$main_results$pair_id), ]
    inst <- f$instruments[f$instruments$exposure == pair$exposure & f$instruments$marker %in% markers, ]
    separate <- v064_batch(f, pairs = pair, instruments = inst, ...)
    expect_identical(separate$main_results$status, "ok")
    expect_identical(separate$main_results$n_loci, length(markers))
    expect_true(v064_same_values(row, separate$main_results, cols), label = paste(row$pair_id, row$subset))
    checked <- checked + 1L
  }
  checked
}

test_that("batch: a subset row equals a separate batch run on that subset's loci", {
  f <- v064_fixture()
  subsets <- v064_subsets[c("F5", "F2", "F2_SF", "zero")]
  expect_identical(v064_batch_equivalence(f, subsets), 12L)
  # With every other fit setting, which the batch no longer takes (0.6.1.9006);
  # sample_policy = "pair_common_complete" was removed.
  expect_identical(qtl2mr_with_settings(v064_batch_equivalence(f, subsets), heterogeneity_method = "simulation",
    heterogeneity_nsim = 25L, heterogeneity_seed = 7L, diagnostic_bootstrap_B = 15L, diagnostic_bootstrap_seed = 3L,
    reported_estimator_policy = "auto_fallback", cross_cov_model = "raw", standardize_phenotypes = FALSE), 12L)
  raw <- qtl2mr_with_settings(v064_batch(f, pairs = f$pairs[1L, ])$main_results, cross_cov_model = "raw",
    standardize_phenotypes = FALSE, reported_estimator_policy = "auto_fallback")
  expect_identical(c(raw$cross_cov_model, raw$reported_estimator_policy), c("raw", "auto_fallback"))
  expect_false(raw$standardize_phenotypes)
})

v064_single <- function(f, inst, exposure, outcome, ...) {
  suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, exposure, outcome),
    output_dir = tempfile("v064-single-"), make_html = FALSE, make_figures = FALSE, progress = FALSE, ...))
}

test_that("single: subset rows equal separate runs, main outputs are unchanged, and batch rows agree", {
  f <- v064_fixture()
  subsets <- v064_subsets[c("F5", "F2", "F2_SF", "zero", "huge")]
  # The Q simulation and bootstrap settings are no longer entry options
  # (0.6.1.9006); every run below uses them in place of the fixed settings.
  with_settings <- function(code) qtl2mr_with_settings(code, heterogeneity_method = "auto", heterogeneity_nsim = 30L,
    heterogeneity_seed = 5L, diagnostic_bootstrap_B = 20L, diagnostic_bootstrap_seed = 2L)
  on <- with_settings(v064_single(f, f$instruments, "x", "y", subsets = subsets))
  off <- with_settings(v064_single(f, f$instruments, "x", "y"))
  expect_identical(on$pilot$summary$Q_p_method, "simulation")
  expect_equal(on$pilot$summary$Q_sim_n_requested, 30)
  expect_identical(on$pilot$summary$diagnostic_bootstrap_B, 20L)
  sub <- on$subset_results
  expect_identical(sub$subset, names(subsets))
  expect_identical(sub$status, c("ok", "ok", "ok", "ok", "too_few_loci"))
  expect_identical(sub$same_loci_as, c(NA, NA, NA, "all", NA))
  expect_identical(names(sub), c(qtl2mr:::qtl2mr_subset_id_columns,
    setdiff(names(on$pilot$summary), qtl2mr:::qtl2mr_subset_dropped_summary)))
  # The file is the table as written by write.csv() (15 significant digits).
  roundtrip <- tempfile(fileext = ".csv")
  utils::write.csv(sub, roundtrip, row.names = FALSE)
  expect_identical(unname(tools::md5sum(file.path(on$output_dir, "data", "subset_results.csv"))),
    unname(tools::md5sum(roundtrip)))
  expect_true("subset_results.csv" %in% basename(on$output_manifest$path_relative))
  # Main outputs, including the existing sensitivity tables, are as without subsets.
  expect_identical(on$pilot$summary, off$pilot$summary)
  for (file in c("main_result.csv", "instrument_sensitivity.csv", "steiger_sensitivity.csv", "loo_mr.csv",
                 "method_comparison.csv", "block_first_stage.csv", "warnings.csv", "q_calibration.csv")) {
    expect_identical(unname(tools::md5sum(file.path(on$output_dir, "data", file))),
      unname(tools::md5sum(file.path(off$output_dir, "data", file))), label = file)
  }
  expect_identical(on$instrument_sensitivity, off$instrument_sensitivity)
  expect_identical(on$steiger_sensitivity, off$steiger_sensitivity)
  # Each ok row equals a separate single analysis on the subset's loci (Q simulation and bootstrap included).
  first <- on$pilot$locus_first_stage
  cols <- setdiff(names(sub), c("subset", "rules", "same_loci_as"))
  for (i in which(sub$status == "ok")) {
    keep <- qtl2mr:::qtl2mr_subset_members(first, subsets[[sub$subset[[i]]]])
    inst <- f$instruments[f$instruments$exposure == "x" & f$instruments$marker %in% first$marker[keep], ]
    separate <- with_settings(v064_single(f, inst, "x", "y"))
    expect_identical(separate$pilot$summary$n_loci, sub$n_loci[[i]])
    expect_true(v064_same_values(sub[i, ], cbind(separate$pilot$summary, status = "ok", error_message = ""), cols),
      label = sub$subset[[i]])
  }
  # The batch gives the same subset rows for the same pair and settings.
  batch <- with_settings(v064_batch(f, pairs = f$pairs[1L, ], subsets = subsets))
  bcols <- setdiff(names(batch$subset_results), "pair_id")
  for (i in seq_len(nrow(sub))) {
    expect_true(v064_same_values(batch$subset_results[i, ], sub[i, ], bcols), label = sub$subset[[i]])
  }
})
