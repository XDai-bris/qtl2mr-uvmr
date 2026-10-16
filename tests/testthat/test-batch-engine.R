make_nonestimable_batch_fixture <- function(n = 70L, seed = 901L) {
  gp_good <- make_sim_genoprobs(n = n, marker = "m_good", seed = seed)
  ids <- dimnames(gp_good[[1L]])[[1L]]
  founders <- dimnames(gp_good[[1L]])[[2L]]
  set.seed(seed + 1L)
  raw_bad <- matrix(stats::rexp(length(ids) * length(founders)),
                    nrow = length(ids), ncol = length(founders),
                    dimnames = list(ids, founders))
  prob_bad <- raw_bad / rowSums(raw_bad)
  prob_bad[, "A"] <- prob_bad[, "A"] + prob_bad[, "H"]
  prob_bad[, "H"] <- 0
  allele_probs <- list("1" = array(
    c(gp_good[[1L]][, , "m_good"], prob_bad),
    dim = c(length(ids), length(founders), 2L),
    dimnames = list(ids, founders, c("m_good", "m_bad"))
  ))
  class(allele_probs) <- c("calc_genoprob", "list")
  contrasts <- make_founder_contrasts(founders)
  z <- allele_probs[["1"]][, , "m_good"] %*% contrasts
  pheno <- data.frame(
    x = as.numeric(0.7 * z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.6 * z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  selected <- data.frame(
    marker = c("m_good", "m_bad"),
    chr = c("1", "1"),
    pos = c(1, 2),
    locus_rank = 1:2,
    stringsAsFactors = FALSE
  )
  selected_file <- tempfile("qtl2mr-nonestimable-loci-", fileext = ".csv")
  utils::write.csv(selected, selected_file, row.names = FALSE)
  marker_map <- selected[, c("marker", "chr", "pos")]
  qtl_scan <- list(x = matrix(
    c(5.5, 5.0),
    ncol = 1L,
    dimnames = list(selected$marker, "lod")
  ))
  resources <- qtl2mr_prepare_resources(
    allele_probs,
    kinship = NULL,
    pheno = pheno,
    covar = covar,
    qtl_scan = qtl_scan,
    marker_map = marker_map
  )
  input <- as_qtl2mr_prepared_input(
    allele_probs,
    pheno,
    covar,
    kinship = NULL,
    marker_map = marker_map,
    selected_loci = selected,
    exposure = "x",
    outcome = "y"
  )
  pairs <- data.frame(
    exposure = "x",
    outcome = "y",
    selected_loci_file = selected_file,
    stringsAsFactors = FALSE
  )
  list(input = input, resources = resources, pairs = pairs, selected = selected)
}

test_that("batch founder LMM matches single-marker LMM", {
  gp1 <- make_sim_genoprobs(n = 70, marker = "m1", seed = 801)
  set.seed(802)
  raw2 <- matrix(stats::rexp(70 * 8), nrow = 70, ncol = 8)
  probs2 <- raw2 / rowSums(raw2)
  arr <- array(NA_real_, dim = c(70, 8, 2),
               dimnames = list(dimnames(gp1[[1]])[[1]], LETTERS[1:8], c("m1", "m2")))
  arr[, , "m1"] <- gp1[[1]][, , "m1"]
  arr[, , "m2"] <- probs2
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  ids <- dimnames(arr)[[1L]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- arr[, , "m1"] %*% A
  y <- stats::setNames(as.numeric(0.8 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.4)), ids)
  covar <- matrix(stats::rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))

  ctx <- qtl2mr:::prepare_lmm_context(y, addcovar = covar, sample_ids = ids,
                                      standardize_pheno = TRUE, trait_name = "x", chr = "1")
  batch <- qtl2mr:::fit_founder_lmm_blocks_batch(gp[["1"]], c("m1", "m2"), ctx, A,
                                                 fallback_genoprobs = gp)
  single <- fit_founder_lmm_block(gp, y, addcovar = covar, contrasts = A,
                                  marker = "m1", sample_ids = ids,
                                  standardize_pheno = TRUE)
  expect_equal(batch$m1$beta_fc, single$beta_fc, tolerance = 1e-10)
  expect_equal(batch$m1$vcov_fc, single$vcov_fc, tolerance = 1e-10)
  expect_equal(batch$m1$sigma, single$sigma, tolerance = 1e-10)
  expect_equal(batch$m1$influence_star, single$influence_star, tolerance = 1e-10)
})

test_that("batch founder LMM matches single-marker LMM with kinship", {
  gp <- make_sim_genoprobs(n = 64, marker = "m1", seed = 806)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  set.seed(807)
  y <- stats::setNames(as.numeric(0.7 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)), ids)
  covar <- matrix(stats::rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  raw_k <- matrix(stats::rnorm(length(ids) * 5), nrow = length(ids))
  kinship <- tcrossprod(raw_k) / 5
  kinship <- kinship / mean(diag(kinship))
  dimnames(kinship) <- list(ids, ids)

  ctx <- qtl2mr:::prepare_lmm_context(
    y, kinship = kinship, addcovar = covar, hsq = 0.35,
    sample_ids = ids, standardize_pheno = TRUE, trait_name = "x", chr = "1"
  )
  batch <- qtl2mr:::fit_founder_lmm_blocks_batch(gp[["1"]], "m1", ctx, A,
                                                 fallback_genoprobs = gp)
  single <- fit_founder_lmm_block(
    gp, y, kinship = kinship, addcovar = covar, contrasts = A,
    hsq = 0.35, marker = "m1", sample_ids = ids,
    standardize_pheno = TRUE
  )
  expect_equal(batch$m1$beta_fc, single$beta_fc, tolerance = 1e-10)
  expect_equal(batch$m1$vcov_fc, single$vcov_fc, tolerance = 1e-10)
  expect_equal(batch$m1$influence_star, single$influence_star, tolerance = 1e-10)
})

test_that("batch MR locus blocks match per-locus overlap fitter", {
  gp <- make_sim_genoprobs(n = 72, marker = "m1", seed = 811)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  x <- stats::setNames(as.numeric(0.7 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.5)), ids)
  y <- stats::setNames(as.numeric(0.5 * x + 0.2 * Z[, "FC3"] + stats::rnorm(length(ids), sd = 0.5)), ids)
  covar <- matrix(stats::rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  selected <- data.frame(marker = "m1", chr = "1", stringsAsFactors = FALSE)

  batch <- qtl2mr:::fit_mr_locus_blocks_batch(
    gp, selected, x, y, addcovar = covar, contrasts = A,
    sample_ids_x = ids[1:68], sample_ids_y = ids[5:72],
    standardize_phenotypes = TRUE
  )$m1
  single <- fit_mr_locus_block_overlap(
    gp, x, y, addcovar_x = covar, addcovar_y = covar, contrasts = A,
    marker = "m1", sample_ids_x = ids[1:68], sample_ids_y = ids[5:72],
    standardize_phenotypes = TRUE
  )
  expect_equal(batch$beta_x, single$beta_x, tolerance = 1e-10)
  expect_equal(batch$beta_y, single$beta_y, tolerance = 1e-10)
  expect_equal(batch$Sigma_x, single$Sigma_x, tolerance = 1e-10)
  expect_equal(batch$Sigma_y, single$Sigma_y, tolerance = 1e-10)
  expect_equal(batch$Sigma_xy, single$Sigma_xy, tolerance = 1e-10)
  expect_equal(batch$O_xy_summary$n_overlap, single$O_xy_summary$n_overlap)
})

test_that("batch MR drops selected loci with non-estimable founder contrasts", {
  gp1 <- make_sim_genoprobs(n = 70, marker = "m_good", seed = 813)
  ids <- dimnames(gp1[[1]])[[1]]
  founders <- dimnames(gp1[[1]])[[2]]
  set.seed(814)
  raw_bad <- matrix(stats::rexp(length(ids) * 8), nrow = length(ids), ncol = 8,
                    dimnames = list(ids, founders))
  prob_bad <- raw_bad / rowSums(raw_bad)
  prob_bad[, "A"] <- prob_bad[, "A"] + prob_bad[, "H"]
  prob_bad[, "H"] <- 0
  arr <- array(NA_real_, dim = c(length(ids), 8, 2),
               dimnames = list(ids, founders, c("m_good", "m_bad")))
  arr[, , "m_good"] <- gp1[[1]][, , "m_good"]
  arr[, , "m_bad"] <- prob_bad
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  A <- make_founder_contrasts(founders)
  Z <- arr[, , "m_good"] %*% A
  x <- stats::setNames(as.numeric(0.7 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)), ids)
  y <- stats::setNames(as.numeric(0.6 * x + stats::rnorm(length(ids), sd = 0.5)), ids)
  selected <- data.frame(marker = c("m_good", "m_bad"), chr = c("1", "1"), stringsAsFactors = FALSE)

  blocks <- qtl2mr:::fit_mr_locus_blocks_batch(
    gp, selected, x, y, contrasts = A, sample_ids_x = ids, sample_ids_y = ids,
    standardize_phenotypes = TRUE, drop_nonestimable_loci = TRUE
  )
  dropped <- attr(blocks, "dropped_loci")
  expect_equal(names(blocks), "m_good")
  expect_equal(attr(blocks, "n_loci_selected"), 2L)
  expect_equal(attr(blocks, "n_loci_retained"), 1L)
  expect_equal(attr(blocks, "n_loci_dropped"), 1L)
  expect_true(is.data.frame(dropped))
  expect_true("m_bad" %in% dropped$marker)
  expect_true(all(dropped$drop_reason == "nonestimable_founder_contrast"))

  expect_error(
    qtl2mr:::fit_mr_locus_blocks_batch(
      gp, selected, x, y, contrasts = A, sample_ids_x = ids, sample_ids_y = ids,
      standardize_phenotypes = TRUE, drop_nonestimable_loci = FALSE
    ),
    "Founder contrast coefficients are not estimable"
  )
})

test_that("high-level and exact workflows fail closed on non-estimable loci", {
  fixture <- make_nonestimable_batch_fixture()
  expect_error(
    run_qtl2mr_analysis(
      fixture$input,
      adapter = "prepared",
      output_dir = tempfile("qtl2mr-confirmatory-fail-"),
      overwrite = TRUE,
      make_html = FALSE,
      make_figures = FALSE,
      verbose = FALSE,
      heterogeneity_method = "none",
      diagnostic_bootstrap_B = 0L,
      run_loo_mr = FALSE,
      run_instrument_sensitivity = FALSE
    ),
    "Founder contrast coefficients are not estimable"
  )

  exact <- run_qtl2mr_batch(
    fixture$resources,
    fixture$pairs,
    fit_options = list(diagnostic_bootstrap_B = 0L),
    progress = FALSE
  )
  expect_equal(exact$main_results$status, "error")
  expect_match(exact$main_results$error_message, "Founder contrast coefficients are not estimable")
  expect_equal(exact$run_manifest$nonestimable_locus_policy, "fail_closed")

  grouped <- run_qtl2mr_batch_exact_grouped(
    fixture$resources,
    fixture$pairs,
    fit_options = list(diagnostic_bootstrap_B = 0L),
    progress = FALSE
  )
  expect_equal(grouped$main_results$status, "error")
  expect_match(grouped$main_results$error_message, "Founder contrast coefficients are not estimable")
  expect_equal(grouped$run_manifest$nonestimable_locus_policy, "fail_closed")
})

test_that("fast screening drops non-estimable loci with a complete audit", {
  fixture <- make_nonestimable_batch_fixture(seed = 911L)
  fast <- run_qtl2mr_batch_fast(
    fixture$resources,
    fixture$pairs,
    fit_options = list(diagnostic_bootstrap_B = 0L),
    progress = FALSE
  )

  expect_equal(fast$main_results$status, "ok")
  expect_equal(fast$main_results$locus_fit_status, "ok_with_dropped_nonestimable_loci")
  expect_equal(fast$main_results$n_loci_selected, 2L)
  expect_equal(fast$main_results$n_loci, 1L)
  expect_equal(fast$main_results$n_loci_dropped_nonestimable, 1L)
  expect_equal(fast$main_results$reported_estimator_policy, "diagnostic_only")
  expect_equal(nrow(fast$selected_loci), 2L)
  expect_equal(sort(fast$selected_loci$locus_status), c("dropped_nonestimable", "retained"))
  dropped <- fast$selected_loci[fast$selected_loci$locus_status == "dropped_nonestimable", , drop = FALSE]
  expect_equal(dropped$marker, "m_bad")
  expect_match(dropped$drop_reason, "nonestimable_founder_contrast")
  expect_true(nzchar(dropped$error_message))
  expect_true("dropped_nonestimable_loci" %in% fast$warnings$warning_id)
  expect_equal(fast$run_manifest$nonestimable_locus_policy, "drop_audited_screening")
  expect_equal(fast$run_manifest$n_loci_dropped_nonestimable, 1L)
})

test_that("batch LMM context cache reuses trait-chromosome contexts", {
  gp <- make_sim_genoprobs(n = 66, marker = "m1", seed = 812)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  x <- stats::setNames(as.numeric(0.6 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.5)), ids)
  y <- stats::setNames(as.numeric(0.4 * x + stats::rnorm(length(ids), sd = 0.5)), ids)
  covar <- matrix(stats::rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  selected <- data.frame(marker = "m1", chr = "1", stringsAsFactors = FALSE)
  cache <- qtl2mr:::qtl2mr_context_cache_new()
  fit_cache <- qtl2mr:::qtl2mr_context_cache_new()

  first <- qtl2mr:::fit_mr_locus_blocks_batch(
    gp, selected, x, y, addcovar = covar, contrasts = A,
    sample_ids_x = ids, sample_ids_y = ids, standardize_phenotypes = TRUE,
    context_cache = cache, fit_cache = fit_cache,
    exposure_name = "x", outcome_name = "y"
  )$m1
  stats_first <- qtl2mr:::qtl2mr_context_cache_stats(cache)
  fit_stats_first <- qtl2mr:::qtl2mr_context_cache_stats(fit_cache)
  second <- qtl2mr:::fit_mr_locus_blocks_batch(
    gp, selected, x, y, addcovar = covar, contrasts = A,
    sample_ids_x = ids, sample_ids_y = ids, standardize_phenotypes = TRUE,
    context_cache = cache, fit_cache = fit_cache,
    exposure_name = "x", outcome_name = "y"
  )$m1
  stats_second <- qtl2mr:::qtl2mr_context_cache_stats(cache)
  fit_stats_second <- qtl2mr:::qtl2mr_context_cache_stats(fit_cache)

  expect_equal(stats_first$misses, 2L)
  expect_equal(stats_first$hits, 0L)
  expect_equal(fit_stats_first$misses, 2L)
  expect_equal(fit_stats_first$hits, 0L)
  expect_equal(stats_second$misses, 2L)
  expect_equal(stats_second$hits, 2L)
  expect_equal(fit_stats_second$misses, 2L)
  expect_equal(fit_stats_second$hits, 2L)
  expect_equal(first$Sigma_x, second$Sigma_x, tolerance = 1e-12)
  expect_equal(first$Sigma_xy, second$Sigma_xy, tolerance = 1e-12)
})

test_that("scan locus selector reproduces contiguous LOD island maxima", {
  markers <- paste0("m", 1:10)
  scan <- matrix(c(0, 4.5, 5, 4.2, 0, 4.1, 4.3, 0, 5.5, 4.8), ncol = 1,
                 dimnames = list(markers, "lod"))
  marker_map <- data.frame(marker = markers, chr = c(rep("1", 8), "X", "2"),
                           pos = seq_along(markers), stringsAsFactors = FALSE)
  selected <- qtl2mr_select_loci_from_scan(list(x = scan), marker_map, "x", lod_thresh = 4, exclude_chr = "X")
  # The chr2 marker must not be swallowed by the adjacent chrX island.
  expect_equal(selected$marker, c("m3", "m10", "m7"))
  expect_equal(selected$locus_rank, 1:3)
})

test_that("multi-trait selected loci subset allele probabilities to union markers", {
  gp1 <- make_sim_genoprobs(n = 30, marker = "m1", seed = 816)
  ids <- dimnames(gp1[[1]])[[1]]
  founders <- dimnames(gp1[[1]])[[2]]
  set.seed(817)
  arr1 <- array(stats::runif(30 * 8 * 3), dim = c(30, 8, 3),
                dimnames = list(ids, founders, c("m1", "m2", "m3")))
  arr1 <- sweep(arr1, c(1, 3), apply(arr1, c(1, 3), sum), "/")
  arr2 <- array(stats::runif(30 * 8 * 2), dim = c(30, 8, 2),
                dimnames = list(ids, founders, c("m4", "m5")))
  arr2 <- sweep(arr2, c(1, 3), apply(arr2, c(1, 3), sum), "/")
  gp <- list("1" = arr1, "2" = arr2)
  class(gp) <- c("calc_genoprob", "list")
  marker_map <- data.frame(
    marker = paste0("m", 1:5),
    chr = c("1", "1", "1", "2", "2"),
    pos = 1:5,
    stringsAsFactors = FALSE
  )
  qtl_scan <- list(
    x = matrix(c(4.5, 0, 0, 5.2, 0), ncol = 1, dimnames = list(marker_map$marker, "lod")),
    y = matrix(c(0, 4.8, 0, 5.1, 0), ncol = 1, dimnames = list(marker_map$marker, "lod"))
  )
  selected <- qtl2mr_select_loci_for_traits(qtl_scan, marker_map, c("x", "y"), lod_thresh = 4)
  expect_true(all(c("instrument_trait", "marker", "chr", "exposure_lod") %in% names(selected)))
  expect_equal(sort(unique(selected$instrument_trait)), c("x", "y"))
  subset_gp <- qtl2mr_subset_allele_probs_to_loci(gp, selected)
  expect_s3_class(subset_gp, "calc_genoprob")
  expect_equal(dimnames(subset_gp[["1"]])[[3]], c("m1", "m2"))
  expect_equal(dimnames(subset_gp[["2"]])[[3]], "m4")
  expect_equal(dim(subset_gp[["1"]])[1:2], c(30, 8))
  expect_equal(attr(subset_gp, "qtl2mr_subset")$n_retained_loci, 3L)
})

test_that("LOD selection validates the requested instrument trait before peak picking", {
  marker_map <- data.frame(
    marker = paste0("m", 1:4), chr = "1", pos = 1:4,
    stringsAsFactors = FALSE
  )
  scan <- matrix(
    c(5, 0, 0, 0,
      0, 0, 0, 6),
    ncol = 2,
    dimnames = list(marker_map$marker, c("traitA", "traitB"))
  )
  expect_equal(unname(qtl2mr:::qtl2mr_lod_vector_from_scan(scan, "traitA")["m1"]), 5)
  expect_equal(unname(qtl2mr:::qtl2mr_lod_vector_from_scan(structure(scan, class = c("scan1", "matrix")), "traitB")["m4"]), 6)
  expect_equal(unname(qtl2mr:::qtl2mr_lod_vector_from_scan(as.data.frame(scan), "traitB")["m4"]), 6)
  expect_error(qtl2mr:::qtl2mr_lod_vector_from_scan(scan, "missing"), "missing instrument trait")
  duplicate <- scan
  colnames(duplicate) <- c("traitA", "traitA")
  expect_error(qtl2mr:::qtl2mr_lod_vector_from_scan(duplicate, "traitA"), "duplicated")
  unnamed_multi <- scan
  colnames(unnamed_multi) <- NULL
  expect_error(qtl2mr:::qtl2mr_lod_vector_from_scan(unnamed_multi, "traitA"), "non-empty trait column names")

  one <- scan[, "traitA", drop = FALSE]
  expect_silent(qtl2mr:::qtl2mr_lod_vector_from_scan(one, "traitA"))
  expect_error(qtl2mr:::qtl2mr_lod_vector_from_scan(one, "traitB"), "does not match")
  unnamed_one <- unname(one)
  rownames(unnamed_one) <- marker_map$marker
  expect_warning(qtl2mr:::qtl2mr_lod_vector_from_scan(unnamed_one, "traitA"), "could not be validated")
  expect_equal(unname(qtl2mr:::qtl2mr_lod_vector_from_scan(list(traitA = one), "traitA")["m1"]), 5)

  selected <- qtl2mr_select_loci_for_traits(scan, marker_map, c("traitA", "traitB"), lod_thresh = 4, exclude_chr = NULL)
  expect_equal(selected$marker, c("m1", "m4"))
  expect_equal(selected$instrument_trait, c("traitA", "traitB"))
  expect_equal(qtl2mr_select_loci_from_scan(scan, marker_map, "traitB", lod_thresh = 4, exclude_chr = NULL)$instrument_trait, "traitB")

  ids <- c("id1", "id2")
  gp <- list("1" = array(
    rep(1 / 8, length(ids) * 8 * nrow(marker_map)),
    dim = c(length(ids), 8L, nrow(marker_map)),
    dimnames = list(ids, LETTERS[1:8], marker_map$marker)
  ))
  resources <- qtl2mr_prepare_resources(
    allele_probs = gp,
    kinship = NULL,
    pheno = data.frame(traitA = 1:2, traitB = 3:4, outcome = 5:6, row.names = ids),
    covar = data.frame(covar = c(0, 1), row.names = ids),
    qtl_scan = scan,
    marker_map = marker_map
  )
  pairs <- data.frame(
    pair_id = c("traitA__to__outcome", "traitB__to__outcome"),
    exposure = c("traitA", "traitB"), outcome = "outcome",
    instrument_trait = c("traitA", "traitB"), stringsAsFactors = FALSE
  )
  selected_by_batch_route <- qtl2mr:::qtl2mr_selected_loci_for_pairs(
    resources, pairs, list(lod_thresh = 4, exclude_chr = NULL, max_loci = NULL)
  )
  expect_equal(selected_by_batch_route$marker, c("m1", "m4"))
  expect_equal(selected_by_batch_route$instrument_trait, c("traitA", "traitB"))
})

test_that("LOD threshold and marker identifiers fail before peak selection", {
  scan <- matrix(c(5, 0), ncol = 1, dimnames = list(c("m1", "m2"), "traitA"))
  marker_map <- data.frame(marker = c("m1", "m2"), chr = "1", pos = 1:2)
  for (bad in list(0, -1, NA_real_, NaN, Inf, c(1, 2))) {
    expect_error(qtl2mr_select_loci_from_scan(scan, marker_map, "traitA", lod_thresh = bad), "finite numeric value > 0")
  }
  bad_markers <- scan
  rownames(bad_markers) <- c("m1", "m1")
  expect_error(qtl2mr:::qtl2mr_lod_vector_from_scan(bad_markers, "traitA"), "unique marker")

  pair <- data.frame(lod_thresh = "not-a-number", stringsAsFactors = FALSE)
  expect_error(qtl2mr:::qtl2mr_pair_lod_thresh(pair, 4), "finite numeric value > 0")
  expect_equal(qtl2mr:::qtl2mr_pair_lod_thresh(data.frame(lod_thresh = NA_real_), 4), 4)
})

test_that("marker drop provenance distinguishes nonestimable and general fitting errors", {
  nonestimable <- qtl2mr:::qtl2mr_lmm_marker_error("m1", "1", "exposure", simpleError("Founder contrast coefficients are not estimable"))
  generic <- qtl2mr:::qtl2mr_lmm_marker_error("m2", "1", "outcome", simpleError("unexpected numerical failure"))
  nonestimable_row <- qtl2mr:::qtl2mr_lmm_marker_error_row(nonestimable)
  generic_row <- qtl2mr:::qtl2mr_lmm_marker_error_row(generic)
  expect_equal(nonestimable_row$drop_reason, "nonestimable_founder_contrast")
  expect_equal(generic_row$drop_reason, "marker_fit_error")
  expect_equal(generic_row$error_class, "simpleError")
  expect_match(generic_row$error_message, "unexpected numerical failure")
  expect_equal(qtl2mr:::qtl2mr_marker_drop_reason(list(error_message = NA_character_)), "marker_fit_error")
  expect_equal(qtl2mr:::qtl2mr_marker_drop_reason(list()), "marker_fit_error")
})

test_that("parallel group failures retain their group key without aggregation errors", {
  groups <- list(
    list(key = "traitA", pairs = data.frame(exposure = "traitA", outcome = "outcomeA")),
    list(key = "traitB", pairs = data.frame(exposure = "traitB", outcome = "outcomeB"))
  )
  failure <- qtl2mr:::qtl2mr_parallel_group_failure(groups[[1L]], simpleError("worker root cause"))
  successful <- list(run_manifest = data.frame(context_cache_hits = 0L))
  collected <- qtl2mr:::qtl2mr_collect_parallel_group_results(list(failure, successful), groups, fail_fast = FALSE)
  expect_equal(nrow(collected$failures), 1L)
  expect_equal(collected$failures$exposure_group, "traitA")
  expect_match(collected$failures$error_message, "worker root cause")
  expect_length(collected$successes, 1L)
  expect_error(
    qtl2mr:::qtl2mr_collect_parallel_group_results(list(failure), groups[1L], fail_fast = TRUE),
    "traitA.*worker root cause"
  )
})

test_that("stripped MR block cache refits first-stage and Steiger filters", {
  make_block <- function(marker, bx_value, by_value, sx = 0.05, sy = 0.10) {
    coords <- paste0("FC", 1:7)
    beta_x <- stats::setNames(rep(bx_value, 7), coords)
    beta_y <- stats::setNames(rep(by_value, 7), coords)
    Sx <- diag(sx^2, 7)
    Sy <- diag(sy^2, 7)
    dimnames(Sx) <- list(coords, coords)
    dimnames(Sy) <- list(coords, coords)
    list(
      marker = marker,
      chr = "1",
      beta_x = beta_x,
      beta_y = beta_y,
      Sigma_x = Sx,
      Sigma_y = Sy,
      Sigma_xy = matrix(0, 7, 7, dimnames = list(coords, coords)),
      n_overlap = 100L,
      O_xy_summary = list(n_x = 100L, n_y = 100L, n_overlap = 100L),
      cross_cov_model = "star",
      fit_x = list(n = 100L, df = 90, sigma2 = 1),
      fit_y = list(n = 100L, df = 90, sigma2 = 1)
    )
  }
  blocks <- list(
    m1 = make_block("m1", bx_value = 0.20, by_value = 0.12),
    m2 = make_block("m2", bx_value = 0.02, by_value = 0.01),
    m3 = make_block("m3", bx_value = 0.15, by_value = 0.50)
  )
  pilot <- list(
    locus_blocks = blocks,
    diagnostic_layer = list(per_locus = data.frame()),
    selected_loci = data.frame(marker = names(blocks), chr = "1", locus_rank = 1:3),
    summary = data.frame(theta = 0.5),
    overlap_summary = data.frame()
  )
  pair <- data.frame(pair_id = "x__to__y", exposure = "x", outcome = "y",
                     instrument_trait = "x", stringsAsFactors = FALSE)
  cache <- qtl2mr:::qtl2mr_make_mr_block_cache(pair, pilot)
  fs <- qtl2mr:::qtl2mr_mr_block_cache_first_stage_table(cache)
  expect_true(all(c("F_X_block", "F_Y_block", "partial_R2_X_block", "partial_R2_Y_block") %in% names(fs)))
  expect_lt(fs$F_X_block[fs$marker == "m2"], 2)
  direct <- fit_multilocus_founder_mr(blocks, heterogeneity_method = "none")
  from_cache <- refit_qtl2mr_mr_block_cache(cache, F_X_min = NA_real_)
  expect_equal(from_cache$theta, direct$theta, tolerance = 1e-10)
  f2 <- refit_qtl2mr_mr_block_cache(cache, F_X_min = 2)
  expect_equal(f2$n_loci_kept, 2L)
  expect_equal(f2$n_loci_dropped_F, 1L)
  steiger <- refit_qtl2mr_mr_block_cache(cache, F_X_min = NA_real_, steiger_filter = TRUE)
  expect_lt(steiger$n_loci_kept, from_cache$n_loci_kept)
})

test_that("prepared batch engine matches locus engine for primary fit", {
  gp <- make_sim_genoprobs(n = 74, marker = "m1", seed = 821)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.7 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.6 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  selected <- data.frame(marker = "m1", chr = "1", locus_rank = 1, stringsAsFactors = FALSE)
  input <- as_qtl2mr_prepared_input(gp, pheno, covar, kinship = NULL,
                                    marker_map = marker_map, selected_loci = selected,
                                    exposure = "x", outcome = "y")
  loc <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", engine = "locus",
                                      heterogeneity_method = "none",
                                      diagnostic_bootstrap_B = 0L)
  bat <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", engine = "batch",
                                      heterogeneity_method = "none",
                                      diagnostic_bootstrap_B = 0L)
  expect_equal(bat$summary$theta, loc$summary$theta, tolerance = 1e-10)
  expect_equal(bat$summary$se, loc$summary$se, tolerance = 1e-10)
  expect_equal(bat$summary$n_loci, 1L)
})

test_that("compact batch runner is fail-soft and writes compact outputs", {
  gp <- make_sim_genoprobs(n = 76, marker = "m1", seed = 831)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.6 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.5 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(x = matrix(5, nrow = 1, dimnames = list("m1", "lod")))
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(exposure = c("x", "missing"), outcome = c("y", "y"),
                      stringsAsFactors = FALSE)
  out_dir <- tempfile("qtl2mr-batch-")
  batch <- run_qtl2mr_batch(resources, pairs, output_dir = out_dir, progress = FALSE)
  expect_equal(nrow(batch$main_results), 2L)
  expect_true("ok" %in% batch$main_results$status)
  expect_true("error" %in% batch$main_results$status)
  expect_identical(
    batch$main_results$lmm_backend[batch$main_results$status == "ok"],
    "ordinary_least_squares"
  )
  expect_false(any(
    batch$main_results$lmm_backend[batch$main_results$status == "error"] == "internal_base_lmm",
    na.rm = TRUE
  ))
  expect_identical(batch$run_manifest$lmm_backend, "ordinary_least_squares")
  expected <- c(
    "batch_main_results.csv", "batch_method_comparison.csv", "batch_selected_loci.csv",
    "batch_overlap_qc.csv", "batch_first_stage_summary.csv", "batch_warnings.csv",
    "batch_errors.csv", "batch_run_manifest.csv", "batch_sample_qc.csv",
    "batch_config.json"
  )
  expect_true(all(file.exists(file.path(out_dir, expected))))
  expect_false(batch$run_manifest$sample_approximation)
  expect_false(file.exists(file.path(out_dir, "report.html")))
  expect_false(dir.exists(file.path(out_dir, "figures")))
})

test_that("batch_common_complete sample rule uses one sample set across pairs", {
  gp <- make_sim_genoprobs(n = 50, marker = "m1", seed = 841)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.6 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.5 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    z = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  pheno[ids[1:2], "y"] <- NA
  pheno[ids[3:4], "z"] <- NA
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(
    x = matrix(5, nrow = 1, dimnames = list("m1", "lod")),
    y = matrix(5, nrow = 1, dimnames = list("m1", "lod"))
  )
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(exposure = c("x", "y"), outcome = c("y", "z"),
                      stringsAsFactors = FALSE)
  batch <- run_qtl2mr_batch(
    resources, pairs,
    sample_rule = list(mode = "batch_common_complete"),
    fit_options = list(cache_contexts = TRUE, diagnostic_bootstrap_B = 0L, min_samples = 20L),
    progress = FALSE
  )
  expect_equal(batch$main_results$n_overlap_median, c(46, 46))
  expect_equal(batch$run_manifest$sample_mode, "batch_common_complete")
  expect_equal(batch$run_manifest$batch_common_n, 46L)
  expect_gt(batch$run_manifest$context_cache_hits, 0)
})

test_that("fast batch wrapper uses batch-common samples and caches", {
  gp <- make_sim_genoprobs(n = 60, marker = "m1", seed = 851)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.6 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.5 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    z = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(x = matrix(5, nrow = 1, dimnames = list("m1", "lod")))
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(exposure = c("x", "x"), outcome = c("y", "z"),
                      stringsAsFactors = FALSE)
  out_dir <- tempfile("qtl2mr-fast-")
  fast <- run_qtl2mr_batch_fast(resources, pairs, output_dir = out_dir, progress = FALSE)
  exact <- run_qtl2mr_batch(resources, pairs, progress = FALSE)

  expect_equal(fast$run_manifest$sample_mode, "batch_common_complete")
  expect_true(fast$run_manifest$sample_approximation)
  expect_equal(fast$sample_qc$sample_retention_ratio, 1)
  expect_equal(fast$run_manifest$sample_retention_gate_status, "pass")
  expect_gt(fast$run_manifest$lmm_fit_cache_hits, 0)
  expect_true(file.exists(file.path(out_dir, "batch_sample_qc.csv")))
  expect_equal(fast$main_results$theta, exact$main_results$theta, tolerance = 1e-10)
})

test_that("fast batch wrapper enforces and can override sample-retention gate", {
  gp <- make_sim_genoprobs(n = 60, marker = "m1", seed = 861)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.6 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    y = as.numeric(0.5 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    z = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5)),
    row.names = ids
  )
  pheno[ids[1:30], "y"] <- NA
  pheno[ids[26:50], "z"] <- NA
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(x = matrix(5, nrow = 1, dimnames = list("m1", "lod")))
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(exposure = c("x", "x"), outcome = c("y", "z"),
                      stringsAsFactors = FALSE)
  out_dir <- tempfile("qtl2mr-fast-stop-")
  expect_error(
    run_qtl2mr_batch_fast(resources, pairs, fit_options = list(min_samples = 9L),
                        output_dir = out_dir, progress = FALSE),
    "sample retention"
  )
  expect_true(file.exists(file.path(out_dir, "batch_sample_qc.csv")))
  stopped_qc <- utils::read.csv(file.path(out_dir, "batch_sample_qc.csv"))
  expect_equal(stopped_qc$sample_retention_gate_status, "fail")

  fast <- run_qtl2mr_batch_fast(
    resources, pairs,
    allow_sample_loss = TRUE,
    fit_options = list(min_samples = 9L),
    progress = FALSE
  )
  expect_equal(fast$sample_qc$sample_retention_gate_status, "override")
  expect_true("batch_common_sample_retention_low" %in% fast$warnings$warning_id)
})

test_that("parallel-by-exposure batch wrapper merges group results", {
  gp <- make_sim_genoprobs(n = 62, marker = "m1", seed = 871)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.7 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    w = as.numeric(0.5 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    y = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    z = as.numeric(0.3 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    row.names = ids
  )
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(
    x = matrix(5.5, nrow = 1, dimnames = list("m1", "lod")),
    w = matrix(5.2, nrow = 1, dimnames = list("m1", "lod"))
  )
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(
    exposure = c("x", "x", "w"),
    outcome = c("y", "z", "z"),
    stringsAsFactors = FALSE
  )
  out_dir <- tempfile("qtl2mr-parallel-")
  parallel_batch <- run_qtl2mr_batch_parallel_by_exposure(
    resources, pairs,
    n_workers = 1L,
    output_dir = out_dir,
    progress = FALSE
  )
  fast_batch <- run_qtl2mr_batch_fast(resources, pairs, progress = FALSE)

  expect_equal(parallel_batch$main_results$pair_id, fast_batch$main_results$pair_id)
  expect_equal(parallel_batch$main_results$theta, fast_batch$main_results$theta, tolerance = 1e-10)
  expect_equal(parallel_batch$run_manifest$parallel_strategy, "by_exposure")
  expect_equal(parallel_batch$run_manifest$parallel_backend, "sequential")
  expect_equal(parallel_batch$run_manifest$n_groups, 2L)
  expect_equal(parallel_batch$run_manifest$n_workers, 1L)
  expect_true(parallel_batch$run_manifest$sample_approximation)
  expect_true(dir.exists(file.path(out_dir, "by_exposure", "x")))
  expect_true(file.exists(file.path(out_dir, "batch_main_results.csv")))

  exact_parallel <- run_qtl2mr_batch_parallel_by_exposure(
    resources, pairs,
    fast = FALSE,
    n_workers = 1L,
    progress = FALSE
  )
  expect_false(exact_parallel$run_manifest$sample_approximation)
  expect_equal(exact_parallel$run_manifest$sample_mode, "common_complete")
})

test_that("exact grouped batch plan groups identical pair sample sets", {
  gp <- make_sim_genoprobs(n = 64, marker = "m1", seed = 881)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.7 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    w = as.numeric(0.5 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    y = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    z = as.numeric(0.3 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    row.names = ids
  )
  pheno[ids[1:4], c("y", "z")] <- NA
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(
    x = matrix(5.5, nrow = 1, dimnames = list("m1", "lod")),
    w = matrix(5.2, nrow = 1, dimnames = list("m1", "lod"))
  )
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(
    exposure = c("x", "x", "w"),
    outcome = c("y", "z", "z"),
    stringsAsFactors = FALSE
  )

  plan <- plan_qtl2mr_batch_exact_groups(resources, pairs)
  expect_s3_class(plan, "qtl2mr_exact_group_plan")
  expect_equal(plan$workload_summary$n_instrument_requests, 2L)
  expect_equal(plan$workload_summary$n_exact_sample_groups, 2L)
  expect_equal(
    plan$sample_groups$n_pairs[plan$sample_groups$exposure == "x"],
    2L
  )
  expect_false(plan$workload_summary$sample_approximation)
  expect_equal(
    unique(plan$pair_plan$sample_group_id[plan$pair_plan$exposure == "x"]),
    plan$pair_plan$sample_group_id[plan$pair_plan$pair_id == "x__to__y"]
  )
})

test_that("exact grouped batch runner matches exact lazy batch", {
  gp <- make_sim_genoprobs(n = 66, marker = "m1", seed = 891)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  pheno <- data.frame(
    x = as.numeric(0.7 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    w = as.numeric(0.5 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    y = as.numeric(0.4 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.45)),
    z = as.numeric(0.3 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.45)),
    row.names = ids
  )
  pheno[ids[1:3], c("y", "z")] <- NA
  covar <- data.frame(c1 = stats::rnorm(length(ids)), row.names = ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  scan <- list(
    x = matrix(5.5, nrow = 1, dimnames = list("m1", "lod")),
    w = matrix(5.2, nrow = 1, dimnames = list("m1", "lod"))
  )
  resources <- qtl2mr_prepare_resources(gp, kinship = NULL, pheno = pheno,
                                        covar = covar, qtl_scan = scan,
                                        marker_map = marker_map)
  pairs <- data.frame(
    exposure = c("x", "x", "w"),
    outcome = c("y", "z", "z"),
    stringsAsFactors = FALSE
  )
  out_dir <- tempfile("qtl2mr-exact-grouped-")
  grouped <- run_qtl2mr_batch_exact_grouped(resources, pairs, output_dir = out_dir, progress = FALSE)
  exact <- run_qtl2mr_batch(resources, pairs, progress = FALSE)

  expect_equal(grouped$main_results$pair_id, exact$main_results$pair_id)
  expect_equal(grouped$main_results$theta, exact$main_results$theta, tolerance = 1e-10)
  expect_equal(grouped$main_results$se, exact$main_results$se, tolerance = 1e-10)
  expect_equal(grouped$main_results$n_overlap_median, exact$main_results$n_overlap_median)
  expect_equal(grouped$run_manifest$execution_strategy, "exact_sample_grouped")
  expect_false(grouped$run_manifest$sample_approximation)
  expect_true(file.exists(file.path(out_dir, "batch_exact_group_plan.csv")))
  expect_true(file.exists(file.path(out_dir, "batch_sample_groups.csv")))
  expect_true(file.exists(file.path(out_dir, "batch_workload_summary.csv")))
})
