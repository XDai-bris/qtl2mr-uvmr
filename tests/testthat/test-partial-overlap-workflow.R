test_that("partial-overlap locus block uses O_xy sample matching", {
  gp <- make_sim_genoprobs(n = 90, marker = "m1", seed = 401)
  ids <- dimnames(gp[[1]])[[1]]
  H <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% H
  set.seed(402)
  x <- as.numeric(0.8 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.5))
  y <- as.numeric(0.6 * x + 0.2 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  ids_x <- ids[1:70]
  ids_y <- ids[21:90]
  overlap <- intersect(ids_x, ids_y)
  block <- fit_mr_locus_block_overlap(
    gp,
    exposure = x,
    outcome = y,
    contrasts = H,
    marker = "m1",
    sample_ids_x = ids_x,
    sample_ids_y = ids_y,
    overlap_ids = overlap,
    cross_cov_model = "star"
  )
  expect_equal(block$n_overlap, length(overlap))
  expect_equal(dim(block$Sigma_xy), c(7L, 7L))
  expect_true(isTRUE(block$joint_psd_valid))

  ids_y_perm <- rev(ids_y)
  O_perm <- make_overlap_matrix(ids_x, ids_y_perm, overlap_ids = overlap)
  block_perm <- fit_mr_locus_block_overlap(
    gp,
    exposure = x,
    outcome = y,
    contrasts = H,
    marker = "m1",
    sample_ids_x = ids_x,
    sample_ids_y = ids_y_perm,
    O_xy = O_perm,
    cross_cov_model = "star"
  )
  expect_equal(block_perm$n_overlap, length(overlap))
  expect_equal(block_perm$Sigma_xy, block$Sigma_xy, tolerance = 1e-8)
})

test_that("prepared partial-overlap with kinship uses star-space sample matching", {
  testthat::skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 90, marker = "m1", seed = 411)
  ids <- dimnames(gp[[1]])[[1]]
  H <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% H
  set.seed(412)
  x <- as.numeric(0.9 * Z[, "FC1"] + stats::rnorm(length(ids), sd = 0.55))
  y <- as.numeric(0.7 * x + 0.25 * Z[, "FC2"] + stats::rnorm(length(ids), sd = 0.55))
  pheno <- data.frame(x = x, y = y, row.names = ids)
  pheno$x[71:90] <- NA_real_
  pheno$y[1:20] <- NA_real_
  covar <- data.frame(age = as.numeric(scale(seq_along(ids))), row.names = ids)
  M <- matrix(stats::rnorm(length(ids) * 6), nrow = length(ids))
  K <- tcrossprod(M) / ncol(M)
  K <- cov2cor(K + diag(0.05, length(ids)))
  dimnames(K) <- list(ids, ids)
  marker_map <- data.frame(marker = "m1", chr = "1", pos = 1, stringsAsFactors = FALSE)
  input <- as_qtl2mr_prepared_input(
    allele_probs = gp,
    pheno = pheno,
    covar = covar,
    kinship = list("1" = K),
    marker_map = marker_map,
    selected_loci = "m1",
    exposure = "x",
    outcome = "y",
    hsq = list(exposure = 0.2, outcome = 0.3)
  )
  fit <- qtl2mr:::run_prepared_qtl2mr(input, exposure = "x", outcome = "y", loci = 1,
                                      write_outputs = FALSE, verbose = FALSE)
  block <- fit$locus_blocks[[1]]
  expect_equal(block$n_overlap, 50L)
  expect_equal(block$O_xy_summary$n_x, 70L)
  expect_equal(block$O_xy_summary$n_y, 70L)
  expect_equal(dim(block$Sigma_xy), c(7L, 7L))
  expect_true(isTRUE(block$joint_psd_valid))

  input_perm <- input
  input_perm$outcome_samples <- rev(ids[21:90])
  fit_perm <- qtl2mr:::run_prepared_qtl2mr(input_perm, exposure = "x", outcome = "y", loci = 1,
                                           write_outputs = FALSE, verbose = FALSE)
  expect_equal(fit_perm$locus_blocks[[1]]$Sigma_xy, block$Sigma_xy, tolerance = 1e-8)
})
