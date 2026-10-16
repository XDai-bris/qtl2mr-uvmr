make_lmm_backend_fixture <- function(n = 72L, seed = 7101L) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  A <- matrix(stats::rnorm(n * 6L), nrow = n)
  K <- tcrossprod(A) / ncol(A) + diag(0.15, n)
  K <- K / mean(diag(K))
  dimnames(K) <- list(ids, ids)
  covar <- cbind(c1 = seq(-1, 1, length.out = n), c2 = stats::rnorm(n))
  rownames(covar) <- ids
  L <- chol(0.45 * (2 * K) + 0.55 * diag(n))
  y <- as.numeric(covar %*% c(0.4, -0.2) + t(L) %*% stats::rnorm(n))
  names(y) <- ids
  list(ids = ids, K = K, covar = covar, y = y)
}

test_that("internal kinship decomposition reconstructs 2K across spectra", {
  fx <- make_lmm_backend_fixture(n = 24L)
  dec <- qtl2mr_decompose_kinship(fx$K, fx$ids)
  reconstructed <- t(dec$vectors) %*% diag(dec$values) %*% dec$vectors
  expect_equal(reconstructed, 2 * fx$K, tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(all(diff(dec$values) >= 0))
  expect_identical(dec$backend, "internal_base_lmm")

  repeated <- diag(c(0, 0.5, 0.5, 1, 1, 1), 6L)
  dimnames(repeated) <- list(letters[1:6], letters[1:6])
  dec_repeated <- qtl2mr_decompose_kinship(repeated)
  expect_equal(
    t(dec_repeated$vectors) %*% diag(dec_repeated$values) %*% dec_repeated$vectors,
    2 * repeated,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("kinship validation clamps only scale-small numerical negatives", {
  ids <- paste0("i", 1:4)
  near_psd <- diag(c(1, 0.5, 0.2, -1e-10))
  dimnames(near_psd) <- list(ids, ids)
  dec <- qtl2mr_decompose_kinship(near_psd)
  expect_equal(min(dec$values), 0)
  expect_equal(dec$clamped_negative_eigenvalues, 1)

  bad <- diag(c(1, 0.5, 0.2, -1e-3))
  dimnames(bad) <- list(ids, ids)
  expect_error(qtl2mr_decompose_kinship(bad), "positive semidefinite")
  asymmetric <- diag(4)
  asymmetric[1, 2] <- 0.1
  dimnames(asymmetric) <- list(ids, ids)
  expect_error(qtl2mr_decompose_kinship(asymmetric), "symmetric")
  nonfinite <- diag(4)
  nonfinite[1, 1] <- Inf
  dimnames(nonfinite) <- list(ids, ids)
  expect_error(qtl2mr_decompose_kinship(nonfinite), "finite")

  extreme <- diag(c(1e120, 3e119, 1e118, 0), 4L)
  dimnames(extreme) <- list(ids, ids)
  dec_extreme <- qtl2mr_decompose_kinship(extreme)
  expect_equal(
    t(dec_extreme$vectors) %*% diag(dec_extreme$values) %*% dec_extreme$vectors,
    2 * extreme,
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("kinship sample identifiers are validated and permutations align", {
  fx <- make_lmm_backend_fixture(n = 20L)
  perm <- rev(fx$ids)
  dec <- qtl2mr_decompose_kinship(fx$K, perm)
  expect_identical(colnames(dec$vectors), perm)
  expect_equal(
    t(dec$vectors) %*% diag(dec$values) %*% dec$vectors,
    2 * fx$K[perm, perm],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  unnamed <- unname(fx$K)
  expect_error(qtl2mr_decompose_kinship(unnamed), "sample IDs")
  duplicated <- fx$K
  rownames(duplicated)[2] <- rownames(duplicated)[1]
  expect_error(qtl2mr_decompose_kinship(duplicated), "unique")
  expect_error(qtl2mr_decompose_kinship(fx$K, c(fx$ids[-1], "missing")), "missing sample IDs")
})

test_that("internal ML and REML hsq agree with optional qtl2", {
  skip_if_not_installed("qtl2")
  fx <- make_lmm_backend_fixture()
  for (reml in c(FALSE, TRUE)) {
    for (covariates in list(NULL, fx$covar)) {
      internal <- qtl2mr_estimate_lmm_hsq(fx$y, fx$K, covariates, reml = reml)
      reference <- qtl2::est_herit(fx$y, fx$K, addcovar = covariates, reml = reml)
      expect_lte(abs(as.numeric(internal) - as.numeric(reference)), 1e-6)
      expect_identical(attr(internal, "backend"), "internal_base_lmm")
      expect_true(attr(internal, "status") %in% c("success", "boundary_zero", "boundary_one"))
    }
  }
})

test_that("hsq estimation handles covariates, missing rows, permutations, and boundaries", {
  fx <- make_lmm_backend_fixture(n = 64L)
  rank_deficient <- cbind(fx$covar, duplicate = fx$covar[, 1])
  fit_rank <- NULL
  expect_warning(
    fit_rank <- qtl2mr_estimate_lmm_hsq(fx$y, fx$K, rank_deficient, reml = TRUE),
    "linearly dependent"
  )
  expect_true(is.finite(fit_rank))

  gp <- make_sim_genoprobs(n = 64L, seed = 7102L)
  ids <- dimnames(gp[[1]])[[1]]
  K <- fx$K
  dimnames(K) <- list(ids, ids)
  covar <- rank_deficient
  rownames(covar) <- ids
  y <- fx$y
  names(y) <- ids
  y[c(3, 8)] <- NA_real_
  covar[5, 2] <- NA_real_
  block <- NULL
  expect_warning(
    block <- fit_founder_lmm_block(gp, y, kinship = K, addcovar = covar, hsq = NULL),
    "linearly dependent"
  )
  expect_equal(block$n, 61L)
  expect_true(is.finite(block$hsq))
  expect_identical(block$lmm_backend, "internal_base_lmm")

  keep <- block$sample_ids
  perm <- rev(keep)
  block_perm <- NULL
  expect_warning(
    block_perm <- fit_founder_lmm_block(
      gp, y, kinship = K, addcovar = covar, hsq = block$hsq, sample_ids = perm
    ),
    "linearly dependent"
  )
  expect_equal(block_perm$beta_fc, block$beta_fc, tolerance = 1e-9)
  expect_equal(block_perm$vcov_fc, block$vcov_fc, tolerance = 1e-9)

  boundary_ids <- paste0("b", seq_len(40L))
  boundary_K <- diag(seq(0.05, 2, length.out = 40L))
  dimnames(boundary_K) <- list(boundary_ids, boundary_ids)
  ascending <- seq(0, 1, length.out = 40L)^8 + 1e-4 * seq_len(40L)
  h_one <- qtl2mr_estimate_lmm_hsq(stats::setNames(ascending, boundary_ids), boundary_K, reml = TRUE)
  h_zero <- qtl2mr_estimate_lmm_hsq(stats::setNames(rev(ascending), boundary_ids), boundary_K, reml = TRUE)
  expect_equal(as.numeric(h_one), 1)
  expect_identical(attr(h_one, "status"), "boundary_one")
  expect_equal(as.numeric(h_zero), 0)
  expect_identical(attr(h_zero, "status"), "boundary_zero")

  rank_deficient_K <- diag(c(rep(1, 35L), rep(0, 5L)))
  dimnames(rank_deficient_K) <- list(boundary_ids, boundary_ids)
  h_rank_deficient <- qtl2mr_estimate_lmm_hsq(
    stats::setNames(sin(seq_len(40L) / 4), boundary_ids), rank_deficient_K, reml = TRUE
  )
  expect_true(is.finite(h_rank_deficient))
  expect_true(h_rank_deficient >= 0 && h_rank_deficient <= 1)
})

test_that("internal backend remains practical for a LOCO-scale sample count", {
  n <- 760L
  ids <- paste0("loco", seq_len(n))
  K <- diag(seq(0.7, 1.3, length.out = n))
  dimnames(K) <- list(ids, ids)
  y <- stats::setNames(sin(seq_len(n) / 13) + stats::rnorm(n, sd = 0.1), ids)
  fit <- qtl2mr_estimate_lmm_hsq(y, K, reml = TRUE)
  expect_true(is.finite(fit))
  expect_true(fit >= 0 && fit <= 1)
})

test_that("internal REML recognizes a deterministic near-one hsq simulation", {
  set.seed(7105L)
  n <- 240L
  ids <- paste0("near", seq_len(n))
  A <- matrix(stats::rnorm(n * 20L), nrow = n)
  K <- tcrossprod(A) / ncol(A) + diag(0.01, n)
  K <- K / mean(diag(K))
  dimnames(K) <- list(ids, ids)
  true_hsq <- 0.98
  covariance <- true_hsq * (2 * K) + (1 - true_hsq) * diag(n)
  y <- stats::setNames(as.numeric(t(chol(covariance)) %*% stats::rnorm(n)), ids)
  estimated <- qtl2mr_estimate_lmm_hsq(y, K, reml = TRUE)
  expect_gt(estimated, 0.9)
  expect_lte(estimated, 1)
})

test_that("estimated-hsq direct fit preserves qtl2 beta SE and vcov parity", {
  skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 84L, seed = 7103L)
  ids <- dimnames(gp[[1]])[[1]]
  fx <- make_lmm_backend_fixture(n = 84L, seed = 7104L)
  K <- fx$K
  dimnames(K) <- list(ids, ids)
  covar <- fx$covar
  rownames(covar) <- ids
  A <- make_founder_contrasts(LETTERS[1:8])
  z <- gp[[1]][, , 1] %*% A
  y <- stats::setNames(as.numeric(0.5 * z[, "FC2"] + 0.3 * covar[, 1] + fx$y), ids)
  h_internal <- qtl2mr_estimate_lmm_hsq(y, K, covar, reml = TRUE)
  h_qtl2 <- qtl2::est_herit(y, K, addcovar = covar, reml = TRUE)
  expect_lte(abs(as.numeric(h_internal) - as.numeric(h_qtl2)), 1e-6)
  internal <- fit_founder_lmm_block(gp, y, kinship = K, addcovar = covar, contrasts = A, hsq = h_internal)
  reference <- compare_founder_lmm_to_scan1coef(
    gp, y, kinship = K, addcovar = covar, contrasts = A, hsq = h_qtl2
  )
  expect_equal(internal$beta_fc, reference$qtl2_beta, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(internal$se_fc, reference$qtl2_se, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(diag(internal$vcov_fc), reference$qtl2_se^2, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("hsq and decomposition invalid inputs fail with package diagnostics", {
  expect_error(qtl2mr_validate_hsq(NA_real_), "one finite numeric")
  expect_error(qtl2mr_validate_hsq(c(0.2, 0.3)), "one finite numeric")
  expect_error(qtl2mr_validate_hsq(TRUE), "one finite numeric")
  fx <- make_lmm_backend_fixture(n = 16L)
  bad_y <- fx$y
  bad_y[1] <- Inf
  expect_error(qtl2mr_estimate_lmm_hsq(bad_y, fx$K), "finite")
})

test_that("qtl2 is optional package metadata and backend provenance is stable", {
  metadata <- utils::packageDescription("qtl2mr")
  imports <- strsplit(metadata$Imports, "[,[:space:]]+")[[1L]]
  suggests <- strsplit(metadata$Suggests, "[,[:space:]]+")[[1L]]
  expect_false("qtl2" %in% imports)
  expect_true("qtl2" %in% suggests)
  expect_false("qtl2" %in% names(getNamespaceImports("qtl2mr")))
  provenance <- qtl2mr_lmm_backend_provenance()
  expect_identical(provenance$backend, "not_available")
  expect_length(provenance$fit_backends, 0L)
  expect_identical(provenance$qtl2_role, "optional_parity_validation")

  ols <- list(a = list(fit_x = list(lmm_backend = "ordinary_least_squares"),
                       fit_y = list(lmm_backend = "ordinary_least_squares")))
  internal <- list(a = list(fit_x = list(lmm_backend = "internal_base_lmm"),
                            fit_y = list(lmm_backend = "internal_base_lmm")))
  mixed <- c(ols, internal)
  expect_identical(qtl2mr_lmm_backend_provenance(ols)$backend, "ordinary_least_squares")
  expect_identical(qtl2mr_lmm_backend_provenance(internal)$backend, "internal_base_lmm")
  expect_identical(qtl2mr_lmm_backend_provenance(mixed)$backend, "mixed")
  expect_identical(qtl2mr_lmm_backend_provenance(mixed)$fit_backends,
                   c("internal_base_lmm", "ordinary_least_squares"))
})

test_that("generic qtl2 validator reports compact explicit-input parity", {
  skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 72L, marker = "m1", seed = 7120L)
  ids <- dimnames(gp[[1L]])[[1L]]
  set.seed(7121L)
  y <- stats::setNames(stats::rnorm(length(ids)), ids)
  out <- validate_qtl2_scan1coef(gp, y, marker = "m1", tolerance = 1e-6)
  expect_identical(out$validation_status, "pass")
  expect_true(out$validation_available)
  expect_identical(out$lmm_backend, "ordinary_least_squares")
  expect_lte(out$max_beta_abs_diff, 1e-6)
  expect_lte(out$max_se_abs_diff, 1e-6)
  expect_true(nzchar(out$qtl2_version))
})

test_that("generic qtl2 validator fails explicitly when optional qtl2 is absent", {
  gp <- make_sim_genoprobs(n = 24L, marker = "m1", seed = 7122L)
  ids <- dimnames(gp[[1L]])[[1L]]
  y <- stats::setNames(seq_along(ids), ids)
  testthat::with_mocked_bindings(
    expect_error(validate_qtl2_scan1coef(gp, y, marker = "m1"), "optional qtl2 unavailable"),
    qtl2mr_require_qtl2 = function() stop("optional qtl2 unavailable", call. = FALSE),
    .package = "qtl2mr"
  )
})
