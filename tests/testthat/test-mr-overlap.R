test_that("overlap matrix preserves matching rather than row order", {
  O <- make_overlap_matrix(c("a", "b", "c"), c("a", "b", "c"))
  expected_identity <- diag(3)
  dimnames(expected_identity) <- list(c("a", "b", "c"), c("a", "b", "c"))
  expect_equal(O, expected_identity)

  nested <- make_overlap_matrix(c("a", "b", "c"), c("c", "a"))
  expected_nested <- matrix(
    c(0, 1,
      0, 0,
      1, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("c", "a"))
  )
  expect_equal(nested, expected_nested)

  reordered <- make_overlap_matrix(c("a", "b", "c"), c("c", "b", "a"))
  expect_equal(reordered["a", "a"], 1)
  expect_equal(reordered["b", "b"], 1)
  expect_equal(reordered["c", "c"], 1)
  expect_equal(sum(reordered), 3)

  disjoint <- make_overlap_matrix(c("a", "b"), c("c", "d"))
  expect_equal(sum(disjoint), 0)
  expect_error(make_overlap_matrix(c("a", "a"), c("a", "b")), "duplicate")
  expect_error(make_overlap_matrix(c("a", "b"), c("a", "b"), overlap_ids = c("a", "z")), "both")
})

test_that("bilinear rho_star matches cor and handles reordered/nested samples", {
  rx <- c(a = -1, b = 0, c = 1, d = 2)
  ry <- c(a = -2, b = 0, c = 2, d = 3)
  O <- make_overlap_matrix(names(rx), names(ry))
  expect_equal(estimate_rho_star(rx, ry, O, min_overlap_warn = 0L), stats::cor(rx, ry), tolerance = 1e-12)

  ry_reordered <- ry[c("d", "b", "a", "c")]
  O_reordered <- make_overlap_matrix(names(rx), names(ry_reordered))
  expect_equal(estimate_rho_star(rx, ry_reordered, O_reordered, min_overlap_warn = 0L), stats::cor(rx, ry), tolerance = 1e-12)

  O_nested <- make_overlap_matrix(names(rx), c("d", "b"))
  expect_equal(
    estimate_rho_star(rx, ry[c("d", "b")], O_nested, min_overlap_warn = 0L),
    stats::cor(rx[c("b", "d")], ry[c("b", "d")]),
    tolerance = 1e-12
  )
  expect_error(
    estimate_rho_star(c(a = 1, b = 1), c(a = 2, b = 3), make_overlap_matrix(c("a", "b"), c("a", "b")), min_overlap_warn = 0L),
    "zero"
  )
})

test_that("rho_star warns for very small nonzero overlap", {
  rx <- stats::setNames(seq_len(6), paste0("id", seq_len(6)))
  ry <- stats::setNames(seq_len(6) + c(0.1, -0.1, 0.05, -0.05, 0.02, -0.02), paste0("id", seq_len(6)))
  O <- make_overlap_matrix(names(rx), names(ry))
  expect_warning(estimate_rho_star(rx, ry, O), "Only 6 overlapping samples")

  rx0 <- stats::setNames(seq_len(6), paste0("x", seq_len(6)))
  ry0 <- stats::setNames(seq_len(6), paste0("y", seq_len(6)))
  O0 <- make_overlap_matrix(names(rx0), names(ry0))
  expect_error(estimate_rho_star(rx0, ry0, O0), "no usable sample overlap")
})

test_that("one-sample locus block estimates nonzero Sigma_xy with correlated residuals", {
  gp <- make_sim_genoprobs(n = 90)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  covar <- matrix(rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  e1 <- rnorm(length(ids), sd = 0.6)
  e2 <- 0.7 * e1 + rnorm(length(ids), sd = 0.3)
  x <- as.numeric(Z[, "FC1"] * 0.8 + covar[, 1] * 0.2 + e1)
  y <- as.numeric(0.5 * x + Z[, "FC3"] * 0.1 + covar[, 1] * 0.1 + e2)
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, addcovar = covar, contrasts = A)
  expect_equal(block$cross_cov_model, "star")
  expect_gt(abs(block$rho_xy), 0.2)
  expect_gt(abs(block$rho_xy_raw), 0.2)
  expect_gt(abs(block$rho_xy_star), 0.2)
  expect_gt(max(abs(block$Sigma_xy)), 0)
  common <- intersect(names(block$fit_x$residual_star), names(block$fit_y$residual_star))
  expected_star <- block$rho_xy_star * block$fit_x$sigma * block$fit_y$sigma *
    (block$fit_x$influence_star[, common, drop = FALSE] %*%
       t(block$fit_y$influence_star[, common, drop = FALSE]))
  expect_equal(block$Sigma_xy, expected_star, tolerance = 1e-12)
  expect_false(isTRUE(block$star_space_psd_adjusted))
  expect_lte(block$star_space_max_cross_singular_value, 1 + 1e-8)
  expect_equal(block$Sigma_xy_raw_space, block$Sigma_xy_star_space, tolerance = 1e-10)

  block0 <- fit_mr_locus_block(gp, x, y, addcovar = covar, contrasts = A, cross_cov_model = "zero")
  expect_equal(block0$rho_xy, 0)
  expect_equal(max(abs(block0$Sigma_xy)), 0)
  expect_gt(abs(block0$rho_xy_star), 0.2)
})

test_that("star-space cross-covariance is sample-order invariant with O_xy", {
  gp <- make_sim_genoprobs(n = 75)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  x <- as.numeric(Z[, "FC1"] * 0.5 + rnorm(length(ids), sd = 0.4))
  y <- as.numeric(Z[, "FC2"] * 0.2 + 0.3 * x + rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  fit_x <- fit_founder_lmm_block(gp, x, contrasts = A)
  fit_y <- fit_founder_lmm_block(gp, y, contrasts = A)
  base <- qtl2mr:::qtl2mr_lmm_cross_covariance(fit_x, fit_y, paste0("FC", 1:7), space = "star")

  perm <- rev(ids)
  fit_y_perm <- fit_y
  fit_y_perm$residual_star <- fit_y$residual_star[perm]
  fit_y_perm$residual_raw <- fit_y$residual_raw[perm]
  fit_y_perm$influence_star <- fit_y$influence_star[, perm, drop = FALSE]
  fit_y_perm$influence_raw <- fit_y$influence_raw[, perm, drop = FALSE]
  fit_y_perm$sample_ids <- perm
  reordered <- qtl2mr:::qtl2mr_lmm_cross_covariance(fit_x, fit_y_perm, paste0("FC", 1:7), space = "star")
  expect_equal(reordered$Sigma_xy, base$Sigma_xy, tolerance = 1e-12)
  expect_equal(reordered$rho_xy, base$rho_xy, tolerance = 1e-12)
})

test_that("overlap-aware multilocus MR equals independence when Sigma_xy is zero", {
  gp <- make_sim_genoprobs(n = 95)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  x <- as.numeric(Z[, "FC1"] * 0.8 + rnorm(length(ids), sd = 0.4))
  y <- as.numeric(0.6 * x + rnorm(length(ids), sd = 0.4))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = A)
  block$Sigma_xy[,] <- 0
  fit1 <- fit_multilocus_founder_mr(list(block), allow_cross_cov = TRUE, compute_independence = FALSE)
  fit0 <- fit_multilocus_founder_mr(list(block), allow_cross_cov = FALSE, compute_independence = FALSE)
  expect_equal(fit1$theta, fit0$theta, tolerance = 1e-10)
  expect_equal(fit1$Q, fit0$Q, tolerance = 1e-10)
})

test_that("large joint PSD violations are flagged rather than silently repaired", {
  Sigma_x <- diag(2)
  Sigma_y <- diag(2)
  Sigma_xy <- diag(c(1.2, 0.4))
  psd <- qtl2mr:::qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, Sigma_xy)
  expect_false(psd$psd_valid)
  expect_false(psd$psd_valid_before)
  expect_false(psd$psd_valid_after)
  expect_false(psd$adjusted)
  expect_equal(psd$status, "invalid_large")
  expect_equal(psd$Sigma_xy, Sigma_xy, tolerance = 1e-12)

  repaired <- qtl2mr:::qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, Sigma_xy, repair = "always")
  expect_true(repaired$psd_valid)
  expect_false(repaired$psd_valid_before)
  expect_true(repaired$psd_valid_after)
  expect_true(repaired$adjusted)
  expect_gt(repaired$smax_before, 1)
  expect_lte(repaired$smax_after, 1 + 1e-8)
  expect_gt(psd$max_cross_singular_value, 1)
  expect_lte(max(svd(diag(2) %*% repaired$Sigma_xy %*% diag(2))$d), 1 + 1e-8)
})

test_that("star-space Sigma_xy is the default LMM overlap model with kinship", {
  skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 80)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  covar <- matrix(rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  e1 <- rnorm(length(ids), sd = 0.5)
  e2 <- 0.65 * e1 + rnorm(length(ids), sd = 0.35)
  x <- as.numeric(Z[, "FC2"] * 0.7 + covar[, 1] * 0.2 + e1)
  y <- as.numeric(0.6 * x + Z[, "FC5"] * 0.15 + covar[, 1] * 0.1 + e2)
  names(x) <- ids
  names(y) <- ids
  kinship <- matrix(0.15, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  diag(kinship) <- 1

  block <- fit_mr_locus_block(
    gp, x, y,
    kinship = kinship,
    addcovar = covar,
    contrasts = A,
    hsq_x = 0.45,
    hsq_y = 0.40
  )
  expect_equal(block$cross_cov_model, "star")
  expect_equal(block$rho_xy, block$rho_xy_star, tolerance = 1e-12)
  expect_gt(abs(block$rho_xy_raw), 0.2)
  expect_false(isTRUE(block$star_space_psd_adjusted))
  expect_lte(block$star_space_max_cross_singular_value, 1 + 1e-8)
  expect_true(is.finite(block$raw_space_max_cross_singular_value))
  common <- intersect(names(block$fit_x$residual_star), names(block$fit_y$residual_star))
  expected_star <- block$rho_xy_star * block$fit_x$sigma * block$fit_y$sigma *
    (block$fit_x$influence_star[, common, drop = FALSE] %*%
       t(block$fit_y$influence_star[, common, drop = FALSE]))
  expect_equal(block$Sigma_xy, expected_star, tolerance = 1e-12)
  expect_gt(max(abs(block$Sigma_xy_raw_space - block$Sigma_xy_star_space)), 1e-8)
})
