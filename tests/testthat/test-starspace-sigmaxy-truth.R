test_that("star-space Sigma_xy equals residual-correlation influence formula", {
  ids <- paste0("id", 1:12)
  fc <- paste0("FC", 1:7)
  set.seed(521)
  Ax <- matrix(rnorm(length(fc) * length(ids), sd = 0.08), nrow = length(fc),
               dimnames = list(fc, ids))
  Ay <- matrix(rnorm(length(fc) * length(ids), sd = 0.07), nrow = length(fc),
               dimnames = list(fc, ids))
  rx <- stats::setNames(seq(-1.2, 1.4, length.out = length(ids)), ids)
  ry <- stats::setNames(seq(-0.8, 1.2, length.out = length(ids)) + rep(c(-0.05, 0.04), length.out = length(ids)), ids)
  fit_x <- list(
    residual_star = rx,
    influence_star = Ax,
    sigma = 1.7,
    sample_ids = ids
  )
  fit_y <- list(
    residual_star = ry,
    influence_star = Ay,
    sigma = 0.9,
    sample_ids = ids
  )
  O <- make_overlap_matrix(ids, ids)

  out <- qtl2mr:::qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc, space = "star", O_xy = O)
  rho <- estimate_rho_star(rx, ry, O, min_overlap_warn = 0L)
  expected <- rho * fit_x$sigma * fit_y$sigma * (Ax %*% O %*% t(Ay))
  dimnames(expected) <- list(fc, fc)

  expect_equal(out$rho_xy, rho, tolerance = 1e-12)
  expect_equal(out$Sigma_xy, expected, tolerance = 1e-12)
})

test_that("public PSD validator accepts valid joint covariance and rejects invalid cross-covariance", {
  Sigma_x <- diag(c(0.2, 0.3, 0.5))
  Sigma_y <- diag(c(0.4, 0.6, 0.8))
  valid_cross <- 0.2 * qtl2mr:::qtl2mr_matrix_power(Sigma_x, 0.5) %*%
    qtl2mr:::qtl2mr_matrix_power(Sigma_y, 0.5)
  invalid_cross <- 6.0 * valid_cross

  valid <- validate_joint_psd(Sigma_x, Sigma_y, valid_cross)
  invalid <- validate_joint_psd(Sigma_x, Sigma_y, invalid_cross)

  expect_true(valid$psd_valid)
  expect_true(valid$canonical_pass)
  expect_false(invalid$psd_valid)
  expect_false(invalid$canonical_pass)
})
