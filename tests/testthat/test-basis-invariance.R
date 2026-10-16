test_that("multi-locus SO estimate is invariant to orthonormal founder-contrast basis rotation", {
  gp1 <- make_sim_genoprobs(n = 96, marker = "m1", seed = 511)
  set.seed(512)
  raw2 <- matrix(stats::rexp(96 * 8), nrow = 96, ncol = 8)
  probs2 <- raw2 / rowSums(raw2)
  arr <- array(
    NA_real_,
    dim = c(96, 8, 2),
    dimnames = list(dimnames(gp1[[1]])[[1]], LETTERS[1:8], c("m1", "m2"))
  )
  arr[, , "m1"] <- gp1[[1]][, , "m1"]
  arr[, , "m2"] <- probs2
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  ids <- dimnames(gp[[1]])[[1]]

  A <- make_founder_contrasts(LETTERS[1:8])
  qr_rot <- qr(matrix(rnorm(49), nrow = 7))
  R <- qr.Q(qr_rot)
  A_rot <- A
  A_rot[, -1L] <- A[, -1L] %*% R
  colnames(A_rot) <- colnames(A)
  rownames(A_rot) <- rownames(A)

  Z1 <- gp[[1]][, , "m1"] %*% A
  Z2 <- gp[[1]][, , "m2"] %*% A
  x <- as.numeric(0.55 * Z1[, "FC1"] - 0.35 * Z2[, "FC3"] + rnorm(length(ids), sd = 0.45))
  y <- as.numeric(0.72 * x + 0.15 * Z1[, "FC4"] + rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids

  blocks <- list(
    m1 = fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m1"),
    m2 = fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m2")
  )
  blocks_rot <- list(
    m1 = fit_mr_locus_block(gp, x, y, contrasts = A_rot, marker = "m1"),
    m2 = fit_mr_locus_block(gp, x, y, contrasts = A_rot, marker = "m2")
  )

  fit <- fit_multilocus_founder_mr(blocks, heterogeneity_method = "chisq")
  fit_rot <- fit_multilocus_founder_mr(blocks_rot, heterogeneity_method = "chisq")

  expect_equal(fit_rot$theta, fit$theta, tolerance = 1e-8)
  expect_equal(fit_rot$Q, fit$Q, tolerance = 1e-7)
  expect_equal(fit_rot$se, fit$se, tolerance = 1e-7)
  expect_equal(fit_rot$profile_ci_lower, fit$profile_ci_lower, tolerance = 1e-6)
  expect_equal(fit_rot$profile_ci_upper, fit$profile_ci_upper, tolerance = 1e-6)
})

make_basis_test_block <- function(seed, marker) {
  set.seed(seed)
  k <- 7L
  Ax <- matrix(rnorm(k * k, sd = 0.15), k, k)
  Ay <- matrix(rnorm(k * k, sd = 0.15), k, k)
  beta_x <- rnorm(k, sd = 3)
  list(
    marker = marker,
    chr = "1",
    beta_x = beta_x,
    beta_y = 0.55 * beta_x + rnorm(k, sd = 0.05),
    Sigma_x = crossprod(Ax) + diag(0.5, k),
    Sigma_y = crossprod(Ay) + diag(0.6, k),
    Sigma_xy = matrix(rnorm(k * k, sd = 0.001), k, k),
    founder_H = rbind(diag(k), -rep(1, k))
  )
}

rotate_basis_test_block <- function(block, Q) {
  block$beta_x <- as.vector(t(Q) %*% block$beta_x)
  block$beta_y <- as.vector(t(Q) %*% block$beta_y)
  block$Sigma_x <- t(Q) %*% block$Sigma_x %*% Q
  block$Sigma_y <- t(Q) %*% block$Sigma_y %*% Q
  block$Sigma_xy <- t(Q) %*% block$Sigma_xy %*% Q
  block$founder_H <- block$founder_H %*% Q
  block
}

test_that("finite-ridge SO fits are invariant under dense, permutation, and sign rotations", {
  blocks <- list(
    make_basis_test_block(901, "m1"),
    make_basis_test_block(902, "m2"),
    make_basis_test_block(904, "m3")
  )
  ridge <- 0.05
  set.seed(903)
  transforms <- list(
    dense_1 = qr.Q(qr(matrix(rnorm(49), 7, 7))),
    dense_2 = qr.Q(qr(matrix(rnorm(49), 7, 7))),
    permutation = diag(7)[c(7, 2, 5, 1, 6, 3, 4), ],
    sign = diag(c(-1, 1, -1, 1, -1, 1, -1))
  )
  base <- fit_multilocus_founder_mr(blocks,
    ridge = ridge, heterogeneity_method = "none", q_sim_n = 0L)
  expect_equal(sum(base$locus_diagnostics$Q_locus), base$Q, tolerance = 1e-12)
  for (Q in transforms) {
    rotated <- lapply(blocks, rotate_basis_test_block, Q = Q)
    fit <- fit_multilocus_founder_mr(rotated,
      ridge = ridge, heterogeneity_method = "none", q_sim_n = 0L)
    expect_lte(abs(fit$theta - base$theta), 1e-8)
    expect_lte(abs(fit$se - base$se), 1e-7)
    expect_lte(abs(fit$profile_ci_lower - base$profile_ci_lower), 1e-7)
    expect_lte(abs(fit$profile_ci_upper - base$profile_ci_upper), 1e-7)
    expect_lte(abs(fit$profile_p - base$profile_p), 1e-7)
    expect_lte(abs(fit$Q - base$Q), 1e-8)
    expect_lte(abs(fit$objective_min - base$objective_min), 1e-8)
    expect_lte(abs(fit$theta_start - base$theta_start), 1e-10)
    expect_identical(fit$bracket_lower, base$bracket_lower)
    expect_identical(fit$bracket_upper, base$bracket_upper)
    expect_identical(fit$optimizer_status, base$optimizer_status)
    expect_identical(fit$bracket_status, base$bracket_status)
    expect_identical(fit$bracket_n_expand, base$bracket_n_expand)
    expect_equal(sum(fit$locus_diagnostics$Q_locus), fit$Q, tolerance = 1e-12)
  }
})

test_that("single-locus SO and regularized Omega are basis invariant", {
  block <- make_basis_test_block(911, "m1")
  set.seed(912)
  Q <- qr.Q(qr(matrix(rnorm(49), 7, 7)))
  rotated <- rotate_basis_test_block(block, Q)
  fit <- second_order_gls_ivw(block$beta_x, block$beta_y, block$Sigma_x, block$Sigma_y,
    allow_cross_cov = TRUE, Sigma_xy = block$Sigma_xy, ridge = 0.05)
  fit_rot <- second_order_gls_ivw(rotated$beta_x, rotated$beta_y, rotated$Sigma_x, rotated$Sigma_y,
    allow_cross_cov = TRUE, Sigma_xy = rotated$Sigma_xy, ridge = 0.05)
  expect_lte(abs(fit_rot$theta - fit$theta), 1e-8)
  expect_lte(abs(fit_rot$se - fit$se), 1e-7)
  expect_lte(abs(fit_rot$profile_ci_lower - fit$profile_ci_lower), 1e-7)
  expect_lte(abs(fit_rot$profile_ci_upper - fit$profile_ci_upper), 1e-7)
  expect_lte(abs(fit_rot$Q - fit$Q), 1e-8)
  expect_lte(abs(fit_rot$objective_min - fit$objective_min), 1e-8)
  expect_lte(abs(fit_rot$theta_start - fit$theta_start), 1e-10)
  expect_identical(fit_rot$optimizer_status, fit$optimizer_status)
  expect_identical(fit_rot$bracket_status, fit$bracket_status)

  theta <- 0.37
  Om <- mr_block_omega(block, theta, ridge = 0.05)
  Om_rot <- mr_block_omega(rotated, theta, ridge = 0.05)
  scale <- max(1, max(abs(Om)))
  expect_lte(max(abs(Om_rot - t(Q) %*% Om %*% Q)), 1e-12 * scale)
  expect_identical(qtl2mr_add_invariant_ridge(Om, 0), Om)
  expect_equal(
    mr_block_q(rotated, theta, ridge = 0.05),
    mr_block_q(block, theta, ridge = 0.05),
    tolerance = 1e-10
  )
})

test_that("invariant ridge handles valid edge scales and rejects invalid inputs", {
  expect_identical(qtl2mr_add_invariant_ridge(matrix(3, 1, 1), 0), matrix(3, 1, 1))
  expect_equal(qtl2mr_add_invariant_ridge(matrix(3, 1, 1), 0.2), matrix(3.6, 1, 1))

  zero <- qtl2mr_add_invariant_ridge(matrix(0, 3, 3), 1e-10)
  expect_true(all(is.finite(zero)))
  expect_true(all(diag(zero) > 0))

  near_singular <- diag(c(1, 1e-16, 0))
  regularized <- qtl2mr_add_invariant_ridge(near_singular, 1e-10)
  expect_true(all(eigen(regularized, symmetric = TRUE, only.values = TRUE)$values > 0))

  extreme <- diag(c(1e300, 5e299))
  extreme_regularized <- qtl2mr_add_invariant_ridge(extreme, 1e-10)
  expect_true(all(is.finite(extreme_regularized)))
  expect_equal(
    diag(extreme_regularized) / 1e300,
    c(1 + 7.5e-11, 0.5 + 7.5e-11),
    tolerance = 1e-14
  )

  expect_error(qtl2mr_add_invariant_ridge(matrix(1, 2, 3), 1e-10), "square")
  expect_error(qtl2mr_add_invariant_ridge(matrix("x", 1, 1), 1e-10), "numeric")
  expect_error(qtl2mr_add_invariant_ridge(matrix(TRUE, 1, 1), 1e-10), "numeric")
  expect_error(qtl2mr_add_invariant_ridge(matrix(1 + 0i, 1, 1), 1e-10), "numeric")
  expect_error(qtl2mr_add_invariant_ridge(matrix(factor("x"), 1, 1), 1e-10), "numeric")
  expect_error(qtl2mr_add_invariant_ridge(matrix(c(1, NA, NA, 1), 2), 1e-10), "finite")
  expect_error(qtl2mr_add_invariant_ridge(diag(2), -1), "non-negative")
  expect_error(qtl2mr_add_invariant_ridge(diag(2), c(0, 1)), "one finite")
  expect_error(qtl2mr_add_invariant_ridge(diag(2), TRUE), "one finite")
  expect_error(qtl2mr_add_invariant_ridge(diag(2), 1 + 0i), "one finite")
  expect_error(qtl2mr_add_invariant_ridge(diag(2), factor("1")), "one finite")
  expect_error(qtl2mr_add_invariant_ridge(matrix(c(1, 0.2, 0, 1), 2), 1e-10), "symmetric")
  expect_error(qtl2mr_add_invariant_ridge(matrix(c(1, 0.2, 0, 1), 2), 0), "symmetric")
  expect_error(qtl2mr_add_invariant_ridge(diag(c(1, -0.1)), 1e-10), "positive semidefinite")
  expect_error(qtl2mr_add_invariant_ridge(diag(c(1, -0.1)), 0), "positive semidefinite")
})

test_that("trace-scaled ridge is equivariant at near-singular scale", {
  set.seed(921)
  Q <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  Omega <- Q %*% diag(c(4, 1, 1e-8, 1e-14, 0)) %*% t(Q)
  set.seed(922)
  U <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  lhs <- qtl2mr_add_invariant_ridge(t(U) %*% Omega %*% U, 1e-8)
  rhs <- t(U) %*% qtl2mr_add_invariant_ridge(Omega, 1e-8) %*% U
  expect_equal(lhs, rhs, tolerance = 1e-13)
})
