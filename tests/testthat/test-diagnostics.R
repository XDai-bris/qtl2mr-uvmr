test_that("block diagnostics match manuscript formulas", {
  gp <- make_sim_genoprobs(n = 95, seed = 11)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  x <- as.numeric(0.8 * Z[, "FC1"] - 0.25 * Z[, "FC4"] + rnorm(length(ids), sd = 0.45))
  y <- as.numeric(0.6 * x + 0.2 * Z[, "FC2"] + rnorm(length(ids), sd = 0.45))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = A)
  theta <- 0.5

  expected_omega <- qtl2mr:::qtl2mr_sym(
    block$Sigma_y + theta^2 * block$Sigma_x - theta * (block$Sigma_xy + t(block$Sigma_xy))
  )
  expect_equal(mr_block_omega(block, theta), expected_omega, tolerance = 1e-12)

  d <- matrix(block$beta_y - theta * block$beta_x, ncol = 1)
  expected_omega_ridged <- qtl2mr:::qtl2mr_add_invariant_ridge(expected_omega, 1e-10)
  expected_q <- as.numeric(t(d) %*% qtl2mr:::qtl2mr_mat_inv(expected_omega_ridged, ridge = 0) %*% d)
  expect_equal(mr_block_q(block, theta), expected_q, tolerance = 1e-10)

  expected_f <- as.numeric(t(matrix(block$beta_x, ncol = 1)) %*%
    qtl2mr:::qtl2mr_mat_inv(block$Sigma_x) %*% matrix(block$beta_x, ncol = 1)) / 7
  expect_equal(mr_block_first_stage_f(block), expected_f, tolerance = 1e-10)
  expect_gt(mr_block_first_stage_f(block), 0)
})

test_that("diagnostics separate raw sampling covariance from ridge working covariance", {
  gp <- make_sim_genoprobs(n = 95, seed = 931)
  ids <- dimnames(gp[[1]])[[1]]
  H <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% H
  set.seed(932)
  x <- as.numeric(0.7 * Z[, "FC1"] - 0.3 * Z[, "FC4"] + rnorm(length(ids), sd = 0.4))
  y <- as.numeric(0.6 * x + 0.2 * Z[, "FC2"] + rnorm(length(ids), sd = 0.4))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = H)
  theta <- 0.4
  ridge <- 0.5

  raw <- mr_block_omega(block, theta, ridge = 0)
  work <- qtl2mr:::qtl2mr_add_invariant_ridge(raw, ridge)
  d <- matrix(block$beta_y - theta * block$beta_x, ncol = 1L)
  expected_q <- as.numeric(t(d) %*% qtl2mr:::qtl2mr_mat_inv(work, ridge = 0) %*% d)
  expect_equal(mr_block_q(block, theta, ridge = ridge), expected_q, tolerance = 1e-12)

  diag_zero <- mr_block_so_diagnostics(block, theta = theta, ridge = 0)
  diag_ridge <- mr_block_so_diagnostics(block, theta = theta, ridge = ridge)
  expect_identical(diag_ridge$omega_min_eigenvalue, diag_zero$omega_min_eigenvalue)
  expect_identical(diag_ridge$omega_condition_number, diag_zero$omega_condition_number)
  expect_false(isTRUE(all.equal(diag_ridge$Q_min, diag_zero$Q_min)))

  residual_cov <- qtl2mr:::mr_block_residual_covariance(block, theta, ridge = ridge)
  expected_vr <- qtl2mr:::qtl2mr_sym(
    residual_cov$P_perp %*% raw %*% t(residual_cov$P_perp)
  )
  work_vr <- qtl2mr:::qtl2mr_sym(
    residual_cov$P_perp %*% work %*% t(residual_cov$P_perp)
  )
  expect_equal(residual_cov$Omega_raw, raw, tolerance = 1e-14)
  expect_equal(residual_cov$Omega_work, work, tolerance = 1e-14)
  expect_equal(residual_cov$V_r, expected_vr, tolerance = 1e-14)
  expect_gt(max(abs(residual_cov$V_r - work_vr)), 1e-4)

})

test_that("utility guardrails preserve dimnames and warn on dropped covariates", {
  covar <- cbind(a = 1:6, b = 2 * (1:6), c = rep(1, 6))
  rownames(covar) <- paste0("id", seq_len(6))
  expect_warning(
    kept <- qtl2mr:::qtl2mr_drop_depcols(covar, add_intercept = TRUE),
    "Dropping linearly dependent covariate"
  )
  expect_equal(colnames(kept), "a")

  A <- diag(2)
  dimnames(A) <- list(c("a1", "a2"), c("a1", "a2"))
  B <- diag(3)
  dimnames(B) <- list(c("b1", "b2", "b3"), c("b1", "b2", "b3"))
  out <- qtl2mr:::block_diag_matrix(list(A, B))
  expect_equal(rownames(out), c("a1", "a2", "b1", "b2", "b3"))
  expect_equal(colnames(out), c("a1", "a2", "b1", "b2", "b3"))
})

test_that("founder residual and ratio diagnostics use stored founder H", {
  gp <- make_sim_genoprobs(n = 90, seed = 12)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  x <- as.numeric(0.7 * Z[, "FC2"] + rnorm(length(ids), sd = 0.4))
  y <- as.numeric(0.55 * x + 0.15 * Z[, "FC5"] + rnorm(length(ids), sd = 0.4))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = A)
  fit <- second_order_gls_ivw(block$beta_x, block$beta_y, block$Sigma_x, block$Sigma_y,
                              allow_cross_cov = TRUE, Sigma_xy = block$Sigma_xy)

  scores <- founder_residual_scores(block, theta = fit$theta)
  expect_equal(nrow(scores), 8L)
  expect_lt(abs(sum(scores$e)), 1e-10)
  expect_true(all(is.finite(scores$t_one_founder_vs_rest)))

  f_ratios <- founder_coordinate_ratios(block, theta_ref = fit$theta)
  expect_equal(nrow(f_ratios), 8L)
  expect_true(all(c("theta_diag", "F_X_diag", "weak_denominator", "cov_xy") %in% names(f_ratios)))
})

test_that("multilocus diagnostic layer returns scalar summaries and Q_between", {
  gp1 <- make_sim_genoprobs(n = 100, marker = "m1", seed = 21)
  set.seed(22)
  raw2 <- matrix(stats::rexp(100 * 8), nrow = 100, ncol = 8)
  probs2 <- raw2 / rowSums(raw2)
  arr <- array(NA_real_, dim = c(100, 8, 2),
               dimnames = list(dimnames(gp1[[1]])[[1]], LETTERS[1:8], c("m1", "m2")))
  arr[, , "m1"] <- gp1[[1]][, , 1]
  arr[, , "m2"] <- probs2
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z1 <- gp[[1]][, , "m1"] %*% A
  Z2 <- gp[[1]][, , "m2"] %*% A
  x <- as.numeric(0.6 * Z1[, "FC1"] + 0.5 * Z2[, "FC2"] + rnorm(length(ids), sd = 0.5))
  y <- as.numeric(0.7 * x + rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  blocks <- list(
    fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m1"),
    fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m2")
  )
  diag <- fit_multilocus_diagnostic_layer(blocks)
  expect_equal(nrow(diag$per_locus), 2L)
  expect_true("Global SO" %in% diag$estimates$estimator)
  expect_true("Two-step IVW" %in% diag$estimates$estimator)
  expect_equal(diag$heterogeneity$df_between, 1L)
  expect_true(is.finite(diag$heterogeneity$Q_between))
})

test_that("weighted median and mode diagnostics get bootstrap-SE Wald CIs without changing RNG", {
  per_locus <- data.frame(
    theta = c(0.42, 0.61, 0.74, 0.88, 1.03),
    se = c(0.11, 0.18, 0.14, 0.16, 0.21),
    stringsAsFactors = FALSE
  )
  point_only <- qtl2mr:::summarize_per_locus_so(per_locus, bootstrap_B = 0L)
  expect_true(is.na(point_only$estimates$ci_lower[point_only$estimates$estimator == "Weighted median"]))
  expect_true(is.na(point_only$estimates$ci_lower[point_only$estimates$estimator == "Weighted mode"]))

  set.seed(704)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  boot <- qtl2mr:::summarize_per_locus_so(per_locus, bootstrap_B = 80L, bootstrap_seed = 17L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed_before)

  median_row <- boot$estimates[boot$estimates$estimator == "Weighted median", , drop = FALSE]
  mode_row <- boot$estimates[boot$estimates$estimator == "Weighted mode", , drop = FALSE]
  expect_equal(median_row$ci_type[[1L]], "bootstrap-SE Wald")
  expect_equal(mode_row$ci_type[[1L]], "bootstrap-SE Wald")
  expect_true(all(is.finite(c(median_row$se, median_row$ci_lower, median_row$ci_upper))))
  expect_true(all(is.finite(c(mode_row$se, mode_row$ci_lower, mode_row$ci_upper))))
  expect_equal(
    boot$estimates$theta[boot$estimates$estimator == "Weighted median"],
    point_only$estimates$theta[point_only$estimates$estimator == "Weighted median"]
  )
  expect_equal(
    boot$estimates$theta[boot$estimates$estimator == "Weighted mode"],
    point_only$estimates$theta[point_only$estimates$estimator == "Weighted mode"]
  )
})

test_that("diagnostic layer tolerates a single locus by omitting weighted mode", {
  gp <- make_sim_genoprobs(n = 85, marker = "m1", seed = 23)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  x <- as.numeric(0.7 * Z[, "FC1"] + rnorm(length(ids), sd = 0.5))
  y <- as.numeric(0.5 * x + rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m1")
  diag <- fit_multilocus_diagnostic_layer(list(m1 = block))
  expect_equal(nrow(diag$per_locus), 1L)
  expect_true(is.na(diag$estimates$theta[diag$estimates$estimator == "Weighted mode"]))
  expect_equal(diag$heterogeneity$df_between, 0L)
})

test_that("base plotting helpers render to an explicit device", {
  gp <- make_sim_genoprobs(n = 70, seed = 41)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  x <- as.numeric(0.8 * Z[, "FC1"] + rnorm(length(ids), sd = 0.5))
  y <- as.numeric(0.5 * x + 0.2 * Z[, "FC4"] + rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  block <- fit_mr_locus_block(gp, x, y, contrasts = A)
  per <- fit_per_locus_so(list(block))
  fit <- fit_multilocus_founder_mr(list(block))
  res <- founder_residual_scores(block, theta = fit$theta)
  ratios <- founder_coordinate_ratios(block, theta_ref = fit$theta)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  expect_silent(plot_per_locus_so_forest(per, global_theta = fit$theta))
  expect_silent(plot_founder_residual_heatmap(res, row = "marker"))
  expect_silent(plot_founder_ratio_forest(ratios))
  grDevices::dev.off()
  expect_true(file.exists(path))
})
