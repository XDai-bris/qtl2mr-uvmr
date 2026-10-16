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

  H_resolved <- qtl2mr:::qtl2mr_resolve_founder_H(block)
  V_e <- qtl2mr:::qtl2mr_sym(H_resolved %*% expected_vr %*% t(H_resolved))
  V_one <- qtl2mr:::qtl2mr_one_vs_rest_matrix(rownames(H_resolved))
  denom <- sqrt(pmax(diag(t(V_one) %*% V_e %*% V_one), .Machine$double.eps))
  nsim <- 200L
  seed <- 932L
  expected_p95 <- qtl2mr:::qtl2mr_with_seed(seed, {
    L <- qtl2mr:::qtl2mr_matrix_sqrt(expected_vr, ridge = 0)
    z <- matrix(stats::rnorm(nsim * nrow(expected_vr)), nrow = nrow(expected_vr), ncol = nsim)
    simulated <- t(V_one) %*% H_resolved %*% L %*% z
    simulated <- sweep(simulated, 1L, denom, "/")
    apply(abs(simulated), 1L, stats::quantile, probs = 0.95, names = FALSE)
  })
  null <- founder_residual_null_simulation(block,
    theta = theta, nsim = nsim, seed = seed, ridge = ridge)
  expect_equal(null$founder$null_p95_abs_t, unname(expected_p95), tolerance = 1e-12)
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
  c_ratios <- founder_contrast_ratios(block, theta_ref = fit$theta)
  expect_equal(nrow(f_ratios), 8L)
  expect_equal(nrow(c_ratios), 7L)
  expect_true(all(c("theta_diag", "F_X_diag", "weak_denominator", "cov_xy") %in% names(f_ratios)))

  alpha <- second_order_nonorthogonality(block, theta = fit$theta)
  expect_true(is.finite(alpha$alpha_nonorthogonality))
  expect_true(is.finite(alpha$overlap_component))

  null <- founder_residual_null_simulation(block, theta = fit$theta, nsim = 80, seed = 1)
  expect_equal(nrow(null$founder), 8L)
  expect_equal(nrow(null$summary), 1L)
  expect_true(all(c("observed_abs_t", "null_p95_abs_t", "observed_percentile") %in% names(null$founder)))
  expect_true(all(null$founder$observed_percentile >= 0 & null$founder$observed_percentile <= 1))

  set.seed(701)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  invisible(founder_residual_null_simulation(block, theta = fit$theta, nsim = 20, seed = 99))
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed_before)
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

test_that("founder score and heteroskedasticity helpers work on named vectors", {
  gp <- make_sim_genoprobs(n = 80, seed = 31)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  beta <- stats::setNames(c(0.5, -0.2, rep(0, 5)), paste0("FC", 1:7))
  score <- founder_local_score(gp, beta_x = beta, contrasts = A)
  manual <- as.numeric((gp[[1]][, , 1] %*% A[, names(beta)]) %*% beta)
  expect_equal(as.numeric(score), manual, tolerance = 1e-12)

  y <- as.numeric(0.4 * score + rnorm(length(score), sd = 0.5 + 0.1 * abs(score)))
  names(y) <- ids
  covar <- matrix(rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  resid <- score_adjusted_residuals(y, score, covariates = covar)
  expect_equal(nrow(resid), length(ids))
  ht <- score_heteroskedasticity_test(resid)
  expect_true(is.finite(ht$joint_p))
  bins <- score_bin_variance(resid, B = 20, seed = 1)
  expect_equal(nrow(bins), 5L)

  set.seed(702)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  invisible(score_bin_variance(resid, B = 20, seed = 12))
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed_before)
})

test_that("score adjusted residual model always includes an intercept", {
  ids <- paste0("id", 1:30)
  score <- stats::setNames(seq(-1, 1, length.out = length(ids)), ids)
  y <- stats::setNames(2 + 0.5 * score + 0.2 * score^2 + rnorm(length(ids), sd = 0.01), ids)
  resid <- score_adjusted_residuals(y, score, covariates = NULL)
  expect_lt(abs(mean(resid$residual)), 1e-10)
  lm0 <- stats::lm(y ~ score + I(score^2))
  expect_equal(as.numeric(resid$residual), as.numeric(stats::resid(lm0)), tolerance = 1e-10)
  expect_match(unique(resid$mean_model_used), "intercept")
})

test_that("multi-locus xPRS diagnostics match summed founder scores", {
  gp1 <- make_sim_genoprobs(n = 85, marker = "m1", seed = 33)
  set.seed(34)
  raw2 <- matrix(stats::rexp(85 * 8), nrow = 85, ncol = 8)
  probs2 <- raw2 / rowSums(raw2)
  arr <- array(NA_real_, dim = c(85, 8, 2),
               dimnames = list(dimnames(gp1[[1]])[[1]], LETTERS[1:8], c("m1", "m2")))
  arr[, , "m1"] <- gp1[[1]][, , 1]
  arr[, , "m2"] <- probs2
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  beta1 <- stats::setNames(c(0.5, -0.2, rep(0, 5)), paste0("FC", 1:7))
  beta2 <- stats::setNames(c(0.1, 0.4, -0.15, rep(0, 4)), paste0("FC", 1:7))
  blocks <- list(
    m1 = list(marker = "m1", beta_x = beta1, contrasts = A),
    m2 = list(marker = "m2", beta_x = beta2, contrasts = A)
  )
  xprs <- multilocus_xprs_scores(gp, blocks, markers = c("m1", "m2"), sample_ids = ids)
  s1 <- founder_local_score(gp, marker = "m1", beta_x = beta1, contrasts = A)
  s2 <- founder_local_score(gp, marker = "m2", beta_x = beta2, contrasts = A)
  expect_equal(xprs$sample_scores$xprs, as.numeric(s1[ids] + s2[ids]), tolerance = 1e-12)
  expect_equal(as.numeric(xprs$locus_score_matrix[, "m1"]), as.numeric(s1[ids]), tolerance = 1e-12)

  loo <- leave_one_locus_xprs_scores(xprs)
  expect_equal(as.numeric(loo[, "m1"]), as.numeric(s2[ids]), tolerance = 1e-12)
  expect_equal(as.numeric(loo[, "m2"]), as.numeric(s1[ids]), tolerance = 1e-12)

  set.seed(35)
  y <- as.numeric(0.6 * xprs$sample_scores$xprs + rnorm(length(ids), sd = 0.5))
  names(y) <- ids
  covar <- matrix(rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  one <- multilocus_xprs_scores(gp, blocks["m1"], markers = "m1", sample_ids = ids)
  global_one <- multilocus_xprs_heteroskedasticity(
    traits = list(exposure_trait = y),
    xprs = one$sample_scores$xprs,
    covariates = covar,
    sample_ids = ids,
    trait_roles = "exposure",
    B = 20,
    seed = 1,
    n_loci = 1L
  )
  s1z <- qtl2mr:::qtl2mr_standardized_named_score(s1[ids])
  direct_resid <- score_adjusted_residuals(y, s1z, covariates = covar, sample_ids = ids)
  direct_ht <- score_heteroskedasticity_test(direct_resid)
  expect_equal(global_one$tests$joint_p, direct_ht$joint_p, tolerance = 1e-12)

  y2 <- as.numeric(0.3 * xprs$sample_scores$xprs + rnorm(length(ids), sd = 0.6))
  names(y2) <- ids
  loo_tests <- leave_one_locus_xprs_heteroskedasticity(
    traits = list(exposure_trait = y, outcome_trait = y2),
    xprs = xprs,
    covariates = covar,
    sample_ids = ids,
    trait_roles = c("exposure", "outcome"),
    marker_info = data.frame(marker = c("m1", "m2"), chr = c("1", "1"), locus_rank = 1:2)
  )
  expect_equal(nrow(loo_tests), 4L)
  expect_true(all(is.finite(loo_tests$joint_p)))

  diag <- fit_xprs_heteroskedasticity_diagnostic(
    traits = list(exposure_trait = y, outcome_trait = y2),
    xprs = xprs,
    covariates = covar,
    sample_ids = ids,
    trait_roles = c("exposure", "outcome"),
    per_locus_scores = xprs$locus_score_matrix,
    marker_info = data.frame(marker = c("m1", "m2"), chr = c("1", "1"), locus_rank = 1:2),
    include_loo = TRUE,
    n_bins = 5,
    B = 20,
    seed = 2,
    n_loci = 2L
  )
  expect_s3_class(diag, "qtl2mr_xprs_heteroskedasticity_diagnostic")
  expect_equal(nrow(diag$global_tests), 2L)
  expect_equal(nrow(diag$loo_tests), 4L)
  expect_equal(nrow(diag$localization_tests), 4L)
  expect_equal(nrow(diag$localization_bins), 20L)
})

test_that("score heteroskedasticity table augmentation handles multi-locus traits", {
  tests <- data.frame(
    marker = rep(paste0("m", 1:4), each = 2),
    chr = rep(1:4, each = 2),
    locus_rank = rep(1:4, each = 2),
    trait_role = rep(c("exposure", "outcome"), 4),
    trait = rep(c("exposure_trait", "outcome_trait"), 4),
    joint_p = c(0.2, 0.001, 0.04, 0.5, 0.006, 0.03, 0.8, 0.9),
    stringsAsFactors = FALSE
  )
  bins <- do.call(rbind, lapply(seq_len(nrow(tests)), function(i) {
    v <- if (tests$marker[[i]] == "m1" && tests$trait[[i]] == "outcome_trait") c(1, 1, 1, 1, 4) else seq(1, 1.4, length.out = 5)
    data.frame(
      marker = tests$marker[[i]],
      chr = tests$chr[[i]],
      locus_rank = tests$locus_rank[[i]],
      trait_role = tests$trait_role[[i]],
      trait = tests$trait[[i]],
      B = 10,
      score_bin = 1:5,
      score_bin_label = paste0("Q", 1:5),
      mean_score = seq(-1, 1, length.out = 5),
      n = 20,
      residual_variance = v,
      ci_lower = v * 0.8,
      ci_upper = v * 1.2,
      stringsAsFactors = FALSE
    )
  }))
  aug <- augment_score_heteroskedasticity_tests(tests, bins)
  expect_equal(nrow(aug), nrow(tests))
  expect_true(all(c("joint_q_bh", "variance_ratio_q5_q1", "variance_ratio_max_min",
                    "bin_variance_slope", "shape_sensitive_qc_warning") %in% names(aug)))
  expect_true(aug$flag_nominal_0.01[aug$marker == "m1" & aug$trait == "outcome_trait"])
  expect_true(aug$shape_sensitive_qc_warning[aug$marker == "m1" & aug$trait == "outcome_trait"])
  expect_true(all(is.finite(aug$joint_q_bh)))
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
  null <- founder_residual_null_simulation(block, theta = fit$theta, nsim = 40, seed = 2)
  bins <- data.frame(score_bin = 1:5, score_bin_label = paste0("Q", 1:5),
                     residual_variance = seq(1, 2, length.out = 5),
                     ci_lower = seq(0.8, 1.8, length.out = 5),
                     ci_upper = seq(1.2, 2.2, length.out = 5))
  tests_all <- data.frame(
    marker = rep(paste0("m", 1:6), each = 2),
    chr = rep(1:6, each = 2),
    locus_rank = rep(1:6, each = 2),
    trait_role = rep(c("exposure", "outcome"), 6),
    trait = rep(c("exposure_trait", "outcome_trait"), 6),
    joint_p = c(0.4, 0.002, 0.03, 0.5, 0.8, 0.04, 0.01, 0.7, 0.2, 0.6, 0.9, 0.8),
    stringsAsFactors = FALSE
  )
  bins_all <- do.call(rbind, lapply(seq_len(nrow(tests_all)), function(i) {
    vv <- seq(0.8, 1.2, length.out = 5) * (1 + 0.05 * tests_all$locus_rank[[i]])
    if (tests_all$joint_p[[i]] < 0.01) vv <- c(0.8, 0.9, 1.0, 1.1, 2.5)
    data.frame(
      marker = tests_all$marker[[i]],
      chr = tests_all$chr[[i]],
      locus_rank = tests_all$locus_rank[[i]],
      trait_role = tests_all$trait_role[[i]],
      trait = tests_all$trait[[i]],
      score_bin = 1:5,
      score_bin_label = paste0("Q", 1:5),
      residual_variance = vv,
      ci_lower = vv * 0.8,
      ci_upper = vv * 1.2,
      stringsAsFactors = FALSE
    )
  }))
  tests_all <- augment_score_heteroskedasticity_tests(tests_all, bins_all)

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  expect_silent(plot_per_locus_so_forest(per, global_theta = fit$theta))
  expect_silent(plot_founder_residual_heatmap(res, row = "marker"))
  expect_silent(plot_founder_ratio_forest(ratios))
  ratios2 <- ratios
  ratios2$marker <- "m2"
  ratio_heat <- rbind(ratios, ratios2)
  ratio_heat$locus_rank <- rep(1:2, each = nrow(ratios))
  expect_silent(plot_founder_coordinate_wald_ratio_heatmap(ratio_heat, marker_order = c("m1", "m2")))
  expect_silent(plot_founder_residual_null(null))
  expect_silent(plot_score_bin_variance(bins))
  expect_silent(plot_score_heteroskedasticity_pvalue_scan(tests_all))
  expect_silent(plot_score_bin_variance_heatmap(bins_all))
  expect_silent(plot_score_heteroskedasticity_dashboard(tests_all, bins_all))
  expect_silent(plot_score_bin_variance_grid(tests_all, bins_all, top_n = 4))
  expect_equal(nrow(plot_score_bin_variance_grid(tests_all, bins_all, top_n = 1)), 1L)
  expect_equal(nrow(plot_score_bin_variance_grid(tests_all, bins_all, top_n = NULL)),
               sum(tests_all$flag_nominal_0.01 | tests_all$flag_fdr_0.05))
  global_tests <- data.frame(
    marker = "multi_locus_xprs",
    chr = "multi",
    locus_rank = NA_integer_,
    trait_role = c("exposure", "outcome"),
    trait = c("exposure_trait", "outcome_trait"),
    n_loci = 6L,
    sample_size = 70L,
    score_sd = 1,
    residual_sd = c(0.9, 1.0),
    joint_test_statistic = c(5.2, 2.1),
    joint_p = c(0.004, 0.12),
    variance_ratio_q5_q1 = c(1.5, 1.2),
    variance_ratio_max_min = c(1.8, 1.4),
    bin_variance_slope = c(0.1, 0.03),
    stringsAsFactors = FALSE
  )
  global_bins <- do.call(rbind, lapply(seq_len(nrow(global_tests)), function(i) {
    vv <- seq(0.8, 1.3, length.out = 5) * (1 + 0.1 * i)
    data.frame(
      marker = "multi_locus_xprs",
      chr = "multi",
      trait_role = global_tests$trait_role[[i]],
      trait = global_tests$trait[[i]],
      score_bin = 1:5,
      score_bin_label = paste0("Q", 1:5),
      residual_variance = vv,
      ci_lower = vv * 0.85,
      ci_upper = vv * 1.15,
      stringsAsFactors = FALSE
    )
  }))
  loo_tests <- do.call(rbind, lapply(seq_len(nrow(tests_all)), function(i) {
    data.frame(
      omitted_marker = tests_all$marker[[i]],
      chr = tests_all$chr[[i]],
      locus_rank = tests_all$locus_rank[[i]],
      trait_role = tests_all$trait_role[[i]],
      trait = tests_all$trait[[i]],
      joint_p = pmin(0.9, tests_all$joint_p[[i]] + 0.03),
      stringsAsFactors = FALSE
    )
  }))
  diag_plot <- structure(
    list(
      global_tests = global_tests,
      global_bins = global_bins,
      loo_tests = loo_tests,
      localization_tests = tests_all,
      localization_bins = bins_all
    ),
    class = c("qtl2mr_xprs_heteroskedasticity_diagnostic", "list")
  )
  expect_silent(plot_multilocus_xprs_heteroskedasticity_global(global_tests))
  expect_silent(plot_multilocus_xprs_heteroskedasticity_global(global_tests, global_bins))
  expect_silent(plot_multilocus_xprs_heteroskedasticity_global(global_tests, global_bins, traits = "outcome_trait"))
  expect_silent(plot_multilocus_xprs_heteroskedasticity_global(diag_plot))
  expect_silent(plot_multilocus_xprs_bin_variance(global_bins))
  expect_silent(plot_multilocus_xprs_bin_variance(diag_plot))
  expect_silent(plot_multilocus_xprs_leave_one_out(loo_tests, global_tests))
  expect_silent(plot_multilocus_xprs_leave_one_out(diag_plot))
  expect_silent(plot_multilocus_xprs_heteroskedasticity_dashboard(global_tests, global_bins, loo_tests, tests_all))
  expect_silent(plot_multilocus_xprs_heteroskedasticity_dashboard(diag_plot))
  expect_output(print(diag_plot), "qtl2mr multi-locus xPRS heteroskedasticity diagnostic")
  expect_silent(plot(diag_plot))
  expect_silent(plot(diag_plot, type = "loo"))
  expect_silent(plot(diag_plot, type = "localization", top_n = 2))
  grDevices::dev.off()
  expect_true(file.exists(path))
})
