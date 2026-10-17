test_that("second-order optimizer reports bracket and convergence metadata", {
  beta_x <- c(0.35, -0.2, 0.15)
  beta_y <- 0.8 * beta_x + c(0.01, -0.005, 0.003)
  Sigma_x <- diag(c(0.015, 0.012, 0.018))
  Sigma_y <- diag(c(0.02, 0.025, 0.022))
  Sigma_xy <- matrix(0, 3, 3)

  fit <- second_order_gls_ivw(
    beta_x = beta_x,
    beta_y = beta_y,
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    allow_cross_cov = TRUE,
    Sigma_xy = Sigma_xy
  )

  expect_true(all(c(
    "theta_start", "bracket_lower", "bracket_upper", "bracket_status",
    "bracket_n_expand", "optimizer_status", "objective_min", "optimizer_warnings"
  ) %in% names(fit)))
  expect_equal(fit$optimizer_status, "success")
  expect_identical(fit$bracket_status, "whole_real_line")
  expect_lte(fit$bracket_lower, fit$theta)
  expect_gte(fit$bracket_upper, fit$theta)
  expect_equal(fit$objective_min, fit$Q)
})

test_that("multi-locus wrapper preserves optimizer metadata", {
  block <- list(
    beta_x = stats::setNames(c(0.35, -0.2, 0.15), paste0("FC", 1:3)),
    beta_y = stats::setNames(c(0.29, -0.16, 0.12), paste0("FC", 1:3)),
    Sigma_x = diag(c(0.015, 0.012, 0.018)),
    Sigma_y = diag(c(0.02, 0.025, 0.022)),
    Sigma_xy = matrix(0, 3, 3),
    marker = "m1",
    chr = "1",
    n_overlap = 50L
  )

  fit <- fit_multilocus_founder_mr(list(m1 = block))

  expect_equal(fit$optimizer_status, "success")
  expect_identical(fit$bracket_status, "whole_real_line")
  expect_lte(fit$bracket_lower, fit$theta)
  expect_gte(fit$bracket_upper, fit$theta)
  expect_equal(fit$objective_min, fit$Q)
  expect_equal(fit$independence$optimizer_status, "success")
})

test_that("weak second-order denominator is flagged as unreliable", {
  fit <- second_order_gls_ivw(
    beta_x = c(0.01, 0.015, 0.02),
    beta_y = c(1, -0.5, 0.8),
    Sigma_x = diag(rep(1, 3)),
    Sigma_y = diag(rep(0.01, 3)),
    allow_cross_cov = TRUE,
    Sigma_xy = matrix(0, 3, 3)
  )

  expect_equal(fit$optimizer_status, "disconnected_profile")
  expect_equal(fit$bracket_status, "whole_real_line")
  expect_true(any(!is.finite(unlist(fit$confidence_set))))
  expect_false(isTRUE(fit$profile_ci_finite))
  expect_false(isTRUE(fit$primary_reliable))
})

test_that("diagnostic_only reports reliable SO and suppresses unreliable SO", {
  diagnostic_layer <- list(estimates = data.frame(
    estimator = "Two-step IVW",
    theta = 0.42,
    se = 0.08,
    ci_lower = 0.2632,
    ci_upper = 0.5768,
    ci_type = "Wald",
    stringsAsFactors = FALSE
  ))
  reliable_fit <- list(
    theta = 0.5,
    se = 0.1,
    profile_ci_lower = 0.3,
    profile_ci_upper = 0.7,
    primary_reliable = TRUE
  )
  reliable <- qtl2mr:::qtl2mr_reported_estimate(reliable_fit, diagnostic_layer)
  expect_equal(reliable$policy, "diagnostic_only")
  expect_equal(reliable$estimator, "Star-space qtl2MR")
  expect_equal(reliable$source, "primary")
  expect_equal(reliable$theta, reliable_fit$theta)
  expect_equal(reliable$se, reliable_fit$se)

  unreliable_fit <- reliable_fit
  unreliable_fit$primary_reliable <- FALSE
  suppressed <- qtl2mr:::qtl2mr_reported_estimate(unreliable_fit, diagnostic_layer)
  expect_equal(suppressed$policy, "diagnostic_only")
  expect_equal(suppressed$estimator, "not_reported_unreliable_so")
  expect_equal(suppressed$source, "none")
  expect_true(all(is.na(unlist(suppressed[c("theta", "se", "ci_lower", "ci_upper")]))))
  expect_equal(unreliable_fit$theta, 0.5)
  expect_equal(diagnostic_layer$estimates$theta, 0.42)
})

test_that("auto_fallback remains an explicit diagnostic compatibility opt-in", {
  fit <- list(
    theta = 0.5,
    se = 0.1,
    profile_ci_lower = NA_real_,
    profile_ci_upper = NA_real_,
    primary_reliable = FALSE
  )
  diagnostic_layer <- list(estimates = data.frame(
    estimator = "Two-step IVW",
    theta = 0.42,
    se = 0.08,
    ci_lower = 0.2632,
    ci_upper = 0.5768,
    ci_type = "Wald",
    stringsAsFactors = FALSE
  ))
  reported <- qtl2mr:::qtl2mr_reported_estimate(
    fit,
    diagnostic_layer,
    policy = "auto_fallback"
  )
  expect_equal(reported$policy, "auto_fallback")
  expect_equal(reported$estimator, "Two-step IVW")
  expect_equal(reported$source, "diagnostic_fallback")
  expect_equal(reported$theta, 0.42)
  expect_equal(reported$se, 0.08)
})

test_that("heterogeneity auto policy chooses chisq for one locus and simulation for multiple loci", {
  block <- list(
    beta_x = stats::setNames(c(0.35, -0.2, 0.15), paste0("FC", 1:3)),
    beta_y = stats::setNames(c(0.29, -0.16, 0.12), paste0("FC", 1:3)),
    Sigma_x = diag(c(0.015, 0.012, 0.018)),
    Sigma_y = diag(c(0.02, 0.025, 0.022)),
    Sigma_xy = matrix(0, 3, 3),
    marker = "m1",
    chr = "1",
    n_overlap = 50L
  )
  one <- fit_multilocus_founder_mr(list(m1 = block), q_sim_n = 10L, q_sim_seed = 1L)
  expect_equal(one$Q_p_method, "chisq")
  expect_equal(one$Q_chisq_df, 2L)
  expect_equal(one$Q_p_reported, one$Q_p_chisq)
  expect_false(is.finite(one$Q_p_sim))

  block2 <- block
  block2$marker <- "m2"
  block2$beta_y <- block2$beta_y + c(0.01, -0.01, 0.005)
  two <- fit_multilocus_founder_mr(list(m1 = block, m2 = block2), q_sim_n = 10L, q_sim_seed = 1L)
  expect_equal(two$Q_p_method, "simulation")
  expect_equal(two$Q_sim_n, 10L)
  expect_equal(two$Q_sim_n_requested, 10L)
  expect_equal(two$Q_p_reported, two$Q_p_sim)
  expect_true(is.finite(two$Q_p_chisq))

  set.seed(703)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  invisible(fit_multilocus_founder_mr(list(m1 = block, m2 = block2), q_sim_n = 5L, q_sim_seed = 44L))
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed_before)

  expect_warning(
    no_sim <- fit_multilocus_founder_mr(list(m1 = block, m2 = block2), q_sim_n = 0L),
    "q_sim_n = 0"
  )
  expect_equal(no_sim$Q_p_method, "simulation")
  expect_false(is.finite(no_sim$Q_p_reported))
  expect_equal(no_sim$Q_sim_n, 0L)
  expect_equal(no_sim$Q_sim_n_requested, 0L)
  expect_true(is.finite(no_sim$Q_p_chisq))
})
