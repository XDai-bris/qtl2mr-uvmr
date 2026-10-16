# All fixtures below are generated from random numbers. They contain no study
# data. Seed 1 reproduces the former unit-scaling failure; seed 2024, draw 23
# reproduces a missed minimum and a three-component confidence set.
v040_scale_block <- function() {
  qtl2mr:::qtl2mr_with_seed(1, {
    k <- 7L
    Ax <- matrix(rnorm(k * k, sd = .15), k)
    Ay <- matrix(rnorm(k * k, sd = .15), k)
    bx <- rnorm(k, sd = 1.5)
    Sx <- crossprod(Ax) + diag(.05, k)
    Sy <- crossprod(Ay) + diag(.06, k)
    list(beta_x = bx, beta_y = .5 * bx + rnorm(k, sd = .2),
         Sigma_x = Sx, Sigma_y = Sy,
         Sigma_xy = .3 * qtl2mr:::qtl2mr_matrix_power(Sx, .5) %*%
           qtl2mr:::qtl2mr_matrix_power(Sy, .5))
  })
}

v040_rescale_block <- function(b, x, y) {
  b$beta_x <- b$beta_x * x
  b$beta_y <- b$beta_y * y
  b$Sigma_x <- b$Sigma_x * x^2
  b$Sigma_y <- b$Sigma_y * y^2
  b$Sigma_xy <- b$Sigma_xy * x * y
  b
}

v040_multiminimum_blocks <- function() {
  qtl2mr:::qtl2mr_with_seed(2024, {
    for (draw in seq_len(23L)) {
      L <- sample(1:3, 1)
      strength <- sample(c(.1, .3, .6, 1.2), 1)
      blocks <- lapply(seq_len(L), function(i) {
        # Keep draw order explicit for reproducibility of this synthetic case.
        k <- 7L
        Ax <- matrix(rnorm(k * k, sd = .3), k)
        Ay <- matrix(rnorm(k * k, sd = .3), k)
        Sx <- crossprod(Ax) + diag(.05, k)
        Sy <- crossprod(Ay) + diag(.05, k)
        rho <- runif(1, -.6, .6)
        Sxy <- rho * qtl2mr:::qtl2mr_matrix_power(Sx, .5) %*%
          qtl2mr:::qtl2mr_matrix_power(Sy, .5)
        ax <- rnorm(k) * strength
        joint <- rbind(cbind(Sx, Sxy), cbind(t(Sxy), Sy))
        e <- qtl2mr:::qtl2mr_rmvnorm1(rep(0, 2 * k), joint)
        list(beta_x = ax + e[seq_len(k)], beta_y = .5 * ax + e[k + seq_len(k)],
             Sigma_x = Sx, Sigma_y = Sy, Sigma_xy = Sxy)
      })
    }
    blocks
  })
}

test_that("second-order estimates and confidence sets respect trait units", {
  b <- v040_scale_block()
  base <- qtl2mr:::qtl2mr_fit_multilocus_once(list(b))
  expect_true(base$primary_reliable)
  expect_identical(base$inference_status, "not_calibrated")
  for (units in list(c(1e4, 1), c(1e-8, 1e5), c(1e8, 1e-5), c(-100, 10))) {
    changed <- v040_rescale_block(b, units[[1L]], units[[2L]])
    ratio <- units[[2L]] / units[[1L]]
    fit <- qtl2mr:::qtl2mr_fit_multilocus_once(list(changed))
    expect_equal(fit$theta / ratio, base$theta, tolerance = 1e-9)
    expect_equal(fit$se / abs(ratio), base$se, tolerance = 1e-9)
    expect_equal(unname(sort(unlist(fit$confidence_set) / ratio)), unname(unlist(base$confidence_set)), tolerance = 1e-8)
    expect_equal(fit$Q, base$Q, tolerance = 1e-10)
    expect_identical(fit$primary_reliable, base$primary_reliable)
    single <- second_order_gls_ivw(changed$beta_x, changed$beta_y,
      changed$Sigma_x, changed$Sigma_y, allow_cross_cov = TRUE, Sigma_xy = changed$Sigma_xy)
    expect_equal(single$theta, fit$theta)
    expect_equal(single$confidence_set, fit$confidence_set)
  }
})

test_that("whole-line search exposes missed minima and disconnected unbounded sets", {
  blocks <- v040_multiminimum_blocks()
  fit <- qtl2mr:::qtl2mr_fit_multilocus_once(blocks)
  expect_lt(fit$Q, fit$objective(-.05824))
  expect_equal(fit$theta, -61.22557, tolerance = 2e-5)
  expect_equal(fit$Q, 20.19209, tolerance = 1e-5)
  expect_identical(fit$confidence_set_type, "disconnected")
  expect_equal(nrow(fit$confidence_set), 3L)
  expect_identical(fit$confidence_set$lower[[1L]], -Inf)
  expect_identical(tail(fit$confidence_set$upper, 1L), Inf)
  expect_true(all(fit$objective_tail < fit$confidence_cutoff))
  expect_false(fit$primary_reliable)
  expect_false(fit$profile_ci_finite)
  expect_true(is.na(fit$profile_ci_lower) && is.na(fit$profile_ci_upper))
  expect_true(fit$numerical_ok)
  expect_false(fit$search_certified)
  # Independent, denser evaluation checks objective minimum and set membership.
  theta <- tan(seq(-pi / 2 + 1e-5, pi / 2 - 1e-5, length.out = 2001L))
  q <- vapply(theta, function(th) sum(vapply(blocks, function(b) {
    S <- b$Sigma_y + th^2 * b$Sigma_x - th * (b$Sigma_xy + t(b$Sigma_xy))
    S <- S + diag(1e-10 * mean(diag(S)), nrow(S))
    r <- b$beta_y - th * b$beta_x
    sum(backsolve(chol(S), r, transpose = TRUE)^2)
  }, numeric(1))), numeric(1))
  expect_lte(fit$Q, min(q) + 1e-8)
  in_set <- vapply(theta, function(th) any(th >= fit$confidence_set$lower & th <= fit$confidence_set$upper), logical(1))
  expect_identical(in_set, q <= fit$confidence_cutoff)
  hinted <- qtl2mr:::qtl2mr_fit_multilocus_once(blocks, theta_grid = c(-1, 0, 1))
  expect_equal(hinted$Q, fit$Q, tolerance = 1e-10)
  expect_equal(hinted$confidence_set, fit$confidence_set, tolerance = 1e-8)
  expect_identical(hinted$theta_grid_role, "supplemental_search_points")
})

test_that("profile-set helper retains bounded, one-sided and all-real topology", {
  grid <- seq(-pi / 2, pi / 2, length.out = 257L)
  profile <- function(fun, qmin = 0) qtl2mr:::qtl2mr_profile_set(
    fun, grid, vapply(grid, fun, numeric(1)), qmin, scale = 1)
  all_real <- profile(function(a) 2, qmin = 2)
  expect_identical(all_real$type, "all_real")
  expect_equal(unname(unlist(all_real$components)), c(-Inf, Inf))
  one_sided <- profile(function(a) 10 * (1 + sin(a)))
  expect_identical(one_sided$type, "unbounded")
  expect_identical(one_sided$components$lower[[1L]], -Inf)
  separated <- profile(function(a) 100 * (sin(a)^2 - .5)^2)
  expect_identical(separated$type, "disconnected")
  expect_equal(nrow(separated$components), 2L)
  expect_true(all(is.finite(unlist(separated$components))))
})

test_that("uninformative exposure does not produce a reliable finite effect", {
  fit <- second_order_gls_ivw(c(0, 0), c(1, -1), diag(2), diag(2))
  expect_identical(fit$confidence_set_type, "all_real")
  expect_false(fit$primary_reliable)
  expect_true(is.na(fit$theta))
  expect_identical(fit$optimizer_status, "minimum_at_infinity")
  wrapped <- fit_multilocus_founder_mr(list(list(beta_x = c(0, 0),
    beta_y = c(1, -1), Sigma_x = diag(2), Sigma_y = diag(2), Sigma_xy = matrix(0, 2, 2))),
    heterogeneity_method = "none")
  expect_false(wrapped$primary_reliable)
  expect_identical(wrapped$confidence_set_type, "all_real")
  expect_true(is.na(wrapped$locus_diagnostics$Q_locus))
  expect_error(second_order_gls_ivw(c(1, 2), c(1, 2), diag(2), diag(2),
    theta_grid = c(0, NA)), "finite numeric")
})

test_that("empty covariates can be normalized more than once", {
  ids <- paste0("s", seq_len(20L))
  empty <- qtl2mr:::qtl2mr_as_addcovar(NULL, ids)
  expect_identical(dim(empty), c(20L, 0L))
  expect_identical(qtl2mr:::qtl2mr_as_addcovar(empty, ids), empty)
})

test_that("unsupported unregularized singular covariance fails closed", {
  expect_error(second_order_gls_ivw(c(1, 0), c(0, 1), diag(c(1, 0)), diag(2),
    ridge = 0), "positive definite.*singular covariance")
  # With positive trace ridge, singular exposure covariance has a defined limit.
  fit <- second_order_gls_ivw(c(1, 0), c(0, 1), diag(c(1, 0)), diag(2), ridge = 1e-10)
  expect_false(fit$primary_reliable)
  expect_true(any(!is.finite(unlist(fit$confidence_set))))
})

test_that("covariate rank checks ignore nonfinite rows without hiding sample loss", {
  covar <- matrix(c(Inf, seq_len(19L)), ncol = 1L,
    dimnames = list(paste0("s", seq_len(20L)), "z"))
  checked <- qtl2mr:::qtl2mr_as_addcovar(covar)
  expect_identical(dim(checked), c(20L, 1L))
  expect_identical(checked[[1L]], Inf)
  # The fitter, not normalization, removes incomplete rows and records the loss.
  covar[,] <- Inf
  expect_error(qtl2mr:::qtl2mr_as_addcovar(covar), "No finite rows")
  expect_error(qtl2mr:::qtl2mr_as_addcovar(matrix("1", 2L, 1L,
    dimnames = list(c("a", "b"), "z"))), "numeric/logical")
})

test_that("nonfinite interior optimizer probes cannot be labelled reliable", {
  gp <- make_sim_genoprobs(n = 60L, marker = "m1", seed = 861L)
  ids <- dimnames(gp[[1L]])[[1L]]
  Z <- gp[[1L]][, , "m1"] %*% make_founder_contrasts(LETTERS[1:8])
  pheno <- data.frame(x = .6 * Z[, "FC1"] + rnorm(60L, sd = .5),
    y = .5 * Z[, "FC1"] + rnorm(60L, sd = .5),
    z = .4 * Z[, "FC1"] + rnorm(60L, sd = .5), row.names = ids)
  covar <- data.frame(c1 = rnorm(60L), row.names = ids)
  # Ten samples, eight founder parameters and one covariate leave one residual
  # degree of freedom. This creates a numerically singular overlap covariance.
  b <- fit_mr_locus_block_overlap(gp, setNames(pheno$x, ids), setNames(pheno$y, ids),
    addcovar_x = covar, addcovar_y = covar, marker = "m1",
    sample_ids_x = ids[51:60], sample_ids_y = ids[51:60], standardize_phenotypes = TRUE)
  expect_silent(fit <- qtl2mr:::qtl2mr_fit_multilocus_once(list(b)))
  expect_true(fit$search_nonfinite_interior)
  expect_false(fit$numerical_ok)
  expect_false(fit$primary_reliable)
  expect_identical(fit$confidence_set_type, "unresolved")
  expect_true(is.na(fit$profile_ci_lower) && is.na(fit$profile_ci_upper))
  expect_true(any(grepl("non-finite interior", fit$optimizer_warnings, fixed = TRUE)))
})
