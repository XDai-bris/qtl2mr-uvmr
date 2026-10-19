v041_numeric_input <- function(kinship = FALSE) {
  set.seed(41017)
  n <- 48L
  ids <- paste0("synthetic", seq_len(n))
  a <- array(rexp(n * 8L * 2L), c(n, 8L, 2L),
    dimnames = list(ids, LETTERS[1:8], c("m1", "m2")))
  for (m in seq_len(2L)) a[, , m] <- a[, , m] / rowSums(a[, , m])
  z <- a[, , 1L] %*% make_founder_contrasts(LETTERS[1:8])
  x <- 3 * z[, 1L] + rnorm(n, sd = .15)
  pheno <- data.frame(x = x, y = .7 * x + rnorm(n, sd = .15), row.names = ids)
  K <- diag(.5, n)
  dimnames(K) <- list(ids, ids)
  list(allele_probs = list(`1` = a), pheno = pheno,
    covar = data.frame(age = seq_len(n), row.names = ids),
    kinship = if (kinship) K else NULL,
    marker_map = data.frame(marker = c("m1", "m2"), chr = "1", pos = 1:2),
    selected_loci = "m1", exposure = "x", outcome = "y")
}

test_that("weighted mode bandwidth and numerical search respect effect units and sign", {
  mode <- qtl2mr:::qtl2mr_weighted_mode
  x <- c(-.0004, -.0001, .00015, .0003, .0007)
  w <- c(4, 2, 7, 3, 1)
  base <- mode(x, w)
  for (scale in c(1e-200, 1e-9, 1e-4, 1e4, 1e9, 1e200, -1e-9, -1e4)) {
    expect_equal(mode(x * scale, w) / scale, base, tolerance = 1e-10)
    expect_equal(mode(x * scale, w, bandwidth = .0001 * abs(scale)) / scale,
      mode(x, w, bandwidth = .0001), tolerance = 1e-10)
  }
  # An exactly tied, separated mixture selects the peak nearest the first
  # component, not always the leftmost peak. No midpoint/valley is reported.
  for (scale in c(1e-8, 1, 1e8, -1e-8, -1, -1e8)) {
    expect_equal(mode(c(-1, 1) * scale, c(1, 1)) / scale,
      mode(c(-1, 1), c(1, 1)), tolerance = 1e-10)
    expect_equal(mode(c(-1, 1) * scale, c(1, 1), bandwidth = .1 * abs(scale)) / scale,
      -1, tolerance = 1e-7)
    expect_equal(mode(rep(.7 * scale, 4), 1:4) / scale, .7)
    expect_identical(mode(rep(0, 4), 1:4), 0)
  }
  expect_error(mode(rep(0, 4), 1:4, bandwidth = -1), "bandwidth")
})

test_that("mode bootstrap keeps fixed-locus Wald construction after unit correction", {
  fun <- qtl2mr:::qtl2mr_bootstrap_scalar_sensitivities
  theta <- c(-.04, .01, .06, .08)
  se <- c(.01, .02, .025, .03)
  baseline <- fun(theta, se, B = 12L, seed = 41L)
  changed <- fun(theta * 1e-6, se * 1e-6, B = 12L, seed = 41L)
  expect_equal(changed$se / 1e-6, baseline$se, tolerance = 1e-8)
  expect_equal(changed$ci_lower / 1e-6, baseline$ci_lower, tolerance = 1e-8)
  expect_equal(changed$ci_upper / 1e-6, baseline$ci_upper, tolerance = 1e-8)
  expect_true(all(baseline$ci_type == "bootstrap-SE Wald"))
  expect_true(all(baseline$bootstrap_complete_B == 12L))
  point <- qtl2mr:::qtl2mr_weighted_mode(theta, 1 / se^2)
  expect_equal(baseline$ci_lower[2], point - qnorm(.975) * baseline$se[2])
  expect_equal(baseline$ci_upper[2], point + qnorm(.975) * baseline$se[2])
})

test_that("hsq formats resolve deliberately without estimation or silent clamping", {
  raw <- v041_numeric_input(kinship = TRUE)
  resolve <- qtl2mr:::qtl2mr_prepared_hsq
  formats <- list(.35, c(exposure = .35, outcome = .35), c(x = .35, y = .35),
    c(`1` = .35), list(`1` = .35), list(exposure = .35, outcome = .35),
    list(`1` = list(x = .35, y = .35)), list(exposure = c(`1` = .35), outcome = list(`1` = .35)))
  for (h in formats) {
    input <- do.call(as_qtl2mr_prepared_input, c(raw, list(hsq = h)))
    expect_identical(resolve(input, "1", "exposure", "x"), .35)
    expect_identical(resolve(input, "1", "outcome", "y"), .35)
  }
  input <- do.call(as_qtl2mr_prepared_input, c(raw, list(hsq = list(exposure = .2))))
  expect_null(resolve(input, "1", "outcome", "y"))
  input$hsq <- .35
  fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none")
  expect_equal(fit$locus_blocks[[1]]$fit_x$hsq, .35)
  expect_equal(fit$locus_blocks[[1]]$fit_y$hsq, .35)
  for (h in list(-.01, 1.01, NA_real_, Inf, "0.3", c(.2, .3),
                 list(unknown = .3), list(exposure = list(x = .3)),
                 c(exposure = .2, exposure = .3))) {
    expect_error(do.call(as_qtl2mr_prepared_input, c(raw, list(hsq = h))), "hsq")
  }
  expect_error(do.call(as_qtl2mr_prepared_input,
    c(raw, list(hsq = list(exposure = .2, x = .3)))), "Conflicting hsq")
  expect_error(do.call(as_qtl2mr_prepared_input,
    c(v041_numeric_input(), list(hsq = .2))), "requires kinship")
  expect_error(qtl2mr:::qtl2mr_validate_hsq(-.001), "not clamped")
  expect_error(qtl2mr:::qtl2mr_validate_hsq(1.001), "not clamped")
  expect_identical(qtl2mr:::qtl2mr_validate_hsq(0), 0)
  expect_identical(qtl2mr:::qtl2mr_validate_hsq(1), 1)
})

test_that("four-to-six decimal probabilities have disclosed bounded normalization", {
  for (digits in 4:6) {
    raw <- v041_numeric_input()
    raw$allele_probs[[1]] <- round(raw$allele_probs[[1]], digits)
    original <- raw$allele_probs
    expect_warning(input <- do.call(as_qtl2mr_prepared_input, raw), "Founder probability rounding")
    expect_identical(input$allele_probs, original)
    record <- input$metadata$probability_normalization
    expect_identical(record$policy, "normalize_rounding")
    expect_equal(record$rounding_sum_tolerance, 4e-4)
    expect_true(record$original_arrays_preserved)
    expect_gt(sum(record$adjustments$n_normalized), 0)
    normalized <- qtl2mr:::qtl2mr_prepared_probability_copy(input, "m1")
    sums <- rowSums(original[[1]][, , "m1"])
    adjust <- abs(sums - 1) > 1e-6 + 16 * .Machine$double.eps
    expect_equal(rowSums(normalized$allele_probs[[1]][adjust, , "m1"]), rep(1, sum(adjust)), ignore_attr = TRUE)
    expect_identical(normalized$allele_probs[[1]][!adjust, , "m1"], original[[1]][!adjust, , "m1"])
    expect_identical(normalized$allele_probs[[1]][, , "m2"], original[[1]][, , "m2"])
    expect_identical(input$allele_probs, original)
    expect_error(do.call(as_qtl2mr_prepared_input, c(raw, list(probability_policy = "strict"))), "sum to 1")
  }
  raw <- v041_numeric_input()
  input <- do.call(as_qtl2mr_prepared_input, raw)
  expect_equal(nrow(input$metadata$probability_normalization$adjustments), 0)
  expect_identical(qtl2mr:::qtl2mr_prepared_probability_copy(input, "m1")$allele_probs, raw$allele_probs)
})

test_that("a disclosed resource rounding policy is retained without per-pair warning spam", {
  raw <- v041_numeric_input()
  raw$allele_probs[[1]] <- round(raw$allele_probs[[1]], 4L)
  raw$selected_loci <- c("m1", "m2")
  expect_warning(all_input <- do.call(as_qtl2mr_prepared_input, raw), "Founder probability rounding")
  raw$selected_loci <- "m1"
  raw$metadata <- all_input$metadata
  expect_no_warning(one_input <- do.call(as_qtl2mr_prepared_input, raw))
  expect_no_warning(local_copy <- qtl2mr:::qtl2mr_prepared_probability_copy(one_input, "m1"))
  expect_gt(sum(local_copy$metadata$probability_normalization$adjustments$n_normalized), 0)
  expect_identical(one_input$allele_probs, raw$allele_probs)
  expect_identical(local_copy$metadata$probability_normalization$adjustments$marker, "m1")
})

test_that("rounding allowance never repairs invalid probabilities or dimensions", {
  for (bad in c(NaN, Inf)) {
    raw <- v041_numeric_input()
    raw$allele_probs[[1]][1, 1, 1] <- bad
    expect_error(do.call(as_qtl2mr_prepared_input, raw), "finite founder allele probabilities")
  }
  for (bad in c(-1e-6, 1 + 1e-6)) {
    raw <- v041_numeric_input()
    raw$allele_probs[[1]][1, 1, 1] <- bad
    expect_error(do.call(as_qtl2mr_prepared_input, raw), "lie in \\[0, 1\\].*observed range")
  }
  # Floating-point overshoot within 1e-8 of the bounds is accepted, as in v0.4.0.
  for (row in list(c(1 + 1e-12, rep(0, 7L)), c(-1e-12, 1 + 1e-12, rep(0, 6L)))) {
    raw <- v041_numeric_input()
    raw$allele_probs[[1]][1, , 1] <- row
    expect_s3_class(do.call(as_qtl2mr_prepared_input, raw), "qtl2mr_prepared_input")
  }
  raw <- v041_numeric_input()
  raw$allele_probs[[1]][1, , 1] <- raw$allele_probs[[1]][1, , 1] * 1.001
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "sum to 1")
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "maximum observed.*allowed tolerance.*normalize_rounding")
  raw <- v041_numeric_input()
  raw$allele_probs[[1]] <- raw$allele_probs[[1]][, -1, , drop = FALSE]
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "sum to 1|founder dimension")
  raw <- v041_numeric_input()
  raw$allele_probs[[1]] <- raw$allele_probs[[1]][, , 1]
  expect_error(do.call(as_qtl2mr_prepared_input, raw), "3D|sample x founder x marker")
})

test_that("unavailable theta and zero exposure diagnostics remain explicit and safe", {
  block <- list(beta_x = setNames(rep(0, 7), paste0("FC", 1:7)),
    beta_y = setNames(rep(1, 7), paste0("FC", 1:7)), Sigma_x = diag(.1, 7),
    Sigma_y = diag(.1, 7), Sigma_xy = matrix(0, 7, 7),
    founder_H = make_founder_contrasts(LETTERS[1:8]), marker = "degenerate", chr = "1")
  diagnostic <- mr_block_so_diagnostics(block)
  expect_true(is.na(diagnostic$theta))
  expect_identical(diagnostic$diagnostic_status, "unavailable")
  expect_true(is.na(diagnostic$omega_min_eigenvalue))
  expect_match(diagnostic$diagnostic_note, "theta is unavailable")
  expect_identical(mr_block_so_diagnostics(block, theta = NA_real_)$diagnostic_status, "unavailable")
  for (theta in list(NULL, NA_real_, 0)) {
    residual <- founder_residual_scores(block, theta = theta)
    expect_equal(nrow(residual), 8L)
    expect_true(all(residual$diagnostic_status == "unavailable"))
    expect_true(all(is.na(residual$z_diag)))
  }
  global <- fit_multilocus_founder_mr(list(block), heterogeneity_method = "none")
  saved <- global
  layer <- fit_multilocus_diagnostic_layer(list(block), global_fit = global, bootstrap_B = 0L)
  expect_identical(global, saved)
  expect_identical(layer$global, saved)
  expect_true(all(is.na(layer$estimates$theta)))
  expect_true(all(layer$per_locus$diagnostic_status == "unavailable"))
  expect_true(is.na(layer$heterogeneity$Q_between))
})

test_that("single prepared covariate exclusions are warned without changing sample policy", {
  raw <- v041_numeric_input()
  raw$covar$age[2] <- Inf
  input <- do.call(as_qtl2mr_prepared_input, raw)
  expect_warning(fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = "none"), "Missing covariate rows or non-finite covariate values")
  expect_true(all(fit$sample_alignment$n_incomplete_covar == 1L))
  expect_false(rownames(raw$covar)[2] %in% fit$locus_blocks[[1]]$fit_x$sample_ids)
  expect_false(rownames(raw$covar)[2] %in% fit$locus_blocks[[1]]$fit_y$sample_ids)
})
