test_that("direct founder LMM returns covariance matching SE", {
  gp <- make_sim_genoprobs(n = 70)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  covar <- matrix(rnorm(length(ids) * 2), ncol = 2, dimnames = list(ids, c("c1", "c2")))
  y <- as.numeric(Z[, "FC1"] * 0.5 - Z[, "FC4"] * 0.25 + covar[, 1] * 0.2 + rnorm(length(ids), sd = 0.5))
  names(y) <- ids
  fit <- fit_founder_lmm_block(gp, y, addcovar = covar, contrasts = A)
  expect_equal(as.numeric(diag(fit$vcov_fc)), as.numeric(fit$se_fc^2), tolerance = 1e-10)
  expect_equal(dim(fit$influence_raw), c(7L, length(ids)))
  expect_equal(names(fit$beta_fc), paste0("FC", 1:7))
})

test_that("direct founder LMM matches qtl2::scan1coef beta and SE without kinship", {
  skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 60)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , 1] %*% A
  covar <- matrix(rnorm(length(ids)), ncol = 1, dimnames = list(ids, "c1"))
  y <- as.numeric(Z[, "FC2"] * 0.7 + covar[, 1] * 0.3 + rnorm(length(ids), sd = 0.4))
  names(y) <- ids
  cmp <- compare_founder_lmm_to_scan1coef(gp, y, addcovar = covar, contrasts = A)
  expect_lt(cmp$max_abs_beta_diff, 1e-8)
  expect_lt(cmp$max_abs_se_diff, 1e-8)
})
