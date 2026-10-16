test_that("direct founder LMM matches qtl2::scan1coef beta and SE with fixed kinship hsq", {
  skip_if_not_installed("qtl2")
  gp <- make_sim_genoprobs(n = 68, marker = "m1", seed = 501)
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z <- gp[[1]][, , "m1"] %*% A
  covar <- matrix(rnorm(length(ids) * 2), ncol = 2, dimnames = list(ids, c("c1", "c2")))
  K0 <- matrix(rnorm(length(ids)^2), nrow = length(ids))
  K <- tcrossprod(K0)
  K <- K / mean(diag(K))
  dimnames(K) <- list(ids, ids)
  y <- as.numeric(0.4 * Z[, "FC2"] - 0.25 * Z[, "FC5"] + 0.3 * covar[, "c1"] + rnorm(length(ids), sd = 0.5))
  names(y) <- ids

  cmp <- compare_founder_lmm_to_scan1coef(
    genoprobs = gp,
    pheno = y,
    kinship = K,
    addcovar = covar,
    contrasts = A,
    hsq = 0.35,
    marker = "m1"
  )

  expect_lt(cmp$max_abs_beta_diff, 1e-8)
  expect_lt(cmp$max_abs_se_diff, 1e-8)
})
