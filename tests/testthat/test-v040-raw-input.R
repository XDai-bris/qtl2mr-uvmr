test_that("nonzero kinship and missing traits retain explicit overlap and units", {
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  d <- qtl2mr_example
  ids <- rownames(d$pheno)
  set.seed(404017)
  background <- matrix(rnorm(length(ids) * 16), nrow = length(ids))
  K <- tcrossprod(background) / ncol(background)
  K <- (K / mean(diag(K)) + diag(length(ids))) / 2
  dimnames(K) <- list(ids, ids)
  masks <- list(full = list(x = integer(), y = integer()),
    partial = list(x = 81:96, y = 1:16),
    none = list(x = 49:96, y = 1:48))
  for (scenario in names(masks)) {
    pheno <- d$pheno
    pheno$x[masks[[scenario]]$x] <- NA_real_
    pheno$y[masks[[scenario]]$y] <- NA_real_
    input <- as_qtl2mr_prepared_input(d$allele_probs, pheno,
      covar = NULL, kinship = list("1" = K), marker_map = d$marker_map,
      selected_loci = d$selected_loci, exposure = "x", outcome = "y")
    fit <- qtl2mr:::run_prepared_qtl2mr(input, "x", "y", verbose = FALSE,
      standardize_phenotypes = FALSE, heterogeneity_method = "none",
      diagnostic_bootstrap_B = 0L)
    expect_true(all(vapply(fit$locus_blocks, function(b) {
      b$n_overlap == sum(is.finite(pheno$x) & is.finite(pheno$y))
    }, logical(1))), info = scenario)
    expect_identical(fit$summary$effect_units, "outcome units per exposure unit")
    expect_equal(fit$phenotype_scaling$n[fit$phenotype_scaling$role == "exposure"],
      rep(sum(is.finite(pheno$x)), 2L))
    expect_equal(fit$sample_alignment$n_dropped[fit$sample_alignment$role == "outcome"],
      rep(sum(!is.finite(pheno$y)), 2L))
    expect_true(all(vapply(fit$locus_blocks, function(b) is.finite(b$fit_x$hsq) && is.finite(b$fit_y$hsq), logical(1))))
    if (scenario == "none") expect_true(all(vapply(fit$locus_blocks,
      function(b) max(abs(b$Sigma_xy)) < 1e-12, logical(1))))
  }
})
