source_all_qtl2mr <- function() {
  root <- normalizePath(file.path(dirname(getwd())), mustWork = FALSE)
  if (!file.exists(file.path(root, "DESCRIPTION"))) root <- getwd()
  for (f in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
}

make_sim_genoprobs <- function(n = 80, marker = "m1", seed = 1) {
  set.seed(seed)
  founders <- LETTERS[1:8]
  raw <- matrix(stats::rexp(n * 8), nrow = n, ncol = 8)
  probs <- raw / rowSums(raw)
  arr <- array(probs, dim = c(n, 8, 1), dimnames = list(paste0("id", seq_len(n)), founders, marker))
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  gp
}

make_fake_qtl2mr_analysis_result <- function() {
  gp1 <- make_sim_genoprobs(n = 72, marker = "m1", seed = 101)
  set.seed(102)
  raw2 <- matrix(stats::rexp(72 * 8), nrow = 72, ncol = 8)
  probs2 <- raw2 / rowSums(raw2)
  arr <- array(
    NA_real_,
    dim = c(72, 8, 2),
    dimnames = list(dimnames(gp1[[1]])[[1]], LETTERS[1:8], c("m1", "m2"))
  )
  arr[, , "m1"] <- gp1[[1]][, , 1]
  arr[, , "m2"] <- probs2
  gp <- list("1" = arr)
  class(gp) <- c("calc_genoprob", "list")
  ids <- dimnames(gp[[1]])[[1]]
  A <- make_founder_contrasts(LETTERS[1:8])
  Z1 <- gp[[1]][, , "m1"] %*% A
  Z2 <- gp[[1]][, , "m2"] %*% A
  x <- as.numeric(0.7 * Z1[, "FC1"] + 0.5 * Z2[, "FC2"] + stats::rnorm(length(ids), sd = 0.5))
  y <- as.numeric(0.6 * x + 0.2 * Z1[, "FC3"] + stats::rnorm(length(ids), sd = 0.5))
  names(x) <- ids
  names(y) <- ids
  blocks <- list(
    m1 = fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m1"),
    m2 = fit_mr_locus_block(gp, x, y, contrasts = A, marker = "m2")
  )
  fit <- fit_multilocus_founder_mr(blocks, q_sim_n = 30L, q_sim_seed = 11L)
  diagnostic_layer <- fit_multilocus_diagnostic_layer(blocks, global_fit = fit, bootstrap_B = 40L, bootstrap_seed = 3L)
  get_est <- function(name, col) {
    rows <- diagnostic_layer$estimates$estimator == name
    if (!any(rows)) return(NA_real_)
    diagnostic_layer$estimates[[col]][which(rows)[1L]]
  }
  summary <- data.frame(
    exposure = "exposure_trait",
    outcome = "outcome_trait",
    n_loci = length(blocks),
    theta = fit$theta,
    se = fit$se,
    profile_ci_lower = fit$profile_ci_lower,
    profile_ci_upper = fit$profile_ci_upper,
    profile_p = fit$profile_p,
    theta_independence = fit$independence$theta,
    theta_overlap_delta = fit$theta_overlap_delta,
    theta_two_step_ivw = get_est("Two-step IVW", "theta"),
    se_two_step_ivw = get_est("Two-step IVW", "se"),
    theta_weighted_median = get_est("Weighted median", "theta"),
    se_weighted_median = get_est("Weighted median", "se"),
    ci_lower_weighted_median = get_est("Weighted median", "ci_lower"),
    ci_upper_weighted_median = get_est("Weighted median", "ci_upper"),
    theta_weighted_mode = get_est("Weighted mode", "theta"),
    se_weighted_mode = get_est("Weighted mode", "se"),
    ci_lower_weighted_mode = get_est("Weighted mode", "ci_lower"),
    ci_upper_weighted_mode = get_est("Weighted mode", "ci_upper"),
    diagnostic_bootstrap_B = 40L,
    diagnostic_bootstrap_seed = 3L,
    Q_between = diagnostic_layer$heterogeneity$Q_between[[1L]],
    Q_between_df = diagnostic_layer$heterogeneity$df_between[[1L]],
    Q_between_p = diagnostic_layer$heterogeneity$Q_between_p[[1L]],
    Q = fit$Q,
    Q_chisq_df = fit$Q_chisq_df %||% NA_integer_,
    Q_p_chisq = fit$Q_p_chisq %||% NA_real_,
    Q_p_sim = fit$Q_p_sim %||% NA_real_,
    Q_sim_n = fit$Q_sim_n %||% NA_integer_,
    Q_sim_seed = fit$Q_sim_seed %||% NA_integer_,
    Q_sim_exceedances = fit$Q_sim_exceedances %||% NA_integer_,
    Q_sim_min_resolvable_p = fit$Q_sim_min_resolvable_p %||% NA_real_,
    Q_sim_null_mean = fit$Q_sim_null_mean %||% NA_real_,
    Q_sim_null_sd = fit$Q_sim_null_sd %||% NA_real_,
    Q_sim_method = fit$Q_sim_method %||% NA_character_,
    Q_p_reported = fit$Q_p_reported %||% NA_real_,
    Q_p_method = fit$Q_p_method %||% NA_character_,
    Q_p_note = fit$Q_p_note %||% NA_character_,
    cross_cov_model = "star",
    n_loci_psd_adjusted = sum(fit$locus_diagnostics$psd_adjusted, na.rm = TRUE),
    n_loci_psd_invalid = sum(!fit$locus_diagnostics$joint_psd_valid, na.rm = TRUE),
    n_loci_raw_space_psd_adjusted = sum(fit$locus_diagnostics$raw_space_psd_adjusted, na.rm = TRUE),
    n_loci_raw_space_psd_invalid = sum(!fit$locus_diagnostics$raw_space_psd_valid, na.rm = TRUE),
    n_loci_star_space_psd_adjusted = sum(fit$locus_diagnostics$star_space_psd_adjusted, na.rm = TRUE),
    n_loci_star_space_psd_invalid = sum(!fit$locus_diagnostics$star_space_psd_valid, na.rm = TRUE),
    standardize_phenotypes = TRUE,
    stringsAsFactors = FALSE
  )
  selected <- data.frame(marker = names(blocks), chr = c("1", "1"), locus_rank = 1:2, stringsAsFactors = FALSE)
  founder_residuals <- do.call(rbind, lapply(blocks, founder_residual_scores, theta = fit$theta))
  founder_ratios <- do.call(rbind, lapply(blocks, founder_coordinate_ratios, theta_ref = fit$theta))

  pilot <- list(
    summary = summary,
    fit = fit,
    locus_blocks = blocks,
    selected_loci = selected,
    locus_diagnostics = fit$locus_diagnostics,
    diagnostic_layer = diagnostic_layer,
    founder_residual_scores_global = founder_residuals,
    founder_coordinate_ratios = founder_ratios,
    sample_ids = ids,
    backend_provenance = qtl2mr:::qtl2mr_lmm_backend_provenance(blocks)
  )
  structure(
    list(
      input = tempdir(),
      exposure = "exposure_trait",
      outcome = "outcome_trait",
      output_dir = tempfile("qtl2mr-result-"),
      created_at = Sys.time(),
      pilot = pilot
    ),
    class = c("qtl2mr_analysis_result", "list")
  )
}

compare_founder_lmm_to_scan1coef <- qtl2mr:::compare_founder_lmm_to_scan1coef
write_qtl2mr_results <- qtl2mr:::write_qtl2mr_results
make_overlap_matrix <- qtl2mr:::make_overlap_matrix
estimate_rho_star <- qtl2mr:::estimate_rho_star
mr_block_omega <- qtl2mr:::mr_block_omega
mr_block_q <- qtl2mr:::mr_block_q
mr_block_first_stage_f <- qtl2mr:::mr_block_first_stage_f
fit_per_locus_so <- qtl2mr:::fit_per_locus_so
