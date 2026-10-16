qtl2mr_matrix_sqrt <- function(S, ridge = 1e-10) {
  qtl2mr_matrix_power(S, 0.5, ridge = ridge)
}

qtl2mr_one_vs_rest_matrix <- function(founders) {
  p <- length(founders)
  V <- matrix(-1 / (p - 1L), nrow = p, ncol = p, dimnames = list(founders, founders))
  diag(V) <- 1
  V
}

founder_residual_null_simulation <- function(block,
                                             theta = NULL,
                                             H = NULL,
                                             nsim = 10000L,
                                             seed = NULL,
                                             allow_cross_cov = TRUE,
                                             ridge = 1e-10) {
  if (is.null(theta)) {
    theta <- second_order_gls_ivw(
      block$beta_x,
      block$beta_y,
      block$Sigma_x,
      block$Sigma_y,
      allow_cross_cov = allow_cross_cov,
      Sigma_xy = qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov),
      ridge = ridge
    )$theta
  }
  H <- qtl2mr_resolve_founder_H(block, H)
  observed <- founder_residual_scores(
    block,
    theta = theta,
    H = H,
    allow_cross_cov = allow_cross_cov,
    ridge = ridge
  )
  founders <- rownames(H)
  bx <- matrix(qtl2mr_block_vector(block, "beta_x"), ncol = 1L)
  residual_cov <- mr_block_residual_covariance(block, theta,
    allow_cross_cov = allow_cross_cov, ridge = ridge)
  V_r <- residual_cov$V_r
  V_e <- qtl2mr_sym(H %*% V_r %*% t(H))
  V_one <- qtl2mr_one_vs_rest_matrix(founders)
  denom <- sqrt(pmax(diag(t(V_one) %*% V_e %*% V_one), .Machine$double.eps))
  ax <- as.numeric(H %*% bx)
  names(ax) <- founders

  qtl2mr_with_seed(seed, {
    L <- qtl2mr_matrix_sqrt(V_r, ridge = 0)
    z <- matrix(stats::rnorm(nsim * nrow(V_r)), nrow = nrow(V_r), ncol = nsim)
    r_sim <- L %*% z
    e_sim <- H %*% r_sim
    t_sim <- t(V_one) %*% e_sim
    t_sim <- sweep(t_sim, 1L, denom, "/")
    rownames(t_sim) <- founders
    abs_t <- abs(t(t_sim))
    colnames(abs_t) <- founders
    top_idx <- max.col(abs_t, ties.method = "random")
    top_founder <- colnames(abs_t)[top_idx]

    rows <- lapply(founders, function(f) {
      null_abs <- abs_t[, f]
      obs <- observed[observed$founder == f, , drop = FALSE]
      obs_abs <- obs$abs_t_one_founder_vs_rest[[1L]]
      p95 <- as.numeric(stats::quantile(null_abs, 0.95, na.rm = TRUE, names = FALSE))
      p99 <- as.numeric(stats::quantile(null_abs, 0.99, na.rm = TRUE, names = FALSE))
      data.frame(
        marker = block$marker %||% NA_character_,
        chr = block$chr %||% NA_character_,
        founder = f,
        abs_aX_f = abs(ax[[f]]),
        null_mean_abs_t = mean(null_abs, na.rm = TRUE),
        null_p95_abs_t = p95,
        null_p99_abs_t = p99,
        null_top_founder_frequency = mean(top_founder == f),
        observed_abs_t = obs_abs,
        observed_percentile = mean(null_abs <= obs_abs, na.rm = TRUE),
        empirical_p_abs_t_ge_observed = mean(null_abs >= obs_abs, na.rm = TRUE),
        interpretation_flag = ifelse(obs_abs > p95, "above_null_p95", "within_null_p95"),
        nsim = nsim,
        theta = theta,
        stringsAsFactors = FALSE
      )
    })
    summary <- data.frame(
      marker = block$marker %||% NA_character_,
      chr = block$chr %||% NA_character_,
      theta = theta,
      nsim = nsim,
      observed_max_abs_t = max(observed$abs_t_one_founder_vs_rest, na.rm = TRUE),
      null_mean_max_abs_t = mean(apply(abs_t, 1L, max), na.rm = TRUE),
      null_p95_max_abs_t = as.numeric(stats::quantile(apply(abs_t, 1L, max), 0.95, names = FALSE)),
      empirical_p_max_ge_observed = mean(apply(abs_t, 1L, max) >= max(observed$abs_t_one_founder_vs_rest, na.rm = TRUE)),
      observed_count_abs_t_gt_2 = sum(observed$abs_t_one_founder_vs_rest > 2, na.rm = TRUE),
      null_mean_count_abs_t_gt_2 = mean(rowSums(abs_t > 2), na.rm = TRUE),
      empirical_p_count_ge_observed = mean(rowSums(abs_t > 2) >= sum(observed$abs_t_one_founder_vs_rest > 2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    list(
      founder = do.call(rbind, rows),
      summary = summary,
      observed = observed,
      null_abs_t = abs_t,
      top_founder = top_founder
    )
  })
}

summarize_founder_residual_null <- function(null_result) {
  if (is.data.frame(null_result)) return(null_result)
  if (!is.null(null_result$founder)) return(null_result$founder)
  stop("null_result must be a data frame or a result from founder_residual_null_simulation()", call. = FALSE)
}
