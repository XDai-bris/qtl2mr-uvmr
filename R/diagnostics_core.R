qtl2mr_block_vector <- function(block, name) {
  x <- block[[name]]
  if (is.null(x)) stop("Locus block is missing ", name, call. = FALSE)
  stats::setNames(as.numeric(x), names(x))
}

qtl2mr_block_cross <- function(block, allow_cross_cov = TRUE) {
  bx <- qtl2mr_block_vector(block, "beta_x")
  if (isTRUE(allow_cross_cov)) {
    Sxy <- block$Sigma_xy
    if (is.null(Sxy)) stop("Locus block is missing Sigma_xy", call. = FALSE)
    return(as.matrix(Sxy))
  }
  matrix(0, nrow = length(bx), ncol = length(bx), dimnames = list(names(bx), names(bx)))
}

qtl2mr_resolve_founder_H <- function(block = NULL, H = NULL) {
  bx_names <- if (!is.null(block)) names(qtl2mr_block_vector(block, "beta_x")) else NULL
  if (is.null(H) && !is.null(block$founder_H)) H <- block$founder_H
  if (is.null(H) && !is.null(block$contrasts)) H <- block$contrasts
  if (is.null(H)) {
    stop("A founder contrast matrix H is required for founder-coordinate diagnostics", call. = FALSE)
  }
  H <- as.matrix(H)
  if (!is.null(bx_names) && !is.null(colnames(H)) && all(bx_names %in% colnames(H))) {
    H <- H[, bx_names, drop = FALSE]
  } else if (!is.null(bx_names) && ncol(H) == length(bx_names) + 1L) {
    H <- H[, -1L, drop = FALSE]
    colnames(H) <- bx_names
  } else if (!is.null(bx_names) && ncol(H) != length(bx_names)) {
    stop("Founder contrast matrix H has incompatible dimensions", call. = FALSE)
  }
  if (is.null(rownames(H))) rownames(H) <- paste0("founder_", seq_len(nrow(H)))
  if (!is.null(bx_names)) colnames(H) <- bx_names
  H
}

mr_block_omega <- function(block,
                           theta,
                           allow_cross_cov = TRUE,
                           ridge = 0) {
  bx <- qtl2mr_block_vector(block, "beta_x")
  Sx <- as.matrix(block$Sigma_x)
  Sy <- as.matrix(block$Sigma_y)
  Sxy <- qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov)
  if (!all(dim(Sx) == length(bx)) || !all(dim(Sy) == length(bx)) || !all(dim(Sxy) == length(bx))) {
    stop("Covariance matrices must match beta_x length", call. = FALSE)
  }
  Om <- qtl2mr_sym(Sy + theta^2 * Sx - theta * (Sxy + t(Sxy)))
  qtl2mr_add_invariant_ridge(Om, ridge)
}

mr_block_residual <- function(block, theta, H = NULL) {
  bx <- qtl2mr_block_vector(block, "beta_x")
  by <- qtl2mr_block_vector(block, "beta_y")
  d <- stats::setNames(by - theta * bx, names(bx))
  out <- list(d = d)
  if (!is.null(H) || !is.null(block$founder_H) || !is.null(block$contrasts)) {
    H <- qtl2mr_resolve_founder_H(block, H)
    e <- as.numeric(H %*% matrix(d, ncol = 1L))
    out$e <- stats::setNames(e, rownames(H))
  }
  out
}

mr_block_q <- function(block,
                       theta,
                       allow_cross_cov = TRUE,
                       ridge = 1e-10) {
  d <- matrix(mr_block_residual(block, theta)$d, ncol = 1L)
  Om <- mr_block_omega(block, theta, allow_cross_cov = allow_cross_cov, ridge = ridge)
  as.numeric(t(d) %*% qtl2mr_mat_inv(Om, ridge = 0) %*% d)
}

mr_block_residual_covariance <- function(block,
                                         theta,
                                         allow_cross_cov = TRUE,
                                         ridge = 1e-10) {
  bx <- matrix(qtl2mr_block_vector(block, "beta_x"), ncol = 1L)
  Omega_raw <- mr_block_omega(block, theta,
    allow_cross_cov = allow_cross_cov, ridge = 0)
  Omega_work <- qtl2mr_add_invariant_ridge(Omega_raw, ridge)
  W <- qtl2mr_mat_inv(Omega_work, ridge = 0)
  denom <- as.numeric(t(bx) %*% W %*% bx)
  if (!is.finite(denom) || denom <= 0) {
    stop("Residual diagnostic unavailable: exposure direction has no positive finite information", call. = FALSE)
  }
  P_perp <- diag(nrow(Omega_raw)) - bx %*% (t(bx) %*% W) / denom
  V_r <- qtl2mr_sym(P_perp %*% Omega_raw %*% t(P_perp))
  list(
    Omega_raw = Omega_raw,
    Omega_work = Omega_work,
    W = W,
    P_perp = P_perp,
    V_r = V_r
  )
}

mr_block_first_stage_f <- function(block, ridge = 1e-10) {
  bx <- matrix(qtl2mr_block_vector(block, "beta_x"), ncol = 1L)
  Sx_inv <- qtl2mr_mat_inv(block$Sigma_x, ridge = ridge)
  as.numeric(t(bx) %*% Sx_inv %*% bx) / length(bx)
}

mr_block_so_diagnostics <- function(block,
                                    theta = NULL,
                                    H = NULL,
                                    allow_cross_cov = TRUE,
                                    ridge = 1e-10) {
  bx <- qtl2mr_block_vector(block, "beta_x")
  fit <- second_order_gls_ivw(
    beta_x = bx,
    beta_y = qtl2mr_block_vector(block, "beta_y"),
    Sigma_x = block$Sigma_x,
    Sigma_y = block$Sigma_y,
    allow_cross_cov = allow_cross_cov,
    Sigma_xy = qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov),
    ridge = ridge
  )
  if (!is.null(theta)) {
    fit$theta <- theta
    fit$Q <- if (length(theta) == 1L && is.finite(theta)) {
      mr_block_q(block, theta, allow_cross_cov = allow_cross_cov, ridge = ridge)
    } else NA_real_
  }
  diagnostic_note <- ""
  geometry <- tryCatch({
    if (length(fit$theta) != 1L || !is.finite(fit$theta)) stop("SO theta is unavailable")
    Om_raw <- mr_block_omega(block, fit$theta, allow_cross_cov = allow_cross_cov, ridge = 0)
    list(ev = eigen(Om_raw, symmetric = TRUE, only.values = TRUE)$values,
      residual = mr_block_residual(block, fit$theta, H = H))
  }, error = function(e) {
    diagnostic_note <<- paste("Diagnostic unavailable:", conditionMessage(e))
    list(ev = NA_real_, residual = NULL)
  })
  ev <- geometry$ev
  list(
    marker = block$marker %||% NA_character_,
    chr = block$chr %||% NA_character_,
    theta = fit$theta,
    se = fit$se,
    profile_ci_lower = fit$profile_ci_lower,
    profile_ci_upper = fit$profile_ci_upper,
    Q_min = fit$Q,
    Q_df = length(bx) - 1L,
    Q_p_chisq = if (nzchar(diagnostic_note)) NA_real_ else stats::pchisq(fit$Q, df = length(bx) - 1L, lower.tail = FALSE),
    F_X_block = mr_block_first_stage_f(block, ridge = ridge),
    omega_min_eigenvalue = min(ev),
    omega_condition_number = max(ev) / min(ev),
    residual = geometry$residual,
    diagnostic_status = if (nzchar(diagnostic_note)) "unavailable" else "available",
    diagnostic_note = diagnostic_note,
    fit = fit
  )
}

founder_residual_scores <- function(block,
                                    theta = NULL,
                                    H = NULL,
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
  unavailable <- function(note) data.frame(
    marker = block$marker %||% NA_character_, chr = block$chr %||% NA_character_,
    theta = if (length(theta) == 1L) theta else NA_real_, founder = rownames(H),
    e = NA_real_, z_diag = NA_real_, t_one_founder_vs_rest = NA_real_,
    abs_t_one_founder_vs_rest = NA_real_, strong_abs_t_gt_2 = NA,
    diagnostic_status = "unavailable", diagnostic_note = paste("Diagnostic unavailable:", note),
    stringsAsFactors = FALSE)
  if (length(theta) != 1L || !is.finite(theta)) return(unavailable("SO theta is unavailable"))
  d <- matrix(mr_block_residual(block, theta)$d, ncol = 1L)
  residual_cov <- tryCatch(mr_block_residual_covariance(block, theta,
    allow_cross_cov = allow_cross_cov, ridge = ridge), error = identity)
  if (inherits(residual_cov, "error")) return(unavailable(conditionMessage(residual_cov)))
  V_r <- residual_cov$V_r
  V_e <- qtl2mr_sym(H %*% V_r %*% t(H))
  e <- as.numeric(H %*% d)
  founders <- rownames(H)
  names(e) <- founders
  t_one <- vapply(seq_along(founders), function(i) {
    v <- rep(-1 / (length(founders) - 1L), length(founders))
    v[[i]] <- 1
    vv <- matrix(v, ncol = 1L)
    den <- sqrt(as.numeric(t(vv) %*% V_e %*% vv))
    as.numeric(t(vv) %*% matrix(e, ncol = 1L)) / den
  }, numeric(1))
  data.frame(
    marker = block$marker %||% NA_character_,
    chr = block$chr %||% NA_character_,
    theta = theta,
    founder = founders,
    e = e,
    z_diag = e / sqrt(pmax(diag(V_e), .Machine$double.eps)),
    t_one_founder_vs_rest = t_one,
    abs_t_one_founder_vs_rest = abs(t_one),
    strong_abs_t_gt_2 = abs(t_one) > 2,
    diagnostic_status = "available",
    diagnostic_note = "",
    stringsAsFactors = FALSE
  )
}

qtl2mr_ratio_table <- function(x,
                               y,
                               Sx,
                               Sy,
                               Sxy,
                               labels,
                               theta_ref = NA_real_,
                               weak_f = 10,
                               level = 0.95,
                               marker = NA_character_,
                               chr = NA_character_,
                               coordinate_type = "coordinate") {
  zcrit <- stats::qnorm((1 + level) / 2)
  rows <- lapply(seq_along(labels), function(i) {
    x_i <- x[[i]]
    y_i <- y[[i]]
    var_x <- Sx[i, i]
    var_y <- Sy[i, i]
    cov_xy <- Sxy[i, i]
    theta_i <- if (is.finite(x_i) && abs(x_i) > .Machine$double.eps) y_i / x_i else NA_real_
    var_theta <- if (is.finite(theta_i)) {
      var_y / x_i^2 + y_i^2 * var_x / x_i^4 - 2 * y_i * cov_xy / x_i^3
    } else {
      NA_real_
    }
    se <- if (is.finite(var_theta) && var_theta >= 0) sqrt(var_theta) else NA_real_
    F_x <- x_i^2 / var_x
    data.frame(
      marker = marker,
      chr = chr,
      coordinate_type = coordinate_type,
      coordinate = labels[[i]],
      theta_ref = theta_ref,
      theta_diag = theta_i,
      se = se,
      ci_lower = theta_i - zcrit * se,
      ci_upper = theta_i + zcrit * se,
      x = x_i,
      y = y_i,
      var_x = var_x,
      var_y = var_y,
      cov_xy = cov_xy,
      F_X_diag = F_x,
      weak_denominator = is.finite(F_x) && F_x < weak_f,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

founder_coordinate_ratios <- function(block,
                                      H = NULL,
                                      theta_ref = NULL,
                                      allow_cross_cov = TRUE,
                                      weak_f = 10,
                                      level = 0.95) {
  H <- qtl2mr_resolve_founder_H(block, H)
  Sxy <- qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov)
  ax <- as.numeric(H %*% matrix(qtl2mr_block_vector(block, "beta_x"), ncol = 1L))
  ay <- as.numeric(H %*% matrix(qtl2mr_block_vector(block, "beta_y"), ncol = 1L))
  Sx <- H %*% as.matrix(block$Sigma_x) %*% t(H)
  Sy <- H %*% as.matrix(block$Sigma_y) %*% t(H)
  Sxy_f <- H %*% Sxy %*% t(H)
  qtl2mr_ratio_table(
    ax,
    ay,
    Sx,
    Sy,
    Sxy_f,
    labels = rownames(H),
    theta_ref = theta_ref %||% NA_real_,
    weak_f = weak_f,
    level = level,
    marker = block$marker %||% NA_character_,
    chr = block$chr %||% NA_character_,
    coordinate_type = "founder"
  )
}

founder_contrast_ratios <- function(block,
                                    theta_ref = NULL,
                                    allow_cross_cov = TRUE,
                                    weak_f = 10,
                                    level = 0.95) {
  bx <- qtl2mr_block_vector(block, "beta_x")
  by <- qtl2mr_block_vector(block, "beta_y")
  qtl2mr_ratio_table(
    as.numeric(bx),
    as.numeric(by),
    as.matrix(block$Sigma_x),
    as.matrix(block$Sigma_y),
    qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov),
    labels = names(bx),
    theta_ref = theta_ref %||% NA_real_,
    weak_f = weak_f,
    level = level,
    marker = block$marker %||% NA_character_,
    chr = block$chr %||% NA_character_,
    coordinate_type = "founder_contrast"
  )
}

second_order_nonorthogonality <- function(block,
                                          theta = NULL,
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
  bx <- matrix(qtl2mr_block_vector(block, "beta_x"), ncol = 1L)
  d <- matrix(mr_block_residual(block, theta)$d, ncol = 1L)
  Sx <- as.matrix(block$Sigma_x)
  Sxy <- qtl2mr_block_cross(block, allow_cross_cov = allow_cross_cov)
  Gamma <- Sxy + t(Sxy)
  Om_raw <- mr_block_omega(block, theta, allow_cross_cov = allow_cross_cov, ridge = 0)
  Om_work <- qtl2mr_add_invariant_ridge(Om_raw, ridge)
  W <- qtl2mr_mat_inv(Om_work, ridge = 0)
  exposure_term <- -theta * as.numeric(t(d) %*% W %*% Sx %*% W %*% d)
  overlap_term <- 0.5 * as.numeric(t(d) %*% W %*% Gamma %*% W %*% d)
  denom <- as.numeric(t(bx) %*% W %*% bx)
  data.frame(
    marker = block$marker %||% NA_character_,
    chr = block$chr %||% NA_character_,
    theta = theta,
    alpha_nonorthogonality = (exposure_term + overlap_term) / denom,
    exposure_uncertainty_component = exposure_term,
    overlap_component = overlap_term,
    denominator = denom,
    bx_w_d = as.numeric(t(bx) %*% W %*% d),
    stringsAsFactors = FALSE
  )
}

founder_local_score <- function(genoprobs,
                                marker = NULL,
                                beta_x,
                                contrasts = NULL) {
  selected <- qtl2mr_select_marker_genoprobs(genoprobs, marker = marker)
  gp_raw <- selected$genoprobs
  if (is.null(contrasts)) contrasts <- make_founder_contrasts(dimnames(gp_raw)[[2L]])
  Z <- qtl2mr_genoprobs_design(gp_raw, contrasts)
  beta_x <- qtl2mr_block_vector(list(beta_x = beta_x), "beta_x")
  if (is.null(names(beta_x)) || any(!nzchar(names(beta_x)))) {
    stop("beta_x must be a named founder-contrast vector", call. = FALSE)
  }
  if (!all(names(beta_x) %in% colnames(Z))) {
    stop("beta_x names must match founder contrast columns in the genotype design", call. = FALSE)
  }
  score <- as.numeric(Z[, names(beta_x), drop = FALSE] %*% matrix(beta_x, ncol = 1L))
  stats::setNames(score, rownames(Z))
}

qtl2mr_standardized_named_score <- function(score) {
  score <- qtl2mr_as_named_numeric(score, "score")
  s <- stats::sd(score, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop("Cannot standardize score with non-positive SD", call. = FALSE)
  stats::setNames(as.numeric((score - mean(score, na.rm = TRUE)) / s), names(score))
}

qtl2mr_xprs_vector <- function(xprs, sample_ids = NULL, column = "xprs") {
  if (is.list(xprs) && !is.null(xprs$sample_scores)) xprs <- xprs$sample_scores
  if (is.data.frame(xprs)) {
    if (!all(c("sample_id", column) %in% names(xprs))) {
      stop("xprs data frame must contain sample_id and ", column, " columns", call. = FALSE)
    }
    return(stats::setNames(as.numeric(xprs[[column]]), as.character(xprs$sample_id)))
  }
  out <- as.numeric(xprs)
  nm <- names(xprs)
  if ((is.null(nm) || any(!nzchar(nm))) && !is.null(sample_ids) && length(sample_ids) == length(out)) {
    nm <- as.character(sample_ids)
  }
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("xprs must be named, or sample_ids must match its length", call. = FALSE)
  }
  stats::setNames(out, nm)
}

qtl2mr_resolve_locus_block <- function(locus_blocks, marker) {
  if (!is.null(names(locus_blocks)) && marker %in% names(locus_blocks)) return(locus_blocks[[marker]])
  block_markers <- vapply(locus_blocks, function(b) as.character(b$marker %||% NA_character_), character(1))
  idx <- match(marker, block_markers)
  if (is.na(idx)) stop("No locus block found for marker: ", marker, call. = FALSE)
  locus_blocks[[idx]]
}

multilocus_xprs_scores <- function(genoprobs,
                                   locus_blocks,
                                   markers = NULL,
                                   sample_ids = NULL,
                                   standardize = TRUE) {
  if (!is.list(locus_blocks) || !length(locus_blocks)) {
    stop("locus_blocks must be a non-empty list", call. = FALSE)
  }
  if (is.null(markers)) {
    markers <- names(locus_blocks)
    if (is.null(markers) || any(!nzchar(markers))) {
      markers <- vapply(locus_blocks, function(b) as.character(b$marker %||% NA_character_), character(1))
    }
  }
  markers <- as.character(markers)
  markers <- markers[nzchar(markers) & !is.na(markers)]
  if (!length(markers)) stop("No markers supplied for multi-locus xPRS", call. = FALSE)

  score_list <- lapply(markers, function(marker) {
    block <- qtl2mr_resolve_locus_block(locus_blocks, marker)
    founder_local_score(
      genoprobs = genoprobs,
      marker = marker,
      beta_x = block$beta_x,
      contrasts = block$contrasts
    )
  })
  names(score_list) <- markers
  ids <- Reduce(intersect, lapply(score_list, names))
  if (!is.null(sample_ids)) ids <- intersect(as.character(sample_ids), ids)
  ids <- ids[vapply(ids, function(id) all(vapply(score_list, function(s) is.finite(s[[id]]), logical(1))), logical(1))]
  if (!length(ids)) stop("No complete samples available for multi-locus xPRS", call. = FALSE)

  score_mat <- do.call(cbind, lapply(score_list, function(s) as.numeric(s[ids])))
  rownames(score_mat) <- ids
  colnames(score_mat) <- markers
  xprs <- stats::setNames(rowSums(score_mat), ids)
  xprs_z <- if (isTRUE(standardize)) qtl2mr_standardized_named_score(xprs) else stats::setNames(rep(NA_real_, length(xprs)), ids)
  sample_scores <- data.frame(
    sample_id = ids,
    xprs = as.numeric(xprs),
    xprs_z = as.numeric(xprs_z),
    n_loci = ncol(score_mat),
    stringsAsFactors = FALSE
  )
  list(
    sample_scores = sample_scores,
    locus_score_matrix = score_mat,
    markers = markers,
    n_loci = length(markers)
  )
}

leave_one_locus_xprs_scores <- function(xprs) {
  mat <- if (is.list(xprs) && !is.null(xprs$locus_score_matrix)) xprs$locus_score_matrix else xprs
  mat <- as.matrix(mat)
  if (!nrow(mat) || !ncol(mat)) stop("xprs must contain a non-empty locus score matrix", call. = FALSE)
  full <- rowSums(mat)
  out <- matrix(full, nrow = nrow(mat), ncol = ncol(mat)) - mat
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

score_adjusted_residuals <- function(y,
                                     score,
                                     covariates = NULL,
                                     include_score2 = TRUE,
                                     sample_ids = NULL,
                                     tol = 1e-12) {
  y <- qtl2mr_as_named_numeric(y, "y")
  score <- qtl2mr_as_named_numeric(score, "score")
  ids <- intersect(names(y), names(score))
  if (!is.null(sample_ids)) ids <- intersect(ids, sample_ids)
  ids <- ids[is.finite(y[ids]) & is.finite(score[ids])]
  if (length(ids) < 4L) stop("Need at least four complete observations", call. = FALSE)
  C <- qtl2mr_as_addcovar(covariates, ids = ids, tol = tol)
  X <- cbind(intercept = 1, C, score = score[ids])
  if (isTRUE(include_score2)) X <- cbind(X, score2 = score[ids]^2)
  X <- qtl2mr_drop_depcols(X, add_intercept = FALSE, tol = tol)
  fit <- stats::lm.fit(X, y[ids])
  data.frame(
    sample_id = ids,
    y = as.numeric(y[ids]),
    score = as.numeric(score[ids]),
    fitted = as.numeric(fit$fitted.values),
    residual = as.numeric(fit$residuals),
    mean_model_used = if (isTRUE(include_score2)) "intercept + C + S_X + S_X^2" else "intercept + C + S_X",
    stringsAsFactors = FALSE
  )
}

score_heteroskedasticity_test <- function(residual_data,
                                          score = NULL,
                                          include_score2 = TRUE) {
  if (is.data.frame(residual_data)) {
    d <- residual_data
    if (!"residual" %in% names(d) || !"score" %in% names(d)) {
      stop("residual_data must contain residual and score columns", call. = FALSE)
    }
  } else {
    residual <- qtl2mr_as_named_numeric(residual_data, "residual")
    score <- qtl2mr_as_named_numeric(score, "score")
    ids <- intersect(names(residual), names(score))
    d <- data.frame(residual = residual[ids], score = score[ids])
  }
  d <- d[stats::complete.cases(d[, c("residual", "score"), drop = FALSE]), , drop = FALSE]
  d$resid2 <- d$residual^2
  d$score2 <- d$score^2
  full_formula <- if (isTRUE(include_score2)) resid2 ~ score + score2 else resid2 ~ score
  full <- stats::lm(full_formula, data = d)
  null <- stats::lm(resid2 ~ 1, data = d)
  av <- stats::anova(null, full)
  co <- summary(full)$coefficients
  data.frame(
    sample_size = nrow(d),
    include_score2 = isTRUE(include_score2),
    gamma1 = unname(co["score", "Estimate"]),
    gamma1_se = unname(co["score", "Std. Error"]),
    gamma2 = if (isTRUE(include_score2)) unname(co["score2", "Estimate"]) else NA_real_,
    gamma2_se = if (isTRUE(include_score2)) unname(co["score2", "Std. Error"]) else NA_real_,
    joint_test_statistic = unname(av$F[2]),
    joint_p = unname(av$`Pr(>F)`[2]),
    stringsAsFactors = FALSE
  )
}

score_bin_variance <- function(residual_data,
                               score = NULL,
                               n_bins = 5L,
                               B = 600L,
                               seed = NULL) {
  if (is.data.frame(residual_data)) {
    d <- residual_data
    if (!"residual" %in% names(d) || !"score" %in% names(d)) {
      stop("residual_data must contain residual and score columns", call. = FALSE)
    }
  } else {
    residual <- qtl2mr_as_named_numeric(residual_data, "residual")
    score <- qtl2mr_as_named_numeric(score, "score")
    ids <- intersect(names(residual), names(score))
    d <- data.frame(sample_id = ids, residual = residual[ids], score = score[ids])
  }
  d <- d[stats::complete.cases(d[, c("residual", "score"), drop = FALSE]), , drop = FALSE]
  d <- d[order(d$score), , drop = FALSE]
  d$score_bin <- pmin(n_bins, ceiling(seq_len(nrow(d)) / nrow(d) * n_bins))
  qtl2mr_with_seed(seed, {
    rows <- lapply(split(d, d$score_bin), function(dd) {
      v_obs <- stats::var(dd$residual, na.rm = TRUE)
      boots <- replicate(B, stats::var(sample(dd$residual, replace = TRUE), na.rm = TRUE))
      data.frame(
        score_bin = as.integer(dd$score_bin[[1L]]),
        score_bin_label = paste0("Q", as.integer(dd$score_bin[[1L]])),
        mean_score = mean(dd$score, na.rm = TRUE),
        n = nrow(dd),
        residual_variance = v_obs,
        ci_lower = stats::quantile(boots, 0.025, na.rm = TRUE, names = FALSE),
        ci_upper = stats::quantile(boots, 0.975, na.rm = TRUE, names = FALSE),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
}

augment_score_heteroskedasticity_tests <- function(tests,
                                                   bins,
                                                   nominal_alpha = 0.01,
                                                   fdr_alpha = 0.05,
                                                   spike_ratio = 1.5) {
  required_tests <- c("marker", "trait_role", "trait", "joint_p")
  required_bins <- c("marker", "trait_role", "trait", "score_bin", "residual_variance")
  miss_tests <- setdiff(required_tests, names(tests))
  miss_bins <- setdiff(required_bins, names(bins))
  if (length(miss_tests)) stop("tests is missing: ", paste(miss_tests, collapse = ", "), call. = FALSE)
  if (length(miss_bins)) stop("bins is missing: ", paste(miss_bins, collapse = ", "), call. = FALSE)

  tests <- as.data.frame(tests, stringsAsFactors = FALSE)
  bins <- as.data.frame(bins, stringsAsFactors = FALSE)
  tests$joint_q_bh <- NA_real_
  for (trait in unique(as.character(tests$trait))) {
    idx <- which(as.character(tests$trait) == trait & is.finite(tests$joint_p))
    if (length(idx)) tests$joint_q_bh[idx] <- stats::p.adjust(tests$joint_p[idx], method = "BH")
  }
  tests$flag_nominal_0.01 <- is.finite(tests$joint_p) & tests$joint_p < nominal_alpha
  tests$flag_fdr_0.05 <- is.finite(tests$joint_q_bh) & tests$joint_q_bh < fdr_alpha
  tests$variance_ratio_q5_q1 <- NA_real_
  tests$variance_ratio_max_min <- NA_real_
  tests$bin_variance_slope <- NA_real_
  tests$max_bin <- NA_integer_
  tests$max_to_second_bin_ratio <- NA_real_
  tests$shape_sensitive_qc_warning <- FALSE

  for (i in seq_len(nrow(tests))) {
    bb <- bins[
      as.character(bins$marker) == as.character(tests$marker[[i]]) &
        as.character(bins$trait_role) == as.character(tests$trait_role[[i]]) &
        as.character(bins$trait) == as.character(tests$trait[[i]]),
      ,
      drop = FALSE
    ]
    bb <- bb[order(bb$score_bin), , drop = FALSE]
    v <- as.numeric(bb$residual_variance)
    ok <- is.finite(v) & v >= 0
    if (!length(v) || sum(ok) < 2L) next
    tests$variance_ratio_q5_q1[[i]] <- v[[length(v)]] / max(v[[1L]], .Machine$double.eps)
    tests$variance_ratio_max_min[[i]] <- max(v[ok]) / max(min(v[ok]), .Machine$double.eps)
    fit <- stats::lm(v ~ bb$score_bin)
    tests$bin_variance_slope[[i]] <- unname(stats::coef(fit)[[2L]])
    max_i <- which.max(v)
    tests$max_bin[[i]] <- as.integer(bb$score_bin[[max_i]])
    remaining <- v[-max_i]
    if (length(remaining) && any(is.finite(remaining))) {
      tests$max_to_second_bin_ratio[[i]] <- v[[max_i]] / max(max(remaining, na.rm = TRUE), .Machine$double.eps)
    }
    tests$shape_sensitive_qc_warning[[i]] <- isTRUE(tests$flag_nominal_0.01[[i]] || tests$flag_fdr_0.05[[i]]) &&
      is.finite(tests$max_to_second_bin_ratio[[i]]) &&
      tests$max_to_second_bin_ratio[[i]] >= spike_ratio
  }
  tests
}

multilocus_xprs_heteroskedasticity <- function(traits,
                                               xprs,
                                               covariates = NULL,
                                               sample_ids = NULL,
                                               trait_roles = NULL,
                                               n_bins = 5L,
                                               B = 600L,
                                               seed = NULL,
                                               include_score2 = TRUE,
                                               n_loci = NA_integer_,
                                               score_label = "multi_locus_xprs") {
  if (!is.list(traits)) traits <- list(trait = traits)
  trait_names <- names(traits)
  if (is.null(trait_names) || any(!nzchar(trait_names))) trait_names <- paste0("trait", seq_along(traits))
  if (is.null(trait_roles)) trait_roles <- trait_names
  if (length(trait_roles) != length(traits)) stop("trait_roles must match traits length", call. = FALSE)
  xprs <- qtl2mr_standardized_named_score(qtl2mr_xprs_vector(xprs, sample_ids = sample_ids))

  test_rows <- vector("list", length(traits))
  bin_rows <- vector("list", length(traits))
  residual_rows <- vector("list", length(traits))
  for (i in seq_along(traits)) {
    y <- qtl2mr_as_named_numeric(traits[[i]], trait_names[[i]])
    ids <- intersect(names(y), names(xprs))
    if (!is.null(sample_ids)) ids <- intersect(as.character(sample_ids), ids)
    ids <- ids[is.finite(y[ids]) & is.finite(xprs[ids])]
    resid <- score_adjusted_residuals(
      y = y,
      score = xprs,
      covariates = covariates,
      include_score2 = include_score2,
      sample_ids = ids
    )
    ht <- score_heteroskedasticity_test(resid, include_score2 = include_score2)
    bins <- score_bin_variance(
      resid,
      n_bins = n_bins,
      B = B,
      seed = if (is.null(seed)) NULL else seed + i
    )
    test_rows[[i]] <- data.frame(
      marker = score_label,
      chr = "multi",
      locus_rank = NA_integer_,
      trait_role = trait_roles[[i]],
      trait = trait_names[[i]],
      score_type = "summed_exposure_xprs",
      n_loci = n_loci,
      score_mean = mean(resid$score, na.rm = TRUE),
      score_sd = stats::sd(resid$score, na.rm = TRUE),
      score_min = min(resid$score, na.rm = TRUE),
      score_max = max(resid$score, na.rm = TRUE),
      residual_sd = stats::sd(resid$residual, na.rm = TRUE),
      B = B,
      ht,
      stringsAsFactors = FALSE
    )
    bin_rows[[i]] <- cbind(
      data.frame(
        marker = score_label,
        chr = "multi",
        locus_rank = NA_integer_,
        trait_role = trait_roles[[i]],
        trait = trait_names[[i]],
        score_type = "summed_exposure_xprs",
        n_loci = n_loci,
        B = B,
        stringsAsFactors = FALSE
      ),
      bins
    )
    residual_rows[[i]] <- cbind(
      data.frame(
        marker = score_label,
        chr = "multi",
        locus_rank = NA_integer_,
        trait_role = trait_roles[[i]],
        trait = trait_names[[i]],
        score_type = "summed_exposure_xprs",
        n_loci = n_loci,
        stringsAsFactors = FALSE
      ),
      resid
    )
  }
  tests <- do.call(rbind, test_rows)
  bins <- do.call(rbind, bin_rows)
  residuals <- do.call(rbind, residual_rows)
  tests <- augment_score_heteroskedasticity_tests(tests, bins)
  list(tests = tests, bins = bins, residuals = residuals)
}

leave_one_locus_xprs_heteroskedasticity <- function(traits,
                                                    xprs,
                                                    covariates = NULL,
                                                    sample_ids = NULL,
                                                    trait_roles = NULL,
                                                    marker_info = NULL,
                                                    include_score2 = TRUE) {
  if (!is.list(traits)) traits <- list(trait = traits)
  trait_names <- names(traits)
  if (is.null(trait_names) || any(!nzchar(trait_names))) trait_names <- paste0("trait", seq_along(traits))
  if (is.null(trait_roles)) trait_roles <- trait_names
  if (length(trait_roles) != length(traits)) stop("trait_roles must match traits length", call. = FALSE)
  mat <- if (is.list(xprs) && !is.null(xprs$locus_score_matrix)) xprs$locus_score_matrix else xprs
  mat <- as.matrix(mat)
  loo <- leave_one_locus_xprs_scores(mat)
  markers <- colnames(loo)
  rows <- vector("list", length(markers) * length(traits))
  k <- 0L
  for (m in seq_along(markers)) {
    score <- qtl2mr_standardized_named_score(stats::setNames(loo[, m], rownames(loo)))
    marker <- markers[[m]]
    chr <- if (!is.null(marker_info) && "chr" %in% names(marker_info)) {
      hit <- match(marker, as.character(marker_info$marker))
      if (is.na(hit)) NA_character_ else as.character(marker_info$chr[[hit]])
    } else {
      NA_character_
    }
    rank <- if (!is.null(marker_info) && "locus_rank" %in% names(marker_info)) {
      hit <- match(marker, as.character(marker_info$marker))
      if (is.na(hit)) m else marker_info$locus_rank[[hit]]
    } else {
      m
    }
    for (i in seq_along(traits)) {
      y <- qtl2mr_as_named_numeric(traits[[i]], trait_names[[i]])
      ids <- intersect(names(y), names(score))
      if (!is.null(sample_ids)) ids <- intersect(as.character(sample_ids), ids)
      ids <- ids[is.finite(y[ids]) & is.finite(score[ids])]
      resid <- score_adjusted_residuals(
        y = y,
        score = score,
        covariates = covariates,
        include_score2 = include_score2,
        sample_ids = ids
      )
      ht <- score_heteroskedasticity_test(resid, include_score2 = include_score2)
      k <- k + 1L
      rows[[k]] <- data.frame(
        omitted_marker = marker,
        chr = chr,
        locus_rank = rank,
        trait_role = trait_roles[[i]],
        trait = trait_names[[i]],
        score_type = "leave_one_locus_out_xprs",
        sample_size = ht$sample_size,
        include_score2 = ht$include_score2,
        score_mean = mean(resid$score, na.rm = TRUE),
        score_sd = stats::sd(resid$score, na.rm = TRUE),
        score_min = min(resid$score, na.rm = TRUE),
        score_max = max(resid$score, na.rm = TRUE),
        residual_sd = stats::sd(resid$residual, na.rm = TRUE),
        gamma1 = ht$gamma1,
        gamma1_se = ht$gamma1_se,
        gamma2 = ht$gamma2,
        gamma2_se = ht$gamma2_se,
        joint_test_statistic = ht$joint_test_statistic,
        joint_p = ht$joint_p,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$joint_q_bh <- NA_real_
  for (trait in unique(as.character(out$trait))) {
    idx <- which(as.character(out$trait) == trait & is.finite(out$joint_p))
    if (length(idx)) out$joint_q_bh[idx] <- stats::p.adjust(out$joint_p[idx], method = "BH")
  }
  out
}

qtl2mr_xprs_localization_heteroskedasticity <- function(traits,
                                                        per_locus_scores,
                                                        covariates = NULL,
                                                        sample_ids = NULL,
                                                        trait_roles = NULL,
                                                        marker_info = NULL,
                                                        n_bins = 5L,
                                                        B = 600L,
                                                        seed = NULL,
                                                        include_score2 = TRUE) {
  if (!is.list(traits)) traits <- list(trait = traits)
  trait_names <- names(traits)
  if (is.null(trait_names) || any(!nzchar(trait_names))) trait_names <- paste0("trait", seq_along(traits))
  if (is.null(trait_roles)) trait_roles <- trait_names
  if (length(trait_roles) != length(traits)) stop("trait_roles must match traits length", call. = FALSE)
  mat <- as.matrix(per_locus_scores)
  if (is.null(rownames(mat)) || any(!nzchar(rownames(mat)))) {
    stop("per_locus_scores must have sample row names", call. = FALSE)
  }
  if (is.null(colnames(mat)) || any(!nzchar(colnames(mat)))) {
    stop("per_locus_scores must have marker column names", call. = FALSE)
  }

  test_rows <- list()
  bin_rows <- list()
  residual_rows <- list()
  row_i <- 0L
  bin_i <- 0L
  resid_i <- 0L
  for (m in seq_len(ncol(mat))) {
    marker <- colnames(mat)[[m]]
    score <- stats::setNames(as.numeric(mat[, m]), rownames(mat))
    chr <- if (!is.null(marker_info) && "chr" %in% names(marker_info)) {
      hit <- match(marker, as.character(marker_info$marker))
      if (is.na(hit)) NA_character_ else as.character(marker_info$chr[[hit]])
    } else {
      NA_character_
    }
    locus_rank <- if (!is.null(marker_info) && "locus_rank" %in% names(marker_info)) {
      hit <- match(marker, as.character(marker_info$marker))
      if (is.na(hit)) m else marker_info$locus_rank[[hit]]
    } else {
      m
    }
    for (j in seq_along(traits)) {
      y <- qtl2mr_as_named_numeric(traits[[j]], trait_names[[j]])
      ids <- intersect(names(y), names(score))
      if (!is.null(sample_ids)) ids <- intersect(as.character(sample_ids), ids)
      ids <- ids[is.finite(y[ids]) & is.finite(score[ids])]
      resid <- score_adjusted_residuals(
        y = y,
        score = score,
        covariates = covariates,
        include_score2 = include_score2,
        sample_ids = ids
      )
      ht <- score_heteroskedasticity_test(resid, include_score2 = include_score2)
      bins <- score_bin_variance(
        resid,
        n_bins = n_bins,
        B = B,
        seed = if (is.null(seed)) NULL else seed + m * 100L + j
      )
      row_i <- row_i + 1L
      test_rows[[row_i]] <- data.frame(
        marker = marker,
        chr = chr,
        locus_rank = locus_rank,
        trait_role = trait_roles[[j]],
        trait = trait_names[[j]],
        score_mean = mean(resid$score, na.rm = TRUE),
        score_sd = stats::sd(resid$score, na.rm = TRUE),
        score_min = min(resid$score, na.rm = TRUE),
        score_max = max(resid$score, na.rm = TRUE),
        residual_sd = stats::sd(resid$residual, na.rm = TRUE),
        B = B,
        ht,
        stringsAsFactors = FALSE
      )
      bin_i <- bin_i + 1L
      bin_rows[[bin_i]] <- cbind(
        data.frame(
          marker = marker,
          chr = chr,
          locus_rank = locus_rank,
          trait_role = trait_roles[[j]],
          trait = trait_names[[j]],
          B = B,
          stringsAsFactors = FALSE
        ),
        bins
      )
      resid_i <- resid_i + 1L
      residual_rows[[resid_i]] <- cbind(
        data.frame(
          marker = marker,
          chr = chr,
          locus_rank = locus_rank,
          trait_role = trait_roles[[j]],
          trait = trait_names[[j]],
          stringsAsFactors = FALSE
        ),
        resid
      )
    }
  }
  tests <- do.call(rbind, test_rows)
  bins <- do.call(rbind, bin_rows)
  residuals <- do.call(rbind, residual_rows)
  tests <- augment_score_heteroskedasticity_tests(tests, bins)
  list(tests = tests, bins = bins, residuals = residuals)
}

qtl2mr_new_xprs_heteroskedasticity_diagnostic <- function(global,
                                                          xprs = NULL,
                                                          loo_tests = NULL,
                                                          localization = NULL,
                                                          marker_info = NULL,
                                                          settings = list()) {
  out <- list(
    global_tests = global$tests,
    global_bins = global$bins,
    global_residuals = global$residuals,
    xprs_scores = if (is.list(xprs) && !is.null(xprs$sample_scores)) xprs$sample_scores else NULL,
    xprs_locus_score_matrix = if (is.list(xprs) && !is.null(xprs$locus_score_matrix)) xprs$locus_score_matrix else NULL,
    loo_tests = loo_tests,
    localization_tests = if (!is.null(localization)) localization$tests else NULL,
    localization_bins = if (!is.null(localization)) localization$bins else NULL,
    localization_residuals = if (!is.null(localization)) localization$residuals else NULL,
    marker_info = marker_info,
    settings = settings
  )
  class(out) <- c("qtl2mr_xprs_heteroskedasticity_diagnostic", "list")
  out
}

print.qtl2mr_xprs_heteroskedasticity_diagnostic <- function(x, ...) {
  cat("qtl2mr multi-locus xPRS heteroskedasticity diagnostic\n")
  if (!is.null(x$settings$n_loci) && is.finite(x$settings$n_loci)) {
    cat("Loci:", x$settings$n_loci, "\n")
  }
  if (!is.null(x$global_tests) && nrow(x$global_tests)) {
    cat("Primary global tests:", nrow(x$global_tests), "trait(s)\n")
    cols <- intersect(c("trait", "sample_size", "n_loci", "joint_p", "variance_ratio_q5_q1", "variance_ratio_max_min"), names(x$global_tests))
    print(x$global_tests[, cols, drop = FALSE], row.names = FALSE)
  }
  if (!is.null(x$loo_tests) && nrow(x$loo_tests)) {
    cat("Optional leave-one-locus-out tests:", nrow(x$loo_tests), "rows\n")
  }
  if (!is.null(x$localization_tests) && nrow(x$localization_tests)) {
    cat("Optional per-locus localization tests:", nrow(x$localization_tests), "rows\n")
  }
  invisible(x)
}

fit_xprs_heteroskedasticity_diagnostic <- function(traits,
                                                   xprs,
                                                   covariates = NULL,
                                                   sample_ids = NULL,
                                                   trait_roles = NULL,
                                                   per_locus_scores = NULL,
                                                   marker_info = NULL,
                                                   include_loo = TRUE,
                                                   include_localization = !is.null(per_locus_scores),
                                                   n_bins = 5L,
                                                   B = 600L,
                                                   seed = NULL,
                                                   include_score2 = TRUE,
                                                   n_loci = NA_integer_) {
  global <- multilocus_xprs_heteroskedasticity(
    traits = traits,
    xprs = xprs,
    covariates = covariates,
    sample_ids = sample_ids,
    trait_roles = trait_roles,
    n_bins = n_bins,
    B = B,
    seed = seed,
    include_score2 = include_score2,
    n_loci = n_loci
  )
  loo_tests <- NULL
  if (isTRUE(include_loo)) {
    mat <- if (is.list(xprs) && !is.null(xprs$locus_score_matrix)) xprs$locus_score_matrix else NULL
    if (!is.null(mat)) {
      loo_tests <- leave_one_locus_xprs_heteroskedasticity(
        traits = traits,
        xprs = xprs,
        covariates = covariates,
        sample_ids = sample_ids,
        trait_roles = trait_roles,
        marker_info = marker_info,
        include_score2 = include_score2
      )
    }
  }
  localization <- NULL
  if (isTRUE(include_localization) && !is.null(per_locus_scores)) {
    localization <- qtl2mr_xprs_localization_heteroskedasticity(
      traits = traits,
      per_locus_scores = per_locus_scores,
      covariates = covariates,
      sample_ids = sample_ids,
      trait_roles = trait_roles,
      marker_info = marker_info,
      n_bins = n_bins,
      B = B,
      seed = seed,
      include_score2 = include_score2
    )
  }
  qtl2mr_new_xprs_heteroskedasticity_diagnostic(
    global = global,
    xprs = xprs,
    loo_tests = loo_tests,
    localization = localization,
    marker_info = marker_info,
    settings = list(
      n_bins = n_bins,
      B = B,
      seed = seed,
      include_score2 = include_score2,
      n_loci = n_loci,
      include_loo = isTRUE(include_loo),
      include_localization = isTRUE(include_localization)
    )
  )
}
