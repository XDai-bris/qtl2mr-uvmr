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

