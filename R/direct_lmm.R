qtl2mr_select_marker_genoprobs <- function(genoprobs, marker = NULL) {
  if (!is.list(genoprobs)) stop("genoprobs must be a qtl2-style list of arrays", call. = FALSE)
  chr_names <- names(genoprobs)
  if (is.null(chr_names)) chr_names <- as.character(seq_along(genoprobs))

  hits <- list()
  for (i in seq_along(genoprobs)) {
    arr <- genoprobs[[i]]
    if (length(dim(arr)) != 3L) next
    markers <- dimnames(arr)[[3L]]
    if (is.null(markers)) markers <- as.character(seq_len(dim(arr)[3L]))
    if (is.null(marker)) {
      if (dim(arr)[3L] == 1L) hits[[length(hits) + 1L]] <- c(chr = chr_names[[i]], index = 1L)
    } else {
      idx <- match(marker, markers)
      if (!is.na(idx)) hits[[length(hits) + 1L]] <- c(chr = chr_names[[i]], index = idx)
    }
  }
  if (length(hits) != 1L) {
    if (is.null(marker)) {
      stop("marker is required unless genoprobs contains exactly one marker", call. = FALSE)
    }
    stop("Expected marker '", marker, "' in exactly one chromosome; found ", length(hits), call. = FALSE)
  }
  chr <- unname(hits[[1L]][["chr"]])
  idx <- as.integer(hits[[1L]][["index"]])
  arr <- genoprobs[[match(chr, chr_names)]][, , idx, drop = FALSE]
  marker_name <- dimnames(arr)[[3L]][[1L]]
  list(genoprobs = arr, marker = marker_name, chr = chr)
}

qtl2mr_genoprobs_design <- function(gp_raw, contrasts) {
  if (length(dim(gp_raw)) != 3L || dim(gp_raw)[3L] != 1L) {
    stop("Selected genoprobs must be an n x founder x 1 array", call. = FALSE)
  }
  G <- gp_raw[, , 1L, drop = TRUE]
  if (is.null(dim(G))) G <- matrix(G, ncol = dim(gp_raw)[2L])
  allele_names <- colnames(G)
  if (is.null(allele_names)) allele_names <- rownames(contrasts)
  if (is.null(allele_names) || any(!allele_names %in% rownames(contrasts))) {
    stop("Founder allele names must match contrast matrix row names", call. = FALSE)
  }
  contrasts <- contrasts[allele_names, , drop = FALSE]
  Z <- G %*% contrasts
  rownames(Z) <- rownames(G)
  colnames(Z) <- colnames(contrasts)
  Z
}

qtl2mr_common_lmm_ids <- function(gp_raw, pheno, addcovar, kinship = NULL, sample_ids = NULL) {
  gp_ids <- dimnames(gp_raw)[[1L]]
  if (is.null(gp_ids)) stop("genoprobs arrays must have sample row names", call. = FALSE)
  base_ids <- if (is.null(sample_ids)) gp_ids else as.character(sample_ids)
  ids <- base_ids[base_ids %in% gp_ids & base_ids %in% names(pheno)]
  if (!is.null(addcovar) && nrow(addcovar)) ids <- ids[ids %in% rownames(addcovar)]
  if (!is.null(kinship)) ids <- ids[ids %in% rownames(as.matrix(kinship))]
  if (!length(ids)) stop("No common sample IDs across genoprobs, phenotype, covariates, and kinship", call. = FALSE)

  G <- gp_raw[ids, , 1L, drop = TRUE]
  if (is.null(dim(G))) G <- matrix(G, nrow = length(ids), dimnames = list(ids, NULL))
  cov_ok <- rep(TRUE, length(ids))
  if (!is.null(addcovar) && nrow(addcovar) && ncol(addcovar)) {
    cov_ok <- qtl2mr_complete_rows(addcovar[ids, , drop = FALSE])
  }
  keep <- is.finite(pheno[ids]) & cov_ok & apply(G, 1L, function(z) all(is.finite(z)))
  ids <- ids[keep]
  if (length(ids) <= 8L) stop("Need more than 8 complete samples; got ", length(ids), call. = FALSE)
  ids
}

estimate_founder_lmm_hsq <- function(pheno, kinship, addcovar = NULL, reml = TRUE) {
  qtl2mr_estimate_lmm_hsq(pheno, kinship, addcovar = addcovar, reml = reml)
}

qtl2mr_lmm_transform <- function(y, X, kinship = NULL, hsq = NULL, reml = TRUE, addcovar_for_hsq = NULL) {
  ids <- names(y)
  if (!is.null(hsq)) {
    validated_hsq <- qtl2mr_validate_hsq(hsq)
    if (is.null(kinship) && validated_hsq != 0) {
      stop("Nonzero hsq requires kinship; it cannot be used in an OLS fit", call. = FALSE)
    }
  }
  if (is.null(kinship)) {
    raw_transform <- diag(1, nrow = length(y), ncol = length(y))
    dimnames(raw_transform) <- list(ids, ids)
    star_rotation <- diag(1, nrow = length(y), ncol = length(y))
    dimnames(star_rotation) <- list(ids, ids)
    return(list(
      ystar = as.numeric(y),
      Xstar = X,
      raw_transform = raw_transform,
      star_rotation = star_rotation,
      inverse_residual = function(resid_star) stats::setNames(as.numeric(resid_star), ids),
      hsq = NA_real_,
      hsq_status = "not_applicable",
      lmm_backend = "ordinary_least_squares",
      kinship_decomp = NULL
    ))
  }

  kinship0 <- as.matrix(kinship)[ids, ids, drop = FALSE]
  if (is.null(hsq)) hsq <- estimate_founder_lmm_hsq(y, kinship0, addcovar_for_hsq, reml = reml)
  hsq_status <- attr(hsq, "status") %||% "fixed"
  hsq <- qtl2mr_validate_hsq(hsq)

  kinship_decomp <- qtl2mr_decompose_kinship(kinship0, sample_ids = ids)
  eigenvec <- kinship_decomp$vectors
  wts <- 1 / sqrt(hsq * kinship_decomp$values + (1 - hsq))
  Xstar <- sweep(eigenvec %*% X, 1L, wts, "*")
  ystar <- as.numeric(wts * (eigenvec %*% as.numeric(y)))
  rownames(Xstar) <- ids
  colnames(Xstar) <- colnames(X)
  raw_transform <- diag(wts, nrow = length(wts)) %*% eigenvec
  colnames(raw_transform) <- ids
  star_rotation <- eigenvec
  colnames(star_rotation) <- ids
  inverse_residual <- function(resid_star) {
    out <- as.numeric(t(eigenvec) %*% (as.numeric(resid_star) / wts))
    stats::setNames(out, ids)
  }
  list(
    ystar = ystar,
    Xstar = Xstar,
    raw_transform = raw_transform,
    star_rotation = star_rotation,
    inverse_residual = inverse_residual,
    hsq = hsq,
    hsq_status = hsq_status,
    lmm_backend = "internal_base_lmm",
    kinship_decomp = kinship_decomp
  )
}

fit_founder_lmm_block <- function(genoprobs,
                                  pheno,
                                  kinship = NULL,
                                  addcovar = NULL,
                                  contrasts = NULL,
                                  hsq = NULL,
                                  marker = NULL,
                                  sample_ids = NULL,
                                  standardize_pheno = FALSE,
                                  reml = TRUE,
                                  fc_cols = paste0("FC", 1:7),
                                  tol = 1e-12) {
  selected <- qtl2mr_select_marker_genoprobs(genoprobs, marker = marker)
  gp_raw <- selected$genoprobs
  pheno <- qtl2mr_as_named_numeric(pheno)

  allele_names <- dimnames(gp_raw)[[2L]]
  if (is.null(contrasts)) contrasts <- make_founder_contrasts(allele_names)
  if (!is.matrix(contrasts)) contrasts <- as.matrix(contrasts)
  validate_founder_contrasts(contrasts, rownames(contrasts))

  ids0 <- qtl2mr_common_lmm_ids(gp_raw, pheno, addcovar, kinship = kinship, sample_ids = sample_ids)
  addcovar0 <- qtl2mr_as_addcovar(addcovar, ids = ids0, tol = tol)
  y <- pheno[ids0]
  if (isTRUE(standardize_pheno)) {
    s <- stats::sd(y)
    if (!is.finite(s) || s <= 0) stop("Cannot standardize phenotype with non-positive SD", call. = FALSE)
    y <- (y - mean(y)) / s
  }
  names(y) <- ids0

  gp_raw <- gp_raw[ids0, , , drop = FALSE]
  qtl_design <- qtl2mr_genoprobs_design(gp_raw, contrasts)
  X <- cbind(qtl_design, addcovar0)
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  colnames(X)[seq_len(ncol(qtl_design))] <- colnames(qtl_design)

  missing_fc <- setdiff(fc_cols, colnames(X))
  if (length(missing_fc)) stop("Missing founder contrast columns: ", paste(missing_fc, collapse = ", "), call. = FALSE)

  tr <- qtl2mr_lmm_transform(
    y = y,
    X = X,
    kinship = kinship,
    hsq = hsq,
    reml = reml,
    addcovar_for_hsq = addcovar0
  )
  Xstar <- tr$Xstar
  ystar <- tr$ystar

  q <- qr(Xstar, tol = tol)
  beta_full <- as.numeric(qr.coef(q, ystar))
  names(beta_full) <- colnames(Xstar)
  fitted_star <- as.numeric(qr.fitted(q, ystar))
  residual_star <- ystar - fitted_star
  df <- length(ystar) - q$rank
  if (df <= 0L) stop("Non-positive residual degrees of freedom", call. = FALSE)
  sigma2 <- as.numeric(crossprod(residual_star)) / df
  sigma <- sqrt(sigma2)

  coef_map_star <- qr.coef(q, diag(nrow(Xstar)))
  rownames(coef_map_star) <- colnames(Xstar)
  colnames(coef_map_star) <- ids0
  if (any(!is.finite(coef_map_star[fc_cols, , drop = FALSE])) ||
      any(!is.finite(beta_full[fc_cols]))) {
    stop("Founder contrast coefficients are not estimable at marker ", selected$marker, call. = FALSE)
  }

  C_full <- coef_map_star %*% t(coef_map_star)
  dimnames(C_full) <- list(colnames(Xstar), colnames(Xstar))
  vcov_full <- sigma2 * C_full
  vcov_fc <- qtl2mr_sym(vcov_full[fc_cols, fc_cols, drop = FALSE])
  beta_fc <- stats::setNames(as.numeric(beta_full[fc_cols]), fc_cols)
  se_fc <- stats::setNames(sqrt(diag(vcov_fc)), fc_cols)

  influence_raw <- coef_map_star[fc_cols, , drop = FALSE] %*% tr$raw_transform
  colnames(influence_raw) <- ids0
  rownames(influence_raw) <- fc_cols
  influence_star <- coef_map_star[fc_cols, , drop = FALSE] %*% tr$star_rotation
  colnames(influence_star) <- ids0
  rownames(influence_star) <- fc_cols
  residual_raw <- tr$inverse_residual(residual_star)
  residual_star_symmetric <- as.numeric(t(tr$star_rotation) %*% as.numeric(residual_star))

  list(
    marker = selected$marker,
    chr = selected$chr,
    founders = allele_names,
    contrasts = contrasts,
    founder_H = contrasts[, fc_cols, drop = FALSE],
    qtl_design = qtl_design,
    beta_fc = beta_fc,
    se_fc = se_fc,
    vcov_fc = vcov_fc,
    beta_full = beta_full,
    vcov_full = vcov_full,
    sigma = sigma,
    sigma2 = sigma2,
    residual_raw = residual_raw,
    residual_star = stats::setNames(residual_star_symmetric, ids0),
    residual_star_fit = stats::setNames(as.numeric(residual_star), ids0),
    influence_raw = influence_raw,
    influence_star = influence_star,
    influence_star_fit = coef_map_star[fc_cols, , drop = FALSE],
    hsq = tr$hsq,
    hsq_status = tr$hsq_status,
    lmm_backend = tr$lmm_backend,
    n = length(ids0),
    df = df,
    rank = q$rank,
    sample_ids = ids0,
    design_colnames = colnames(Xstar),
    fc_cols = fc_cols,
    standardize_pheno = isTRUE(standardize_pheno)
  )
}

qtl2mr_lmm_cross_covariance <- function(fit_x,
                                        fit_y,
                                        fc_cols,
                                        space = c("star", "raw", "zero"),
                                        O_xy = NULL,
                                        overlap_ids = NULL) {
  space <- match.arg(space)
  residual_slot <- if (identical(space, "star")) "residual_star" else "residual_raw"
  influence_slot <- if (identical(space, "star")) "influence_star" else "influence_raw"
  sample_ids_x <- if (identical(space, "zero")) fit_x$sample_ids else names(fit_x[[residual_slot]])
  sample_ids_y <- if (identical(space, "zero")) fit_y$sample_ids else names(fit_y[[residual_slot]])
  if (is.null(O_xy)) {
    O_xy <- make_overlap_matrix(sample_ids_x, sample_ids_y, overlap_ids = overlap_ids)
  } else {
    O_xy <- as.matrix(O_xy)
    if (nrow(O_xy) != length(sample_ids_x) || ncol(O_xy) != length(sample_ids_y)) {
      stop("O_xy dimensions do not match exposure/outcome sample IDs", call. = FALSE)
    }
    if (is.null(rownames(O_xy))) rownames(O_xy) <- sample_ids_x
    if (is.null(colnames(O_xy))) colnames(O_xy) <- sample_ids_y
  }
  n_overlap <- sum(O_xy != 0)
  O_summary <- list(
    n_x = nrow(O_xy),
    n_y = ncol(O_xy),
    n_overlap = n_overlap,
    exposure_overlap_fraction = if (nrow(O_xy)) n_overlap / nrow(O_xy) else NA_real_,
    outcome_overlap_fraction = if (ncol(O_xy)) n_overlap / ncol(O_xy) else NA_real_,
    jaccard_overlap_fraction = if ((nrow(O_xy) + ncol(O_xy) - n_overlap) > 0) {
      n_overlap / (nrow(O_xy) + ncol(O_xy) - n_overlap)
    } else {
      NA_real_
    },
    complete_overlap = identical(rownames(O_xy), colnames(O_xy)) &&
      n_overlap == nrow(O_xy) && all(diag(O_xy) == 1),
    identical_order = identical(rownames(O_xy), colnames(O_xy))
  )
  if (identical(space, "zero")) {
    Sigma_xy <- matrix(0, nrow = length(fc_cols), ncol = length(fc_cols),
                       dimnames = list(fc_cols, fc_cols))
    return(list(Sigma_xy = Sigma_xy, rho_xy = 0, rho_star = 0,
                n_overlap = n_overlap, O_xy_summary = O_summary))
  }
  if (n_overlap == 0L) {
    Sigma_xy <- matrix(0, nrow = length(fc_cols), ncol = length(fc_cols),
                       dimnames = list(fc_cols, fc_cols))
    return(list(Sigma_xy = Sigma_xy, rho_xy = NA_real_, rho_star = NA_real_,
                n_overlap = n_overlap, O_xy_summary = O_summary))
  }

  rho_xy <- estimate_rho_star(fit_x[[residual_slot]], fit_y[[residual_slot]], O_xy)
  if (!is.finite(rho_xy)) stop("Residual cross-trait correlation is not finite in ", space, " space", call. = FALSE)
  if (!all(rownames(O_xy) %in% colnames(fit_x[[influence_slot]]))) {
    stop("Exposure influence matrix does not contain all O_xy row IDs", call. = FALSE)
  }
  if (!all(colnames(O_xy) %in% colnames(fit_y[[influence_slot]]))) {
    stop("Outcome influence matrix does not contain all O_xy column IDs", call. = FALSE)
  }
  Ax <- fit_x[[influence_slot]][, rownames(O_xy), drop = FALSE]
  Ay <- fit_y[[influence_slot]][, colnames(O_xy), drop = FALSE]
  Sigma_xy <- rho_xy * fit_x$sigma * fit_y$sigma * (Ax %*% O_xy %*% t(Ay))
  rownames(Sigma_xy) <- fc_cols
  colnames(Sigma_xy) <- fc_cols
  list(Sigma_xy = Sigma_xy, rho_xy = rho_xy, rho_star = rho_xy,
       n_overlap = n_overlap, O_xy_summary = O_summary)
}

fit_mr_locus_block <- function(genoprobs,
                               exposure,
                               outcome,
                               kinship = NULL,
                               addcovar = NULL,
                               contrasts = NULL,
                               hsq_x = NULL,
                               hsq_y = NULL,
                               marker = NULL,
                               sample_ids = NULL,
                               standardize_phenotypes = FALSE,
                               reml = TRUE,
                               fc_cols = paste0("FC", 1:7),
                               ridge = 1e-10,
                               cross_cov_model = c("star", "raw", "zero"),
                               enforce_psd = TRUE,
                               stop_on_invalid_psd = TRUE,
                               tol = 1e-12) {
  cross_cov_model <- match.arg(cross_cov_model)
  selected <- qtl2mr_select_marker_genoprobs(genoprobs, marker = marker)
  gp_raw <- selected$genoprobs
  x <- qtl2mr_as_named_numeric(exposure, "exposure")
  y <- qtl2mr_as_named_numeric(outcome, "outcome")
  joint_ids <- qtl2mr_common_lmm_ids(
    gp_raw = gp_raw,
    pheno = x[names(x) %in% names(y)],
    addcovar = addcovar,
    kinship = kinship,
    sample_ids = sample_ids
  )
  joint_ids <- joint_ids[joint_ids %in% names(y)]
  if (!is.null(addcovar) && nrow(as.matrix(addcovar))) {
    cov0 <- qtl2mr_as_addcovar(addcovar, ids = joint_ids, tol = tol)
    joint_ids <- joint_ids[stats::complete.cases(x[joint_ids], y[joint_ids], cov0)]
  } else {
    joint_ids <- joint_ids[is.finite(x[joint_ids]) & is.finite(y[joint_ids])]
  }
  if (length(joint_ids) <= 8L) stop("Need more than 8 complete one-sample observations", call. = FALSE)

  fit_x <- fit_founder_lmm_block(
    genoprobs = genoprobs,
    pheno = x,
    kinship = kinship,
    addcovar = addcovar,
    contrasts = contrasts,
    hsq = hsq_x,
    marker = selected$marker,
    sample_ids = joint_ids,
    standardize_pheno = standardize_phenotypes,
    reml = reml,
    fc_cols = fc_cols,
    tol = tol
  )
  fit_y <- fit_founder_lmm_block(
    genoprobs = genoprobs,
    pheno = y,
    kinship = kinship,
    addcovar = addcovar,
    contrasts = contrasts,
    hsq = hsq_y,
    marker = selected$marker,
    sample_ids = joint_ids,
    standardize_pheno = standardize_phenotypes,
    reml = reml,
    fc_cols = fc_cols,
    tol = tol
  )

  Sigma_x <- fit_x$vcov_fc
  Sigma_y <- fit_y$vcov_fc
  O_xy <- make_overlap_matrix(fit_x$sample_ids, fit_y$sample_ids)
  cross_raw <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "raw", O_xy = O_xy)
  cross_star <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "star", O_xy = O_xy)
  cross_zero <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "zero", O_xy = O_xy)
  selected_cross <- switch(
    cross_cov_model,
    raw = cross_raw,
    star = cross_star,
    zero = cross_zero
  )

  psd_raw <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_raw$Sigma_xy, ridge = ridge)
  psd_star <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_star$Sigma_xy, ridge = ridge)
  psd_zero <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_zero$Sigma_xy, ridge = ridge)
  psd <- switch(
    cross_cov_model,
    raw = psd_raw,
    star = psd_star,
    zero = psd_zero
  )

  Sigma_x_raw <- Sigma_x
  Sigma_y_raw <- Sigma_y
  Sigma_xy_raw <- selected_cross$Sigma_xy
  Sigma_xy <- selected_cross$Sigma_xy
  if (isTRUE(enforce_psd)) {
    if (!isTRUE(psd$psd_valid) && isTRUE(stop_on_invalid_psd)) {
      stop(
        "Selected ", cross_cov_model, " Sigma_xy fails joint PSD validation at marker ",
        selected$marker, " (status=", psd$status,
        ", max singular value=", signif(psd$max_cross_singular_value, 6), ")",
        call. = FALSE
      )
    }
    Sigma_x <- psd$Sigma_x
    Sigma_y <- psd$Sigma_y
    Sigma_xy <- psd$Sigma_xy
  }

  list(
    marker = selected$marker,
    chr = selected$chr,
    beta_x = fit_x$beta_fc,
    beta_y = fit_y$beta_fc,
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    Sigma_xy = Sigma_xy,
    Sigma_x_raw = Sigma_x_raw,
    Sigma_y_raw = Sigma_y_raw,
    Sigma_xy_raw = Sigma_xy_raw,
    Sigma_xy_raw_space = cross_raw$Sigma_xy,
    Sigma_xy_star_space = cross_star$Sigma_xy,
    Sigma_xy_zero = cross_zero$Sigma_xy,
    rho_xy = selected_cross$rho_xy,
    rho_star = selected_cross$rho_star,
    rho_xy_raw = cross_raw$rho_xy,
    rho_xy_star = cross_star$rho_xy,
    n_overlap = selected_cross$n_overlap,
    n_overlap_raw = cross_raw$n_overlap,
    n_overlap_star = cross_star$n_overlap,
    O_xy_summary = selected_cross$O_xy_summary,
    fit_x = fit_x,
    fit_y = fit_y,
    founders = fit_x$founders,
    contrasts = fit_x$contrasts,
    founder_H = fit_x$founder_H,
    hsq_x = fit_x$hsq,
    hsq_y = fit_y$hsq,
    cross_cov_model = cross_cov_model,
    psd_adjusted = isTRUE(enforce_psd) && isTRUE(psd$adjusted),
    joint_psd_candidate_adjusted = isTRUE(psd$adjusted),
    joint_psd_valid = isTRUE(psd$psd_valid),
    joint_psd_status = psd$status,
    joint_psd_large_violation = isTRUE(psd$large_violation),
    joint_min_eigenvalue = psd$min_eigenvalue,
    joint_max_cross_singular_value = psd$max_cross_singular_value %||% NA_real_,
    raw_space_psd_adjusted = isTRUE(psd_raw$adjusted),
    raw_space_psd_valid = isTRUE(psd_raw$psd_valid),
    raw_space_psd_status = psd_raw$status,
    raw_space_psd_large_violation = isTRUE(psd_raw$large_violation),
    raw_space_joint_min_eigenvalue = psd_raw$min_eigenvalue,
    raw_space_max_cross_singular_value = psd_raw$max_cross_singular_value %||% NA_real_,
    star_space_psd_adjusted = isTRUE(psd_star$adjusted),
    star_space_psd_valid = isTRUE(psd_star$psd_valid),
    star_space_psd_status = psd_star$status,
    star_space_psd_large_violation = isTRUE(psd_star$large_violation),
    star_space_joint_min_eigenvalue = psd_star$min_eigenvalue,
    star_space_max_cross_singular_value = psd_star$max_cross_singular_value %||% NA_real_,
    standardize_phenotypes = isTRUE(standardize_phenotypes)
  )
}

compare_founder_lmm_to_scan1coef <- function(genoprobs,
                                             pheno,
                                             kinship = NULL,
                                             addcovar = NULL,
                                             contrasts = NULL,
                                             hsq = NULL,
                                             marker = NULL,
                                             sample_ids = NULL,
                                             standardize_pheno = FALSE,
                                             reml = TRUE,
                                             fc_cols = paste0("FC", 1:7),
                                             tol = 1e-12) {
  qtl2mr_require_qtl2()
  fit <- fit_founder_lmm_block(
    genoprobs = genoprobs,
    pheno = pheno,
    kinship = kinship,
    addcovar = addcovar,
    contrasts = contrasts,
    hsq = hsq,
    marker = marker,
    sample_ids = sample_ids,
    standardize_pheno = standardize_pheno,
    reml = reml,
    fc_cols = fc_cols,
    tol = tol
  )
  y <- qtl2mr_as_named_numeric(pheno)[fit$sample_ids]
  if (isTRUE(standardize_pheno)) y <- (y - mean(y)) / stats::sd(y)
  addcovar0 <- qtl2mr_as_addcovar(addcovar, ids = fit$sample_ids, tol = tol)
  if (!ncol(addcovar0)) addcovar0 <- NULL
  selected <- qtl2mr_select_marker_genoprobs(genoprobs, marker = fit$marker)
  gp_one <- list(selected$genoprobs)
  names(gp_one) <- selected$chr
  class(gp_one) <- class(genoprobs)
  coef_obj <- qtl2::scan1coef(
    genoprobs = gp_one,
    pheno = y,
    kinship = if (is.null(kinship)) NULL else as.matrix(kinship)[fit$sample_ids, fit$sample_ids, drop = FALSE],
    addcovar = addcovar0,
    contrasts = contrasts %||% make_founder_contrasts(dimnames(selected$genoprobs)[[2L]]),
    se = TRUE,
    hsq = if (is.na(fit$hsq)) NULL else fit$hsq,
    reml = reml
  )
  se_obj <- attr(coef_obj, "SE")
  qtl2_beta <- stats::setNames(as.numeric(coef_obj[1L, fc_cols, drop = TRUE]), fc_cols)
  qtl2_se <- stats::setNames(as.numeric(se_obj[1L, fc_cols, drop = TRUE]), fc_cols)
  list(
    direct = fit,
    qtl2_beta = qtl2_beta,
    qtl2_se = qtl2_se,
    max_abs_beta_diff = max(abs(fit$beta_fc - qtl2_beta)),
    max_abs_se_diff = max(abs(fit$se_fc - qtl2_se))
  )
}
