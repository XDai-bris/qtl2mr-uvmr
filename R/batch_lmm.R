qtl2mr_chr_option <- function(x, chr) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.data.frame(x)) return(x[[as.character(chr)]] %||% x[[1L]])
  x
}

# Exact key of an LMM context: every input that can change the context or the
# warnings it signals. The kinship enters through the fingerprint of its
# sample-set submatrix, and covariate_source holds the covariates as supplied,
# whose dependent columns are named in warnings.
qtl2mr_lmm_context_key <- function(trait_name,
                                   chr,
                                   sample_ids,
                                   covariates,
                                   standardize_pheno,
                                   hsq,
                                   kinship,
                                   phenotype = NULL,
                                   reml = TRUE,
                                   tol = 1e-12,
                                   covariate_source = NULL) {
  qtl2mr_fingerprint_object(list(schema = "lmm-context-v062-1", trait_name = trait_name,
    chr = chr, sample_ids = sample_ids, covariates = covariates,
    phenotype = phenotype, standardize_pheno = standardize_pheno,
    hsq = hsq, kinship = qtl2mr_kinship_fingerprint(kinship, sample_ids),
    reml = reml, tol = tol, covariate_source = covariate_source))
}

# A cache scope counts hits and misses as one batch run (or exposure-group
# worker) sees them: a key is a miss the first time the scope asks for it.
# The values themselves live in the session registry, shared by every scope.
qtl2mr_context_cache_new <- function() {
  cache <- new.env(parent = emptyenv())
  cache$.stats <- list(hits = 0L, misses = 0L)
  cache
}

qtl2mr_context_cache_stats <- function(cache) {
  if (is.null(cache) || !is.environment(cache)) {
    return(list(hits = 0L, misses = 0L, size = 0L))
  }
  stats <- cache$.stats %||% list(hits = 0L, misses = 0L)
  list(
    hits = as.integer(stats$hits %||% 0L),
    misses = as.integer(stats$misses %||% 0L),
    size = length(setdiff(ls(cache, all.names = TRUE), ".stats"))
  )
}

qtl2mr_context_cache_bump <- function(cache, field) {
  if (is.null(cache) || !is.environment(cache)) return(invisible(NULL))
  stats <- cache$.stats %||% list(hits = 0L, misses = 0L)
  stats[[field]] <- as.integer(stats[[field]] %||% 0L) + 1L
  cache$.stats <- stats
  invisible(NULL)
}

qtl2mr_cache_seen <- function(cache, key) {
  is.environment(cache) && exists(key, envir = cache, inherits = FALSE)
}

qtl2mr_cache_mark <- function(cache, key) {
  if (is.environment(cache)) assign(key, TRUE, envir = cache)
  invisible(key)
}

qtl2mr_context_lookup_key <- function(pheno, kinship, addcovar, hsq, sample_ids, standardize_pheno,
                                      reml, trait_name, chr, tol) {
  pheno0 <- qtl2mr_as_named_numeric(pheno, trait_name %||% "phenotype")
  ids <- if (is.null(sample_ids)) names(pheno0) else as.character(sample_ids)
  ids <- ids[ids %in% names(pheno0)]
  if (!is.null(addcovar) && nrow(as.matrix(addcovar))) {
    ids <- ids[ids %in% rownames(as.matrix(addcovar))]
  }
  if (!is.null(kinship)) ids <- ids[ids %in% rownames(as.matrix(kinship))]
  addcovar_key <- qtl2mr_as_addcovar(addcovar, ids = ids, tol = tol)
  qtl2mr_lmm_context_key(
    trait_name = trait_name,
    chr = chr,
    sample_ids = ids,
    covariates = addcovar_key,
    standardize_pheno = standardize_pheno,
    hsq = hsq,
    kinship = kinship,
    phenotype = pheno0[ids], reml = reml, tol = tol,
    covariate_source = if (is.null(addcovar)) NULL else addcovar[ids, , drop = FALSE]
  )
}

# LMM context of one trait on one chromosome. With a cache scope this behaves
# as the batch cache always did (the key is computed, warnings included, and
# hits and misses are counted); without one, as a direct prepare_lmm_context()
# call. Either way an identical context computed earlier in the session is
# reused, and the warnings its computation signalled are signalled again
# wherever the computation would have run. A reused context is returned
# without its n x n matrices; fit_founder_lmm_blocks_batch() rebuilds them only
# when a fit is missing.
qtl2mr_prepare_lmm_context_cached <- function(pheno,
                                             kinship = NULL,
                                             addcovar = NULL,
                                             hsq = NULL,
                                             sample_ids = NULL,
                                             standardize_pheno = FALSE,
                                             reml = TRUE,
                                             trait_name = NULL,
                                             chr = NULL,
                                             tol = 1e-12,
                                             context_cache = NULL) {
  compute <- function() prepare_lmm_context(
    pheno = pheno, kinship = kinship, addcovar = addcovar, hsq = hsq,
    sample_ids = sample_ids, standardize_pheno = standardize_pheno,
    reml = reml, trait_name = trait_name, chr = chr, tol = tol
  )
  scoped <- is.environment(context_cache)
  lookup <- function() qtl2mr_context_lookup_key(pheno, kinship, addcovar, hsq, sample_ids,
    standardize_pheno, reml, trait_name, chr, tol)
  # Without a scope the earlier code computed no key, so the key adds no
  # warnings or errors here; if it cannot be formed the context is computed
  # directly and fails exactly as before.
  key <- if (scoped) lookup() else tryCatch(qtl2mr_quietly(lookup()), error = function(e) NULL)
  if (is.null(key)) return(compute())
  reuse <- function(ctx) {
    ctx$materialise <- function() qtl2mr_quietly(compute())
    ctx
  }
  if (scoped && qtl2mr_cache_seen(context_cache, key)) {
    qtl2mr_context_cache_bump(context_cache, "hits")
    entry <- qtl2mr_registry_get("context", key)
    if (!is.null(entry)) return(reuse(entry$value))
    # A hit whose stored value was released: recompute it silently, since a hit
    # signals nothing.
    rec <- qtl2mr_recording(compute(), muffle = TRUE)
    qtl2mr_count_computed("context")
    qtl2mr_registry_put("context", key, qtl2mr_compact_context(rec$value), conditions = rec$conditions,
      info = list(trait = trait_name, chr = as.character(chr %||% NA_character_), marker = NA_character_))
    return(rec$value)
  }
  entry <- qtl2mr_registry_get("context", key)
  if (!is.null(entry)) {
    qtl2mr_replay_conditions(entry$conditions)
    if (scoped) {
      qtl2mr_cache_mark(context_cache, key)
      qtl2mr_context_cache_bump(context_cache, "misses")
    }
    return(reuse(entry$value))
  }
  rec <- qtl2mr_recording(compute())
  qtl2mr_count_computed("context")
  qtl2mr_registry_put("context", key, qtl2mr_compact_context(rec$value), conditions = rec$conditions,
    info = list(trait = trait_name, chr = as.character(chr %||% NA_character_), marker = NA_character_))
  if (scoped) {
    qtl2mr_cache_mark(context_cache, key)
    qtl2mr_context_cache_bump(context_cache, "misses")
  }
  rec$value
}

prepare_lmm_context <- function(pheno,
                                kinship = NULL,
                                addcovar = NULL,
                                hsq = NULL,
                                sample_ids = NULL,
                                standardize_pheno = FALSE,
                                reml = TRUE,
                                trait_name = NULL,
                                chr = NULL,
                                tol = 1e-12) {
  pheno0 <- qtl2mr_as_named_numeric(pheno, trait_name %||% "phenotype")
  base_ids <- if (is.null(sample_ids)) names(pheno0) else as.character(sample_ids)
  ids <- base_ids[base_ids %in% names(pheno0)]
  if (!is.null(addcovar) && nrow(as.matrix(addcovar))) {
    ids <- ids[ids %in% rownames(as.matrix(addcovar))]
  }
  if (!is.null(kinship)) ids <- ids[ids %in% rownames(as.matrix(kinship))]
  if (!length(ids)) stop("No common sample IDs for LMM context", call. = FALSE)
  addcovar0 <- qtl2mr_as_addcovar(addcovar, ids = ids, tol = tol)
  keep <- is.finite(pheno0[ids])
  if (ncol(addcovar0)) keep <- keep & qtl2mr_complete_rows(addcovar0)
  ids <- ids[keep]
  if (length(ids) <= 8L) stop("Need more than 8 complete samples; got ", length(ids), call. = FALSE)
  addcovar0 <- qtl2mr_as_addcovar(addcovar, ids = ids, tol = tol)
  y <- pheno0[ids]
  y_mean <- mean(y)
  y_sd <- stats::sd(y)
  if (isTRUE(standardize_pheno)) {
    if (!is.finite(y_sd) || y_sd <= 0) stop("Cannot standardize phenotype with non-positive SD", call. = FALSE)
    y <- (y - y_mean) / y_sd
  }
  names(y) <- ids
  C <- addcovar0
  colnames(C) <- make.names(colnames(C), unique = TRUE)

  if (!is.null(hsq)) {
    validated_hsq <- qtl2mr_validate_hsq(hsq)
    if (is.null(kinship) && validated_hsq != 0) {
      stop("Nonzero hsq requires kinship; it cannot be used in an OLS fit", call. = FALSE)
    }
  }
  if (is.null(kinship)) {
    raw_transform <- diag(1, nrow = length(ids), ncol = length(ids))
    dimnames(raw_transform) <- list(ids, ids)
    star_rotation <- raw_transform
    ystar <- as.numeric(y)
    Cstar <- C
    inverse_residual <- function(resid_star) stats::setNames(as.numeric(resid_star), ids)
    kinship_decomp <- NULL
    hsq_fit <- NA_real_
  } else {
    kinship0 <- as.matrix(kinship)[ids, ids, drop = FALSE]
    # The hsq fit and the transform decompose the same matrix; the
    # decomposition is computed once and shared with traits on the same samples.
    decompose <- function() qtl2mr_decompose_kinship_cached(kinship, ids, kinship0)
    if (is.null(hsq)) hsq <- qtl2mr_estimate_lmm_hsq(y, kinship0, addcovar0, reml = reml, decompose = decompose)
    hsq_status <- attr(hsq, "status") %||% "fixed"
    hsq_fit <- qtl2mr_validate_hsq(hsq)
    kinship_decomp <- decompose()
    eigenvec <- kinship_decomp$vectors
    wts <- 1 / sqrt(hsq_fit * kinship_decomp$values + (1 - hsq_fit))
    # Row scaling, bit for bit the product diag(wts) %*% eigenvec: every entry
    # of that product is the single term wts[i] * eigenvec[i, j] plus exact
    # zeros, and "+ 0" turns an exact -0 into the +0 the product gives.
    raw_transform <- wts * eigenvec + 0
    colnames(raw_transform) <- ids
    rownames(raw_transform) <- ids
    star_rotation <- eigenvec
    colnames(star_rotation) <- ids
    rownames(star_rotation) <- ids
    ystar <- as.numeric(raw_transform %*% as.numeric(y))
    Cstar <- raw_transform %*% C
    rownames(Cstar) <- ids
    colnames(Cstar) <- colnames(C)
    inverse_residual <- function(resid_star) {
      out <- as.numeric(t(eigenvec) %*% (as.numeric(resid_star) / wts))
      stats::setNames(out, ids)
    }
  }

  if (is.null(rownames(Cstar))) rownames(Cstar) <- ids
  list(
    trait_name = trait_name %||% "trait",
    chr = as.character(chr %||% NA_character_),
    sample_ids = ids,
    y = y,
    y_raw = pheno0[ids],
    y_mean = y_mean,
    y_sd = y_sd,
    addcovar = addcovar0,
    C = C,
    ystar = ystar,
    Cstar = Cstar,
    raw_transform = raw_transform,
    star_rotation = star_rotation,
    inverse_residual = inverse_residual,
    hsq = hsq_fit,
    hsq_status = if (is.null(kinship)) "not_applicable" else hsq_status,
    lmm_backend = if (is.null(kinship)) "ordinary_least_squares" else "internal_base_lmm",
    kinship_decomp = kinship_decomp,
    kinship = if (is.null(kinship)) NULL else kinship0,
    standardize_pheno = isTRUE(standardize_pheno),
    key = qtl2mr_lmm_context_key(trait_name, chr, ids, addcovar0, standardize_pheno, hsq, kinship,
                               phenotype = pheno0[ids], reml = reml, tol = tol)
  )
}

qtl2mr_batch_marker_design <- function(genoprobs_chr,
                                       marker,
                                       sample_ids,
                                       contrasts) {
  markers <- dimnames(genoprobs_chr)[[3L]]
  idx <- match(marker, markers)
  if (is.na(idx)) stop("Marker not found in chromosome allele probabilities: ", marker, call. = FALSE)
  gp_raw <- genoprobs_chr[sample_ids, , idx, drop = FALSE]
  list(
    gp_raw = gp_raw,
    qtl_design = qtl2mr_genoprobs_design(gp_raw, contrasts)
  )
}

qtl2mr_lmm_marker_error <- function(marker,
                                    chr,
                                    role,
                                    error) {
  structure(
    list(
      marker = as.character(marker),
      chr = as.character(chr),
      role = as.character(role %||% NA_character_),
      error_class = class(error)[[1L]],
      error_message = conditionMessage(error)
    ),
    class = c("qtl2mr_lmm_marker_error", "list")
  )
}

qtl2mr_is_lmm_marker_error <- function(x) {
  inherits(x, "qtl2mr_lmm_marker_error")
}

qtl2mr_marker_drop_reason <- function(x) {
  message <- as.character(x$error_message %||% "")
  if (isTRUE(grepl("Founder contrast coefficients are not estimable|rank-deficient founder contrast", message, ignore.case = TRUE))) {
    "nonestimable_founder_contrast"
  } else {
    "marker_fit_error"
  }
}

qtl2mr_lmm_marker_error_row <- function(x) {
  data.frame(
    marker = x$marker %||% NA_character_,
    chr = x$chr %||% NA_character_,
    role = x$role %||% NA_character_,
    drop_reason = qtl2mr_marker_drop_reason(x),
    error_class = x$error_class %||% NA_character_,
    error_message = x$error_message %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# Founder LMM fits of one context at the given markers. Fits are keyed by the
# context, the marker and its probabilities; with a cache scope hits and misses
# are counted as before. A fit computed earlier in the session is reused, and
# the full context is rebuilt only if some fit has to be computed.
fit_founder_lmm_blocks_batch <- function(genoprobs_chr,
                                         markers,
                                         context,
                                         contrasts,
                                         fallback_genoprobs = NULL,
                                         fc_cols = paste0("FC", 1:7),
                                         tol = 1e-12,
                                         fit_cache = NULL,
                                         allow_marker_errors = FALSE,
                                         role = NA_character_) {
  if (!length(markers)) return(list())
  sample_ids <- context$sample_ids
  allele_names <- dimnames(genoprobs_chr)[[2L]]
  if (is.null(allele_names)) allele_names <- rownames(contrasts)
  contrasts <- contrasts[allele_names, , drop = FALSE]
  validate_founder_contrasts(contrasts, rownames(contrasts))
  scoped <- is.environment(fit_cache)
  full <- NULL
  compute <- function(design, marker) {
    if (is.null(full)) full <<- qtl2mr_full_context(context)
    qtl2mr_count_computed("fit_attempt")
    fit <- qtl2mr_fit_founder_lmm_marker(design, marker, full, contrasts, allele_names, fallback_genoprobs, fc_cols, tol)
    qtl2mr_count_computed("fit")
    fit
  }
  info <- function(marker) list(trait = context$trait_name, chr = as.character(context$chr), marker = marker)
  out <- vector("list", length(markers))
  names(out) <- markers
  for (i in seq_along(markers)) {
    marker <- markers[[i]]
    out[[i]] <- tryCatch({
      design <- qtl2mr_batch_marker_design(genoprobs_chr, marker, sample_ids, contrasts)
      fit_key <- qtl2mr_fingerprint_object(list(schema = "lmm-fit-v040-1", context = context$key,
        marker = marker, fc_cols = fc_cols, probabilities = design$gp_raw, contrasts = contrasts, tol = tol))
      entry <- qtl2mr_registry_get("fit", fit_key)
      if (scoped && qtl2mr_cache_seen(fit_cache, fit_key)) {
        qtl2mr_context_cache_bump(fit_cache, "hits")
        if (!is.null(entry)) entry$value else {
          rec <- qtl2mr_recording(compute(design, marker), muffle = TRUE)
          qtl2mr_registry_put("fit", fit_key, rec$value, conditions = rec$conditions, info = info(marker))
          rec$value
        }
      } else if (!is.null(entry)) {
        qtl2mr_replay_conditions(entry$conditions)
        if (scoped) {
          qtl2mr_cache_mark(fit_cache, fit_key)
          qtl2mr_context_cache_bump(fit_cache, "misses")
        }
        entry$value
      } else {
        rec <- qtl2mr_recording(compute(design, marker))
        qtl2mr_registry_put("fit", fit_key, rec$value, conditions = rec$conditions, info = info(marker))
        if (scoped) {
          qtl2mr_cache_mark(fit_cache, fit_key)
          qtl2mr_context_cache_bump(fit_cache, "misses")
        }
        rec$value
      }
    }, error = function(e) {
      if (!isTRUE(allow_marker_errors)) stop(e)
      qtl2mr_lmm_marker_error(marker, context$chr, role = role %||% context$trait_name, error = e)
    })
  }
  out
}


# One founder LMM fit at one marker, given the full context.
qtl2mr_fit_founder_lmm_marker <- function(design, marker, context, contrasts, allele_names, fallback_genoprobs,
                                          fc_cols, tol) {
  sample_ids <- context$sample_ids
  if (!all(is.finite(design$gp_raw))) {
    if (is.null(fallback_genoprobs)) {
      stop("Non-finite founder probabilities at marker ", marker, "; fallback_genoprobs is required", call. = FALSE)
    }
    fit <- fit_founder_lmm_block(
      genoprobs = fallback_genoprobs,
      pheno = context$y_raw,
      kinship = context$kinship,
      addcovar = context$addcovar,
      contrasts = contrasts,
      hsq = if (is.na(context$hsq)) NULL else context$hsq,
      marker = marker,
      sample_ids = sample_ids,
      standardize_pheno = context$standardize_pheno,
      fc_cols = fc_cols,
      tol = tol
    )
    fit
  } else {
    qtl_design <- design$qtl_design
    Zstar <- context$raw_transform %*% qtl_design
    rownames(Zstar) <- sample_ids
    colnames(Zstar) <- colnames(qtl_design)
    Xstar <- cbind(Zstar, context$Cstar)
    colnames(Xstar) <- make.names(colnames(Xstar), unique = TRUE)
    colnames(Xstar)[seq_len(ncol(qtl_design))] <- colnames(qtl_design)
    missing_fc <- setdiff(fc_cols, colnames(Xstar))
    if (length(missing_fc)) stop("Missing founder contrast columns: ", paste(missing_fc, collapse = ", "), call. = FALSE)

    q <- qr(Xstar, tol = tol)
    beta_full <- as.numeric(qr.coef(q, context$ystar))
    names(beta_full) <- colnames(Xstar)
    fitted_star <- as.numeric(qr.fitted(q, context$ystar))
    residual_star_fit <- context$ystar - fitted_star
    df <- length(context$ystar) - q$rank
    if (df <= 0L) stop("Non-positive residual degrees of freedom at marker ", marker, call. = FALSE)
    sigma2 <- as.numeric(crossprod(residual_star_fit)) / df
    sigma <- sqrt(sigma2)
    coef_map_star <- qr.coef(q, diag(nrow(Xstar)))
    rownames(coef_map_star) <- colnames(Xstar)
    colnames(coef_map_star) <- sample_ids
    if (any(!is.finite(coef_map_star[fc_cols, , drop = FALSE])) ||
        any(!is.finite(beta_full[fc_cols]))) {
      stop("Founder contrast coefficients are not estimable at marker ", marker, call. = FALSE)
    }
    C_full <- coef_map_star %*% t(coef_map_star)
    dimnames(C_full) <- list(colnames(Xstar), colnames(Xstar))
    vcov_full <- sigma2 * C_full
    vcov_fc <- qtl2mr_sym(vcov_full[fc_cols, fc_cols, drop = FALSE])
    beta_fc <- stats::setNames(as.numeric(beta_full[fc_cols]), fc_cols)
    se_fc <- stats::setNames(sqrt(diag(vcov_fc)), fc_cols)
    influence_raw <- coef_map_star[fc_cols, , drop = FALSE] %*% context$raw_transform
    colnames(influence_raw) <- sample_ids
    rownames(influence_raw) <- fc_cols
    influence_star <- coef_map_star[fc_cols, , drop = FALSE] %*% context$star_rotation
    colnames(influence_star) <- sample_ids
    rownames(influence_star) <- fc_cols
    residual_raw <- context$inverse_residual(residual_star_fit)
    residual_star <- as.numeric(t(context$star_rotation) %*% as.numeric(residual_star_fit))
    fit <- list(
      marker = marker,
      chr = as.character(context$chr),
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
      residual_star = stats::setNames(residual_star, sample_ids),
      residual_star_fit = stats::setNames(as.numeric(residual_star_fit), sample_ids),
      influence_raw = influence_raw,
      influence_star = influence_star,
      influence_star_fit = coef_map_star[fc_cols, , drop = FALSE],
      hsq = context$hsq,
      hsq_status = context$hsq_status,
      lmm_backend = context$lmm_backend,
      n = length(sample_ids),
      df = df,
      rank = q$rank,
      sample_ids = sample_ids,
      design_colnames = colnames(Xstar),
      fc_cols = fc_cols,
      standardize_pheno = context$standardize_pheno
    )
    fit
  }
}

qtl2mr_assemble_mr_locus_block <- function(marker,
                                           chr,
                                           fit_x,
                                           fit_y,
                                           cross_cov_model = c("star", "raw", "zero"),
                                           O_xy = NULL,
                                           overlap_ids = NULL,
                                           ridge = 1e-10,
                                           fc_cols = paste0("FC", 1:7),
                                           enforce_psd = TRUE,
                                           stop_on_invalid_psd = TRUE,
                                           overlap_products = NULL) {
  cross_cov_model <- match.arg(cross_cov_model)
  if (is.null(O_xy)) {
    O_xy <- make_overlap_matrix(fit_x$sample_ids, fit_y$sample_ids, overlap_ids = overlap_ids)
  }
  cross_raw <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "raw", O_xy = O_xy,
    overlap_products = overlap_products)
  cross_star <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "star", O_xy = O_xy,
    overlap_products = overlap_products)
  cross_zero <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "zero", O_xy = O_xy)
  selected_cross <- switch(cross_cov_model, raw = cross_raw, star = cross_star, zero = cross_zero)
  Sigma_x <- fit_x$vcov_fc
  Sigma_y <- fit_y$vcov_fc
  psd_raw <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_raw$Sigma_xy, ridge = ridge)
  psd_star <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_star$Sigma_xy, ridge = ridge)
  psd_zero <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_zero$Sigma_xy, ridge = ridge)
  psd <- switch(cross_cov_model, raw = psd_raw, star = psd_star, zero = psd_zero)
  Sigma_xy <- selected_cross$Sigma_xy
  if (isTRUE(enforce_psd)) {
    if (!isTRUE(psd$psd_valid) && isTRUE(stop_on_invalid_psd)) {
      stop("Selected ", cross_cov_model, " Sigma_xy fails joint PSD validation at marker ", marker,
           " (status=", psd$status, ")", call. = FALSE)
    }
    Sigma_x <- psd$Sigma_x
    Sigma_y <- psd$Sigma_y
    Sigma_xy <- psd$Sigma_xy
  }
  list(
    marker = marker,
    chr = as.character(chr),
    beta_x = fit_x$beta_fc,
    beta_y = fit_y$beta_fc,
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    Sigma_xy = Sigma_xy,
    Sigma_x_raw = fit_x$vcov_fc,
    Sigma_y_raw = fit_y$vcov_fc,
    Sigma_xy_raw = selected_cross$Sigma_xy,
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
    standardize_phenotypes = isTRUE(fit_x$standardize_pheno)
  )
}

fit_mr_locus_blocks_batch <- function(genoprobs,
                                      selected_loci,
                                      exposure,
                                      outcome,
                                      kinship = NULL,
                                      addcovar = NULL,
                                      contrasts = NULL,
                                      hsq_x = NULL,
                                      hsq_y = NULL,
                                      sample_ids_x = NULL,
                                      sample_ids_y = NULL,
                                      standardize_phenotypes = FALSE,
                                      cross_cov_model = c("star", "raw", "zero"),
                                      ridge = 1e-10,
                                      reml = TRUE,
                                      fc_cols = paste0("FC", 1:7),
                                      tol = 1e-12,
                                      progress_callback = NULL,
                                      context_cache = NULL,
                                      fit_cache = NULL,
                                      exposure_name = "exposure",
                                      outcome_name = "outcome",
                                      drop_nonestimable_loci = FALSE,
                                      min_retained_loci = 1L) {
  cross_cov_model <- match.arg(cross_cov_model)
  if (!is.data.frame(selected_loci) || !all(c("marker", "chr") %in% names(selected_loci))) {
    stop("selected_loci must be a data frame with marker and chr columns", call. = FALSE)
  }
  if (anyNA(selected_loci$marker) || anyDuplicated(selected_loci$marker)) stop("selected_loci markers must be non-missing and unique", call. = FALSE)
  x <- qtl2mr_as_named_numeric(exposure, "exposure")
  y <- qtl2mr_as_named_numeric(outcome, "outcome")
  founders <- dimnames(genoprobs[[1L]])[[2L]]
  if (is.null(contrasts)) contrasts <- make_founder_contrasts(founders)
  out <- list()
  dropped <- list()
  # The overlap matrix and its products D_x = O O', D_y = O'O depend only on the
  # two fitted sample sets, which are the same at every locus of a chromosome
  # and usually of a pair; they are formed once per distinct pair of sets.
  overlaps <- list()
  overlap_for <- function(ids_x, ids_y) {
    for (o in overlaps) if (identical(o$ids_x, ids_x) && identical(o$ids_y, ids_y)) return(o)
    O <- make_overlap_matrix(ids_x, ids_y)
    o <- list(ids_x = ids_x, ids_y = ids_y, O_xy = O, products = qtl2mr_overlap_products(O))
    overlaps[[length(overlaps) + 1L]] <<- o
    o
  }
  for (chr in unique(as.character(selected_loci$chr))) {
    rows <- which(as.character(selected_loci$chr) == chr)
    markers <- selected_loci$marker[rows]
    gp_chr <- genoprobs[[chr]]
    if (is.null(gp_chr)) stop("genoprobs is missing selected chromosome: ", chr, call. = FALSE)
    K <- if (is.list(kinship) && !is.matrix(kinship) && !is.data.frame(kinship)) kinship[[chr]] else kinship
    qtl2mr_progress(progress_callback, 2L, "Fit selected-locus founder LMM/MR blocks",
                    sprintf("chr%s batch: %d loci", chr, length(markers)))
    ctx_x <- qtl2mr_prepare_lmm_context_cached(
      pheno = x,
      kinship = K,
      addcovar = addcovar,
      hsq = qtl2mr_chr_option(hsq_x, chr),
      sample_ids = sample_ids_x,
      standardize_pheno = standardize_phenotypes,
      reml = reml,
      trait_name = exposure_name,
      chr = chr,
      tol = tol,
      context_cache = context_cache
    )
    ctx_y <- qtl2mr_prepare_lmm_context_cached(
      pheno = y,
      kinship = K,
      addcovar = addcovar,
      hsq = qtl2mr_chr_option(hsq_y, chr),
      sample_ids = sample_ids_y,
      standardize_pheno = standardize_phenotypes,
      reml = reml,
      trait_name = outcome_name,
      chr = chr,
      tol = tol,
      context_cache = context_cache
    )
    fits_x <- fit_founder_lmm_blocks_batch(gp_chr, markers, ctx_x, contrasts,
                                           fallback_genoprobs = genoprobs,
                                           fc_cols = fc_cols, tol = tol,
                                           fit_cache = fit_cache,
                                           allow_marker_errors = drop_nonestimable_loci,
                                           role = "exposure")
    fits_y <- fit_founder_lmm_blocks_batch(gp_chr, markers, ctx_y, contrasts,
                                           fallback_genoprobs = genoprobs,
                                           fc_cols = fc_cols, tol = tol,
                                           fit_cache = fit_cache,
                                           allow_marker_errors = drop_nonestimable_loci,
                                           role = "outcome")
    for (marker in markers) {
      err_rows <- list()
      if (qtl2mr_is_lmm_marker_error(fits_x[[marker]])) {
        err_rows[[length(err_rows) + 1L]] <- qtl2mr_lmm_marker_error_row(fits_x[[marker]])
      }
      if (qtl2mr_is_lmm_marker_error(fits_y[[marker]])) {
        err_rows[[length(err_rows) + 1L]] <- qtl2mr_lmm_marker_error_row(fits_y[[marker]])
      }
      if (length(err_rows)) {
        dropped[[length(dropped) + 1L]] <- do.call(rbind, err_rows)
        next
      }
      overlap <- overlap_for(fits_x[[marker]]$sample_ids, fits_y[[marker]]$sample_ids)
      block <- qtl2mr_assemble_mr_locus_block(
        marker = marker,
        chr = chr,
        fit_x = fits_x[[marker]],
        fit_y = fits_y[[marker]],
        cross_cov_model = cross_cov_model,
        O_xy = overlap$O_xy,
        ridge = ridge,
        fc_cols = fc_cols,
        overlap_products = overlap$products
      )
      out[[marker]] <- block
    }
  }
  retained <- out[selected_loci$marker[selected_loci$marker %in% names(out)]]
  dropped_df <- if (length(dropped)) do.call(rbind, dropped) else data.frame(
    marker = character(),
    chr = character(),
    role = character(),
    drop_reason = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  rownames(dropped_df) <- NULL
  n_selected <- nrow(selected_loci)
  n_dropped_markers <- length(unique(dropped_df$marker))
  attr(retained, "dropped_loci") <- dropped_df
  attr(retained, "n_loci_selected") <- n_selected
  attr(retained, "n_loci_retained") <- length(retained)
  attr(retained, "n_loci_dropped") <- n_dropped_markers
  if (length(retained) < min_retained_loci) {
    stop(
      "Too few retained loci after dropping non-estimable selected loci: retained ",
      length(retained), " of ", n_selected, " selected loci.",
      call. = FALSE
    )
  }
  retained
}
