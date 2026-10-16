qtl2mr_validate_hsq <- function(hsq) {
  if (length(hsq) != 1L || !is.numeric(hsq) || is.object(hsq) || !is.finite(hsq)) {
    stop("hsq must be one finite numeric value", call. = FALSE)
  }
  min(max(as.numeric(hsq), 0), 1)
}

qtl2mr_lmm_backend_provenance <- function(blocks = NULL) {
  fit_backends <- character()
  if (is.list(blocks) && length(blocks)) {
    fit_backends <- unlist(lapply(blocks, function(block) {
      c(block$fit_x$lmm_backend %||% character(), block$fit_y$lmm_backend %||% character())
    }), use.names = FALSE)
  }
  fit_backends <- sort(unique(fit_backends[!is.na(fit_backends) & nzchar(fit_backends)]))
  backend <- if (!length(fit_backends)) {
    "not_available"
  } else if (length(fit_backends) == 1L) {
    fit_backends[[1L]]
  } else {
    "mixed"
  }
  list(
    backend = backend,
    fit_backends = fit_backends,
    qtl2_role = "optional_parity_validation"
  )
}

qtl2mr_backend_label <- function(backends) {
  backends <- sort(unique(as.character(backends)))
  backends <- backends[!is.na(backends) & nzchar(backends)]
  if (!length(backends)) return("not_available")
  if (length(backends) == 1L) return(backends[[1L]])
  "mixed"
}

qtl2mr_decompose_kinship <- function(kinship,
                                     sample_ids = NULL,
                                     symmetry_tol = sqrt(.Machine$double.eps),
                                     psd_tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(kinship)) kinship <- as.matrix(kinship)
  if (!is.numeric(kinship) || is.object(kinship) || length(dim(kinship)) != 2L ||
      nrow(kinship) < 1L || nrow(kinship) != ncol(kinship)) {
    stop("kinship must be a non-empty square numeric matrix", call. = FALSE)
  }
  if (any(!is.finite(kinship))) stop("kinship must contain only finite values", call. = FALSE)
  rn <- rownames(kinship)
  cn <- colnames(kinship)
  if (is.null(rn) || is.null(cn) || anyNA(rn) || anyNA(cn) ||
      any(!nzchar(rn)) || any(!nzchar(cn)) || anyDuplicated(rn) || anyDuplicated(cn)) {
    stop("kinship must have unique, non-blank row and column sample IDs", call. = FALSE)
  }
  if (!setequal(rn, cn)) stop("kinship row and column sample IDs must match", call. = FALSE)
  if (!is.null(sample_ids)) {
    sample_ids <- as.character(sample_ids)
    if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
      stop("sample_ids must be unique and non-blank", call. = FALSE)
    }
    missing_ids <- setdiff(sample_ids, intersect(rn, cn))
    if (length(missing_ids)) {
      stop("kinship is missing sample IDs: ", paste(missing_ids, collapse = ", "), call. = FALSE)
    }
    kinship <- kinship[sample_ids, sample_ids, drop = FALSE]
    rn <- cn <- sample_ids
  } else if (!identical(rn, cn)) {
    kinship <- kinship[rn, rn, drop = FALSE]
    cn <- rn
  }
  scale <- max(abs(kinship), .Machine$double.eps)
  asymmetry <- max(abs(kinship - t(kinship)))
  if (asymmetry > symmetry_tol * scale) {
    stop("kinship must be symmetric within numerical tolerance", call. = FALSE)
  }
  kinship <- qtl2mr_sym(kinship)
  kinship_model <- 2 * kinship
  eig <- eigen(kinship_model, symmetric = TRUE)
  eigen_scale <- max(abs(eig$values), .Machine$double.eps)
  negative_tol <- psd_tol * eigen_scale
  if (min(eig$values) < -negative_tol) {
    stop(
      "kinship must be positive semidefinite; minimum eigenvalue is ",
      signif(min(eig$values), 8),
      call. = FALSE
    )
  }
  values <- pmax(eig$values, 0)
  order_ascending <- order(values)
  values <- values[order_ascending]
  vectors <- t(eig$vectors[, order_ascending, drop = FALSE])
  colnames(vectors) <- rn
  rownames(vectors) <- rn
  structure(
    list(
      values = values,
      vectors = vectors,
      backend = "internal_base_lmm",
      doubled_kinship = TRUE,
      min_eigenvalue_raw = min(eig$values),
      negative_eigenvalue_tolerance = negative_tol,
      clamped_negative_eigenvalues = sum(eig$values < 0)
    ),
    class = c("qtl2mr_kinship_decomposition", "list")
  )
}

qtl2mr_hsq_profile <- function(h,
                               y_rotated,
                               design_rotated,
                               eigenvalues,
                               reml,
                               logdet_design_crossprod,
                               tol = 1e-12) {
  variance <- (1 - h) + h * eigenvalues
  variance_tol <- .Machine$double.eps * max(1, max(abs(variance)))
  if (any(!is.finite(variance)) || any(variance <= variance_tol)) return(Inf)
  weights <- 1 / sqrt(variance)
  y_weighted <- weights * y_rotated
  design_weighted <- design_rotated * weights
  qr_design <- qr(design_weighted, tol = tol)
  p <- ncol(design_weighted)
  n <- length(y_weighted)
  if (qr_design$rank != p || n <= p) return(Inf)
  residual <- qr.resid(qr_design, y_weighted)
  rss <- as.numeric(crossprod(residual))
  if (!is.finite(rss) || rss <= .Machine$double.eps) return(Inf)
  if (isTRUE(reml)) {
    logdet_weighted <- as.numeric(determinant(crossprod(design_weighted), logarithm = TRUE)$modulus)
    if (!is.finite(logdet_weighted)) return(Inf)
    sum(log(variance)) + logdet_weighted - logdet_design_crossprod +
      (n - p) * log(rss / (n - p))
  } else {
    sum(log(variance)) + n * log(rss / n)
  }
}

qtl2mr_estimate_lmm_hsq <- function(pheno,
                                    kinship,
                                    addcovar = NULL,
                                    reml = TRUE,
                                    tol = 1e-12) {
  pheno <- qtl2mr_as_named_numeric(pheno)
  ids <- names(pheno)
  if (any(!is.finite(pheno))) stop("pheno must contain only finite values for hsq estimation", call. = FALSE)
  addcovar <- qtl2mr_as_addcovar(addcovar, ids = ids, tol = tol)
  if (ncol(addcovar) && any(!is.finite(addcovar))) {
    stop("addcovar must contain only finite values for hsq estimation", call. = FALSE)
  }
  decomposition <- qtl2mr_decompose_kinship(kinship, sample_ids = ids)
  design <- cbind(`(Intercept)` = 1, addcovar)
  if (qr(design, tol = tol)$rank != ncol(design)) {
    stop("hsq fixed-effect design is rank deficient after covariate filtering", call. = FALSE)
  }
  logdet_design <- as.numeric(determinant(crossprod(design), logarithm = TRUE)$modulus)
  y_rotated <- as.numeric(decomposition$vectors %*% as.numeric(pheno))
  design_rotated <- decomposition$vectors %*% design
  objective <- function(h) qtl2mr_hsq_profile(
    h = h,
    y_rotated = y_rotated,
    design_rotated = design_rotated,
    eigenvalues = decomposition$values,
    reml = reml,
    logdet_design_crossprod = logdet_design,
    tol = tol
  )
  interior <- stats::optimize(objective, interval = c(0, 1), tol = max(tol, 1e-10))
  candidates <- c(0, interior$minimum, 1)
  objective_values <- vapply(candidates, objective, numeric(1))
  if (!any(is.finite(objective_values))) {
    stop("Internal hsq profile likelihood is not finite on [0, 1]", call. = FALSE)
  }
  best <- which.min(objective_values)
  hsq <- candidates[[best]]
  fit_variance <- (1 - hsq) + hsq * decomposition$values
  weights <- 1 / sqrt(fit_variance)
  qr_design <- qr(design_rotated * weights, tol = tol)
  rss <- as.numeric(crossprod(qr.resid(qr_design, weights * y_rotated)))
  denominator <- if (isTRUE(reml)) length(pheno) - qr_design$rank else length(pheno)
  attr(hsq, "backend") <- "internal_base_lmm"
  attr(hsq, "status") <- if (best == 1L) "boundary_zero" else if (best == 3L) "boundary_one" else "success"
  attr(hsq, "objective") <- objective_values[[best]]
  attr(hsq, "sigmasq") <- rss / denominator
  attr(hsq, "reml") <- isTRUE(reml)
  attr(hsq, "n") <- length(pheno)
  attr(hsq, "rank") <- qr_design$rank
  hsq
}
