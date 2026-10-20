qtl2mr_sym <- function(x) {
  x <- as.matrix(x)
  (x + t(x)) / 2
}

qtl2mr_add_invariant_ridge <- function(Omega, ridge) {
  Omega <- as.matrix(Omega)
  if (length(dim(Omega)) != 2L || nrow(Omega) != ncol(Omega) || nrow(Omega) < 1L) {
    stop("Omega must be a non-empty square matrix", call. = FALSE)
  }
  if (!typeof(Omega) %in% c("double", "integer") || is.object(Omega)) {
    stop("Omega must be numeric", call. = FALSE)
  }
  ridge_is_numeric_scalar <- length(ridge) == 1L &&
    typeof(ridge) %in% c("double", "integer") && !is.object(ridge)
  if (!ridge_is_numeric_scalar || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be one finite non-negative number", call. = FALSE)
  }
  if (any(!is.finite(Omega))) stop("Omega must contain only finite values", call. = FALSE)
  matrix_scale <- max(abs(Omega), .Machine$double.eps)
  if (max(abs(Omega - t(Omega))) > sqrt(.Machine$double.eps) * matrix_scale) {
    stop("Omega must be symmetric", call. = FALSE)
  }
  Omega <- qtl2mr_sym(Omega)
  eigenvalues <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  eigen_scale <- max(abs(eigenvalues), .Machine$double.eps)
  if (min(eigenvalues) < -sqrt(.Machine$double.eps) * eigen_scale) {
    stop("Omega must be positive semidefinite", call. = FALSE)
  }
  if (ridge == 0) return(Omega)
  scale <- mean(diag(Omega))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  Omega + diag(ridge * scale, nrow(Omega))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

qtl2mr_with_seed <- function(seed = NULL, expr) {
  if (is.null(seed)) return(force(expr))
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (old_seed_exists) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

qtl2mr_require_qtl2 <- function() {
  if (!requireNamespace("qtl2", quietly = TRUE)) {
    stop("Package 'qtl2' is required only for the requested optional qtl2 parity validation", call. = FALSE)
  }
  invisible(TRUE)
}

qtl2mr_as_named_numeric <- function(x, label = "phenotype") {
  if (is.data.frame(x) || is.matrix(x)) {
    if (ncol(x) != 1L) stop(label, " must have exactly one column", call. = FALSE)
    rn <- rownames(x)
    x <- x[, 1L]
    names(x) <- rn
  }
  out <- as.numeric(x)
  names(out) <- names(x)
  if (is.null(names(out)) || any(!nzchar(names(out)))) {
    stop(label, " must be a named vector or one-column object with row names", call. = FALSE)
  }
  out
}

make_overlap_matrix <- function(sample_ids_x, sample_ids_y, overlap_ids = NULL) {
  sample_ids_x <- as.character(sample_ids_x)
  sample_ids_y <- as.character(sample_ids_y)
  if (anyNA(sample_ids_x) || any(!nzchar(sample_ids_x))) {
    stop("sample_ids_x must not contain missing or blank IDs", call. = FALSE)
  }
  if (anyNA(sample_ids_y) || any(!nzchar(sample_ids_y))) {
    stop("sample_ids_y must not contain missing or blank IDs", call. = FALSE)
  }
  dup_x <- unique(sample_ids_x[duplicated(sample_ids_x)])
  dup_y <- unique(sample_ids_y[duplicated(sample_ids_y)])
  if (length(dup_x)) {
    stop("sample_ids_x contains duplicate IDs: ", paste(dup_x, collapse = ", "), call. = FALSE)
  }
  if (length(dup_y)) {
    stop("sample_ids_y contains duplicate IDs: ", paste(dup_y, collapse = ", "), call. = FALSE)
  }
  if (is.null(overlap_ids)) {
    overlap_ids <- intersect(sample_ids_x, sample_ids_y)
  } else {
    overlap_ids <- as.character(overlap_ids)
    if (anyNA(overlap_ids) || any(!nzchar(overlap_ids))) {
      stop("overlap_ids must not contain missing or blank IDs", call. = FALSE)
    }
    dup_o <- unique(overlap_ids[duplicated(overlap_ids)])
    if (length(dup_o)) {
      stop("overlap_ids contains duplicate IDs: ", paste(dup_o, collapse = ", "), call. = FALSE)
    }
    missing_x <- setdiff(overlap_ids, sample_ids_x)
    missing_y <- setdiff(overlap_ids, sample_ids_y)
    if (length(missing_x) || length(missing_y)) {
      stop("overlap_ids must be present in both sample ID vectors", call. = FALSE)
    }
  }
  O_xy <- matrix(0, nrow = length(sample_ids_x), ncol = length(sample_ids_y),
                 dimnames = list(sample_ids_x, sample_ids_y))
  if (length(overlap_ids)) {
    O_xy[cbind(match(overlap_ids, sample_ids_x), match(overlap_ids, sample_ids_y))] <- 1
  }
  O_xy
}

estimate_rho_star <- function(resid_x_star,
                              resid_y_star,
                              O_xy,
                              center = TRUE,
                              min_overlap_warn = 10L,
                              eps = 1e-12) {
  O_xy <- as.matrix(O_xy)
  rx <- as.numeric(resid_x_star)
  ry <- as.numeric(resid_y_star)
  names(rx) <- names(resid_x_star)
  names(ry) <- names(resid_y_star)

  if (!is.null(rownames(O_xy)) && !is.null(names(rx))) {
    missing_x <- setdiff(rownames(O_xy), names(rx))
    if (length(missing_x)) {
      stop("resid_x_star is missing IDs required by O_xy rows", call. = FALSE)
    }
    rx <- rx[rownames(O_xy)]
  } else if (length(rx) != nrow(O_xy)) {
    stop("length(resid_x_star) must match nrow(O_xy)", call. = FALSE)
  }
  if (!is.null(colnames(O_xy)) && !is.null(names(ry))) {
    missing_y <- setdiff(colnames(O_xy), names(ry))
    if (length(missing_y)) {
      stop("resid_y_star is missing IDs required by O_xy columns", call. = FALSE)
    }
    ry <- ry[colnames(O_xy)]
  } else if (length(ry) != ncol(O_xy)) {
    stop("length(resid_y_star) must match ncol(O_xy)", call. = FALSE)
  }
  if (!all(is.finite(rx)) || !all(is.finite(ry))) {
    stop("Residual vectors must be finite", call. = FALSE)
  }

  D_x <- O_xy %*% t(O_xy)
  D_y <- t(O_xy) %*% O_xy
  one_x <- rep(1, nrow(O_xy))
  one_y <- rep(1, ncol(O_xy))
  denom_mean_x <- as.numeric(t(one_x) %*% D_x %*% one_x)
  denom_mean_y <- as.numeric(t(one_y) %*% D_y %*% one_y)
  if (!is.finite(denom_mean_x) || !is.finite(denom_mean_y) ||
      denom_mean_x <= eps || denom_mean_y <= eps) {
    stop("O_xy contains no usable sample overlap", call. = FALSE)
  }
  n_overlap <- sum(O_xy != 0)
  if (n_overlap > 0L && n_overlap < min_overlap_warn) {
    warning(
      "Only ", n_overlap, " overlapping samples are available for rho_star; ",
      "Sigma_xy may be unstable.",
      call. = FALSE
    )
  }

  ss_x <- as.numeric(t(rx) %*% D_x %*% rx)
  ss_y <- as.numeric(t(ry) %*% D_y %*% ry)
  if (isTRUE(center)) {
    mx <- as.numeric(t(one_x) %*% D_x %*% rx) / denom_mean_x
    my <- as.numeric(t(one_y) %*% D_y %*% ry) / denom_mean_y
    rx <- rx - mx
    ry <- ry - my
  }

  vx <- as.numeric(t(rx) %*% D_x %*% rx)
  vy <- as.numeric(t(ry) %*% D_y %*% ry)
  # Relative to the uncentred sum of squares, so the check does not depend on
  # the units of the traits.
  if (!is.finite(vx) || !is.finite(vy) || vx <= eps * ss_x || vy <= eps * ss_y) {
    stop("Residual variance on the overlap is zero or non-finite", call. = FALSE)
  }
  rho <- as.numeric(t(rx) %*% O_xy %*% ry) / sqrt(vx * vy)
  if (!is.finite(rho)) stop("rho_star is not finite", call. = FALSE)
  if (abs(rho) <= 1 + eps) {
    rho <- max(min(rho, 1), -1)
  } else {
    stop("rho_star is outside [-1, 1] beyond numerical tolerance", call. = FALSE)
  }
  rho
}

qtl2mr_as_addcovar <- function(addcovar, ids = NULL, tol = 1e-12) {
  if (is.null(addcovar)) {
    if (is.null(ids)) return(matrix(numeric(), nrow = 0L, ncol = 0L))
    return(matrix(numeric(), nrow = length(ids), ncol = 0L, dimnames = list(ids, character())))
  }
  if (is.data.frame(addcovar)) {
    bad <- !vapply(addcovar, function(x) is.numeric(x) || is.integer(x) || is.logical(x), logical(1))
    if (any(bad)) {
      stop(
        "addcovar must be numeric/logical after preprocessing; non-numeric columns: ",
        paste(names(addcovar)[bad], collapse = ", "),
        call. = FALSE
      )
    }
    addcovar <- data.matrix(addcovar)
  } else {
    addcovar <- as.matrix(addcovar)
    if (!typeof(addcovar) %in% c("double", "integer", "logical")) {
      stop("addcovar must be numeric/logical after preprocessing", call. = FALSE)
    }
    storage.mode(addcovar) <- "double"
  }
  if (is.null(rownames(addcovar)) || any(!nzchar(rownames(addcovar)))) {
    stop("addcovar must have row names matching sample IDs", call. = FALSE)
  }
  if (ncol(addcovar) && is.null(colnames(addcovar))) {
    colnames(addcovar) <- paste0("covar", seq_len(ncol(addcovar)))
  }
  if (!is.null(ids)) addcovar <- addcovar[ids, , drop = FALSE]
  if (!ncol(addcovar)) return(addcovar)
  qtl2mr_drop_depcols(addcovar, add_intercept = TRUE, tol = tol)
}

qtl2mr_drop_depcols <- function(covar, add_intercept = TRUE, tol = 1e-12) {
  if (is.null(covar) || !ncol(covar)) return(covar)
  covar <- as.matrix(covar)
  if (!typeof(covar) %in% c("double", "integer", "logical")) {
    stop("covar must be numeric/logical for the rank check", call. = FALSE)
  }
  original_colnames <- colnames(covar) %||% paste0("covar", seq_len(ncol(covar)))
  X <- covar[rowSums(!is.finite(covar)) == 0L, , drop = FALSE]
  if (!nrow(X)) stop("No finite rows available for covariate rank check", call. = FALSE)
  base <- if (isTRUE(add_intercept)) matrix(1, nrow = nrow(X), ncol = 1L) else matrix(nrow = nrow(X), ncol = 0L)
  keep <- logical(ncol(X))
  current <- base
  current_rank <- if (ncol(current)) qr(current, tol = tol)$rank else 0L
  for (j in seq_len(ncol(X))) {
    candidate <- cbind(current, X[, j, drop = FALSE])
    candidate_rank <- qr(candidate, tol = tol)$rank
    if (candidate_rank > current_rank) {
      keep[[j]] <- TRUE
      current <- candidate
      current_rank <- candidate_rank
    }
  }
  dropped <- original_colnames[!keep]
  if (length(dropped)) {
    warning(
      "Dropping linearly dependent covariate column(s): ",
      paste(dropped, collapse = ", "),
      call. = FALSE
    )
  }
  covar[, keep, drop = FALSE]
}

block_diag_matrix <- function(mats) {
  if (length(mats) == 0L) stop("Need at least one matrix", call. = FALSE)
  nr <- vapply(mats, nrow, integer(1))
  nc <- vapply(mats, ncol, integer(1))
  row_names <- unlist(lapply(seq_along(mats), function(i) {
    rn <- rownames(mats[[i]])
    if (is.null(rn)) rn <- paste0("block", i, "_r", seq_len(nr[[i]]))
    rn
  }), use.names = FALSE)
  col_names <- unlist(lapply(seq_along(mats), function(i) {
    cn <- colnames(mats[[i]])
    if (is.null(cn)) cn <- paste0("block", i, "_c", seq_len(nc[[i]]))
    cn
  }), use.names = FALSE)
  out <- matrix(0, sum(nr), sum(nc), dimnames = list(row_names, col_names))
  r0 <- 0L
  c0 <- 0L
  for (i in seq_along(mats)) {
    rr <- seq.int(r0 + 1L, r0 + nr[[i]])
    cc <- seq.int(c0 + 1L, c0 + nc[[i]])
    out[rr, cc] <- mats[[i]]
    r0 <- r0 + nr[[i]]
    c0 <- c0 + nc[[i]]
  }
  out
}

# Floors for eigenvalues are relative to the largest one, so results do not
# depend on the units of the traits.
qtl2mr_eigen_scale <- function(values) {
  scale <- max(abs(values))
  if (!is.finite(scale) || scale <= 0) 1 else scale
}

qtl2mr_mat_inv <- function(S, ridge = 1e-10) {
  S <- qtl2mr_sym(S)
  ee <- eigen(S, symmetric = TRUE)
  scale <- qtl2mr_eigen_scale(ee$values)
  lam <- pmax(ee$values, max(ridge, .Machine$double.eps) * scale)
  ee$vectors %*% (diag(1 / lam, nrow = length(lam)) %*% t(ee$vectors))
}

qtl2mr_psd_clip <- function(S, ridge = 1e-10) {
  S <- qtl2mr_sym(S)
  ee <- eigen(S, symmetric = TRUE)
  scale <- qtl2mr_eigen_scale(ee$values)
  lam <- pmax(ee$values, max(ridge, .Machine$double.eps) * scale)
  out <- ee$vectors %*% (diag(lam, nrow = length(lam)) %*% t(ee$vectors))
  qtl2mr_sym(out)
}

qtl2mr_matrix_power <- function(S, power, ridge = 1e-10) {
  S <- qtl2mr_sym(S)
  ee <- eigen(S, symmetric = TRUE)
  scale <- qtl2mr_eigen_scale(ee$values)
  lam <- pmax(ee$values, max(ridge, .Machine$double.eps) * scale)
  out <- ee$vectors %*% (diag(lam^power, nrow = length(lam)) %*% t(ee$vectors))
  qtl2mr_sym(out)
}

qtl2mr_joint_psd_diagnostics <- function(Sigma_x,
                                         Sigma_y,
                                         Sigma_xy,
                                         ridge = 1e-10,
                                         tol = 1e-8,
                                         pass_tol = 1e-8) {
  joint <- rbind(cbind(Sigma_x, Sigma_xy), cbind(t(Sigma_xy), Sigma_y))
  joint <- qtl2mr_sym(joint)
  vals <- eigen(joint, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(vals), .Machine$double.eps)
  Sx_inv_half <- qtl2mr_matrix_power(Sigma_x, -0.5, ridge = ridge)
  Sy_inv_half <- qtl2mr_matrix_power(Sigma_y, -0.5, ridge = ridge)
  whitened_cross <- Sx_inv_half %*% Sigma_xy %*% Sy_inv_half
  sv <- svd(whitened_cross)
  smax <- if (length(sv$d)) max(sv$d) else 0
  eigen_pass <- min(vals) >= -tol * scale
  canonical_pass <- smax <= 1 + pass_tol
  list(
    min_eigen = min(vals),
    eigen_scale = scale,
    eigen_pass = eigen_pass,
    smax = smax,
    canonical_pass = canonical_pass,
    psd_valid = isTRUE(eigen_pass) && isTRUE(canonical_pass),
    singular_values = sv$d
  )
}

validate_joint_psd <- function(Sigma_x,
                               Sigma_y,
                               Sigma_xy,
                               ridge = 1e-10,
                               tol = 1e-8,
                               pass_tol = 1e-8) {
  out <- qtl2mr_joint_psd_diagnostics(
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    Sigma_xy = Sigma_xy,
    ridge = ridge,
    tol = tol,
    pass_tol = pass_tol
  )
  list(
    psd_valid = isTRUE(out$psd_valid),
    min_eigenvalue = out$min_eigen,
    eigen_scale = out$eigen_scale,
    eigen_pass = isTRUE(out$eigen_pass),
    max_cross_singular_value = out$smax,
    canonical_pass = isTRUE(out$canonical_pass),
    singular_values = out$singular_values
  )
}

qtl2mr_enforce_joint_psd <- function(Sigma_x,
                                     Sigma_y,
                                     Sigma_xy,
                                     ridge = 1e-10,
                                     tol = 1e-8,
                                     pass_tol = 1e-8,
                                     clip_tol = 1e-6,
                                     repair = c("tiny", "always", "none")) {
  repair <- match.arg(repair)
  p <- nrow(Sigma_x)
  Sigma_x <- qtl2mr_sym(Sigma_x)
  Sigma_y <- qtl2mr_sym(Sigma_y)
  Sigma_xy <- as.matrix(Sigma_xy)
  Sigma_xy_input <- Sigma_xy
  before <- qtl2mr_joint_psd_diagnostics(Sigma_x, Sigma_y, Sigma_xy,
                                         ridge = ridge, tol = tol, pass_tol = pass_tol)
  Sx_inv_half <- qtl2mr_matrix_power(Sigma_x, -0.5, ridge = ridge)
  Sy_inv_half <- qtl2mr_matrix_power(Sigma_y, -0.5, ridge = ridge)
  Sx_half <- qtl2mr_matrix_power(Sigma_x, 0.5, ridge = ridge)
  Sy_half <- qtl2mr_matrix_power(Sigma_y, 0.5, ridge = ridge)
  whitened_cross <- Sx_inv_half %*% Sigma_xy %*% Sy_inv_half
  sv <- svd(whitened_cross)
  max_allowed <- 1 - max(ridge, .Machine$double.eps)
  max_sv_before <- before$smax
  cross_excess <- max_sv_before - 1
  tiny_violation <- cross_excess > pass_tol && cross_excess <= clip_tol
  large_violation <- cross_excess > clip_tol
  candidate_psd_valid <- before$psd_valid
  can_clip <- identical(repair, "always") ||
    (identical(repair, "tiny") && !isTRUE(large_violation) && !isTRUE(candidate_psd_valid))
  adjusted <- FALSE
  if (isTRUE(can_clip) && !isTRUE(candidate_psd_valid)) {
    clipped <- pmin(sv$d, max_allowed)
    Sigma_xy <- Sx_half %*% (sv$u %*% (diag(clipped, nrow = length(clipped)) %*% t(sv$v))) %*% Sy_half
    dimnames(Sigma_xy) <- dimnames(Sigma_xy_input)
    adjusted <- TRUE
  }
  after <- qtl2mr_joint_psd_diagnostics(Sigma_x, Sigma_y, Sigma_xy,
                                        ridge = ridge, tol = tol, pass_tol = pass_tol)
  status <- if (isTRUE(candidate_psd_valid)) {
    "pass"
  } else if (isTRUE(adjusted) && isTRUE(after$psd_valid) && isTRUE(tiny_violation)) {
    "clipped_tiny"
  } else if (isTRUE(adjusted) && isTRUE(after$psd_valid)) {
    "clipped"
  } else if (isTRUE(adjusted)) {
    "clip_failed"
  } else if (isTRUE(large_violation)) {
    "invalid_large"
  } else {
    "invalid"
  }
  list(
    Sigma_x = Sigma_x,
    Sigma_xy = Sigma_xy,
    Sigma_xy_input = Sigma_xy_input,
    Sigma_y = Sigma_y,
    adjusted = adjusted,
    candidate_psd_valid = candidate_psd_valid,
    psd_valid = isTRUE(after$psd_valid),
    status = status,
    large_violation = large_violation,
    tiny_violation = tiny_violation,
    min_eigenvalue = after$min_eigen,
    max_cross_singular_value = after$smax,
    min_eigen_joint_before = before$min_eigen,
    smax_before = before$smax,
    psd_valid_before = before$psd_valid,
    min_eigen_joint_after = after$min_eigen,
    smax_after = after$smax,
    psd_valid_after = after$psd_valid
  )
}

qtl2mr_profile_set <- function(evaluate, anchors, values, q_min, scale, level = 0.95) {
  cutoff <- q_min + stats::qchisq(level, df = 1)
  residual <- values - cutoff
  roots <- anchors[is.finite(residual) & residual == 0]
  f <- function(alpha) evaluate(alpha) - cutoff
  failed <- FALSE
  for (i in seq_len(length(anchors) - 1L)) {
    ends <- residual[c(i, i + 1L)]
    if (anyNA(ends) || !all(is.finite(ends))) {
      # Infinite tail objectives are legitimate when exposure variance is zero.
      # uniroot accepts an infinite endpoint, but not an undefined objective.
      if (anyNA(ends)) { failed <- TRUE; next }
    }
    if (sign(ends[[1L]]) == sign(ends[[2L]])) next
    if (any(ends == 0)) next
    root <- try(stats::uniroot(f, anchors[c(i, i + 1L)],
      f.lower = ends[[1L]], f.upper = ends[[2L]], tol = 1e-11)$root, silent = TRUE)
    if (inherits(root, "try-error")) failed <- TRUE else roots <- c(roots, root)
  }
  cuts <- sort(unique(c(-pi / 2, roots, pi / 2)))
  intervals <- lapply(seq_len(length(cuts) - 1L), function(i) {
    lo <- cuts[[i]]
    hi <- cuts[[i + 1L]]
    mid <- f((lo + hi) / 2)
    if (!is.finite(mid)) return(NULL)
    if (mid > 0) return(NULL)
    c(lower = if (lo == -pi / 2) -Inf else scale * tan(lo),
      upper = if (hi == pi / 2) Inf else scale * tan(hi))
  })
  intervals <- Filter(Negate(is.null), intervals)
  # A threshold-touching extremum may split two adjacent accepted intervals,
  # or contribute an isolated point. Preserve the closed-set interpretation.
  for (root in roots) {
    theta <- scale * tan(root)
    covered <- any(vapply(intervals, function(z) theta >= z[[1L]] && theta <= z[[2L]], logical(1)))
    if (!covered) intervals <- c(intervals, list(c(lower = theta, upper = theta)))
  }
  components <- if (length(intervals)) as.data.frame(do.call(rbind, intervals)) else {
    data.frame(lower = numeric(), upper = numeric())
  }
  if (nrow(components) > 1L) {
    components <- components[order(components$lower, components$upper), , drop = FALSE]
    merged <- components[1L, , drop = FALSE]
    for (i in seq.int(2L, nrow(components))) {
      last <- nrow(merged)
      if (components$lower[[i]] <= merged$upper[[last]]) {
        merged$upper[[last]] <- max(merged$upper[[last]], components$upper[[i]])
      } else merged <- rbind(merged, components[i, , drop = FALSE])
    }
    rownames(merged) <- NULL
    components <- merged
  }
  type <- if (failed || !nrow(components)) "unresolved" else if (nrow(components) > 1L) {
    "disconnected"
  } else if (identical(unname(unlist(components)), c(-Inf, Inf))) {
    "all_real"
  } else if (any(!is.finite(unlist(components)))) "unbounded" else "bounded"
  list(components = components, type = type, cutoff = cutoff)
}
