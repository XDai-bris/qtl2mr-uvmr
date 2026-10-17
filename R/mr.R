gls_ivw <- function(beta_x, beta_y, Sigma_y, ridge = 0) {
  beta_x <- matrix(as.numeric(beta_x), ncol = 1L)
  beta_y <- matrix(as.numeric(beta_y), ncol = 1L)
  Sigma_y <- qtl2mr_sym(Sigma_y)
  ok <- is.finite(beta_x[, 1L]) & is.finite(beta_y[, 1L]) &
    is.finite(diag(Sigma_y)) & diag(Sigma_y) > 0
  beta_x <- beta_x[ok, , drop = FALSE]
  beta_y <- beta_y[ok, , drop = FALSE]
  Sigma_y <- Sigma_y[ok, ok, drop = FALSE]
  if (length(beta_x) < 2L) stop("Need at least two valid founder contrasts", call. = FALSE)
  if (ridge > 0) Sigma_y <- qtl2mr_psd_clip(Sigma_y, ridge = ridge)
  W <- solve(Sigma_y)
  denom <- as.numeric(t(beta_x) %*% W %*% beta_x)
  numer <- as.numeric(t(beta_x) %*% W %*% beta_y)
  theta <- numer / denom
  se <- sqrt(1 / denom)
  z <- theta / se
  resid <- beta_y - theta * beta_x
  Q <- as.numeric(t(resid) %*% W %*% resid)
  Q_df <- length(beta_x) - 1L
  Q_p <- stats::pchisq(Q, df = Q_df, lower.tail = FALSE)
  list(
    beta = theta,
    theta = theta,
    se = se,
    z = z,
    p = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    Q = Q,
    Q_df = Q_df,
    Q_chisq_df = Q_df,
    Q_p = Q_p,
    Q_p_chisq = Q_p
  )
}

# Trait-wise scaling and the angular parameterization are computational changes
# only: Q(theta) is the same second-order quadratic objective. The scale uses
# Euclidean norms and traces, so orthonormal founder-basis changes do not affect it.
qtl2mr_so_coordinates <- function(blocks, allow_cross_cov, ridge) {
  if (!is.list(blocks) || !length(blocks)) stop("locus_blocks must be a non-empty list", call. = FALSE)
  if (length(ridge) != 1L || !typeof(ridge) %in% c("double", "integer") ||
      is.object(ridge) || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be one finite non-negative number", call. = FALSE)
  }
  blocks <- lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    bx <- b$beta_x
    by <- b$beta_y
    if (!is.numeric(bx) || !is.numeric(by) || length(bx) < 2L ||
        length(by) != length(bx) || any(!is.finite(c(bx, by)))) {
      stop("Each locus block needs at least two finite numeric beta_x/beta_y contrasts", call. = FALSE)
    }
    k <- length(bx)
    get_cov <- function(S, label) {
      if (!is.matrix(S) || !is.numeric(S) || !identical(dim(S), c(k, k)) ||
          any(!is.finite(S))) {
        stop(label, " must be a finite numeric matrix matching beta lengths", call. = FALSE)
      }
      S
    }
    Sx <- get_cov(as.matrix(b$Sigma_x), "Sigma_x")
    Sy <- get_cov(as.matrix(b$Sigma_y), "Sigma_y")
    Sxy <- if (isTRUE(allow_cross_cov)) get_cov(as.matrix(b$Sigma_xy), "Sigma_xy") else matrix(0, k, k)
    # Validate once, before the repeated objective evaluations.
    qtl2mr_add_invariant_ridge(Sx, 0)
    qtl2mr_add_invariant_ridge(Sy, 0)
    if (ridge == 0 && (inherits(try(chol(Sx), silent = TRUE), "try-error") ||
                      inherits(try(chol(Sy), silent = TRUE), "try-error"))) {
      stop("With ridge = 0, Sigma_x and Sigma_y must be positive definite; singular covariance limits are not supported", call. = FALSE)
    }
    list(beta_x = as.numeric(bx), beta_y = as.numeric(by),
         Sigma_x = Sx, Sigma_y = Sy, Sigma_xy = Sxy)
  })
  trait_scale <- function(beta_name, sigma_name) {
    terms <- unlist(lapply(blocks, function(b) c(b[[beta_name]], sqrt(diag(b[[sigma_name]])))), use.names = FALSE)
    peak <- max(abs(terms))
    if (peak == 0) return(1)
    peak * sqrt(sum((terms / peak)^2) / length(terms))
  }
  sx <- trait_scale("beta_x", "Sigma_x")
  sy <- trait_scale("beta_y", "Sigma_y")
  scale <- sy / sx
  if (!is.finite(scale) || scale <= 0) stop("Trait scales exceed the supported floating-point range", call. = FALSE)
  normalized <- lapply(blocks, function(b) {
    b$beta_x <- b$beta_x / sx
    b$beta_y <- b$beta_y / sy
    b$Sigma_x <- (b$Sigma_x / sx) / sx
    b$Sigma_y <- (b$Sigma_y / sy) / sy
    b$Sigma_xy <- (b$Sigma_xy / sx) / sy
    b$cross <- b$Sigma_xy + t(b$Sigma_xy)
    joint <- rbind(cbind(b$Sigma_x, b$Sigma_xy), cbind(t(b$Sigma_xy), b$Sigma_y))
    qtl2mr_add_invariant_ridge(joint, 0)
    b
  })
  evaluate <- function(alpha, derivatives = FALSE) {
    at_tail <- abs(alpha) >= pi / 2
    cc <- if (at_tail) 0 else cos(alpha)
    ss <- if (at_tail) sign(alpha) else sin(alpha)
    total <- c(Q = 0, first = 0, second = 0)
    for (b in normalized) {
      # If exposure sampling variance is identically zero, the leading
      # homogeneous term vanishes. Evaluate its actual limit separately.
      if (at_tail && all(b$Sigma_x == 0)) {
        if (any(b$beta_x != 0)) return(if (derivatives) c(Q = Inf, first = NA, second = NA) else Inf)
        cf <- chol(qtl2mr_add_invariant_ridge(b$Sigma_y, ridge))
        z <- backsolve(cf, b$beta_y, transpose = TRUE)
        total[["Q"]] <- total[["Q"]] + sum(z^2)
        next
      }
      S <- cc^2 * b$Sigma_y + ss^2 * b$Sigma_x - ss * cc * b$cross
      S <- qtl2mr_sym(S)
      # Normalized matrices have ordinary numerical scale. This trace ridge is
      # identical to the existing invariant ridge wherever trace(S) is positive.
      S <- S + diag(ridge * mean(diag(S)), nrow(S))
      cf <- try(chol(S), silent = TRUE)
      if (inherits(cf, "try-error")) {
        return(if (derivatives) c(Q = Inf, first = NA, second = NA) else Inf)
      }
      solve_S <- function(v) backsolve(cf, backsolve(cf, v, transpose = TRUE))
      r <- cc * b$beta_y - ss * b$beta_x
      z <- solve_S(r)
      total[["Q"]] <- total[["Q"]] + sum(r * z)
      if (derivatives) {
        rp <- -ss * b$beta_y - cc * b$beta_x
        Sp <- 2 * ss * cc * (b$Sigma_x - b$Sigma_y) - (cc^2 - ss^2) * b$cross
        Spp <- 2 * (cc^2 - ss^2) * (b$Sigma_x - b$Sigma_y) + 4 * ss * cc * b$cross
        Sp <- Sp + diag(ridge * mean(diag(Sp)), nrow(Sp))
        Spp <- Spp + diag(ridge * mean(diag(Spp)), nrow(Spp))
        v <- rp - as.numeric(Sp %*% z)
        total[["first"]] <- total[["first"]] + 2 * sum(rp * z) - sum(z * (Sp %*% z))
        total[["second"]] <- total[["second"]] - 2 * sum(r * z) +
          2 * sum(v * solve_S(v)) - sum(z * (Spp %*% z))
      }
    }
    if (derivatives) total else unname(total[["Q"]])
  }
  list(blocks = normalized, scale = scale, evaluate = evaluate)
}

qtl2mr_multilocus_objective <- function(blocks, allow_cross_cov = TRUE, ridge = 1e-10) {
  coordinates <- qtl2mr_so_coordinates(blocks, allow_cross_cov, ridge)
  function(theta) coordinates$evaluate(atan(theta / coordinates$scale))
}

qtl2mr_so_fit <- function(locus_blocks, allow_cross_cov = TRUE, ridge = 1e-10,
                         theta_grid = NULL) {
  coordinates <- qtl2mr_so_coordinates(locus_blocks, allow_cross_cov, ridge)
  scale <- coordinates$scale
  nonfinite_interior <- FALSE
  evaluate <- function(alpha, derivatives = FALSE) {
    value <- coordinates$evaluate(alpha, derivatives)
    q <- if (derivatives) value[["Q"]] else value
    if (abs(alpha) < pi / 2 && !is.finite(q)) nonfinite_interior <<- TRUE
    value
  }
  # optimize() cannot consume Inf without substituting its own value and warning.
  # A penalty is only an exploration device; all reported Q values are evaluated
  # with the true objective and any interior failure makes the result unresolved.
  optimizer_value <- function(alpha, direction = 1) {
    q <- evaluate(alpha)
    if (is.finite(q)) direction * q else .Machine$double.xmax / 100
  }
  objective <- function(theta) evaluate(atan(theta / scale))
  blocks <- coordinates$blocks
  bx <- unlist(lapply(blocks, `[[`, "beta_x"), use.names = FALSE)
  by <- unlist(lapply(blocks, `[[`, "beta_y"), use.names = FALSE)
  Sy <- block_diag_matrix(lapply(blocks, `[[`, "Sigma_y"))
  first <- try(gls_ivw(bx, by, Sy, ridge = ridge), silent = TRUE)
  if (inherits(first, "try-error")) first <- list(theta = NA_real_, se = NA_real_)
  first$theta <- first$theta * scale
  first$se <- first$se * scale
  if (!is.null(theta_grid) && (!is.numeric(theta_grid) || !length(theta_grid) ||
                              any(!is.finite(theta_grid)))) {
    stop("theta_grid must contain finite numeric supplemental search points", call. = FALSE)
  }
  # theta_grid supplements the whole-real-line search; it is never a bound.
  hints <- c(theta_grid, first$theta,
             first$theta + first$se * c(-4, -1, 1, 4))
  hints <- atan(hints[is.finite(hints)] / scale)
  grid <- sort(unique(c(seq(-pi / 2, pi / 2, length.out = 257L), hints)))
  values <- vapply(grid, evaluate, numeric(1))
  if (!any(is.finite(values))) stop("The second-order objective is non-finite over the search grid", call. = FALSE)
  minima <- maxima <- numeric()
  for (i in seq.int(2L, length(grid) - 1L)) {
    yi <- values[[i]]
    near <- values[c(i - 1L, i + 1L)]
    if (!is.finite(yi)) next
    if (yi <= min(near) && yi < max(near)) {
      opt <- stats::optimize(optimizer_value, grid[c(i - 1L, i + 1L)], tol = 1e-11)
      a <- opt$minimum
      # An analytic Newton correction avoids unit- and basis-dependent
      # finite-difference curvature and optimizer stopping effects.
      for (j in seq_len(3L)) {
        d <- evaluate(a, derivatives = TRUE)
        step <- d[["first"]] / d[["second"]]
        candidate <- a - step
        if (!is.finite(step) || d[["second"]] <= 0 ||
            candidate <= grid[[i - 1L]] || candidate >= grid[[i + 1L]]) break
        if (evaluate(candidate) > d[["Q"]] + 1e-10 * max(1, abs(d[["Q"]]))) break
        a <- candidate
      }
      minima <- c(minima, a)
    }
    if (yi >= max(near) && yi > min(near) && all(is.finite(near))) {
      maxima <- c(maxima, stats::optimize(function(a) optimizer_value(a, -1),
        grid[c(i - 1L, i + 1L)], tol = 1e-11)$minimum)
    }
  }
  # Tail-adjacent minima can lie between the endpoint and the first grid point.
  for (ii in list(c(1L, 2L), c(length(grid) - 1L, length(grid)))) {
    opt <- stats::optimize(optimizer_value, grid[ii], tol = 1e-11)
    if (opt$objective < min(values[ii])) minima <- c(minima, opt$minimum)
  }
  anchors <- sort(unique(c(grid, minima, maxima)))
  q_anchor <- vapply(anchors, evaluate, numeric(1))
  best <- which.min(q_anchor)
  alpha_hat <- anchors[[best]]
  q_min <- q_anchor[[best]]
  tails <- c(negative = values[[1L]], positive = values[[length(values)]])
  at_infinity <- best %in% c(1L, length(anchors)) ||
    min(tails) <= q_min + 1e-10 * max(1, abs(q_min))
  theta_hat <- if (at_infinity) NA_real_ else scale * tan(alpha_hat)
  curvature <- if (at_infinity) NA_real_ else evaluate(alpha_hat, derivatives = TRUE)[["second"]]
  se <- if (is.finite(curvature) && curvature > 0) {
    scale / cos(alpha_hat)^2 * sqrt(2 / curvature)
  } else NA_real_
  confidence <- qtl2mr_profile_set(evaluate, anchors, q_anchor, q_min, scale)
  evaluation_ok <- !nonfinite_interior && all(is.finite(q_anchor[-c(1L, length(q_anchor))]))
  if (!evaluation_ok) confidence$type <- "unresolved"
  ci <- confidence$components
  single <- nrow(ci) == 1L && confidence$type != "unresolved"
  finite_ci <- single && all(is.finite(unlist(ci)))
  numerical_ok <- is.finite(q_min) && confidence$type != "unresolved"
  reliable <- numerical_ok && !at_infinity && finite_ci &&
    is.finite(theta_hat) && is.finite(se)
  status <- if (!numerical_ok) "failed" else if (at_infinity) {
    "minimum_at_infinity"
  } else if (confidence$type == "disconnected") {
    "disconnected_profile"
  } else if (!finite_ci) {
    "unbounded_profile"
  } else if (!is.finite(se)) "nonpositive_curvature" else "success"
  warnings <- "Finite angular-grid search with local refinement is not a proof of global optimality."
  if (!evaluation_ok) warnings <- c(warnings, "A non-finite interior objective was encountered; the confidence set is unresolved.")
  if (!finite_ci) warnings <- c(warnings, paste0("Profile confidence set is ", confidence$type, "; do not report a finite interval."))
  if (at_infinity) warnings <- c(warnings, "A finite optimum is not numerically distinguishable from an infinite-effect limit; no finite theta is reported.")
  df <- max(length(bx) - 1L, 1L)
  # Compatibility bracket fields describe the full searched domain.
  list(
    theta = theta_hat, beta = theta_hat, se = se,
    z = theta_hat / se,
    p = 2 * stats::pnorm(abs(theta_hat / se), lower.tail = FALSE),
    profile_p = stats::pchisq(max(objective(0) - q_min, 0), df = 1, lower.tail = FALSE),
    profile_ci_lower = if (single) ci$lower[[1L]] else NA_real_,
    profile_ci_upper = if (single) ci$upper[[1L]] else NA_real_,
    confidence_set = ci, confidence_set_type = confidence$type,
    confidence_cutoff = confidence$cutoff,
    Q = q_min, Q_df = df, Q_chisq_df = df,
    Q_p_chisq = stats::pchisq(q_min, df = df, lower.tail = FALSE),
    first_order_theta = first$theta, first_order_se = first$se,
    theta_start = first$theta,
    bracket_lower = -Inf, bracket_upper = Inf,
    bracket_status = "whole_real_line", bracket_n_expand = 0L,
    optimizer_status = status, profile_ci_finite = finite_ci,
    numerical_ok = numerical_ok,
    primary_reliable = reliable,
    reliability_scope = "numerically_bounded_single_component_only",
    inference_status = "not_calibrated",
    objective_min = q_min, objective_tail = tails,
    optimizer_warnings = warnings,
    search_method = "scaled_arctan_grid_with_local_extrema_refinement",
    search_grid_n = length(grid), search_local_minima_n = length(minima),
    search_nonfinite_interior = nonfinite_interior,
    search_certified = FALSE, theta_grid_role = "supplemental_search_points",
    lower = -Inf, upper = Inf,
    objective = objective, converged = numerical_ok && !at_infinity,
    allow_cross_cov = isTRUE(allow_cross_cov)
  )
}

second_order_gls_ivw <- function(beta_x, beta_y, Sigma_x, Sigma_y,
                                 theta_grid = NULL, allow_cross_cov = FALSE,
                                 Sigma_xy = NULL, ridge = 1e-10) {
  if (allow_cross_cov && is.null(Sigma_xy)) {
    stop("Sigma_xy is required when allow_cross_cov = TRUE", call. = FALSE)
  }
  block <- list(beta_x = beta_x, beta_y = beta_y,
                Sigma_x = Sigma_x, Sigma_y = Sigma_y, Sigma_xy = Sigma_xy)
  qtl2mr_so_fit(list(block), allow_cross_cov, ridge, theta_grid)
}

qtl2mr_fit_multilocus_once <- function(locus_blocks, allow_cross_cov = TRUE,
                                       ridge = 1e-10, theta_grid = NULL) {
  qtl2mr_so_fit(locus_blocks, allow_cross_cov, ridge, theta_grid)
}

qtl2mr_rmvnorm1 <- function(mu, Sigma, ridge = 1e-10) {
  mu <- as.numeric(mu)
  Sigma <- qtl2mr_psd_clip(Sigma, ridge = ridge)
  ee <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  as.numeric(mu + ee$vectors %*% (sqrt(vals) * stats::rnorm(length(mu))))
}

calibrate_multilocus_q <- function(locus_blocks,
                                   fit = NULL,
                                   nsim = 1000L,
                                   seed = 1L,
                                   allow_cross_cov = TRUE,
                                   ridge = 1e-10,
                                   theta_grid = NULL,
                                   keep_simulated = FALSE,
                                   progress_callback = NULL) {
  nsim <- as.integer(nsim %||% 0L)
  if (!is.finite(nsim) || nsim < 1L) {
    return(list(
      Q_p_sim = NA_real_,
      Q_sim_n = 0L,
      Q_sim_n_requested = 0L,
      Q_sim_seed = NA_integer_,
      Q_sim_exceedances = NA_integer_,
      Q_sim_min_resolvable_p = NA_real_,
      Q_sim_observed = if (!is.null(fit)) fit$Q else NA_real_,
      Q_sim_method = "not run",
      Q_sim_values = NULL
    ))
  }
  if (is.null(fit)) {
    fit <- qtl2mr_fit_multilocus_once(
      locus_blocks = locus_blocks,
      allow_cross_cov = allow_cross_cov,
      ridge = ridge,
      theta_grid = theta_grid
    )
  }
  theta0 <- fit$theta
  q_obs <- fit$Q
  if (!is.finite(theta0) || !is.finite(q_obs)) {
    return(list(
      Q_p_sim = NA_real_,
      Q_sim_n = 0L,
      Q_sim_n_requested = nsim,
      Q_sim_seed = as.integer(seed),
      Q_sim_exceedances = NA_integer_,
      Q_sim_min_resolvable_p = NA_real_,
      Q_sim_observed = q_obs,
      Q_sim_method = "not available because theta or Q is non-finite",
      Q_sim_values = NULL
    ))
  }
  qtl2mr_with_seed(seed, {
    sim_info <- lapply(locus_blocks, function(b) {
      bx <- as.numeric(b$beta_x)
      Sxy <- if (isTRUE(allow_cross_cov)) b$Sigma_xy else matrix(0, length(bx), length(bx))
      Om <- qtl2mr_add_invariant_ridge(
        qtl2mr_sym(b$Sigma_y + theta0^2 * b$Sigma_x - theta0 * (Sxy + t(Sxy))),
        ridge = 0
      )
      ee <- eigen(Om, symmetric = TRUE)
      list(
        block = b,
        beta_x = bx,
        beta_y_names = names(b$beta_y),
        vectors = ee$vectors,
        scales = sqrt(pmax(ee$values, 0))
      )
    })
    q_sim <- numeric(nsim)
    started <- proc.time()[["elapsed"]]
    report_every <- max(1L, ceiling(nsim / 10L))
    qtl2mr_progress(progress_callback, 2L, "Simulate heterogeneity Q",
      sprintf("0/%d draws; each draw refits the common-effect model", nsim))
    for (s in seq_len(nsim)) {
      sim_blocks <- lapply(sim_info, function(info) {
        b <- info$block
        by <- theta0 * info$beta_x +
          as.numeric(info$vectors %*% (info$scales * stats::rnorm(length(info$beta_x))))
        names(by) <- info$beta_y_names
        b$beta_y <- by
        b
      })
      q_sim[[s]] <- qtl2mr_fit_multilocus_once(
        locus_blocks = sim_blocks,
        allow_cross_cov = allow_cross_cov,
        ridge = ridge,
        theta_grid = theta_grid
      )$Q
      if (s %% report_every == 0L || s == nsim) {
        qtl2mr_progress(progress_callback, 2L, "Simulate heterogeneity Q",
          sprintf("%d/%d draws; %.1f s elapsed", s, nsim,
            proc.time()[["elapsed"]] - started))
      }
    }
    exceed <- sum(q_sim >= q_obs, na.rm = TRUE)
    list(
      Q_p_sim = (exceed + 1) / (sum(is.finite(q_sim)) + 1),
      Q_sim_n = sum(is.finite(q_sim)),
      Q_sim_n_requested = nsim,
      Q_sim_seed = as.integer(seed),
      Q_sim_exceedances = as.integer(exceed),
      Q_sim_min_resolvable_p = 1 / (sum(is.finite(q_sim)) + 1),
      Q_sim_observed = q_obs,
      Q_sim_null_mean = mean(q_sim, na.rm = TRUE),
      Q_sim_null_sd = stats::sd(q_sim, na.rm = TRUE),
      Q_sim_method = "conditional residual simulation under the fitted common-effect null",
      Q_sim_values = if (isTRUE(keep_simulated)) q_sim else NULL
    )
  })
}

fit_multilocus_founder_mr <- function(locus_blocks,
                                      allow_cross_cov = TRUE,
                                      ridge = 1e-10,
                                      theta_grid = NULL,
                                      compute_independence = TRUE,
                                      heterogeneity_method = c("auto", "simulation", "chisq", "none"),
                                      q_sim_n = 0L,
                                      q_sim_seed = 1L,
                                      keep_q_sim = FALSE,
                                      progress_callback = NULL) {
  heterogeneity_method <- match.arg(heterogeneity_method)
  if (!is.list(locus_blocks) || !length(locus_blocks)) stop("locus_blocks must be a non-empty list", call. = FALSE)
  needed <- c("beta_x", "beta_y", "Sigma_x", "Sigma_y", "Sigma_xy")
  for (i in seq_along(locus_blocks)) {
    miss <- setdiff(needed, names(locus_blocks[[i]]))
    if (length(miss)) stop("Locus block ", i, " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  fit <- qtl2mr_fit_multilocus_once(
    locus_blocks = locus_blocks,
    allow_cross_cov = allow_cross_cov,
    ridge = ridge,
    theta_grid = theta_grid
  )
  per_locus <- do.call(rbind, lapply(seq_along(locus_blocks), function(i) {
    b <- locus_blocks[[i]]
    data.frame(
      marker = b$marker %||% paste0("locus_", i),
      chr = b$chr %||% NA_character_,
      n_overlap = b$n_overlap %||% NA_integer_,
      n_x = b$O_xy_summary$n_x %||% b$fit_x$n %||% NA_integer_,
      n_y = b$O_xy_summary$n_y %||% b$fit_y$n %||% NA_integer_,
      exposure_overlap_fraction = b$O_xy_summary$exposure_overlap_fraction %||% NA_real_,
      outcome_overlap_fraction = b$O_xy_summary$outcome_overlap_fraction %||% NA_real_,
      jaccard_overlap_fraction = b$O_xy_summary$jaccard_overlap_fraction %||% NA_real_,
      complete_overlap = b$O_xy_summary$complete_overlap %||% NA,
      independent_samples = isTRUE((b$n_overlap %||% NA_integer_) == 0L),
      rho_xy = b$rho_xy %||% NA_real_,
      rho_xy_raw = b$rho_xy_raw %||% NA_real_,
      rho_xy_star = b$rho_xy_star %||% NA_real_,
      Q_locus = if (is.finite(fit$theta)) mr_block_q(b, fit$theta,
        allow_cross_cov = allow_cross_cov, ridge = ridge) else NA_real_,
      cross_cov_model = b$cross_cov_model %||% NA_character_,
      psd_adjusted = b$psd_adjusted %||% FALSE,
      joint_psd_candidate_adjusted = b$joint_psd_candidate_adjusted %||% NA,
      joint_psd_valid = b$joint_psd_valid %||% NA,
      joint_psd_status = b$joint_psd_status %||% NA_character_,
      joint_psd_large_violation = b$joint_psd_large_violation %||% NA,
      joint_min_eigenvalue = b$joint_min_eigenvalue %||% NA_real_,
      joint_max_cross_singular_value = b$joint_max_cross_singular_value %||% NA_real_,
      raw_space_psd_adjusted = b$raw_space_psd_adjusted %||% NA,
      raw_space_psd_valid = b$raw_space_psd_valid %||% NA,
      raw_space_psd_status = b$raw_space_psd_status %||% NA_character_,
      raw_space_psd_large_violation = b$raw_space_psd_large_violation %||% NA,
      raw_space_max_cross_singular_value = b$raw_space_max_cross_singular_value %||% NA_real_,
      star_space_psd_adjusted = b$star_space_psd_adjusted %||% NA,
      star_space_psd_valid = b$star_space_psd_valid %||% NA,
      star_space_psd_status = b$star_space_psd_status %||% NA_character_,
      star_space_psd_large_violation = b$star_space_psd_large_violation %||% NA,
      star_space_max_cross_singular_value = b$star_space_max_cross_singular_value %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  fit$locus_blocks <- locus_blocks
  fit$locus_diagnostics <- per_locus
  if (isTRUE(compute_independence)) {
    fit$independence <- qtl2mr_fit_multilocus_once(
      locus_blocks = locus_blocks,
      allow_cross_cov = FALSE,
      ridge = ridge,
      theta_grid = theta_grid
    )
    fit$theta_overlap_delta <- fit$theta - fit$independence$theta
  }
  reported_method <- if (identical(heterogeneity_method, "auto")) {
    if (length(locus_blocks) == 1L) "chisq" else "simulation"
  } else {
    heterogeneity_method
  }
  q_sim_run_n <- if (identical(reported_method, "simulation")) q_sim_n else 0L
  if (identical(reported_method, "simulation") && as.integer(q_sim_run_n) < 1L) {
    warning(
      "heterogeneity_method selects simulation but q_sim_n = 0; ",
      "Q_p_sim and Q_p_reported will be NA. Set q_sim_n > 0, ",
      "or use heterogeneity_method = 'chisq' or 'none'.",
      call. = FALSE
    )
  }
  q_cal <- calibrate_multilocus_q(
    locus_blocks = locus_blocks,
    fit = fit,
    nsim = q_sim_run_n,
    seed = q_sim_seed,
    allow_cross_cov = allow_cross_cov,
    ridge = ridge,
    theta_grid = theta_grid,
    keep_simulated = keep_q_sim,
    progress_callback = progress_callback
  )
  fit$Q_p_sim <- q_cal$Q_p_sim
  fit$Q_sim_n <- q_cal$Q_sim_n
  fit$Q_sim_n_requested <- q_cal$Q_sim_n_requested
  fit$Q_sim_seed <- q_cal$Q_sim_seed
  fit$Q_sim_exceedances <- q_cal$Q_sim_exceedances
  fit$Q_sim_min_resolvable_p <- q_cal$Q_sim_min_resolvable_p
  fit$Q_sim_observed <- q_cal$Q_sim_observed
  fit$Q_sim_null_mean <- q_cal$Q_sim_null_mean %||% NA_real_
  fit$Q_sim_null_sd <- q_cal$Q_sim_null_sd %||% NA_real_
  fit$Q_sim_method <- q_cal$Q_sim_method
  if (isTRUE(keep_q_sim)) fit$Q_sim_values <- q_cal$Q_sim_values
  fit$Q_p_method <- reported_method
  fit$Q_p_reported <- switch(
    reported_method,
    chisq = fit$Q_p_chisq,
    simulation = fit$Q_p_sim,
    none = NA_real_
  )
  fit$Q_p_note <- switch(
    reported_method,
    chisq = paste0("chi-square reference with df=", fit$Q_chisq_df),
    simulation = if (is.finite(fit$Q_p_sim)) {
      "conditional residual simulation under the fitted common-effect null"
    } else {
      "simulation selected but not run or not available"
    },
    none = "heterogeneity p-value not reported"
  )
  fit
}
