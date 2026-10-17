qtl2mr_weighted_median <- function(x, w) {
  if (length(x) != length(w)) stop("x and w must have equal lengths", call. = FALSE)
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) stop("No valid weighted median inputs", call. = FALSE)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / max(w)
  w <- w / sum(w)
  x[which(cumsum(w) >= 0.5)[1L]]
}

qtl2mr_weighted_mode <- function(x, w, bandwidth = NULL, grid_step = 1 / 8) {
  if (length(x) != length(w)) stop("x and w must have equal lengths", call. = FALSE)
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) < 2L) stop("Need at least two valid weighted mode inputs", call. = FALSE)
  w <- w / max(w)
  w <- w / sum(w)
  if ((!is.null(bandwidth) && (length(bandwidth) != 1L || !is.numeric(bandwidth) ||
      !is.finite(bandwidth) || bandwidth <= 0)) ||
      length(grid_step) != 1L || !is.finite(grid_step) || grid_step <= 0 || grid_step > 1 / 8) {
    stop("bandwidth must be positive and grid_step in (0, 1/8]", call. = FALSE)
  }
  if (all(x == x[[1L]])) return(x[[1L]])
  # Search in dimensionless coordinates: bandwidth and optimizer tolerance must
  # not introduce an absolute effect-unit floor. Scaling before centering also
  # avoids overflow when valid effects have very large magnitudes.
  unit <- max(abs(x))
  scaled <- x / unit
  origin <- scaled[[1L]]
  spread <- max(abs(scaled - origin))
  x <- (scaled - origin) / spread
  bandwidth <- if (is.null(bandwidth)) stats::bw.nrd0(x) else (bandwidth / unit) / spread
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("Mode bandwidth is not representable at the input scale", call. = FALSE)
  # Search near every component, not on a fixed-size grid spanning the extremes.
  # Evaluate the full mixture, then refine every detected local maximum.
  grid <- sort(unique(as.numeric(outer(x, seq(-8, 8, by = grid_step) * bandwidth, `+`))))
  if (length(grid) < 3L || any(!is.finite(grid))) stop("Mode search grid is not finite", call. = FALSE)
  density <- function(at) {
    as.vector(exp(-0.5 * (outer(at, x, `-`) / bandwidth)^2) %*% w) /
      (sqrt(2 * pi) * bandwidth)
  }
  dens <- density(grid)
  interior <- seq.int(2L, length(grid) - 1L)
  peaks <- interior[dens[interior] >= dens[interior - 1L] &
                     dens[interior] >= dens[interior + 1L] &
                     (dens[interior] > dens[interior - 1L] | dens[interior] > dens[interior + 1L])]
  if (!length(peaks)) stop("No finite mode peak found", call. = FALSE)
  slope <- function(at) {
    delta <- (x - at) / bandwidth
    sum(delta * exp(-0.5 * delta^2) * w)
  }
  candidates <- vapply(peaks, function(i) {
    bracket <- grid[c(i - 1L, i + 1L)]
    ends <- vapply(bracket, slope, numeric(1))
    if (all(is.finite(ends)) && ends[[1L]] >= 0 && ends[[2L]] <= 0) {
      # Root refinement avoids a flat density maximum limiting optimize() to
      # square-root-machine-precision accuracy after a change of effect units.
      stats::uniroot(slope, bracket, f.lower = ends[[1L]], f.upper = ends[[2L]],
        tol = 1e-12 * bandwidth)$root
    } else stats::optimize(function(at) -density(at), bracket,
      tol = 1e-10 * bandwidth)$minimum
  }, numeric(1))
  heights <- density(candidates)
  tied <- which(heights >= max(heights) * (1 - 1e-10))
  # Equal-height modes have no unique mathematical answer. Resolve numerical
  # ties by proximity to input components in their supplied order, preserving
  # reflection equivariance rather than always selecting the leftmost peak.
  for (anchor in x) {
    distance <- abs(candidates[tied] - anchor)
    tied <- tied[distance <= min(distance) + 1e-10 * bandwidth]
    if (length(tied) == 1L) break
  }
  (origin + candidates[tied[[1L]]] * spread) * unit
}

qtl2mr_bootstrap_scalar_sensitivities <- function(theta,
                                                  se,
                                                  B = 500L,
                                                  seed = 1L,
                                                  level = 0.95,
                                                  mode_bandwidth = NULL) {
  B <- B %||% 0L
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) ||
      B < 0 || B > .Machine$integer.max || B != trunc(B) || B == 1) {
    stop("B must be 0 or an integer of at least 2", call. = FALSE)
  }
  B <- as.integer(B)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be a finite probability between 0 and 1", call. = FALSE)
  }
  if (length(theta) != length(se) || any(!is.finite(theta)) ||
      any(!is.finite(se) | se <= 0)) {
    stop("theta and se must have matching finite inputs, with positive se", call. = FALSE)
  }
  w <- 1 / se^2
  if (any(!is.finite(w) | w <= 0)) stop("Inverse-variance weights must be finite and positive", call. = FALSE)
  out <- data.frame(
    estimator = c("Weighted median", "Weighted mode"),
    se = c(NA_real_, NA_real_),
    ci_lower = c(NA_real_, NA_real_),
    ci_upper = c(NA_real_, NA_real_),
    ci_type = c(NA_character_, NA_character_),
    ci_level = rep(level, 2L),
    bootstrap_B = c(max(B, 0L), max(B, 0L)),
    bootstrap_seed = c(seed %||% NA_integer_, seed %||% NA_integer_),
    bootstrap_complete_B = c(0L, 0L),
    bootstrap_failed_B = c(0L, 0L),
    note = c(
      "Diagnostic sensitivity estimate; bootstrap CI not requested or unavailable.",
      "Diagnostic sensitivity estimate; bootstrap CI not requested or unavailable."
    ),
    stringsAsFactors = FALSE
  )
  if (B <= 0L || length(theta) < 2L) {
    out$bootstrap_B <- 0L
    if (length(theta) < 2L) {
      out$note <- "Diagnostic sensitivity estimate; bootstrap CI requires at least two valid loci."
    }
    return(out)
  }
  point <- c(qtl2mr_weighted_median(theta, w),
             qtl2mr_weighted_mode(theta, w, bandwidth = mode_bandwidth))
  z <- stats::qnorm((1 + level) / 2)
  boot <- qtl2mr_with_seed(seed, {
    median <- rep(NA_real_, B)
    mode <- rep(NA_real_, B)
    n <- length(theta)
    for (b in seq_len(B)) {
      tb <- theta + se * stats::rnorm(n)
      if (any(!is.finite(tb))) next
      median[[b]] <- qtl2mr_weighted_median(tb, w)
      mode[[b]] <- tryCatch(
        qtl2mr_weighted_mode(tb, w, bandwidth = mode_bandwidth),
        error = function(e) NA_real_
      )
    }
    data.frame(median = median, mode = mode)
  })
  fill <- function(row, values) {
    finite <- is.finite(values)
    out$bootstrap_complete_B[[row]] <<- sum(finite)
    out$bootstrap_failed_B[[row]] <<- sum(!finite)
    if (!all(finite)) {
      out$note[[row]] <<- sprintf(
        "Diagnostic CI unavailable: %d/%d parametric-bootstrap draws failed; no draws replaced or omitted from uncertainty calculation.",
        sum(!finite), B)
      return(invisible(NULL))
    }
    scale <- max(abs(values))
    sd_boot <- if (scale == 0) 0 else stats::sd(values / scale) * scale
    lower <- point[[row]] - z * sd_boot
    upper <- point[[row]] + z * sd_boot
    if (any(!is.finite(c(sd_boot, lower, upper)))) {
      out$note[[row]] <<- sprintf(
        "Diagnostic CI unavailable: non-finite numerical SD or bounds despite %d/%d finite draws.", B, B)
      return(invisible(NULL))
    }
    out$se[[row]] <<- sd_boot
    out$ci_lower[[row]] <<- lower
    out$ci_upper[[row]] <<- upper
    out$ci_type[[row]] <<- "bootstrap-SE Wald"
    out$note[[row]] <<- sprintf(
      paste0("Approximate %.0f%% normal CI: estimate +/- %.6g * parametric-bootstrap SD (%d draws). ",
        "Fixed loci and inverse-variance weights; independent normal draws of per-locus effects using their SEs. Not a profile CI or coverage guarantee."),
      100 * level, z, B
    )
    out[row, , drop = FALSE]
  }
  invisible(fill(1L, boot$median))
  invisible(fill(2L, boot$mode))
  out
}

fit_per_locus_so <- function(locus_blocks,
                             allow_cross_cov = TRUE,
                             ridge = 1e-10) {
  if (!is.list(locus_blocks) || !length(locus_blocks)) {
    stop("locus_blocks must be a non-empty list", call. = FALSE)
  }
  rows <- lapply(seq_along(locus_blocks), function(i) {
    b <- locus_blocks[[i]]
    d <- mr_block_so_diagnostics(b, allow_cross_cov = allow_cross_cov, ridge = ridge)
    data.frame(
      locus_index = i,
      marker = d$marker,
      chr = d$chr,
      theta = d$theta,
      se = d$se,
      profile_ci_lower = d$profile_ci_lower,
      profile_ci_upper = d$profile_ci_upper,
      Q_min = d$Q_min,
      Q_df = d$Q_df,
      Q_p_chisq = d$Q_p_chisq,
      F_X_block = d$F_X_block,
      omega_min_eigenvalue = d$omega_min_eigenvalue,
      diagnostic_status = d$diagnostic_status,
      diagnostic_note = d$diagnostic_note,
      cross_cov_model = b$cross_cov_model %||% NA_character_,
      rho_xy = b$rho_xy %||% NA_real_,
      joint_psd_valid = b$joint_psd_valid %||% NA,
      joint_max_cross_singular_value = b$joint_max_cross_singular_value %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$weight <- 1 / out$se^2
  out
}

summarize_per_locus_so <- function(per_locus,
                                   level = 0.95,
                                   mode_bandwidth = NULL,
                                   bootstrap_B = 500L,
                                   bootstrap_seed = 1L) {
  required <- c("theta", "se")
  miss <- setdiff(required, names(per_locus))
  if (length(miss)) stop("per_locus is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  ok <- is.finite(per_locus$theta) & is.finite(per_locus$se) & per_locus$se > 0
  if (!any(ok)) {
    estimates <- data.frame(estimator = c("Two-step IVW", "Weighted median", "Weighted mode"),
      theta = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      ci_type = NA_character_, ci_level = level, bootstrap_B = 0L,
      bootstrap_seed = bootstrap_seed %||% NA_integer_, bootstrap_complete_B = 0L,
      bootstrap_failed_B = 0L, role = c("diagnostic", "sensitivity", "sensitivity"),
      note = "Diagnostic unavailable: no finite per-locus SO estimate with positive SE.",
      stringsAsFactors = FALSE)
    return(list(estimates = estimates, heterogeneity = data.frame(Q_between = NA_real_,
      df_between = NA_integer_, Q_between_p = NA_real_, n_loci = 0L), per_locus = per_locus))
  }
  theta <- per_locus$theta[ok]
  se <- per_locus$se[ok]
  w <- 1 / se^2
  z <- stats::qnorm((1 + level) / 2)
  theta_ivw <- sum(w * theta) / sum(w)
  se_ivw <- sqrt(1 / sum(w))
  q_between <- sum(w * (theta - theta_ivw)^2)
  df_between <- length(theta) - 1L
  boot <- qtl2mr_bootstrap_scalar_sensitivities(
    theta,
    se,
    B = bootstrap_B,
    seed = bootstrap_seed,
    level = level,
    mode_bandwidth = mode_bandwidth
  )
  theta_mode <- if (length(theta) >= 2L) qtl2mr_weighted_mode(theta, w, bandwidth = mode_bandwidth) else NA_real_
  median_boot <- boot[boot$estimator == "Weighted median", , drop = FALSE]
  mode_boot <- boot[boot$estimator == "Weighted mode", , drop = FALSE]
  estimates <- data.frame(
    estimator = c("Two-step IVW", "Weighted median", "Weighted mode"),
    theta = c(
      theta_ivw,
      qtl2mr_weighted_median(theta, w),
      theta_mode
    ),
    se = c(se_ivw, median_boot$se[[1L]], mode_boot$se[[1L]]),
    ci_lower = c(theta_ivw - z * se_ivw, median_boot$ci_lower[[1L]], mode_boot$ci_lower[[1L]]),
    ci_upper = c(theta_ivw + z * se_ivw, median_boot$ci_upper[[1L]], mode_boot$ci_upper[[1L]]),
    ci_type = c("Wald", median_boot$ci_type[[1L]], mode_boot$ci_type[[1L]]),
    ci_level = rep(level, 3L),
    bootstrap_B = c(0L, median_boot$bootstrap_B[[1L]], mode_boot$bootstrap_B[[1L]]),
    bootstrap_seed = c(NA_integer_, median_boot$bootstrap_seed[[1L]], mode_boot$bootstrap_seed[[1L]]),
    bootstrap_complete_B = c(0L, median_boot$bootstrap_complete_B[[1L]], mode_boot$bootstrap_complete_B[[1L]]),
    bootstrap_failed_B = c(0L, median_boot$bootstrap_failed_B[[1L]], mode_boot$bootstrap_failed_B[[1L]]),
    role = c("diagnostic", "sensitivity", "sensitivity"),
    note = c(
      "Two-step IVW diagnostic comparator; Wald CI uses inverse-variance weights across per-locus SO estimates.",
      median_boot$note[[1L]],
      mode_boot$note[[1L]]
    ),
    stringsAsFactors = FALSE
  )
  heterogeneity <- data.frame(
    Q_between = q_between,
    df_between = df_between,
    Q_between_p = stats::pchisq(q_between, df = df_between, lower.tail = FALSE),
    n_loci = length(theta),
    stringsAsFactors = FALSE
  )
  list(estimates = estimates, heterogeneity = heterogeneity, per_locus = per_locus)
}

fit_multilocus_diagnostic_layer <- function(locus_blocks,
                                            global_fit = NULL,
                                            allow_cross_cov = TRUE,
                                            ridge = 1e-10,
                                            level = 0.95,
                                            bootstrap_B = 500L,
                                            bootstrap_seed = 1L) {
  if (is.null(global_fit)) {
    global_fit <- fit_multilocus_founder_mr(
      locus_blocks,
      allow_cross_cov = allow_cross_cov,
      ridge = ridge,
      heterogeneity_method = "chisq"
    )
  }
  per_locus <- fit_per_locus_so(
    locus_blocks,
    allow_cross_cov = allow_cross_cov,
    ridge = ridge
  )
  scalar <- summarize_per_locus_so(
    per_locus,
    level = level,
    bootstrap_B = bootstrap_B,
    bootstrap_seed = bootstrap_seed
  )
  global_row <- data.frame(
    estimator = "Global SO",
    theta = global_fit$theta,
    se = global_fit$se,
    ci_lower = global_fit$profile_ci_lower,
    ci_upper = global_fit$profile_ci_upper,
    ci_type = "profile",
    ci_level = 0.95,
    bootstrap_B = 0L,
    bootstrap_seed = NA_integer_,
    bootstrap_complete_B = 0L,
    bootstrap_failed_B = 0L,
    role = "primary",
    note = "Primary overlap-aware SO Founder GLS-IVW estimate; CI is profile-based.",
    stringsAsFactors = FALSE
  )
  scalar$estimates <- rbind(global_row, scalar$estimates)
  total_q <- data.frame(
    Q_global_min = global_fit$Q,
    Q_global_p_reported = global_fit$Q_p_reported %||% NA_real_,
    Q_global_p_method = global_fit$Q_p_method %||% NA_character_,
    Q_global_p_note = global_fit$Q_p_note %||% NA_character_,
    Q_global_p_chisq = global_fit$Q_p_chisq %||% NA_real_,
    Q_global_chisq_df = global_fit$Q_chisq_df %||% NA_integer_,
    Q_global_p_sim = global_fit$Q_p_sim %||% NA_real_,
    Q_sim_n = global_fit$Q_sim_n %||% NA_integer_,
    Q_sim_exceedances = global_fit$Q_sim_exceedances %||% NA_integer_,
    Q_sim_seed = global_fit$Q_sim_seed %||% NA_integer_,
    interpretation = "reported p-value follows heterogeneity_method: auto uses chi-square for single-locus and simulation for multi-locus",
    stringsAsFactors = FALSE
  )
  list(
    global = global_fit,
    per_locus = per_locus,
    estimates = scalar$estimates,
    heterogeneity = scalar$heterogeneity,
    global_lack_of_fit = total_q
  )
}
