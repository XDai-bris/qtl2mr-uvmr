# 0.6.1.9002: speed-ups must leave every result unchanged. Each test compares
# with identical(): against the earlier computation, between worker counts,
# between a fresh and a warm session, and between batch and single analyses.

perf_fixture <- function(round_digits = NULL) {
  set.seed(62031)
  n <- 90L
  ids <- paste0("s", seq_len(n))
  make_chr <- function(markers) {
    a <- array(rexp(n * 8L * length(markers)), c(n, 8L, length(markers)), list(ids, LETTERS[1:8], markers))
    for (j in seq_along(markers)) a[, , j] <- a[, , j] / rowSums(a[, , j])
    if (!is.null(round_digits)) a <- round(a, round_digits)
    a
  }
  gp <- list(`1` = make_chr(c("a1", "a2", "a3")), `2` = make_chr(c("b1", "b2")), `3` = make_chr(c("c1", "c2")))
  G <- matrix(rnorm(n * 25L), n)
  K <- tcrossprod(G) / 25 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  kinship <- list(`1` = K, `2` = K * 0.98, `3` = K * 1.01)
  sex <- rep(0:1, n / 2L)
  x <- 3 * gp[["1"]][, 1L, "a1"] + 2 * gp[["2"]][, 2L, "b2"] + rnorm(n, sd = 0.3)
  w <- 2 * gp[["1"]][, 3L, "a3"] + 2 * gp[["3"]][, 4L, "c1"] + rnorm(n, sd = 0.3)
  y <- 0.6 * x + 0.3 * w + rnorm(n, sd = 0.3)
  q <- 0.4 * x + rnorm(n, sd = 0.4)
  pheno <- data.frame(x = x, w = w, y = y, q = q, row.names = ids)
  pheno$y[1:6] <- NA_real_
  pheno$x[7:8] <- NA_real_
  # q is measured in one sex only, so sex is constant on its complete samples
  # and every q context warns about a dependent covariate.
  pheno$q[sex == 1L] <- NA_real_
  covar <- data.frame(sex = sex, age = round(runif(n, 50, 90)), row.names = ids)
  map <- data.frame(marker = c("a1", "a2", "a3", "b1", "b2", "c1", "c2"), chr = c("1", "1", "1", "2", "2", "3", "3"),
    pos = c(1, 2, 3, 1, 2, 1, 2))
  scan <- matrix(c(8, 4.5, 1, 1, 6, 1, 1,
                   1, 1, 7, 4.2, 1, 6, 1,
                   5, 1, 4.5, 1, 4.4, 1, 1,
                   4.6, 1, 1, 1, 5.1, 1, 1), ncol = 4L, dimnames = list(map$marker, c("x", "w", "y", "q")))
  resources <- suppressWarnings(qtl2mr_prepare_resources(gp, kinship, pheno, covar, map))
  list(gp = gp, kinship = kinship, pheno = pheno, covar = covar, map = map, scan = scan, resources = resources,
    pairs = data.frame(exposure = c("x", "w", "x", "y", "q", "x"), outcome = c("y", "y", "w", "q", "x", "q")))
}

# qtl2mr_allele_prob_markers() of 0.6.1.9001 (c5a0e52), kept verbatim as the reference.
old_allele_prob_markers <- function(allele_probs, markers_to_validate = NULL,
                                    probability_policy = c("strict", "normalize_rounding")) {
  probability_policy <- match.arg(probability_policy)
  qtl2mr_validate_ids <- qtl2mr:::qtl2mr_validate_ids
  `%||%` <- function(x, y) if (is.null(x)) y else x
  if (!is.list(allele_probs) || !length(allele_probs)) stop("allele_probs must be a non-empty qtl2-style list of 3D arrays", call. = FALSE)
  qtl2mr_validate_ids(names(allele_probs), "allele_probs chromosome names")
  chr_names <- names(allele_probs) %||% as.character(seq_along(allele_probs))
  rows <- list()
  adjustments <- list()
  for (i in seq_along(allele_probs)) {
    arr <- allele_probs[[i]]
    if (length(dim(arr)) != 3L) stop("Each allele_probs element must be a sample x founder x marker array", call. = FALSE)
    dn <- dimnames(arr)
    for (d in seq_len(3L)) qtl2mr_validate_ids(dn[[d]], paste0("allele_probs dimension ", d))
    if (!is.numeric(arr)) stop("allele_probs must be numeric", call. = FALSE)
    first_ids <- dimnames(allele_probs[[1L]])[[1L]]
    if (i > 1L && !setequal(dn[[1L]], first_ids)) stop("allele_probs chromosomes must hold the same samples", call. = FALSE)
    checked <- if (is.null(markers_to_validate)) seq_len(dim(arr)[3L]) else which(dn[[3L]] %in% markers_to_validate)
    for (first in seq.int(1L, max(length(checked), 1L), by = 128L)) {
      if (!length(checked)) break
      chunk <- arr[, , checked[seq.int(first, min(first + 127L, length(checked)))], drop = FALSE]
      if (any(!is.finite(chunk))) stop("allele_probs must contain finite founder allele probabilities", call. = FALSE)
      if (any(chunk < -1e-8 | chunk > 1 + 1e-8)) stop("range", call. = FALSE)
      deviations <- abs(apply(chunk, c(1L, 3L), sum) - 1)
      sum_tolerance <- if (probability_policy == "strict") 1e-6 else dim(arr)[2L] * 5e-5
      if (any(deviations > sum_tolerance + 16 * .Machine$double.eps)) stop("sum", call. = FALSE)
      adjustments[[length(adjustments) + 1L]] <- data.frame(chr = chr_names[[i]],
        marker = dimnames(chunk)[[3L]], n_normalized = colSums(deviations > 1e-6 + 16 * .Machine$double.eps),
        max_abs_sum_deviation = apply(deviations, 2L, max), stringsAsFactors = FALSE)
    }
    if (i > 1L && !identical(dn[[2L]], dimnames(allele_probs[[1L]])[[2L]])) stop("founders", call. = FALSE)
    markers <- dn[[3L]]
    rows[[length(rows) + 1L]] <- data.frame(marker = markers, chr = chr_names[[i]], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  qtl2mr_validate_ids(out$marker, "allele_probs marker names")
  attr(out, "probability_normalization") <- do.call(rbind, adjustments)
  out
}

perf_random_probs <- function(n, markers, digits = NULL, seed = 1L, names_dims = FALSE) {
  set.seed(seed)
  a <- array(rexp(n * 8L * markers), c(n, 8L, markers),
    list(paste0("i", seq_len(n)), LETTERS[1:8], paste0("m", seed, "_", seq_len(markers))))
  for (j in seq_len(markers)) {
    m <- matrix(a[, , j], nrow = n)
    a[, , j] <- m / rowSums(m)
  }
  if (!is.null(digits)) a <- round(a, digits)
  if (names_dims) names(dimnames(a)) <- c("ind", "gen", "mar")
  a
}

# Results without timestamps and output paths, which differ between runs.
perf_strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    x <- x[setdiff(names(x) %||% character(), c("created_at", "output_dir", "paths", "output_manifest", "figure_manifest"))]
    x[] <- lapply(x, perf_strip)
  } else if (is.data.frame(x)) {
    # 0.6.1.9005 manifest counters describe the cache state of the run, by design.
    x <- x[, setdiff(names(x), c("runtime_sec", "wall_sec", "mean_wall_sec_per_pair", "run_id", "created_at",
      "context_computed", "lmm_fit_computed", "cache_released_entries", "cache_released_bytes", "cache_budget_bytes")), drop = FALSE]
  }
  x
}
`%||%` <- function(x, y) if (is.null(x)) y else x

perf_capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

perf_same <- function(a, b) identical(perf_strip(a), perf_strip(b), ignore.environment = TRUE)

test_that("the probability check gives the 0.6.1.9001 marker table and record", {
  cases <- list(
    list(`1` = perf_random_probs(30, 3), `2` = perf_random_probs(30, 1, seed = 2)),
    list(`1` = perf_random_probs(25, 300, digits = 4, seed = 3)),
    list(`1` = perf_random_probs(12, 129, digits = 5, seed = 4), `5` = perf_random_probs(12, 2, digits = 6, seed = 5)),
    list(`1` = perf_random_probs(1, 4, digits = 4, seed = 6)),
    list(`1` = perf_random_probs(20, 5, digits = 4, seed = 7, names_dims = TRUE))
  )
  cases[[3]][["5"]] <- cases[[3]][["5"]][dimnames(cases[[3]][["1"]])[[1L]], , , drop = FALSE]
  for (probs in cases) {
    for (policy in c("strict", "normalize_rounding")) {
      for (markers in list(NULL, dimnames(probs[[1L]])[[3L]][1L])) {
        old <- tryCatch(old_allele_prob_markers(probs, markers, policy), error = function(e) conditionMessage(e))
        new <- tryCatch(qtl2mr:::qtl2mr_allele_prob_markers(probs, markers, policy), error = function(e) conditionMessage(e))
        if (is.character(old)) expect_true(is.character(new)) else expect_identical(new, old)
      }
    }
  }
  f <- perf_fixture(round_digits = 4)
  w_new <- perf_capture(qtl2mr_prepare_resources(f$gp, f$kinship, f$pheno, f$covar, f$map))
  ref <- list(allele_probs = f$gp, kinship = f$kinship, pheno = as.data.frame(f$pheno), covar = as.data.frame(f$covar),
    marker_map = qtl2mr:::qtl2mr_marker_map_table(f$map), sample_ids = dimnames(f$gp[[1L]])[[1L]],
    metadata = list(), probability_policy = "normalize_rounding")
  ref$marker_map$chr <- as.character(ref$marker_map$chr)
  class(ref) <- c("qtl2mr_resources", "list")
  validate_qtl2mr_input(ref)
  w_ref <- perf_capture(qtl2mr:::qtl2mr_record_probability_policy(ref, warn = TRUE,
    marker_info = old_allele_prob_markers(ref$allele_probs, NULL, "normalize_rounding")))
  ref <- w_ref$value
  ref$sample_alignment <- qtl2mr:::qtl2mr_input_alignment(ref)
  expect_identical(w_new$value, ref)
  expect_identical(w_new$warnings, w_ref$warnings)
  expect_gt(nrow(ref$metadata$probability_normalization$adjustments), 0L)
})

test_that("a kinship that passed validation does not excuse a different one", {
  f <- perf_fixture()
  inst <- data.frame(exposure = "x", marker = "a1")
  expect_silent(validate_qtl2mr_input(qtl2mr_pair_input(f$resources, inst, "x", "y")))
  bad <- f$resources
  bad$kinship[["1"]][1L, 2L] <- bad$kinship[["1"]][1L, 2L] + 1
  expect_error(validate_qtl2mr_input(qtl2mr:::qtl2mr_resources_input(bad, "a1", "x", "y")), "symmetric")
})

test_that("only the requested markers are copied when probabilities are normalized", {
  f <- perf_fixture(round_digits = 4)
  input <- suppressWarnings(qtl2mr_pair_input(f$resources, data.frame(exposure = "x", marker = c("a1", "b2")), "x", "y"))
  copy <- suppressWarnings(qtl2mr:::qtl2mr_prepared_probability_copy(input, c("a1", "b2")))
  adjusted <- unique(copy$metadata$probability_normalization$adjustments$chr)
  expect_true(length(adjusted) > 0L)
  for (chr in adjusted) {
    expect_identical(dimnames(copy$allele_probs[[chr]])[[3L]], intersect(dimnames(f$gp[[chr]])[[3L]], c("a1", "b2")))
  }
  expect_identical(input$allele_probs, f$resources$allele_probs)
})


# qtl2mr_so_fit() of 0.6.1.9001 (c5a0e52), kept verbatim (apart from the name and
# the internal helpers it reaches) as the reference for the anchor reuse.
old_so_fit <- function(locus_blocks, allow_cross_cov = TRUE, ridge = 1e-10,
                       theta_grid = NULL) {
  qtl2mr_so_coordinates <- qtl2mr:::qtl2mr_so_coordinates
  gls_ivw <- qtl2mr:::gls_ivw
  block_diag_matrix <- qtl2mr:::block_diag_matrix
  qtl2mr_profile_set <- qtl2mr:::qtl2mr_profile_set
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


perf_fit_same <- function(a, b, probe = c(-3, -0.4, 0, 0.25, 2)) {
  identical(a[setdiff(names(a), c("objective", "independence"))], b[setdiff(names(b), c("objective", "independence"))]) &&
    identical(vapply(probe, a$objective, numeric(1)), vapply(probe, b$objective, numeric(1)))
}

test_that("second-order fits equal the 0.6.1.9001 search, and Q-only refits its minimum", {
  f <- perf_fixture()
  inst <- data.frame(exposure = "x", marker = c("a1", "a2", "b2", "c1"))
  pilot <- qtl2mr:::qtl2mr_fit_prepared(qtl2mr_pair_input(f$resources, inst, "x", "y"), "x", "y",
    qtl2mr:::qtl2mr_fit_settings("fail", "diagnostic_only", "chisq", 0L, 1L, 0L, 1L, TRUE, "star"))
  blocks <- pilot$locus_blocks
  for (k in seq_along(blocks)) {
    b <- blocks[seq_len(k)]
    for (cross in c(TRUE, FALSE)) {
      full <- qtl2mr:::qtl2mr_so_fit(b, cross, 1e-10)
      expect_true(perf_fit_same(full, old_so_fit(b, cross, 1e-10)))
      expect_true(perf_fit_same(qtl2mr:::qtl2mr_so_fit(b, cross, 1e-10, theta_grid = c(-1, 0.5)),
        old_so_fit(b, cross, 1e-10, theta_grid = c(-1, 0.5))))
      expect_identical(qtl2mr:::qtl2mr_so_fit(b, cross, 1e-10, q_only = TRUE), list(Q = full$Q))
    }
  }
  # Zero exposure variance at one locus (infinite tail objective).
  z <- blocks[1:2]
  z[[2L]]$Sigma_x[] <- 0
  z[[2L]]$Sigma_xy[] <- 0
  z[[2L]]$beta_x[] <- 0
  expect_true(perf_fit_same(qtl2mr:::qtl2mr_so_fit(z, TRUE, 1e-10), old_so_fit(z, TRUE, 1e-10)))
})

test_that("ordered worker calls signal warnings, messages and the first error as sequential calls do", {
  skip_on_os("windows")
  fun <- function(i) {
    warning("w", i, call. = FALSE)
    if (i == 4L) message("m", i)
    if (i == 5L) stop("stop at ", i, call. = FALSE)
    i^2
  }
  run <- function(workers) {
    old <- qtl2mr:::qtl2mr_set_workers(workers)
    on.exit(qtl2mr:::qtl2mr_set_workers(old), add = TRUE)
    seen <- list()
    value <- tryCatch(withCallingHandlers(qtl2mr:::qtl2mr_ordered_lapply(as.list(1:7), fun),
      warning = function(w) { seen[[length(seen) + 1L]] <<- w; invokeRestart("muffleWarning") },
      message = function(m) { seen[[length(seen) + 1L]] <<- m; invokeRestart("muffleMessage") }),
      error = function(e) e)
    list(value = value, seen = seen)
  }
  # Same conditions, classes and calls, including those raised by FUN itself.
  expect_identical(run(2L), run(1L))
  fun <- function(i) {
    warning("w", i)
    if (i == 3L) stop("raised by FUN at ", i)
    i
  }
  expect_identical(run(2L), run(1L))
  expect_identical(deparse(conditionCall(run(2L)$value)), "FUN(X[[i]], ...)")
  # Under options(warn = 2) refits stay in this process, where a caller's handler can muffle warnings.
  old_warn <- options(warn = 2)
  on.exit(options(old_warn), add = TRUE)
  expect_false(qtl2mr:::qtl2mr_can_fork(2L, 7L))
  muffled <- function(workers) {
    old <- qtl2mr:::qtl2mr_set_workers(workers)
    on.exit(qtl2mr:::qtl2mr_set_workers(old), add = TRUE)
    suppressWarnings(qtl2mr:::qtl2mr_ordered_lapply(as.list(1:2), function(i) { warning("w"); i }))
  }
  expect_identical(muffled(2L), muffled(1L))
  options(old_warn)
  ok <- function(workers) {
    old <- qtl2mr:::qtl2mr_set_workers(workers)
    on.exit(qtl2mr:::qtl2mr_set_workers(old), add = TRUE)
    qtl2mr:::qtl2mr_ordered_lapply(stats::setNames(as.list(1:5), letters[1:5]), function(i) i + 0.5)
  }
  expect_identical(ok(3L), ok(1L))
})

test_that("analysis results are identical with one and two workers", {
  skip_on_os("windows")
  f <- perf_fixture()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  out <- tempfile("v062-workers-")
  # The 0.6.1.9005 single-pair defaults ("auto" Q, which simulates for several
  # loci) with 40 draws and bootstrap replicates, which the entry no longer takes
  # (0.6.1.9006), so that the workers run random refits.
  run <- function(exposure, outcome, workers) {
    qtl2mr:::qtl2mr_cache_reset()
    set.seed(11)
    before <- .Random.seed
    r <- perf_capture(qtl2mr_with_settings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, exposure, outcome),
      output_dir = out, overwrite = TRUE, make_html = FALSE, make_figures = FALSE, progress = FALSE,
      n_workers = workers), heterogeneity_method = "auto", heterogeneity_nsim = 40L, diagnostic_bootstrap_B = 40L))
    r$seed_kept <- identical(.Random.seed, before)
    r
  }
  for (pair in list(c("x", "y"), c("x", "q"))) {
    one <- run(pair[[1L]], pair[[2L]], 1L)
    two <- run(pair[[1L]], pair[[2L]], 2L)
    expect_true(perf_same(one$value, two$value))
    expect_identical(one$warnings, two$warnings)
    expect_true(one$seed_kept && two$seed_kept)
    expect_identical(one$value$pilot$fit$objective(0.37), two$value$pilot$fit$objective(0.37))
    expect_true(is.finite(one$value$pilot$summary$Q_p_sim))
    expect_true(nrow(one$value$loo_mr) > 1L)
  }
  expect_true(any(grepl("Dropping linearly dependent covariate", one$warnings)))
  expect_error(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, "x", "y"), n_workers = 0L), "n_workers")
})

test_that("single analyses do not depend on what the session has cached", {
  f <- perf_fixture()
  # The 0.6.1.9005 single-pair defaults with 20 draws and replicates (no longer
  # entry options since 0.6.1.9006); the instrument-set sensitivity now always runs.
  run <- function(inst) perf_capture(qtl2mr_with_settings(run_qtl2mr_analysis(
    qtl2mr_pair_input(f$resources, inst, "x", "q"),
    output_dir = tempfile("v062-warm-"), make_html = FALSE, make_figures = FALSE, progress = FALSE),
    heterogeneity_method = "auto", heterogeneity_nsim = 20L, diagnostic_bootstrap_B = 20L))
  qtl2mr:::qtl2mr_cache_reset()
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  cold <- run(inst)
  expect_identical(cold$value$pilot$summary$Q_p_method, "simulation")  # the mock took effect
  qtl2mr:::qtl2mr_cache_reset()
  sel <- perf_capture(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE))
  expect_identical(sel$value, inst)
  after_selection <- run(inst)
  again <- run(inst)
  expect_true(perf_same(after_selection$value, cold$value))
  expect_true(perf_same(again$value, cold$value))
  expect_identical(after_selection$warnings, cold$warnings)
  expect_identical(again$warnings, cold$warnings)
  # The second selection in a warm session reuses every fit and warns the same.
  expect_identical(perf_capture(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), sel)
})

perf_batch <- function(f, inst, ...) {
  perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE, ...))
}

test_that("batch results do not depend on the session cache, worker count or wave size", {
  skip_on_os("windows")
  f <- perf_fixture(round_digits = 5)
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  qtl2mr:::qtl2mr_cache_reset()
  cold <- perf_batch(f, inst)
  warm <- perf_batch(f, inst)
  expect_true(perf_same(warm$value, cold$value))
  expect_identical(warm$warnings, cold$warnings)
  expect_identical(warm$value$run_manifest[c("context_cache_hits", "context_cache_misses", "lmm_fit_cache_hits",
    "lmm_fit_cache_misses")], cold$value$run_manifest[c("context_cache_hits", "context_cache_misses",
    "lmm_fit_cache_hits", "lmm_fit_cache_misses")])
  expect_true(any(grepl("Dropping linearly dependent covariate", cold$warnings)))
  qtl2mr:::qtl2mr_cache_reset()
  # max_loci (removed in 0.6.1.9006) acted after fitting every candidate, so this
  # selection leaves the session cache as before.
  invisible(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)))
  after_selection <- perf_batch(f, inst)
  expect_true(perf_same(after_selection$value, cold$value))
  expect_identical(after_selection$warnings, cold$warnings)
  # Small waves: planning in several steps, with releases in between.
  assign("budget", 2e6, envir = qtl2mr:::qtl2mr_state)
  small <- tryCatch(perf_batch(f, inst), finally = assign("budget", NULL, envir = qtl2mr:::qtl2mr_state))
  expect_true(perf_same(small$value, cold$value))
  expect_identical(small$warnings, cold$warnings)
  # Two workers: groups run in workers (whose warnings are not returned, as before).
  qtl2mr:::qtl2mr_cache_reset()
  forked_old_layout <- perf_batch(f, inst, n_workers = 2L)
  serial <- cold$value
  expect_true(perf_same(forked_old_layout$value[setdiff(names(serial), "run_manifest")],
    serial[setdiff(names(serial), "run_manifest")]))
  # One exposure group with two workers: trait fits are planned in parallel, pairs fitted in order.
  one_group <- f$pairs[f$pairs$exposure == "x", ]
  seq_one <- perf_capture(run_qtl2mr_batch(f$resources, one_group, inst, min_samples = 20L, progress = FALSE))
  qtl2mr:::qtl2mr_cache_reset()
  par_one <- perf_capture(run_qtl2mr_batch(f$resources, one_group, inst, min_samples = 20L, progress = FALSE, n_workers = 2L))
  expect_true(perf_same(par_one$value[setdiff(names(par_one$value), "run_manifest")],
    seq_one$value[setdiff(names(seq_one$value), "run_manifest")]))
  keep <- c("parallel_backend", "context_cache_hits", "context_cache_misses", "lmm_fit_cache_hits", "lmm_fit_cache_misses")
  expect_identical(par_one$value$run_manifest[keep], seq_one$value$run_manifest[keep])
  expect_identical(par_one$warnings, seq_one$warnings)
  # After the batch its entries are the first to be released, but a second call
  # in the session reuses them (same results).
  reg <- qtl2mr:::qtl2mr_registry()
  traits <- unique(c(one_group$exposure, one_group$outcome))
  left <- Filter(function(id) reg$entries[[id]]$info$trait %in% traits, ls(reg$entries))
  expect_gt(length(left), 0L)
  expect_true(all(vapply(left, function(id) reg$entries[[id]]$tick, numeric(1)) == 0))
  again <- perf_capture(run_qtl2mr_batch(f$resources, one_group, inst, min_samples = 20L, progress = FALSE))
  expect_true(perf_same(again$value[setdiff(names(again$value), "run_manifest")],
    seq_one$value[setdiff(names(seq_one$value), "run_manifest")]))
  expect_identical(again$value$run_manifest[keep], seq_one$value$run_manifest[keep])
  expect_identical(again$warnings, seq_one$warnings)
})

test_that("a batch in a fresh R process gives the same results", {
  skip_on_cran()
  skip_on_os("windows")
  f <- perf_fixture()
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  here <- perf_batch(f, inst)$value
  here_again <- perf_batch(f, inst)$value
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  saveRDS(list(resources = f$resources, pairs = f$pairs, instruments = inst), input)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    "suppressPackageStartupMessages(library(qtl2mr))",
    "d <- readRDS(args[[1L]])",
    "b <- suppressWarnings(run_qtl2mr_batch(d$resources, d$pairs, d$instruments, min_samples = 20L, progress = FALSE))",
    "saveRDS(b, args[[2L]])"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), shQuote(input), shQuote(output)),
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)), stdout = FALSE, stderr = FALSE)
  skip_if(status != 0L || !file.exists(output), "child R process unavailable")
  fresh <- readRDS(output)
  expect_true(perf_same(here, fresh))
  expect_true(perf_same(here_again, fresh))
})

test_that("batch rows equal single analyses on each trait's own complete samples", {
  f <- perf_fixture()
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  cols <- qtl2mr:::qtl2mr_batch_summary_columns
  results <- list()
  # sample_policy = "pair_common_complete" was removed in 0.6.1.9006; the run is
  # repeated to check that it is reproducible.
  for (repeat_run in 1:2) {
    batch <- perf_batch(f, inst)$value
    expect_identical(batch$main_results$status, rep("ok", nrow(f$pairs)))
    for (i in seq_len(nrow(f$pairs))) {
      single <- suppressWarnings(run_qtl2mr_analysis(qtl2mr_pair_input(f$resources, inst, f$pairs$exposure[[i]],
        f$pairs$outcome[[i]]), output_dir = tempfile("v062-policy-"), make_html = FALSE,
        make_figures = FALSE, progress = FALSE))
      # rbind() of batch rows can promote an integer median to double (as in 0.6.1.9001); compare those by value.
      row <- as.list(batch$main_results[i, cols])
      ref <- as.list(single$pilot$summary[1L, cols])
      num <- vapply(row, is.numeric, logical(1)) & vapply(ref, is.numeric, logical(1))
      row[num] <- lapply(row[num], as.double)
      ref[num] <- lapply(ref[num], as.double)
      expect_identical(row, ref)
    }
    results[[length(results) + 1L]] <- batch
  }
  expect_true(perf_same(results[[2L]], results[[1L]]))
  expect_null(results[[1L]]$config$sample_policy)
  expect_null(qtl2mr_pair_input(f$resources, inst, "x", "q")$exposure_samples)
  expect_error(qtl2mr_pair_input(f$resources, inst, "x", "q", sample_policy = "pair_common_complete"), "unused argument")
  # q has 45 complete samples of its own.
  enough <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 45L, progress = FALSE))$value
  expect_identical(enough$main_results$status, rep("ok", nrow(f$pairs)))
})

# The main_results columns that depend on the stored blocks alone, rebuilt from them
# exactly as run_prepared_qtl2mr() builds them, with the run's settings.
perf_row_from_blocks <- function(blocks, method, nsim, seed, boot, boot_seed, policy) {
  meta <- attr(blocks, "qtl2mr_locus_blocks")
  fit <- fit_multilocus_founder_mr(blocks, allow_cross_cov = TRUE, ridge = 1e-10, heterogeneity_method = method,
    q_sim_n = nsim, q_sim_seed = seed)
  layer <- qtl2mr:::fit_multilocus_diagnostic_layer(blocks, global_fit = fit, ridge = 1e-10, bootstrap_B = boot,
    bootstrap_seed = boot_seed)
  reported <- qtl2mr:::qtl2mr_reported_estimate(fit, layer, policy)
  est <- function(name) {
    rows <- layer$estimates$estimator == name
    if (!any(rows)) NA_real_ else layer$estimates$theta[which(rows)[1L]]
  }
  q <- qtl2mr:::qtl2mr_add_q_fields(data.frame(Q = fit$Q, Q_chisq_df = fit$Q_chisq_df %||% NA_integer_,
    Q_p_chisq = fit$Q_p_chisq %||% NA_real_, Q_p_reported = fit$Q_p_reported %||% NA_real_,
    Q_p_method = fit$Q_p_method %||% NA_character_, stringsAsFactors = FALSE), fit = fit)
  overlap <- qtl2mr:::qtl2mr_locus_overlap_table(blocks)
  list(
    exact = list(pair_id = meta$pair_id, exposure = meta$exposure, outcome = meta$outcome,
      n_loci = length(blocks),
      n_loci_steiger_fail = sum(vapply(blocks, `[[`, logical(1), "steiger_pass") %in% FALSE),
      theta_reported = reported$theta, se_reported = reported$se, reported_ci_lower = reported$ci_lower,
      reported_ci_upper = reported$ci_upper, reported_ci_type = reported$ci_type,
      reported_estimator = reported$estimator, reported_estimator_source = reported$source,
      reported_estimator_policy = reported$policy, reported_estimator_reason = reported$reason,
      theta = fit$theta, se = fit$se, profile_ci_lower = fit$profile_ci_lower, profile_ci_upper = fit$profile_ci_upper,
      profile_p = fit$profile_p, confidence_set_type = fit$confidence_set_type %||% "unresolved",
      confidence_set_bounded = qtl2mr:::qtl2mr_confidence_set_bounded(fit),
      primary_optimizer_status = fit$optimizer_status %||% NA_character_,
      primary_bracket_status = fit$bracket_status %||% NA_character_,
      primary_bracket_n_expand = fit$bracket_n_expand %||% NA_integer_,
      primary_profile_ci_finite = isTRUE(fit$profile_ci_finite), primary_reliable = isTRUE(fit$primary_reliable),
      numerical_ok = isTRUE(fit$numerical_ok), inference_status = fit$inference_status %||% "not_calibrated",
      effect_units = if (isTRUE(meta$standardize_phenotypes)) "SD_y per SD_x" else "outcome units per exposure unit",
      first_order_theta = fit$first_order_theta %||% NA_real_, first_order_se = fit$first_order_se %||% NA_real_,
      theta_independence = fit$independence$theta, theta_overlap_delta = fit$theta_overlap_delta,
      theta_two_step_ivw = est("Two-step IVW"), theta_weighted_median = est("Weighted median"),
      theta_weighted_mode = est("Weighted mode"),
      Q = q$Q[[1L]], Q_chisq_df = q$Q_chisq_df[[1L]], Q_p_chisq = q$Q_p_chisq[[1L]], Q_p_reported = q$Q_p_reported[[1L]],
      Q_p_method = q$Q_p_method[[1L]], cross_cov_model = meta$cross_cov_model,
      standardize_phenotypes = meta$standardize_phenotypes),
    # rbind() of batch rows may store these as double (as in 0.6.1.9001); compared by value.
    by_value = list(n_x_median = stats::median(overlap$n_x), n_y_median = stats::median(overlap$n_y),
      n_overlap_median = stats::median(overlap$n_overlap))
  )
}

test_that("per-locus blocks are written on request and refit to the batch estimates", {
  f <- perf_fixture()
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  plain_dir <- tempfile("v062-plain-")
  plain <- perf_batch(f, inst, output_dir = plain_dir)$value
  expect_false(dir.exists(file.path(plain_dir, "batch_locus_blocks")))
  expect_false(file.exists(file.path(plain_dir, "batch_locus_blocks_manifest.csv")))
  expect_null(plain$locus_blocks_manifest)
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, write_locus_blocks = TRUE), "requires output_dir")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, write_locus_blocks = NA), "TRUE or FALSE")

  dir <- tempfile("v062-blocks-")
  # The fit settings of each run. They are fixed at the batch entry since
  # 0.6.1.9006 (the first run uses the fixed settings); the others replace them
  # internally. sample_policy = "pair_common_complete" was removed, so the
  # blocks are always of the own-complete samples.
  runs <- list(
    list(),
    list(heterogeneity_method = "simulation", heterogeneity_nsim = 25L,
      heterogeneity_seed = 3L, diagnostic_bootstrap_B = 20L, diagnostic_bootstrap_seed = 5L,
      reported_estimator_policy = "auto_fallback", cross_cov_model = "raw", standardize_phenotypes = FALSE),
    list(heterogeneity_method = "none", cross_cov_model = "zero", reported_estimator_policy = "strict_so"))
  # Columns of main_results that the blocks cannot give: they need fields not stored (LMM backend, loci
  # dropped before fitting, trait SDs). Listed in IMPLEMENTATION_NOTES, O1. With the founder QC on by
  # default since 0.6.1.9006, n_loci_dropped_founder_qc (loci dropped before fitting) is one of them.
  not_from_blocks <- c("lmm_backend", "n_loci_selected", "n_loci_dropped_nonestimable", "exposure_sd", "outcome_sd",
    "raw_conversion_available", "theta_raw_units", "se_raw_units", "status", "error_message", "locus_fit_status",
    "runtime_sec", "n_loci_dropped_founder_qc")
  for (run in runs) {
    policy <- "own_complete"
    run_settings <- do.call(qtl2mr_test_settings, run)
    batch <- with_mocked_bindings(perf_batch(f, inst, output_dir = dir, write_locus_blocks = TRUE, overwrite = TRUE),
      qtl2mr_fixed_fit_settings = function() run_settings)$value
    settings <- batch$config
    if (identical(run, runs[[1L]])) {
      expect_true(perf_same(batch[setdiff(names(batch), c("config", "locus_blocks_manifest"))],
        plain[setdiff(names(plain), "config")]))
    }
    manifest <- utils::read.csv(file.path(dir, "batch_locus_blocks_manifest.csv"), stringsAsFactors = FALSE)
    expect_identical(manifest$pair_id, batch$main_results$pair_id)
    expect_identical(batch$locus_blocks_manifest$file, manifest$file)
    main <- batch$main_results
    loci <- batch$locus_results
    for (i in seq_len(nrow(manifest))) {
      blocks <- readRDS(file.path(dir, manifest$file[[i]]))
      expect_identical(attr(blocks, "qtl2mr_locus_blocks")$sample_policy, policy)
      expect_true(all(c("marker", "chr", "pos", "beta_x", "beta_y", "Sigma_x", "Sigma_y", "Sigma_xy", "F_X_block",
        "F_Y_block", "partial_R2_X_block", "partial_R2_Y_block", "steiger_pass", "df_x", "df_y", "n_x", "n_y",
        "n_overlap", "hsq_x", "hsq_y") %in% names(blocks[[1L]])))
      j <- match(manifest$pair_id[[i]], main$pair_id)
      rebuilt <- perf_row_from_blocks(blocks, settings$heterogeneity_method, settings$heterogeneity_nsim,
        settings$heterogeneity_seed, settings$diagnostic_bootstrap_B, settings$diagnostic_bootstrap_seed,
        settings$reported_estimator_policy)
      for (column in names(rebuilt$exact)) {
        expect_identical(rebuilt$exact[[column]], main[[column]][[j]], label = paste(policy, column))
      }
      for (column in names(rebuilt$by_value)) {
        expect_identical(as.double(rebuilt$by_value[[column]]), as.double(main[[column]][[j]]), label = column)
      }
      # Every main_results column is either rebuilt above or listed as not given by the blocks.
      expect_setequal(setdiff(names(main), c(names(rebuilt$exact), names(rebuilt$by_value))), not_from_blocks)
      l <- loci[loci$pair_id == manifest$pair_id[[i]], ]
      expect_identical(unname(vapply(blocks, `[[`, numeric(1), "F_X_block")), l$F_X_block)
      expect_identical(unname(vapply(blocks, `[[`, logical(1), "steiger_pass")), l$steiger_pass)
      expect_identical(unname(vapply(blocks, `[[`, integer(1), "n_y")), l$n_y)
      expect_identical(as.numeric(blocks[[1L]]$fit_x$df), blocks[[1L]]$df_x)
    }
  }
  # Resume completes the same run and rewrites the block files from the pair records.
  partial <- f$pairs[1:3, ]
  rdir <- tempfile("v062-resume-")
  first <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, write_locus_blocks = TRUE))$value
  unlink(file.path(rdir, "batch_locus_blocks", "pair_000002.rds"))
  unlink(file.path(rdir, "batch_pair_records", "pair_000004.rds"))
  resumed <- perf_capture(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, write_locus_blocks = TRUE, resume = TRUE))$value
  expect_identical(resumed$run_manifest$n_pairs_resumed, nrow(f$pairs) - 1L)
  expect_true(perf_same(resumed[setdiff(names(resumed), "run_manifest")], first[setdiff(names(first), "run_manifest")]))
  expect_true(file.exists(file.path(rdir, "batch_locus_blocks", "pair_000002.rds")))
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, resume = TRUE), "fingerprint changed")
  expect_error(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE, output_dir = rdir),
    "already exists")
  cleared <- perf_capture(run_qtl2mr_batch(f$resources, partial, inst, min_samples = 20L, progress = FALSE,
    output_dir = rdir, overwrite = TRUE))$value
  expect_false(dir.exists(file.path(rdir, "batch_locus_blocks")))
  expect_false(file.exists(file.path(rdir, "batch_locus_blocks_manifest.csv")))
  expect_identical(cleared$main_results$theta, first$main_results$theta[1:3])
})

perf_capture_all <- function(expr) {
  conditions <- character()
  value <- withCallingHandlers(expr,
    warning = function(w) { conditions <<- c(conditions, paste("W", conditionMessage(w))); invokeRestart("muffleWarning") },
    message = function(m) { conditions <<- c(conditions, paste("M", conditionMessage(m))); invokeRestart("muffleMessage") })
  list(value = value, conditions = conditions)
}

test_that("several exposure groups with two workers run in workers, in every wave", {
  skip_on_os("windows")
  f <- perf_fixture()
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  run <- function(workers, budget = NULL) {
    qtl2mr:::qtl2mr_cache_reset()
    assign("budget", budget, envir = qtl2mr:::qtl2mr_state)
    on.exit(assign("budget", NULL, envir = qtl2mr:::qtl2mr_state), add = TRUE)
    perf_capture_all(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = TRUE,
      n_workers = workers))
  }
  sequential <- run(1L)
  forked <- run(2L)
  small <- run(2L, budget = 2e5)
  # The planner makes several waves of at least two groups for the small budget.
  # The plan of the batch: since 0.6.1.9006 without a sample policy, and with the
  # founder QC the batch applies by default.
  plan <- qtl2mr:::qtl2mr_batch_plan(f$resources, small$value$pairs,
    qtl2mr:::qtl2mr_normalize_instruments(inst, f$resources), min_founder_copies = 5)
  groups <- unname(split(seq_len(nrow(f$pairs)), factor(f$pairs$exposure, levels = unique(f$pairs$exposure))))
  assign("budget", 2e5, envir = qtl2mr:::qtl2mr_state)
  waves <- qtl2mr:::qtl2mr_batch_waves(plan, groups, 0.5 * qtl2mr:::qtl2mr_cache_budget(),
    qtl2mr:::qtl2mr_batch_entry_sizes(f$resources), min_units = 2L)
  assign("budget", NULL, envir = qtl2mr:::qtl2mr_state)
  expect_gt(length(waves), 1L)
  expect_true(all(vapply(waves, function(w) sum(vapply(groups, function(g) all(g %in% w), logical(1))),
    integer(1)) >= 2L))
  # As in 0.6.1.9001: groups run in forked workers, so their progress lines and
  # warnings are not signalled here, while the sequential run signals them.
  expect_true(any(startsWith(sequential$conditions, "M [qtl2mr batch]")))
  expect_true(any(grepl("Dropping linearly dependent covariate", sequential$conditions)))
  expect_identical(forked$conditions, character())
  expect_identical(small$conditions, character())
  expect_identical(small$value$run_manifest$parallel_backend, "parallel::mclapply")
  expect_true(perf_same(small$value, forked$value))
  drop <- "run_manifest"
  expect_true(perf_same(small$value[setdiff(names(small$value), drop)],
    sequential$value[setdiff(names(sequential$value), drop)]))
})

test_that("a planned wave fits the cache budget, so no fit is computed twice", {
  f <- perf_fixture()
  # The table max_loci = 2L gave (max_loci was removed in 0.6.1.9006).
  inst <- qtl2mr_test_top_loci(suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, progress = FALSE)), 2L)
  sizes <- qtl2mr:::qtl2mr_batch_entry_sizes(f$resources)
  run <- function(budget) {
    qtl2mr:::qtl2mr_cache_reset()
    assign("budget", budget, envir = qtl2mr:::qtl2mr_state)
    on.exit(assign("budget", NULL, envir = qtl2mr:::qtl2mr_state), add = TRUE)
    batch <- suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE))
    reg <- qtl2mr:::qtl2mr_registry()
    fits <- Filter(function(id) startsWith(id, "fit:"), ls(reg$entries))
    # The batch's plan (no sample policy since 0.6.1.9006; the default founder QC).
    plan <- qtl2mr:::qtl2mr_batch_plan(f$resources, batch$pairs,
      qtl2mr:::qtl2mr_normalize_instruments(inst, f$resources), min_founder_copies = 5)
    list(batch = batch, computed = qtl2mr:::qtl2mr_state$computed[["fit"]], distinct = length(fits),
      fit_bytes = max(vapply(fits, function(id) reg$entries[[id]]$bytes, numeric(1))),
      waves = length(qtl2mr:::qtl2mr_batch_waves(plan, as.list(seq_len(nrow(batch$pairs))),
        0.5 * qtl2mr:::qtl2mr_cache_budget(), sizes)))
  }
  ref <- run(NULL)
  expect_identical(ref$waves, 1L)
  expect_identical(ref$computed, as.numeric(ref$distinct))
  # The planner charges entries as the registry does.
  expect_gte(sizes$fit, ref$fit_bytes)
  expect_lte(sizes$fit, 1.05 * ref$fit_bytes)
  for (budget in c(4e6, 2e6, 1.6e6, 1.2e6)) {
    r <- run(budget)
    expect_true(perf_same(r$batch, ref$batch))
    expect_identical(r$computed, ref$computed)
  }
  expect_gt(run(1.2e6)$waves, 2L)
})

test_that("the session cache stays within its budget and keeps pinned entries", {
  qtl2mr:::qtl2mr_cache_reset()
  on.exit(qtl2mr:::qtl2mr_cache_reset(), add = TRUE)
  assign("budget", 5e5, envir = qtl2mr:::qtl2mr_state)
  value <- function(i) rep(as.numeric(i), 1e4)
  qtl2mr:::qtl2mr_registry_put("fit", "k1", value(1), info = list(trait = "a", chr = "1", marker = "m1"))
  qtl2mr:::qtl2mr_registry_pin("fit:k1")
  for (i in 2:12) qtl2mr:::qtl2mr_registry_put("fit", paste0("k", i), value(i), info = list(trait = "a", chr = "1", marker = paste0("m", i)))
  reg <- qtl2mr:::qtl2mr_registry()
  expect_lte(reg$bytes, 5e5)
  expect_false(is.null(qtl2mr:::qtl2mr_registry_get("fit", "k1")))
  expect_false(is.null(qtl2mr:::qtl2mr_registry_get("fit", "k12")))
  expect_null(qtl2mr:::qtl2mr_registry_get("fit", "k2"))
  expect_equal(reg$bytes, sum(vapply(ls(reg$entries), function(id) reg$entries[[id]]$bytes, numeric(1))))
  # Demoted entries go first, even when recently used.
  qtl2mr:::qtl2mr_registry_demote("fit:k12")
  qtl2mr:::qtl2mr_registry_put("fit", "k13", value(13), info = list(trait = "a", chr = "1", marker = "m13"))
  expect_null(qtl2mr:::qtl2mr_registry_get("fit", "k12"))
  expect_false(is.null(qtl2mr:::qtl2mr_registry_get("fit", "k11")))
})
