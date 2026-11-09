# 0.6.1.9005: performance changes that must leave every value bit for bit unchanged.

# qtl2mr_so_coordinates() of 0.6.1.9004 (5a6616f), kept verbatim (renamed) as the reference.
v0904_so_coordinates <- function(blocks, allow_cross_cov, ridge) {
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

environment(v0904_so_coordinates) <- asNamespace("qtl2mr")

# The 0.6.2 performance fixture (test-v062-performance.R), copied.
v065_fixture <- function(round_digits = NULL) {
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


v065_blocks <- function(K, seed, cross = TRUE, zero_x = FALSE, scale_x = 1, perfect = FALSE, k = 7L) {
  set.seed(seed)
  fc <- paste0("FC", seq_len(k))
  lapply(seq_len(K), function(i) {
    A <- matrix(rnorm(4 * k * k), 2 * k)
    J <- crossprod(A) / (2 * k) * 0.015 + diag(0.001, 2 * k)
    Sx <- J[1:k, 1:k] * scale_x
    Sy <- J[k + 1:k, k + 1:k]
    Sxy <- J[1:k, k + 1:k] * sqrt(scale_x)
    if (perfect) Sx <- Sy <- Sxy <- diag(0.01, k)
    bx <- stats::setNames(rnorm(k, 0, 0.2), fc)
    if (zero_x && i == 1L) {
      Sx[] <- 0
      Sxy[] <- 0
      bx[] <- 0
    }
    dimnames(Sx) <- dimnames(Sy) <- dimnames(Sxy) <- list(fc, fc)
    list(beta_x = bx, beta_y = stats::setNames(0.3 * bx + rnorm(k, 0, 0.1), fc),
      Sigma_x = Sx, Sigma_y = Sy, Sigma_xy = if (cross) Sxy else NULL)
  })
}

test_that("the SO objective and its derivatives are bit for bit those of 0.6.1.9004", {
  ns <- asNamespace("qtl2mr")
  cases <- list(
    list(K = 1L, cross = TRUE, ridge = 1e-10),
    list(K = 3L, cross = FALSE, ridge = 1e-10),
    list(K = 12L, cross = TRUE, ridge = 0),
    list(K = 12L, cross = TRUE, ridge = 1e-10, zero_x = TRUE),
    list(K = 5L, cross = TRUE, ridge = 1e-10, scale_x = 1e-6),
    list(K = 2L, cross = TRUE, ridge = 1e-10, perfect = TRUE)
  )
  alphas <- c(seq(-pi / 2, pi / 2, length.out = 101L), -pi / 2 - 0.1, pi / 2 + 0.1, 0, 1e-12, pi / 4, -pi / 4)
  for (cs in cases) {
    bl <- v065_blocks(cs$K, seed = 65L + cs$K, cross = cs$cross, zero_x = isTRUE(cs$zero_x),
      scale_x = cs$scale_x %||% 1, perfect = isTRUE(cs$perfect))
    old <- v0904_so_coordinates(bl, cs$cross, cs$ridge)
    new <- ns$qtl2mr_so_coordinates(bl, cs$cross, cs$ridge)
    expect_true(identical(new$scale, old$scale))
    expect_true(identical(new$blocks, old$blocks))
    for (a in alphas) {
      expect_true(identical(new$evaluate(a), old$evaluate(a), num.eq = FALSE))
      expect_true(identical(new$evaluate(a, TRUE), old$evaluate(a, TRUE), num.eq = FALSE))
    }
  }
})

test_that("qtl2mr_so_fit() gives the 0.6.1.9004 result bit for bit", {
  ns <- asNamespace("qtl2mr")
  old_fit <- ns$qtl2mr_so_fit
  env <- new.env(parent = ns)
  env$qtl2mr_so_coordinates <- v0904_so_coordinates
  environment(old_fit) <- env
  strip <- function(f) {
    f$objective <- NULL
    f
  }
  for (K in c(1L, 4L, 15L)) for (cross in c(TRUE, FALSE)) {
    bl <- v065_blocks(K, seed = 650L + K, cross = cross)
    a <- ns$qtl2mr_so_fit(bl, cross)
    b <- old_fit(bl, cross)
    expect_true(identical(strip(a), strip(b), num.eq = FALSE))
    expect_true(identical(a$objective(0.3), b$objective(0.3)))
    expect_true(identical(ns$qtl2mr_so_fit(bl, cross, q_only = TRUE), old_fit(bl, cross, q_only = TRUE)))
  }
})

test_that("the context transform equals the diagonal-matrix product bit for bit", {
  ns <- asNamespace("qtl2mr")
  w <- c(2, 3, 0.5)
  E <- matrix(c(-0, 1, 0, -2, 0.25, -0, 3, -0, 1), 3)
  expect_true(identical(w * E + 0, diag(w, nrow = 3L) %*% E, num.eq = FALSE))
  set.seed(651)
  n <- 60L
  ids <- paste0("s", seq_len(n))
  G <- matrix(rnorm(n * 10L), n)
  K <- tcrossprod(G) / 10 + diag(0.5, n)
  dimnames(K) <- list(ids, ids)
  y <- stats::setNames(rnorm(n), ids)
  covar <- data.frame(sex = rep(0:1, n / 2L), row.names = ids)
  for (h in list(NULL, 0, 0.37, 1)) {
    ctx <- ns$prepare_lmm_context(y, kinship = K, addcovar = covar, hsq = h, trait_name = "t", chr = "1")
    d <- ns$qtl2mr_decompose_kinship(K[ctx$sample_ids, ctx$sample_ids], sample_ids = ctx$sample_ids)
    wts <- 1 / sqrt(as.numeric(ctx$hsq) * d$values + (1 - as.numeric(ctx$hsq)))
    ref <- diag(wts, nrow = length(wts)) %*% d$vectors
    dimnames(ref) <- list(ctx$sample_ids, ctx$sample_ids)
    expect_true(identical(ctx$raw_transform, ref, num.eq = FALSE))
  }
})

test_that("kept eigendecompositions are reused across many chromosomes and bounded in bytes", {
  ns <- asNamespace("qtl2mr")
  ns$qtl2mr_cache_reset()
  old_opt <- options(qtl2mr.decomposition_budget = NULL)
  on.exit({
    options(old_opt)
    suppressMessages(untrace("qtl2mr_decompose_kinship", where = ns))
    ns$qtl2mr_cache_reset()
  }, add = TRUE)
  set.seed(652)
  n <- 30L
  ids <- paste0("s", seq_len(n))
  Ks <- lapply(1:12, function(i) {
    G <- matrix(rnorm(n * 8L), n)
    K <- tcrossprod(G) / 8 + diag(0.5, n)
    dimnames(K) <- list(ids, ids)
    K
  })
  cnt <- new.env()
  cnt$n <- 0L
  suppressMessages(trace("qtl2mr_decompose_kinship", tracer = bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
    where = ns, print = FALSE))
  run <- function() lapply(Ks, function(K) ns$qtl2mr_decompose_kinship_cached(K, ids, K))
  first <- run()
  second <- run()
  expect_identical(cnt$n, 12L)
  expect_true(identical(first, second))
  # A budget below one decomposition keeps the minimum of 8, so twelve chromosomes visited in turn are all
  # recomputed (the 0.6.1.9004 behaviour), with the same values.
  ns$qtl2mr_cache_reset()
  cnt$n <- 0L
  options(qtl2mr.decomposition_budget = 0)
  a <- run()
  b <- run()
  expect_identical(cnt$n, 24L)
  expect_true(identical(a, first))
  expect_true(identical(b, first))
  options(qtl2mr.decomposition_budget = -1)
  expect_error(ns$qtl2mr_decomposition_budget(), "non-negative")
})

test_that("the cache budget option is honoured", {
  ns <- asNamespace("qtl2mr")
  old_opt <- options(qtl2mr.cache_budget = NULL)
  on.exit(options(old_opt), add = TRUE)
  expect_identical(ns$qtl2mr_cache_budget(), 4 * 1024^3)
  options(qtl2mr.cache_budget = 12345)
  expect_identical(ns$qtl2mr_cache_budget(), 12345)
  assign("budget", 99, envir = ns$qtl2mr_state)
  expect_identical(ns$qtl2mr_cache_budget(), 99)
  assign("budget", NULL, envir = ns$qtl2mr_state)
  options(qtl2mr.cache_budget = "x")
  expect_error(ns$qtl2mr_cache_budget(), "non-negative")
})

test_that("the run manifest counts computations and releases, in workers too, without changing results", {
  skip_on_os("windows")
  ns <- asNamespace("qtl2mr")
  old_opt <- options(qtl2mr.cache_budget = NULL)
  on.exit({
    options(old_opt)
    ns$qtl2mr_cache_reset()
  }, add = TRUE)
  f <- v065_fixture(round_digits = 5)
  inst <- suppressWarnings(qtl2mr_select_instruments(f$resources, f$scan, max_loci = 2L, progress = FALSE))
  batch <- function(...) suppressWarnings(run_qtl2mr_batch(f$resources, f$pairs, inst, min_samples = 20L, progress = FALSE, ...))
  cols <- c("context_computed", "lmm_fit_computed", "cache_released_entries", "cache_released_bytes", "cache_budget_bytes")
  ns$qtl2mr_cache_reset()
  cold <- batch()
  m <- cold$run_manifest
  expect_true(all(cols %in% names(m)))
  expect_true(m$context_computed > 0 && m$lmm_fit_computed > 0)
  expect_identical(m$cache_released_entries, 0)
  expect_identical(m$cache_budget_bytes, 4 * 1024^3)
  warm <- batch()
  expect_identical(warm$run_manifest$context_computed, 0)
  expect_identical(warm$run_manifest$lmm_fit_computed, 0)
  # The hit and miss columns keep their meaning: a warm run counts as the cold one did.
  hm <- c("context_cache_hits", "context_cache_misses", "lmm_fit_cache_hits", "lmm_fit_cache_misses")
  expect_identical(warm$run_manifest[hm], m[hm])
  # Workers return their computations: the same counts as sequential fitting.
  ns$qtl2mr_cache_reset()
  par <- batch(n_workers = 2L)
  expect_identical(par$run_manifest[c("context_computed", "lmm_fit_computed")], m[c("context_computed", "lmm_fit_computed")])
  # A budget too small to keep anything: the run releases entries, and a second call recomputes although its
  # hit counts are those of a warm run (backlog QB-7); the results are unchanged.
  options(qtl2mr.cache_budget = 1)
  ns$qtl2mr_cache_reset()
  tiny1 <- batch()
  expect_true(tiny1$run_manifest$cache_released_entries > 0 && tiny1$run_manifest$cache_released_bytes > 0)
  expect_true(tiny1$run_manifest$lmm_fit_computed >= m$lmm_fit_computed)
  expect_identical(tiny1$run_manifest$cache_budget_bytes, 1)
  tiny2 <- batch()
  expect_true(tiny2$run_manifest$lmm_fit_computed > 0)
  expect_identical(tiny2$run_manifest[hm], m[hm])
  skip_cols <- c(cols, hm, "run_id", "created_at", "wall_sec", "mean_wall_sec_per_pair", "parallel_backend", "n_workers")
  timeless <- function(d) if (is.data.frame(d)) d[, setdiff(names(d), "runtime_sec"), drop = FALSE] else d
  for (x in list(warm, par, tiny1, tiny2)) {
    for (tab in setdiff(names(cold), c("run_manifest", "config"))) {
      expect_true(identical(timeless(x[[tab]]), timeless(cold[[tab]])), info = tab)
    }
    keep <- setdiff(names(m), skip_cols)
    expect_true(identical(x$run_manifest[keep], m[keep]))
  }
})
