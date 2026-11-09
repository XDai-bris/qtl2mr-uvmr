# Session state for reusing exact intermediate results: trait-level LMM
# contexts and per-locus fits, kinship decompositions, key fingerprints and
# kinship checks. Every entry is keyed by the full content of its inputs, so a
# result is reused only where recomputing it would give the same value. The
# state lives in the package namespace, never in resources or results, so
# saveRDS() of user objects is unchanged and nothing here outlives the session.
qtl2mr_state <- new.env(parent = emptyenv())

# Bound on cached contexts and fits (bytes, as estimated by object.size()):
# option qtl2mr.cache_budget, default 4 GiB (2 GiB before 0.6.1.9005).
# Least recently used entries are released first; a running batch pins the
# entries of the exposure groups it is fitting. The budget changes only what
# is kept for reuse, never a result.
qtl2mr_cache_budget_bytes <- 4 * 1024^3

qtl2mr_cache_budget <- function() {
  if (!is.null(qtl2mr_state$budget)) return(qtl2mr_state$budget)
  b <- getOption("qtl2mr.cache_budget", qtl2mr_cache_budget_bytes)
  if (!is.numeric(b) || length(b) != 1L || is.na(b) || b < 0) {
    stop("option qtl2mr.cache_budget must be one non-negative number of bytes", call. = FALSE)
  }
  as.numeric(b)
}

# Counts kept in this session: contexts and fits actually computed (including
# recomputations of released entries) and registry entries released for the
# budget. A batch reports the change over its run; forked workers return the
# computations they made, which are added here.
qtl2mr_cache_counters <- function() {
  computed <- qtl2mr_state$computed %||% c(context = 0, fit = 0, fit_attempt = 0)
  reg <- qtl2mr_registry()
  c(context = computed[["context"]], fit_attempt = computed[["fit_attempt"]],
    released_entries = reg$released_entries %||% 0, released_bytes = reg$released_bytes %||% 0)
}

qtl2mr_cache_counters_add_computed <- function(delta) {
  if (!is.numeric(delta) || !all(c("context", "fit_attempt") %in% names(delta))) return(invisible(NULL))
  computed <- qtl2mr_state$computed %||% c(context = 0, fit = 0, fit_attempt = 0)
  computed[["context"]] <- computed[["context"]] + delta[["context"]]
  computed[["fit_attempt"]] <- computed[["fit_attempt"]] + delta[["fit_attempt"]]
  assign("computed", computed, envir = qtl2mr_state)
  invisible(NULL)
}

# Forget everything cached in this session (as in a new R process). Internal;
# used by tests.
qtl2mr_cache_reset <- function() {
  rm(list = ls(qtl2mr_state, all.names = TRUE), envir = qtl2mr_state)
  invisible(NULL)
}

qtl2mr_memo_list <- function(name) qtl2mr_state[[name]] %||% list()

qtl2mr_memo_push <- function(name, entry, max_entries) {
  memo <- c(list(entry), qtl2mr_memo_list(name))
  if (length(memo) > max_entries) memo <- memo[seq_len(max_entries)]
  assign(name, memo, envir = qtl2mr_state)
  invisible(entry)
}

# Kinship matrices that passed the matrix-only checks of validation. identical()
# returns at once for the same object, so a resources kinship is checked once.
qtl2mr_kinship_checked <- function(K) {
  for (k in qtl2mr_memo_list("kinship_checked")) if (identical(k, K)) return(TRUE)
  FALSE
}

qtl2mr_kinship_checked_add <- function(K) {
  qtl2mr_memo_push("kinship_checked", K, 64L)
}

# Fingerprint of a kinship matrix restricted to a sample set. The context keys
# include it in place of the matrix itself, so the matrix is serialized once per
# (kinship, sample set) instead of at every key.
qtl2mr_kinship_fingerprint <- function(kinship, sample_ids) {
  if (is.null(kinship)) return(NULL)
  for (e in qtl2mr_memo_list("kinship_fp")) {
    if (identical(e$ids, sample_ids) && identical(e$kinship, kinship)) return(e$fp)
  }
  fp <- qtl2mr_fingerprint_object(as.matrix(kinship)[sample_ids, sample_ids, drop = FALSE])
  qtl2mr_memo_push("kinship_fp", list(kinship = kinship, ids = sample_ids, fp = fp), 512L)
  fp
}

# Eigendecomposition of 2K on a sample set, shared by the hsq fit and the
# transform of one context and by traits fitted on the same samples. The
# decomposition is a deterministic function of the matrix and the samples, so a
# kept decomposition is exactly the one that would be recomputed. Kept
# decompositions are bounded in bytes (option qtl2mr.decomposition_budget,
# default 1 GiB, at least 8 kept); a reused one moves to the front, so a
# selection that visits every chromosome of each exposure in turn reuses them.
qtl2mr_decomposition_budget_bytes <- 1024^3

qtl2mr_decomposition_budget <- function() {
  b <- getOption("qtl2mr.decomposition_budget", qtl2mr_decomposition_budget_bytes)
  if (!is.numeric(b) || length(b) != 1L || is.na(b) || b < 0) {
    stop("option qtl2mr.decomposition_budget must be one non-negative number of bytes", call. = FALSE)
  }
  as.numeric(b)
}

qtl2mr_decompose_kinship_cached <- function(kinship, sample_ids, kinship0) {
  memo <- qtl2mr_memo_list("kinship_decomp")
  for (i in seq_along(memo)) {
    e <- memo[[i]]
    if (identical(e$ids, sample_ids) && identical(e$kinship, kinship)) {
      if (i > 1L) assign("kinship_decomp", c(memo[i], memo[-i]), envir = qtl2mr_state)
      return(e$decomp)
    }
  }
  decomp <- qtl2mr_decompose_kinship(kinship0, sample_ids = sample_ids)
  entry <- list(kinship = kinship, ids = sample_ids, decomp = decomp,
    bytes = as.numeric(utils::object.size(decomp)))
  memo <- c(list(entry), memo)
  bytes <- cumsum(vapply(memo, function(e) e$bytes %||% 0, numeric(1)))
  keep <- max(8L, sum(bytes <= qtl2mr_decomposition_budget()))
  if (length(memo) > keep) memo <- memo[seq_len(keep)]
  assign("kinship_decomp", memo, envir = qtl2mr_state)
  decomp
}

qtl2mr_registry <- function() {
  reg <- qtl2mr_state$registry
  if (is.null(reg)) {
    reg <- new.env(parent = emptyenv())
    reg$entries <- new.env(parent = emptyenv())
    reg$pins <- new.env(parent = emptyenv())
    reg$bytes <- 0
    reg$tick <- 0
    reg$journal <- NULL
    assign("registry", reg, envir = qtl2mr_state)
  }
  reg
}

qtl2mr_registry_key <- function(kind, key) paste0(kind, ":", key)

qtl2mr_registry_get <- function(kind, key) {
  reg <- qtl2mr_registry()
  id <- qtl2mr_registry_key(kind, key)
  entry <- reg$entries[[id]]
  if (is.null(entry)) return(NULL)
  reg$tick <- reg$tick + 1
  entry$tick <- reg$tick
  assign(id, entry, envir = reg$entries)
  entry
}

qtl2mr_registry_put <- function(kind, key, value, conditions = list(), info = list()) {
  reg <- qtl2mr_registry()
  id <- qtl2mr_registry_key(kind, key)
  old <- reg$entries[[id]]
  if (!is.null(old)) reg$bytes <- reg$bytes - old$bytes
  reg$tick <- reg$tick + 1
  entry <- list(value = value, conditions = conditions, info = info,
    bytes = as.numeric(utils::object.size(value)), tick = reg$tick)
  assign(id, entry, envir = reg$entries)
  reg$bytes <- reg$bytes + entry$bytes
  if (is.environment(reg$journal)) assign(id, TRUE, envir = reg$journal)
  if (isTRUE(reg$pin_new)) assign(id, TRUE, envir = reg$pins)
  qtl2mr_registry_evict(reg)
  invisible(entry)
}

qtl2mr_registry_import <- function(entries) {
  reg <- qtl2mr_registry()
  for (id in names(entries)) {
    old <- reg$entries[[id]]
    if (!is.null(old)) reg$bytes <- reg$bytes - old$bytes
    reg$tick <- reg$tick + 1
    entry <- entries[[id]]
    entry$tick <- reg$tick
    assign(id, entry, envir = reg$entries)
    reg$bytes <- reg$bytes + entry$bytes
    if (is.environment(reg$journal)) assign(id, TRUE, envir = reg$journal)
    if (isTRUE(reg$pin_new)) assign(id, TRUE, envir = reg$pins)
  }
  qtl2mr_registry_evict(reg)
  invisible(names(entries))
}

# Internal counts of computed (not reused) LMM contexts and fits, for tests.
qtl2mr_count_computed <- function(kind) {
  counts <- qtl2mr_state$computed %||% c(context = 0, fit = 0, fit_attempt = 0)
  counts[[kind]] <- counts[[kind]] + 1
  assign("computed", counts, envir = qtl2mr_state)
  invisible(NULL)
}

qtl2mr_registry_remove <- function(ids) {
  reg <- qtl2mr_registry()
  ids <- ids[vapply(ids, exists, logical(1), envir = reg$entries, inherits = FALSE)]
  for (id in ids) reg$bytes <- reg$bytes - reg$entries[[id]]$bytes
  if (length(ids)) rm(list = ids, envir = reg$entries)
  invisible(ids)
}

# Over budget, release least recently used unpinned entries down to 90%.
qtl2mr_registry_evict <- function(reg = qtl2mr_registry()) {
  if (reg$bytes <= qtl2mr_cache_budget()) return(invisible(character()))
  ids <- setdiff(ls(reg$entries, all.names = TRUE), ls(reg$pins, all.names = TRUE))
  if (!length(ids)) return(invisible(character()))
  ticks <- vapply(ids, function(id) reg$entries[[id]]$tick, numeric(1))
  bytes <- vapply(ids, function(id) reg$entries[[id]]$bytes, numeric(1))
  ord <- order(ticks)
  excess <- reg$bytes - 0.9 * qtl2mr_cache_budget()
  k <- which(cumsum(bytes[ord]) >= excess)[1L]
  if (is.na(k)) k <- length(ord)
  drop <- ids[ord][seq_len(k)]
  reg$released_entries <- (reg$released_entries %||% 0) + length(drop)
  reg$released_bytes <- (reg$released_bytes %||% 0) + sum(bytes[ord][seq_len(k)])
  qtl2mr_registry_remove(drop)
}

# Entries no longer needed are made the first to be released when the budget is
# exceeded; until then a later call in the session can still reuse them.
qtl2mr_registry_demote <- function(ids) {
  reg <- qtl2mr_registry()
  for (id in intersect(ids, ls(reg$entries, all.names = TRUE))) {
    entry <- reg$entries[[id]]
    entry$tick <- 0
    assign(id, entry, envir = reg$entries)
  }
  invisible(ids)
}

qtl2mr_registry_pin <- function(ids) {
  reg <- qtl2mr_registry()
  for (id in ids) assign(id, TRUE, envir = reg$pins)
  invisible(ids)
}

qtl2mr_registry_unpin <- function(ids) {
  reg <- qtl2mr_registry()
  ids <- intersect(ids, ls(reg$pins, all.names = TRUE))
  if (length(ids)) rm(list = ids, envir = reg$pins)
  invisible(ids)
}

# Entries whose info matches the given (trait, chr, marker) items.
qtl2mr_registry_ids_matching <- function(keep) {
  reg <- qtl2mr_registry()
  ids <- ls(reg$entries, all.names = TRUE)
  ids[vapply(ids, function(id) keep(reg$entries[[id]]$info), logical(1))]
}

# Evaluate expr and record the warnings and messages it signals, in order.
# With muffle = FALSE they still reach the caller's handlers as usual.
qtl2mr_recording <- function(expr, muffle = FALSE) {
  conditions <- list()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      conditions[[length(conditions) + 1L]] <<- w
      if (isTRUE(muffle)) invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (inherits(m, "qtl2mr_progress_message")) return(invisible(NULL))
      conditions[[length(conditions) + 1L]] <<- m
      if (isTRUE(muffle)) invokeRestart("muffleMessage")
    }
  )
  list(value = value, conditions = conditions)
}

# Signal recorded conditions again, as the computation that produced them would.
qtl2mr_replay_conditions <- function(conditions) {
  for (cond in conditions) {
    if (inherits(cond, "warning")) warning(cond) else if (inherits(cond, "message")) message(cond)
  }
  invisible(NULL)
}

qtl2mr_quietly <- function(expr) suppressMessages(suppressWarnings(expr))

# Workers for the refits of one run_qtl2mr_analysis() call; 1 elsewhere and
# inside workers.
qtl2mr_current_workers <- function() as.integer(qtl2mr_state$n_workers %||% 1L)

qtl2mr_set_workers <- function(n_workers) {
  old <- qtl2mr_state$n_workers
  assign("n_workers", n_workers, envir = qtl2mr_state)
  invisible(old)
}

# Under options(warn >= 2) a warning is an error that a caller's handler (for
# example suppressWarnings()) may still muffle; that handler does not exist in
# a forked worker, so refits then run in this process.
qtl2mr_can_fork <- function(n_workers, n_tasks) {
  n_workers > 1L && n_tasks > 1L && identical(.Platform$OS.type, "unix") && getOption("warn", 0L) < 2L
}

# lapply() over X whose calls may run in forked workers. FUN must not draw
# random numbers. Values are returned in order; the warnings and messages of
# each call, and the first error, are signalled in the parent in element order,
# exactly as a sequential lapply() would signal them. after(i) runs in the
# parent once element i has been replayed.
qtl2mr_ordered_lapply <- function(X, FUN, after = NULL) {
  n_workers <- qtl2mr_current_workers()
  if (!qtl2mr_can_fork(n_workers, length(X))) {
    if (!is.function(after)) return(lapply(X, FUN))
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
      out[i] <- lapply(X[i], FUN)
      after(i)
    }
    names(out) <- names(X)
    return(out)
  }
  index <- seq_along(X)
  chunks <- split(index, (index - 1L) %% min(n_workers, length(X)))
  run_chunk <- function(idx) {
    assign("n_workers", 1L, envir = qtl2mr_state)
    lapply(idx, function(i) {
      conditions <- list()
      error <- NULL
      # Through lapply(), so that conditions raised by FUN carry the same call
      # as in a sequential lapply().
      value <- tryCatch(withCallingHandlers(
        lapply(X[i], FUN)[[1L]],
        warning = function(w) {
          conditions[[length(conditions) + 1L]] <<- w
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          conditions[[length(conditions) + 1L]] <<- m
          invokeRestart("muffleMessage")
        }
      ), error = function(e) {
        error <<- e
        NULL
      })
      list(value = value, conditions = conditions, error = error)
    })
  }
  raw <- parallel::mclapply(chunks, run_chunk, mc.cores = length(chunks),
    mc.preschedule = FALSE, mc.set.seed = FALSE)
  failed <- !vapply(raw, is.list, logical(1)) | vapply(raw, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("qtl2mr worker failed: ", paste(as.character(raw[failed][[1L]]), collapse = " "), call. = FALSE)
  }
  results <- vector("list", length(X))
  for (k in seq_along(chunks)) results[chunks[[k]]] <- raw[[k]]
  out <- vector("list", length(X))
  for (i in index) {
    r <- results[[i]]
    qtl2mr_replay_conditions(r$conditions)
    if (!is.null(r$error)) stop(r$error)
    out[i] <- list(r$value)
    if (is.function(after)) after(i)
  }
  names(out) <- names(X)
  out
}

# Contexts are stored without their n x n matrices; a context rebuilt from the
# same arguments is identical, so the matrices are recomputed only when a fit
# is missing.
qtl2mr_context_large_fields <- c("raw_transform", "star_rotation", "inverse_residual", "kinship_decomp", "kinship")

qtl2mr_compact_context <- function(ctx) {
  ctx[intersect(qtl2mr_context_large_fields, names(ctx))] <- NULL
  ctx
}

qtl2mr_full_context <- function(context) {
  if (!is.null(context$raw_transform)) return(context)
  if (!is.function(context$materialise)) stop("Internal error: LMM context cannot be rebuilt", call. = FALSE)
  context$materialise()
}
