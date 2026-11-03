# Session state for reusing exact intermediate results. The state lives in the
# package namespace, never in resources or results, so saveRDS() of user
# objects is unchanged and nothing here outlives the session.
qtl2mr_state <- new.env(parent = emptyenv())

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

qtl2mr_can_fork <- function(n_workers, n_tasks) {
  n_workers > 1L && n_tasks > 1L && identical(.Platform$OS.type, "unix")
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
      out[i] <- list(FUN(X[[i]]))
      after(i)
    }
    names(out) <- names(X)
    return(out)
  }
  index <- seq_along(X)
  chunks <- split(index, (index - 1L) %% min(n_workers, length(X)))
  # Under options(warn >= 2) a warning is an error inside FUN, as in sequential
  # code, so it is not intercepted.
  capture_warnings <- getOption("warn", 0L) < 2L
  run_chunk <- function(idx) {
    assign("n_workers", 1L, envir = qtl2mr_state)
    lapply(idx, function(i) {
      conditions <- list()
      error <- NULL
      value <- tryCatch(withCallingHandlers(
        FUN(X[[i]]),
        warning = function(w) {
          if (!capture_warnings) return(invisible(NULL))
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

