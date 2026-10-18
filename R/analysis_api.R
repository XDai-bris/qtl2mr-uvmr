qtl2mr_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "-", as.character(x))
  x <- gsub("^-|-$", "", x)
  if (!nzchar(x)) "value" else x
}

qtl2mr_analysis_output_dir <- function(exposure, outcome, output_dir = NULL) {
  if (!is.null(output_dir)) return(normalizePath(output_dir, mustWork = FALSE))
  nonce <- basename(tempfile("qtl2mr-"))
  file.path(
    qtl2mr_default_output_root(),
    "results",
    sprintf("result_exp-%s_outcome-%s_%s", qtl2mr_slug(exposure), qtl2mr_slug(outcome), nonce)
  )
}

qtl2mr_fmt_num <- function(x, digits = 4) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg", flag = "#"), "NA")
}

qtl2mr_loci_label <- function(loci) {
  if (is.null(loci)) return("all")
  if (is.numeric(loci) && length(loci) == 1L) return(paste0("top", as.integer(loci)))
  paste0(length(loci), " requested")
}

qtl2mr_add_q_fields <- function(summary, fit = NULL) {
  if (is.null(summary) || !nrow(summary)) return(summary)
  old_q_suffix <- paste0("ap", "prox")
  deprecated <- intersect(c(
    paste0("Q_df_", old_q_suffix),
    paste0("Q_p_", old_q_suffix),
    paste0("Q_p_", old_q_suffix, "_note")
  ), names(summary))
  if (length(deprecated)) summary[deprecated] <- NULL
  fields <- list(
    Q_chisq_df = NA_integer_,
    Q_p_chisq = NA_real_,
    Q_p_sim = NA_real_,
    Q_sim_n = NA_integer_,
    Q_sim_n_requested = NA_integer_,
    Q_sim_seed = NA_integer_,
    Q_sim_exceedances = NA_integer_,
    Q_sim_min_resolvable_p = NA_real_,
    Q_sim_null_mean = NA_real_,
    Q_sim_null_sd = NA_real_,
    Q_sim_method = NA_character_,
    Q_p_reported = NA_real_,
    Q_p_method = NA_character_,
    Q_p_note = NA_character_
  )
  for (nm in names(fields)) {
    if (!nm %in% names(summary)) summary[[nm]] <- fields[[nm]]
  }
  if (!is.null(fit)) {
    for (nm in names(fields)) {
      if (!is.null(fit[[nm]])) summary[[nm]][1L] <- fit[[nm]]
    }
  }
  if ("Q" %in% names(summary) && "Q_chisq_df" %in% names(summary)) {
    can_chisq <- is.finite(summary$Q) & is.finite(summary$Q_chisq_df) & summary$Q_chisq_df > 0
    fill <- is.na(summary$Q_p_chisq) & can_chisq
    summary$Q_p_chisq[fill] <- stats::pchisq(summary$Q[fill], df = summary$Q_chisq_df[fill], lower.tail = FALSE)
  }
  if ("Q_p_method" %in% names(summary) && "Q_p_reported" %in% names(summary)) {
    use_chisq <- identical(as.character(summary$Q_p_method[[1L]]), "chisq") &&
      !is.finite(summary$Q_p_reported[[1L]]) &&
      is.finite(summary$Q_p_chisq[[1L]])
    if (isTRUE(use_chisq)) {
      summary$Q_p_reported[1L] <- summary$Q_p_chisq[[1L]]
      summary$Q_p_note[1L] <- paste0("chi-square reference with df=", summary$Q_chisq_df[[1L]])
    }
    use_sim <- identical(as.character(summary$Q_p_method[[1L]]), "simulation") &&
      !is.finite(summary$Q_p_reported[[1L]]) &&
      is.finite(summary$Q_p_sim[[1L]])
    if (isTRUE(use_sim)) {
      summary$Q_p_reported[1L] <- summary$Q_p_sim[[1L]]
      summary$Q_p_note[1L] <- "conditional residual simulation under the fitted common-effect null"
    }
  }
  summary
}

qtl2mr_normalize_progress <- function(progress, verbose = TRUE) {
  if (is.null(progress)) progress <- isTRUE(verbose)
  if (is.logical(progress) && length(progress) == 1L && !is.na(progress)) {
    return(list(console = isTRUE(progress), log = isTRUE(progress)))
  }
  if (is.character(progress) && length(progress) == 1L && identical(tolower(progress), "log")) {
    return(list(console = FALSE, log = TRUE))
  }
  stop("progress must be TRUE, FALSE, \"log\", or NULL", call. = FALSE)
}

qtl2mr_make_progress_reporter <- function(progress = NULL,
                                          verbose = TRUE,
                                          log_file = NULL) {
  settings <- qtl2mr_normalize_progress(progress, verbose = verbose)
  if (!isTRUE(settings$console) && !isTRUE(settings$log)) return(NULL)
  if (isTRUE(settings$log)) {
    if (is.null(log_file)) stop("log_file is required when progress logging is enabled", call. = FALSE)
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(sprintf("# qtl2mr progress log\n# started\t%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), log_file)
  }
  function(stage, label, detail = NULL) {
    line <- sprintf(
      "[qtl2mr %s] [%s/5] %s%s",
      format(Sys.time(), "%H:%M:%S"),
      as.integer(stage),
      label,
      if (!is.null(detail) && nzchar(detail)) paste0(": ", detail) else ""
    )
    if (isTRUE(settings$console)) {
      condition <- simpleMessage(line)
      class(condition) <- c("qtl2mr_progress_message", class(condition))
      message(condition)
    }
    if (isTRUE(settings$log)) cat(line, "\n", file = log_file, append = TRUE)
    invisible(line)
  }
}

qtl2mr_progress <- function(callback, stage, label, detail = NULL) {
  if (is.function(callback)) {
    # Preserve visible RNG state defensively. Callbacks must still be RNG-free:
    # some generators (e.g. Box-Muller) keep a cache outside .Random.seed.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    saved_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    saved_kind <- RNGkind()
    on.exit({
      if (!identical(RNGkind(), saved_kind)) do.call(RNGkind, as.list(saved_kind))
      if (had_seed) {
        assign(".Random.seed", saved_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    callback(stage = stage, label = label, detail = detail)
  }
  invisible(NULL)
}

qtl2mr_capture_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      if (inherits(m, "qtl2mr_progress_message")) return(invisible(NULL))
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

qtl2mr_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

qtl2mr_q_sim_label <- function(p, exceedances = NA_integer_, min_p = NA_real_, digits = 4) {
  if (!is.finite(p)) return("NA")
  qtl2mr_fmt_num(p, digits = digits)
}

qtl2mr_q_reported_label <- function(summary_row, digits = 4) {
  method <- summary_row$Q_p_method[[1L]] %||% NA_character_
  if (identical(method, "simulation")) {
    return(qtl2mr_q_sim_label(
      summary_row$Q_p_reported[[1L]],
      summary_row$Q_sim_exceedances[[1L]] %||% NA_integer_,
      summary_row$Q_sim_min_resolvable_p[[1L]] %||% NA_real_,
      digits = digits
    ))
  }
  qtl2mr_fmt_num(summary_row$Q_p_reported[[1L]], digits = digits)
}

qtl2mr_q_reference_label <- function(summary_row) {
  method <- qtl2mr_summary_scalar(summary_row, "Q_p_method", "not recorded")
  if (identical(method, "chisq")) {
    return(paste0("chi-square; df=", qtl2mr_summary_scalar(summary_row, "Q_chisq_df")))
  }
  if (identical(method, "simulation")) {
    return(paste0("conditional residual simulation; ",
      qtl2mr_summary_scalar(summary_row, "Q_sim_n"), "/",
      qtl2mr_summary_scalar(summary_row, "Q_sim_n_requested"), " completed draws"))
  }
  if (identical(method, "none")) "not requested" else as.character(method)
}

qtl2mr_summary_scalar <- function(summary_row, name, default = NA_real_) {
  if (is.null(summary_row) || !nrow(summary_row) || !name %in% names(summary_row)) return(default)
  summary_row[[name]][[1L]]
}

qtl2mr_summary_reported <- function(summary_row) {
  list(
    theta = qtl2mr_summary_scalar(summary_row, "theta_reported", qtl2mr_summary_scalar(summary_row, "theta")),
    se = qtl2mr_summary_scalar(summary_row, "se_reported", qtl2mr_summary_scalar(summary_row, "se")),
    ci_lower = qtl2mr_summary_scalar(summary_row, "reported_ci_lower", qtl2mr_summary_scalar(summary_row, "profile_ci_lower")),
    ci_upper = qtl2mr_summary_scalar(summary_row, "reported_ci_upper", qtl2mr_summary_scalar(summary_row, "profile_ci_upper")),
    ci_type = qtl2mr_summary_scalar(summary_row, "reported_ci_type", "profile"),
    estimator = qtl2mr_summary_scalar(summary_row, "reported_estimator", "Star-space qtl2mr"),
    source = qtl2mr_summary_scalar(summary_row, "reported_estimator_source", "primary"),
    policy = qtl2mr_summary_scalar(summary_row, "reported_estimator_policy", "strict_so"),
    reason = qtl2mr_summary_scalar(summary_row, "reported_estimator_reason", "")
  )
}

qtl2mr_html_table <- function(x, digits = 4, scientific = FALSE) {
  if (is.null(x) || !nrow(x)) return("<p>No rows.</p>")
  x[] <- lapply(x, function(col) {
    if (!is.numeric(col)) return(as.character(col))
    if (!isTRUE(scientific)) return(qtl2mr_fmt_num(col, digits = digits))
    out <- rep(NA_character_, length(col))
    finite <- is.finite(col)
    out[finite] <- formatC(col[finite], format = "e", digits = digits)
    out[!finite & !is.na(col)] <- as.character(col[!finite & !is.na(col)])
    out
  })
  header <- paste0("<tr>", paste0("<th>", qtl2mr_html_escape(names(x)), "</th>", collapse = ""), "</tr>")
  rows <- apply(x, 1L, function(z) paste0("<tr>", paste0("<td>", qtl2mr_html_escape(z), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead>", header, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

qtl2mr_html_output_table <- function(x, digits = 4) {
  if (is.null(x) || !nrow(x)) return("<p>No rows.</p>")
  x0 <- x
  x0[] <- lapply(x0, function(col) if (is.numeric(col)) qtl2mr_fmt_num(col, digits = digits) else as.character(col))
  header <- paste0("<tr>", paste0("<th>", qtl2mr_html_escape(names(x0)), "</th>", collapse = ""), "</tr>")
  rows <- apply(x0, 1L, function(z) {
    names(z) <- names(x0)
    cells <- vapply(names(z), function(nm) {
      val <- z[[nm]]
      if (identical(nm, "path_relative") && nzchar(val)) {
        sprintf("<td><a href=\"%s\">%s</a></td>", qtl2mr_html_escape(val), qtl2mr_html_escape(val))
      } else {
        paste0("<td>", qtl2mr_html_escape(val), "</td>")
      }
    }, character(1))
    paste0("<tr>", paste(cells, collapse = ""), "</tr>")
  })
  paste0("<table><thead>", header, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

qtl2mr_write_output_manifest_html <- function(manifest, file) {
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\"><title>qtl2mr Output Manifest</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;color:#1d2528;line-height:1.45;}main{max-width:1180px;margin:auto;}table{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0 22px;}th,td{border:1px solid #d9dfdd;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#eef3f1;}code{background:#eef3f1;padding:1px 4px;}</style>",
    "</head><body><main>",
    "<h1>qtl2mr Output Manifest</h1>",
    "<p>This page lists every file written by the qtl2mr destination workflow. The main report shows only a compact summary.</p>",
    qtl2mr_html_output_table(manifest[, intersect(c("role", "path_relative", "format", "size_kb", "included_in_lite", "notes"), names(manifest)), drop = FALSE], digits = 5),
    "</main></body></html>"
  )
  writeLines(html, file)
  invisible(file)
}

qtl2mr_relative_path <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  root <- sub("/+$", "", root)
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}

qtl2mr_file_size_kb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  round(file.info(path)$size / 1024, 1)
}

qtl2mr_lite_inclusion <- function(path, output_dir, max_file_mb = 10, include_rds = FALSE) {
  rel <- qtl2mr_relative_path(path, output_dir)
  ext <- tolower(tools::file_ext(path))
  size_mb <- qtl2mr_file_size_kb(path) / 1024
  if (identical(ext, "rds") && !isTRUE(include_rds)) return(FALSE)
  if (is.finite(size_mb) && size_mb > max_file_mb) return(FALSE)
  if (startsWith(rel, "logs/")) return(FALSE)
  if (basename(path) %in% c("output_manifest.csv", "output_manifest.html")) return(FALSE)
  TRUE
}

qtl2mr_path_node_info <- function(path) fs::file_info(path, follow = FALSE)

qtl2mr_path_parent_entries <- function(path) {
  fs::dir_ls(path, all = TRUE, recurse = FALSE, type = "any", fail = TRUE)
}

qtl2mr_trim_node_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) return(NA_character_)
  trimmed <- sub("[/\\\\]+$", "", path)
  if (!nzchar(trimmed) || grepl("^[A-Za-z]:$", trimmed)) path else trimmed
}

qtl2mr_path_node_state <- function(path) {
  path <- qtl2mr_trim_node_path(path)
  if (is.na(path)) return("unknown")
  info <- tryCatch(qtl2mr_path_node_info(path), error = function(e) NULL)
  if (is.null(info) || !is.data.frame(info) || nrow(info) != 1L || !"type" %in% names(info)) return("unknown")
  type <- as.character(info$type[[1L]])
  if (!is.na(type)) {
    if (identical(type, "symlink")) return("link")
    if (identical(type, "directory")) return("directory")
    if (identical(type, "file")) return("file")
    return("unknown")
  }

  parent <- dirname(path)
  parent_info <- tryCatch(qtl2mr_path_node_info(parent), error = function(e) NULL)
  if (is.null(parent_info) || !is.data.frame(parent_info) || nrow(parent_info) != 1L ||
      !"type" %in% names(parent_info) || !identical(as.character(parent_info$type[[1L]]), "directory")) {
    return("unknown")
  }
  entries <- tryCatch(qtl2mr_path_parent_entries(parent), error = function(e) NULL)
  if (is.null(entries)) return("unknown")
  entries <- as.character(entries)
  target <- basename(path)
  names_equal <- if (identical(.Platform$OS.type, "windows")) {
    tolower(basename(entries)) == tolower(target)
  } else {
    basename(entries) == target
  }
  if (any(names_equal, na.rm = TRUE)) "unknown" else "absent"
}

qtl2mr_assert_known_path_node <- function(path, purpose = "output path") {
  state <- qtl2mr_path_node_state(path)
  if (identical(state, "unknown")) {
    stop("Unable to verify ", purpose, " safety; refusing filesystem mutation", call. = FALSE)
  }
  state
}

qtl2mr_ensure_real_directory <- function(path, purpose = "output path") {
  path <- qtl2mr_trim_node_path(path)
  if (is.na(path)) stop("Unable to verify ", purpose, " safety; refusing filesystem mutation", call. = FALSE)
  state <- qtl2mr_path_node_state(path)
  parent <- dirname(path)
  if (state %in% c("directory", "link")) {
    canonical <- tryCatch(normalizePath(path, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
    if (is.na(canonical) || !identical(qtl2mr_assert_known_path_node(canonical, purpose), "directory")) {
      stop(purpose, " must resolve to a real writable directory", call. = FALSE)
    }
    return(invisible(canonical))
  }
  if (identical(parent, path)) {
    stop("Unable to verify ", purpose, " safety; refusing filesystem mutation", call. = FALSE)
  }
  canonical_parent <- qtl2mr_ensure_real_directory(parent, purpose)
  canonical <- file.path(canonical_parent, basename(path))
  state <- qtl2mr_assert_known_path_node(canonical, purpose)
  if (identical(state, "directory")) return(invisible(canonical))
  if (!identical(state, "absent") || !dir.create(canonical, recursive = FALSE, showWarnings = FALSE)) {
    state <- qtl2mr_assert_known_path_node(canonical, purpose)
    if (identical(state, "directory")) return(invisible(canonical))
    stop("Could not create ", purpose, call. = FALSE)
  }
  if (!identical(qtl2mr_assert_known_path_node(canonical, purpose), "directory")) {
    stop("Could not create ", purpose, call. = FALSE)
  }
  invisible(canonical)
}

qtl2mr_tree_has_symlink <- function(path) {
  state <- qtl2mr_assert_known_path_node(path, "result tree")
  if (identical(state, "link")) return(TRUE)
  if (!identical(state, "directory")) return(FALSE)
  entries <- tryCatch(qtl2mr_path_parent_entries(path), error = function(e) NULL)
  if (is.null(entries)) stop("Unable to verify result tree safety; refusing filesystem traversal", call. = FALSE)
  for (entry in as.character(entries)) {
    entry_state <- qtl2mr_assert_known_path_node(entry, "result tree")
    if (identical(entry_state, "link")) return(TRUE)
    if (identical(entry_state, "directory") && qtl2mr_tree_has_symlink(entry)) return(TRUE)
  }
  FALSE
}

qtl2mr_path_occupied <- function(path) {
  !identical(qtl2mr_assert_known_path_node(path), "absent")
}

qtl2mr_is_regular_file <- function(path) {
  identical(qtl2mr_path_node_state(path), "file")
}

qtl2mr_is_prior_result_dir <- function(path) {
  required <- c("result_object.rds", "summary.md", "data/output_manifest.csv", "data/main_result.csv", "data/warnings.csv")
  if (!identical(qtl2mr_assert_known_path_node(path, "prior result"), "directory") ||
      !all(vapply(file.path(path, required), qtl2mr_is_regular_file, logical(1))) || qtl2mr_tree_has_symlink(path)) return(FALSE)
  manifest <- tryCatch(
    utils::read.csv(file.path(path, "data", "output_manifest.csv"), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  !is.null(manifest) && all(c("role", "path_relative", "path_absolute") %in% names(manifest)) &&
    all(vapply(manifest$path_relative, qtl2mr_is_safe_relative_path, logical(1))) &&
    all(c("result_object.rds", "summary.md", "data/main_result.csv", "data/warnings.csv") %in% manifest$path_relative)
}

qtl2mr_write_output_transaction_token <- function(stage_dir, token) {
  writeLines(token, file.path(stage_dir, ".qtl2mr-transaction-token"))
  invisible(NULL)
}

qtl2mr_remove_output_transaction_backup <- function(backup_dir) {
  backup_state <- qtl2mr_assert_known_path_node(backup_dir, "transaction backup")
  if (identical(backup_state, "directory") && unlink(backup_dir, recursive = TRUE, force = TRUE) != 0L) {
    stop("Could not remove qtl2mr backup after final-path rewrites", call. = FALSE)
  }
  invisible(NULL)
}

qtl2mr_remove_output_transaction_token <- function(token_file) {
  token_state <- qtl2mr_assert_known_path_node(token_file, "transaction token")
  identical(token_state, "absent") || (identical(token_state, "file") && unlink(token_file, force = TRUE) == 0L)
}

qtl2mr_transaction_path <- function(parent, base, kind) {
  path <- tempfile(paste0(".", base, ".qtl2mr-", kind, "-"), tmpdir = parent)
  if (!identical(qtl2mr_assert_known_path_node(path, "transaction path"), "absent")) {
    stop("Unable to allocate a unique qtl2mr transaction path", call. = FALSE)
  }
  path
}

qtl2mr_begin_output_transaction <- function(output_dir, overwrite = FALSE) {
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty path", call. = FALSE)
  }
  requested_destination <- qtl2mr_trim_node_path(output_dir)
  if (identical(qtl2mr_path_node_state(requested_destination), "link")) {
    stop("output_dir must not be a symlink", call. = FALSE)
  }
  parent <- qtl2mr_ensure_real_directory(dirname(requested_destination), "output_dir parent")
  destination <- file.path(parent, basename(requested_destination))
  destination_state <- qtl2mr_assert_known_path_node(destination)
  if (identical(destination_state, "link")) stop("output_dir must not be a symlink", call. = FALSE)
  if (identical(destination_state, "file")) {
    stop("output_dir exists but is not a directory", call. = FALSE)
  }
  destination <- normalizePath(destination, winslash = "/", mustWork = FALSE)
  existing <- identical(destination_state, "directory")
  nonempty <- existing && length(list.files(destination, all.files = TRUE, no.. = TRUE)) > 0L
  if (nonempty && !isTRUE(overwrite)) {
    stop("output_dir exists and is not empty; use overwrite = TRUE or choose a new output_dir", call. = FALSE)
  }
  if (nonempty && !qtl2mr_is_prior_result_dir(destination)) {
    stop("overwrite = TRUE only permits a recognized prior qtl2mr result directory", call. = FALSE)
  }
  base <- basename(destination)
  stage_dir <- qtl2mr_transaction_path(parent, base, "stage")
  backup_dir <- NULL
  if (existing) {
    backup_dir <- qtl2mr_transaction_path(parent, base, "backup")
    if (!file.rename(destination, backup_dir)) {
      stop("Could not preserve the prior qtl2mr result for replacement", call. = FALSE)
    }
  }
  if (!dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)) {
    if (!is.null(backup_dir) && identical(qtl2mr_assert_known_path_node(backup_dir, "transaction backup"), "directory") && !file.rename(backup_dir, destination)) {
      stop("Could not restore prior result after staging setup failure", call. = FALSE)
    }
    stop("Could not create qtl2mr output staging directory", call. = FALSE)
  }
  token <- paste0("qtl2mr-transaction-", basename(stage_dir))
  token_ok <- tryCatch({
    qtl2mr_write_output_transaction_token(stage_dir, token)
    file.exists(file.path(stage_dir, ".qtl2mr-transaction-token"))
  }, error = function(e) FALSE)
  if (!isTRUE(token_ok)) {
    if (!is.null(backup_dir) && identical(qtl2mr_assert_known_path_node(backup_dir, "transaction backup"), "directory")) {
      if (!file.rename(backup_dir, destination)) {
        stop("Could not create qtl2mr transaction token; prior result remains recoverable at ", backup_dir,
             " and staging remains at ", stage_dir, call. = FALSE)
      }
    }
    stop("Could not create qtl2mr output transaction token; staging remains at ", stage_dir, call. = FALSE)
  }
  transaction <- new.env(parent = emptyenv())
  transaction$destination <- destination
  transaction$stage_dir <- stage_dir
  transaction$backup_dir <- backup_dir
  transaction$token <- token
  transaction$state <- "staging"
  class(transaction) <- "qtl2mr_output_transaction"
  transaction
}

qtl2mr_validate_transaction <- function(transaction) {
  if (!inherits(transaction, "qtl2mr_output_transaction") || !is.environment(transaction)) {
    stop("Invalid internal qtl2mr output transaction", call. = FALSE)
  }
  required_paths <- c("destination", "stage_dir")
  valid_state <- is.character(transaction$state) && length(transaction$state) == 1L &&
    !is.na(transaction$state) && transaction$state %in% c("staging", "published", "finalized", "aborted")
  valid_token <- is.character(transaction$token) && length(transaction$token) == 1L &&
    !is.na(transaction$token) && nzchar(transaction$token)
  valid_paths <- all(vapply(required_paths, function(name) {
    value <- transaction[[name]]
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  }, logical(1)))
  if (!valid_state || !valid_token || !valid_paths) {
    stop("Invalid internal qtl2mr output transaction", call. = FALSE)
  }
  invisible(transaction)
}

qtl2mr_transaction_owned <- function(path, transaction) {
  token_file <- file.path(path, ".qtl2mr-transaction-token")
  is.character(path) && length(path) == 1L && identical(qtl2mr_path_node_state(path), "directory") &&
    identical(qtl2mr_path_node_state(token_file), "file") && identical(
      tryCatch(readLines(token_file, warn = FALSE), error = function(e) character()),
      transaction$token
    )
}

qtl2mr_stop_output_transaction_conflict <- function(transaction) {
  if (qtl2mr_transaction_owned(transaction$stage_dir, transaction)) {
    unlink(transaction$stage_dir, recursive = TRUE, force = TRUE)
  }
  transaction$state <- "aborted"
  has_backup <- !is.null(transaction$backup_dir) && identical(qtl2mr_assert_known_path_node(transaction$backup_dir, "transaction backup"), "directory")
  recovery <- if (has_backup) paste0("; prior result remains recoverable at ", transaction$backup_dir) else ""
  condition <- structure(
    list(
      message = paste0("Output transaction conflict: external destination preserved", recovery),
      call = NULL,
      destination = transaction$destination,
      backup_dir = transaction$backup_dir
    ),
    class = c("qtl2mr_output_transaction_conflict", "error", "condition")
  )
  stop(condition)
}

qtl2mr_abort_output_transaction <- function(transaction) {
  qtl2mr_validate_transaction(transaction)
  if (identical(transaction$state, "finalized") || identical(transaction$state, "aborted")) return(invisible(NULL))
  if (qtl2mr_transaction_owned(transaction$stage_dir, transaction)) {
    if (unlink(transaction$stage_dir, recursive = TRUE, force = TRUE) != 0L) {
      stop("Could not remove owned qtl2mr staging directory", call. = FALSE)
    }
  }
  has_backup <- !is.null(transaction$backup_dir) && identical(qtl2mr_assert_known_path_node(transaction$backup_dir, "transaction backup"), "directory")
  if (qtl2mr_path_occupied(transaction$destination)) {
    owned_publication <- identical(transaction$state, "published") &&
      qtl2mr_transaction_owned(transaction$destination, transaction)
    if (!owned_publication) {
      qtl2mr_stop_output_transaction_conflict(transaction)
    }
    if (unlink(transaction$destination, recursive = TRUE, force = TRUE) != 0L) {
      stop("Could not remove owned published qtl2mr result during rollback", call. = FALSE)
    }
  }
  if (has_backup) {
    if (!file.rename(transaction$backup_dir, transaction$destination)) {
      stop("Could not restore prior qtl2mr result from backup: ", transaction$backup_dir, call. = FALSE)
    }
  }
  transaction$state <- "aborted"
  invisible(NULL)
}

qtl2mr_commit_output_transaction <- function(transaction, keep_backup = FALSE) {
  qtl2mr_validate_transaction(transaction)
  if (!identical(transaction$state, "staging") || !qtl2mr_transaction_owned(transaction$stage_dir, transaction)) {
    stop("qtl2mr staging directory is unavailable for publication", call. = FALSE)
  }
  if (qtl2mr_path_occupied(transaction$destination)) {
    qtl2mr_stop_output_transaction_conflict(transaction)
  }
  if (!file.rename(transaction$stage_dir, transaction$destination)) {
    stop("Could not publish staged qtl2mr result", call. = FALSE)
  }
  transaction$state <- "published"
  # Retain any backup until all final-path rewrites have completed. The
  # keep_backup argument is retained for internal-call compatibility only.
  invisible(transaction$destination)
}

qtl2mr_finalize_output_transaction <- function(transaction) {
  qtl2mr_validate_transaction(transaction)
  if (!identical(transaction$state, "published") || !qtl2mr_transaction_owned(transaction$destination, transaction)) {
    stop("Cannot finalize an unowned qtl2mr publication", call. = FALSE)
  }
  if (!is.null(transaction$backup_dir) && identical(qtl2mr_assert_known_path_node(transaction$backup_dir, "transaction backup"), "directory")) {
    qtl2mr_remove_output_transaction_backup(transaction$backup_dir)
  }
  token_file <- file.path(transaction$destination, ".qtl2mr-transaction-token")
  if (!qtl2mr_remove_output_transaction_token(token_file)) {
    transaction$state <- "finalized"
    warning("qtl2mr result was published, but transaction-token cleanup failed; the result remains valid", call. = FALSE)
    return(invisible(NULL))
  }
  transaction$state <- "finalized"
  invisible(NULL)
}

qtl2mr_is_safe_relative_path <- function(path) {
  is.character(path) && length(path) == 1L && nzchar(path) &&
    !grepl("^(?:[A-Za-z]:[\\\\/]|/|file:)", path, perl = TRUE, ignore.case = TRUE) &&
    !any(strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1L]] %in% c("", ".", ".."))
}

qtl2mr_portable_absolute_path_pattern <- function() {
  "(?<![[:alnum:].:/_-])/(?:Users|private/tmp|tmp|var/folders|home|root|var/tmp|Volumes|mnt|srv|opt|data|scratch)(?:/|$)"
}

qtl2mr_scrub_portable_text <- function(lines, result_dir) {
  lines <- gsub(result_dir, ".", lines, fixed = TRUE)
  lines <- gsub("(?i)file://[^[:space:]\\\"'<>]+", "[redacted-absolute-path]", lines, perl = TRUE)
  lines <- gsub("(?i)(?<![[:alnum:]])[A-Z]:[\\\\/][^[:space:]\\\"'<>]+", "[redacted-absolute-path]", lines, perl = TRUE)
  lines <- gsub("\\\\\\\\[^[:space:]\\\"'<>]+", "[redacted-absolute-path]", lines, perl = TRUE)
  gsub(qtl2mr_portable_absolute_path_pattern(), "[redacted-absolute-path]/", lines, perl = TRUE)
}

qtl2mr_has_portable_path_leak <- function(text, result_dir, lite_root) {
  grepl(result_dir, text, fixed = TRUE) || grepl(lite_root, text, fixed = TRUE) ||
    grepl("(?i)(file://|(?<![[:alnum:]])[A-Z]:[\\\\/]|\\\\\\\\[^\\\\/]+[\\\\/])", text, perl = TRUE) ||
    grepl(qtl2mr_portable_absolute_path_pattern(), text, perl = TRUE)
}

qtl2mr_portable_manifest <- function(output_dir, ...) {
  manifest <- qtl2mr_build_output_manifest(output_dir, ...)
  if (nrow(manifest) && !all(vapply(manifest$path_relative, qtl2mr_is_safe_relative_path, logical(1)))) {
    stop("Portable manifest contains an unsafe path", call. = FALSE)
  }
  manifest[, setdiff(names(manifest), "path_absolute"), drop = FALSE]
}

qtl2mr_build_output_manifest <- function(output_dir,
                                         roles = NULL,
                                         notes = NULL,
                                         max_file_mb = 10,
                                         include_rds = FALSE) {
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- normalizePath(files, mustWork = FALSE)
  rel <- vapply(files, qtl2mr_relative_path, character(1), root = output_dir)
  role_from_rel <- function(x) {
    if (identical(x, "summary.md")) return("summary")
    if (identical(x, "report.html")) return("html_report")
    if (identical(x, "result_object.rds")) return("model_object")
    if (startsWith(x, "figures/")) return("figure")
    if (startsWith(x, "logs/")) return("log")
    if (startsWith(x, "data/")) {
      if (basename(x) == "output_manifest.csv") return("manifest")
      if (basename(x) == "session_info.txt") return("session_info")
      return("data")
    }
    "auxiliary_file"
  }
  default_roles <- vapply(rel, role_from_rel, character(1))
  key <- tools::file_path_sans_ext(basename(rel))
  role <- default_roles
  if (!is.null(roles)) {
    role_lookup <- unname(roles[key])
    role[!is.na(role_lookup)] <- role_lookup[!is.na(role_lookup)]
  }
  note <- rep("", length(rel))
  if (!is.null(notes)) {
    note_lookup <- unname(notes[key])
    note[!is.na(note_lookup)] <- note_lookup[!is.na(note_lookup)]
  }
  data.frame(
    role = unname(role),
    path_relative = rel,
    path_absolute = files,
    format = tolower(tools::file_ext(files)),
    size_kb = vapply(files, qtl2mr_file_size_kb, numeric(1)),
    included_in_lite = vapply(files, qtl2mr_lite_inclusion, logical(1),
                              output_dir = output_dir, max_file_mb = max_file_mb, include_rds = include_rds),
    notes = unname(note),
    stringsAsFactors = FALSE
  )
}

qtl2mr_warning_catalog <- function() {
  catalog <- data.frame(
    warning_id = c(
      "prepared_overlap_engine", "raw_space_psd_invalid",
      "star_space_psd_invalid", "weak_block_first_stage", "primary_so_unbounded",
      "loo_mr_error", "instrument_sensitivity_error", "dropped_nonestimable_loci"
    ),
    default_severity = c(
      "info", "technical_qc",
      "warning", "warning", "warning",
      "warning", "warning", "warning"
    ),
    shown_in_main_report = c(
      FALSE, FALSE,
      TRUE, TRUE, TRUE,
      TRUE, TRUE, TRUE
    ),
    consequence = c(
      "Exposure and outcome use their own complete samples; overlap is derived from harmonized sample IDs.",
      "Technical QC: raw-space candidate covariance is not used by the public estimator.",
      "The selected star-space covariance block may be numerically invalid and should be investigated before interpretation.",
      "Founder-block denominator strength may be limited; this is a heuristic warning for multi-founder blocks.",
      "The theta-dependent second-order objective was not well bracketed or has an unbounded profile interval; the primary theta may be a weak-identification artifact.",
      "Leave-one-locus-out MR sensitivity did not complete for this run.",
      "Instrument-set sensitivity did not complete for this run.",
      "Some selected loci were not estimable in the current sample set because a founder contrast was rank-deficient."
    ),
    recommended_action = c(
      "Check data/overlap_qc.csv and logs/run.log to confirm overlap rates.",
      "Technical QC: raw-space candidate invalid; not used by primary estimator.",
      "Inspect PSD diagnostics and do not use invalid star-space blocks silently.",
      "Weak denominator guide: all retained loci have block F below 10. This is a heuristic warning for founder-block denominators, not a Stock-Yogo calibrated failure. Interpret per-locus ratios and weak-instrument-sensitive diagnostics cautiously.",
      "Treat the primary second-order theta as unstable. Inspect Q(theta), profile CI, first-stage strength, and diagnostic comparators such as Two-step IVW before interpretation.",
      "Inspect logs/loo_mr_error.log; primary MR estimate may still be available.",
      "Inspect logs/instrument_sensitivity_error.log; primary MR estimate may still be available.",
      "Inspect dropped_loci output. This usually means one founder has no information at that marker within the pair-specific sample set; do not impute artificial founder probability."
    ),
    stringsAsFactors = FALSE
  )
  rbind(catalog, data.frame(
    warning_id = c("incomplete_covariates", "probability_rounding_normalized", "diagnostic_unavailable",
      "hsq_boundary", "no_kinship_ols", "batch_shared_samples"),
    default_severity = c(rep("warning", 3L), rep("info", 3L)),
    shown_in_main_report = c(rep(TRUE, 5L), FALSE),
    consequence = c(
      "Samples without a covariate row or with non-finite additive covariates were excluded before fitting.",
      "Small founder-probability sum errors were normalized in an analysis-local copy.",
      "One or more per-locus diagnostic estimates could not be computed.",
      "A null-model heritability estimate lies on the boundary of [0, 1]; at hsq = 0 the kinship term has no effect for that trait and chromosome.",
      "No kinship was supplied, so founder effects were fitted by ordinary least squares without a polygenic random effect.",
      "Exposure and outcome were fitted on the same batch sample set, so sample overlap is complete by construction."),
    recommended_action = c(
      "Inspect data/sample_alignment.csv (single analysis) or batch_sample_exclusions.csv (batch); exclusion causes can overlap.",
      "Inspect data/probability_normalization.csv and input precision; use probability_policy='strict' to reject rounding adjustments.",
      "Inspect diagnostic_status and diagnostic_note in per_locus_results.csv; unavailable diagnostics are not zero estimates.",
      "Inspect data/hsq.csv. A boundary estimate is not an error; check the kinship input and sample size, and supply hsq only with external justification.",
      "For related multiparent samples, supply a chromosome-named LOCO kinship list such as qtl2::calc_kinship(probs, type = 'loco').",
      "Review the batch sample policy and batch_sample_exclusions.csv."),
    stringsAsFactors = FALSE
  ))
}

qtl2mr_enrich_warnings <- function(warnings) {
  if (is.null(warnings) || !nrow(warnings)) {
    return(data.frame(
      warning_id = character(), severity = character(), message = character(),
      consequence = character(), recommended_action = character(),
      shown_in_main_report = logical()
    ))
  }
  cat <- qtl2mr_warning_catalog()
  if (!"severity" %in% names(warnings)) warnings$severity <- NA_character_
  if (!"message" %in% names(warnings)) warnings$message <- warnings$warning_id
  warnings$severity[warnings$severity == "method-limitation"] <- "method_limit"
  idx <- match(warnings$warning_id, cat$warning_id)
  missing_sev <- is.na(warnings$severity) | !nzchar(warnings$severity)
  warnings$severity[missing_sev] <- cat$default_severity[idx[missing_sev]]
  warnings$severity[is.na(warnings$severity)] <- "warning"
  warnings$consequence <- cat$consequence[idx]
  warnings$recommended_action <- cat$recommended_action[idx]
  warnings$shown_in_main_report <- cat$shown_in_main_report[idx]
  warnings$consequence[is.na(warnings$consequence)] <- "No catalog entry available."
  warnings$recommended_action[is.na(warnings$recommended_action)] <- "Review the relevant output and logs."
  warnings$shown_in_main_report[is.na(warnings$shown_in_main_report)] <-
    warnings$severity[is.na(warnings$shown_in_main_report)] %in% c("warning", "error", "method_limit")
  warnings[, c("warning_id", "severity", "message", "consequence", "recommended_action", "shown_in_main_report"), drop = FALSE]
}

qtl2mr_wald_ci <- function(beta, se, level = 0.95) {
  z <- stats::qnorm((1 + level) / 2)
  c(lower = beta - z * se, upper = beta + z * se)
}

qtl2mr_method_comparison_table <- function(result) {
  s <- result$summary
  fit <- result$fit
  diag_est <- result$diagnostic_layer$estimates %||% data.frame()
  row_value <- function(row, col, default) {
    if (col %in% names(row)) row[[col]][[1L]] else default
  }
  rows <- list()
  add <- function(method, estimate, se, ci_lower, ci_upper, ci_type, source, note,
                  ci_level = 0.95, bootstrap_B = 0L, bootstrap_seed = NA_integer_,
                  bootstrap_complete_B = 0L, bootstrap_failed_B = 0L) {
    rows[[length(rows) + 1L]] <<- data.frame(
      method = method,
      estimate = estimate,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_type = ci_type,
      source = source,
      note = note,
      ci_level = ci_level,
      bootstrap_B = bootstrap_B,
      bootstrap_seed = bootstrap_seed,
      bootstrap_complete_B = bootstrap_complete_B,
      bootstrap_failed_B = bootstrap_failed_B,
      stringsAsFactors = FALSE
    )
  }
  add(
    "Star-space qtl2mr",
    s$theta[[1L]], s$se[[1L]], s$profile_ci_lower[[1L]], s$profile_ci_upper[[1L]],
    "profile", "current_package", "Overlap-aware SO Founder GLS-IVW; reported when optimizer/profile diagnostics are reliable"
  )
  if (!is.null(fit$independence)) {
    ci <- qtl2mr_wald_ci(fit$independence$theta, fit$independence$se)
    add("Sigma_xy=0 qtl2mr", fit$independence$theta, fit$independence$se, ci[["lower"]], ci[["upper"]],
        "Wald", "current_package", "No overlap term")
  }
  for (label in c("Two-step IVW", "Weighted median", "Weighted mode")) {
    row <- if (nrow(diag_est) && "estimator" %in% names(diag_est)) {
      diag_est[diag_est$estimator == label, , drop = FALSE][1L, , drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(row) && is.finite(row$theta[[1L]])) {
      add(
        label,
        row$theta[[1L]],
        row_value(row, "se", NA_real_),
        row_value(row, "ci_lower", NA_real_),
        row_value(row, "ci_upper", NA_real_),
        row_value(row, "ci_type", NA_character_),
        "current_package_diagnostic",
        row_value(row, "note", "Multi-locus diagnostic comparator"),
        ci_level = row_value(row, "ci_level", NA_real_),
        bootstrap_B = row_value(row, "bootstrap_B", 0L),
        bootstrap_seed = row_value(row, "bootstrap_seed", NA_integer_),
        bootstrap_complete_B = row_value(row, "bootstrap_complete_B", NA_integer_),
        bootstrap_failed_B = row_value(row, "bootstrap_failed_B", NA_integer_)
      )
    } else {
      nm <- switch(label,
                   "Two-step IVW" = "theta_two_step_ivw",
                   "Weighted median" = "theta_weighted_median",
                   "Weighted mode" = "theta_weighted_mode")
      if (nm %in% names(s) && is.finite(s[[nm]][[1L]])) {
        se_field <- switch(label,
                           "Two-step IVW" = "se_two_step_ivw",
                           "Weighted median" = "se_weighted_median",
                           "Weighted mode" = "se_weighted_mode")
        lower_field <- switch(label,
                              "Two-step IVW" = NA_character_,
                              "Weighted median" = "ci_lower_weighted_median",
                              "Weighted mode" = "ci_lower_weighted_mode")
        upper_field <- switch(label,
                              "Two-step IVW" = NA_character_,
                              "Weighted median" = "ci_upper_weighted_median",
                              "Weighted mode" = "ci_upper_weighted_mode")
        ci_type_field <- switch(label,
                                "Two-step IVW" = NA_character_,
                                "Weighted median" = "ci_type_weighted_median",
                                "Weighted mode" = "ci_type_weighted_mode")
        se <- if (!is.na(se_field) && se_field %in% names(s)) s[[se_field]][[1L]] else NA_real_
        lower <- if (!is.na(lower_field) && lower_field %in% names(s)) s[[lower_field]][[1L]] else NA_real_
        upper <- if (!is.na(upper_field) && upper_field %in% names(s)) s[[upper_field]][[1L]] else NA_real_
        ci_type <- if (!is.na(ci_type_field) && ci_type_field %in% names(s)) s[[ci_type_field]][[1L]] else NA_character_
        if (!is.finite(lower) || !is.finite(upper)) {
          ci <- if (is.finite(se)) qtl2mr_wald_ci(s[[nm]][[1L]], se) else c(lower = NA_real_, upper = NA_real_)
          lower <- ci[["lower"]]
          upper <- ci[["upper"]]
          if (is.finite(se) && !nzchar(ci_type %||% "")) ci_type <- "Wald"
        }
        add(label, s[[nm]][[1L]], se, lower, upper, ci_type, "current_package_diagnostic", "Multi-locus diagnostic comparator")
      }
    }
  }
  do.call(rbind, rows)
}

qtl2mr_diagnostic_ci_caption <- function(methods, detailed = FALSE) {
  needed <- c("method", "ci_type", "ci_lower", "ci_upper")
  if (!all(needed %in% names(methods))) return("")
  rows <- methods[methods$method %in% c("Weighted median", "Weighted mode") &
    methods$ci_type %in% "bootstrap-SE Wald" &
    is.finite(methods$ci_lower) & is.finite(methods$ci_upper), , drop = FALSE]
  if (!nrow(rows)) return("")
  label <- if (length(unique(rows$method)) == 2L) "Weighted median/mode" else rows$method[[1L]]
  levels <- if ("ci_level" %in% names(rows)) unique(rows$ci_level) else NA_real_
  known <- length(levels) == 1L && is.finite(levels) && levels > 0 && levels < 1
  interval <- if (known) sprintf("approximate %g%% CIs", levels * 100) else "normal-approximation CIs"
  z <- if (known) format(round(stats::qnorm((1 + levels) / 2), 3), trim = TRUE) else "z"
  note <- sprintf("%s: %s = estimate +/- %s x parametric-bootstrap SD.", label, interval, z)
  if (isTRUE(detailed)) note <- paste(note,
    "Loci and inverse-variance weights are fixed; per-locus effects are sampled independently using their estimated effects and SEs. These are diagnostic intervals, not profile-likelihood intervals or guaranteed coverage.")
  note
}

plot_method_comparison <- function(methods,
                                   main = "qtl2mr method comparison",
                                   xlab = "MR effect",
                                   show_value_labels = TRUE,
                                   primary_method = NULL,
                                   ...) {
  required <- c("method", "estimate", "ci_lower", "ci_upper")
  miss <- setdiff(required, names(methods))
  if (length(miss)) stop("methods is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- methods[is.finite(methods$estimate), , drop = FALSE]
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  y <- seq_len(nrow(d))
  xlim <- range(c(0, d$ci_lower, d$ci_upper, d$estimate), finite = TRUE)
  pad <- diff(xlim) * 0.12
  if (!is.finite(pad) || pad <= 0) pad <- 1
  xlim <- xlim + c(-pad, pad)
  cols <- ifelse(d$source == "current_package", "#00796b",
                 ifelse(grepl("^external", d$source), "grey45", "#4c63d2"))
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  ci_caption <- qtl2mr_diagnostic_ci_caption(methods)
  graphics::par(mar = c(if (nzchar(ci_caption)) 6.5 else 5, 10, 5, 4))
  graphics::plot(d$estimate, y, xlim = xlim, ylim = c(0.5, length(y) + 0.5),
                 yaxt = "n", xlab = xlab, ylab = "", pch = 16, col = cols, main = main, ...)
  graphics::abline(v = 0, lty = 3, col = "grey70")
  graphics::segments(d$ci_lower, y, d$ci_upper, y, col = cols, lwd = 2)
  graphics::points(d$estimate, y, pch = 16, col = cols)
  graphics::axis(2, at = y, labels = d$method, las = 1, tick = FALSE, cex.axis = 0.8)
  if (nzchar(ci_caption)) {
    graphics::mtext(ci_caption, side = 1, line = 4.5, cex = 0.64, col = "grey30")
  }
  if (isTRUE(show_value_labels)) {
    label <- ifelse(
      is.finite(d$ci_lower) & is.finite(d$ci_upper),
      sprintf(" %.3f [%.3f, %.3f]", d$estimate, d$ci_lower, d$ci_upper),
      sprintf(" %.3f", d$estimate)
    )
    graphics::text(d$estimate, y, labels = label, pos = 4, cex = 0.72, col = cols, xpd = TRUE)
  }
  if (is.null(primary_method)) {
    primary <- methods[methods$method == "Star-space qtl2mr", , drop = FALSE]
    if (nrow(primary)) {
      primary_method <- list(
        estimator = "Star-space qtl2mr",
        theta = primary$estimate[[1L]],
        ci_lower = primary$ci_lower[[1L]],
        ci_upper = primary$ci_upper[[1L]],
        ci_type = primary$ci_type[[1L]]
      )
    }
  }
  if (!is.null(primary_method) && is.finite(primary_method$theta)) {
    ci_type <- as.character(primary_method$ci_type %||% "")
    if (is.na(ci_type)) ci_type <- ""
    graphics::mtext(
      sprintf(
        "Reported estimate (%s): %.4f, %s CI [%.4f, %.4f]",
        primary_method$estimator %||% "estimate",
        primary_method$theta,
        ci_type,
        primary_method$ci_lower,
        primary_method$ci_upper
      ),
      side = 3, line = 0.2, cex = 0.82, col = "#00695c"
    )
  }
  invisible(d)
}

qtl2mr_warning_table <- function(result, sample_policy = NULL) {
  s <- result$summary
  rows <- list()
  add <- function(id, severity, message) {
    rows[[length(rows) + 1L]] <<- data.frame(
      warning_id = id,
      severity = severity,
      message = message,
      stringsAsFactors = FALSE
    )
  }
  if (is.null(sample_policy) || identical(sample_policy, "role_specific_complete")) {
    add("prepared_overlap_engine", "info", "Exposure and outcome use their own complete samples; overlap is derived from harmonized sample IDs.")
  } else {
    add("batch_shared_samples", "info", sprintf("Exposure and outcome were fitted on the same %s sample set.", sample_policy))
  }
  alignment <- result$sample_alignment %||% data.frame()
  covar_counts <- c(alignment$n_incomplete_covar %||% integer(), alignment$n_missing_covar_id %||% integer())
  if (nrow(alignment) && any(covar_counts > 0L, na.rm = TRUE)) {
    add("incomplete_covariates", "warning",
      "Samples without a covariate row or with non-finite covariate values were excluded. Per-locus, per-role counts and sample IDs are recorded in data/sample_alignment.csv; counts across loci are not additive.")
  }
  hsq <- result[["hsq"]]
  if (is.data.frame(hsq) && nrow(hsq) && all(c("hsq_status", "hsq_source") %in% names(hsq))) {
    boundary <- hsq[hsq$hsq_status %in% c("boundary_zero", "boundary_one"), , drop = FALSE]
    if (nrow(boundary)) add("hsq_boundary", "info", sprintf(
      "Null-model hsq is on the boundary for %s.",
      paste(sprintf("%s chr %s (hsq = %g)", boundary$role, boundary$chr, boundary$hsq), collapse = "; ")))
    if (any(hsq$hsq_source == "no_kinship_ols")) add("no_kinship_ols", "info",
      "No kinship was supplied; founder effects were fitted by ordinary least squares.")
  }
  normalization <- result$prepared_metadata$probability_normalization$adjustments %||% data.frame()
  if (is.data.frame(normalization) && nrow(normalization)) {
    add("probability_rounding_normalized", "warning",
      "Founder probabilities with small rounding errors were normalized for analysis; original input arrays were not changed. See data/probability_normalization.csv.")
  }
  if ("n_loci_raw_space_psd_invalid" %in% names(s) && s$n_loci_raw_space_psd_invalid[[1L]] > 0) {
    add("raw_space_psd_invalid", "technical_qc", sprintf("Technical QC: raw-space diagnostic covariance invalid for %s loci; raw-space is not used by the primary estimator.", s$n_loci_raw_space_psd_invalid[[1L]]))
  }
  if ("n_loci_star_space_psd_invalid" %in% names(s) && s$n_loci_star_space_psd_invalid[[1L]] > 0) {
    add("star_space_psd_invalid", "warning", sprintf("Star-space PSD invalid for %s loci.", s$n_loci_star_space_psd_invalid[[1L]]))
  }
  if ("primary_reliable" %in% names(s) && !isTRUE(s$primary_reliable[[1L]])) {
    status <- s$primary_optimizer_status[[1L]] %||% "unknown"
    bracket <- s$primary_bracket_status[[1L]] %||% "unknown"
    add(
      "primary_so_unbounded",
      "warning",
      sprintf(
        "Primary second-order estimate is not well bracketed or has an unbounded profile CI (optimizer=%s; bracket=%s). Treat the primary theta as unstable and review diagnostic comparators.",
        status,
        bracket
      )
    )
  }
  dropped <- result$dropped_loci %||% data.frame()
  if (nrow(dropped)) {
    n_dropped <- length(unique(dropped$marker))
    n_retained <- if ("n_loci" %in% names(s)) s$n_loci[[1L]] else NA_integer_
    add(
      "dropped_nonestimable_loci",
      "warning",
      sprintf(
        "Dropped %d selected loci with non-estimable founder contrasts before MR; retained %s loci.",
        n_dropped,
        n_retained
      )
    )
  }
  per_locus <- result$diagnostic_layer$per_locus
  if (!is.null(per_locus) && "diagnostic_status" %in% names(per_locus) &&
      any(!per_locus$diagnostic_status %in% c("ok", "available"), na.rm = TRUE)) {
    add("diagnostic_unavailable", "warning",
      "Some per-locus diagnostics are unavailable; inspect their status and note. Primary inference is not replaced by missing diagnostic estimates.")
  }
  if (!is.null(per_locus) && "F_X_block" %in% names(per_locus)) {
    n_weak <- sum(per_locus$F_X_block < 10, na.rm = TRUE)
    if (n_weak) add("weak_block_first_stage", "warning", sprintf("%d loci have block F_X < 10 heuristic threshold.", n_weak))
  }
  if (!length(rows)) return(qtl2mr_enrich_warnings(NULL))
  qtl2mr_enrich_warnings(do.call(rbind, rows))
}

qtl2mr_block_first_stage_table <- function(result) {
  per_locus <- result$diagnostic_layer$per_locus
  if (is.null(per_locus) || !"F_X_block" %in% names(per_locus)) {
    return(data.frame(marker = result$selected_loci$marker, chr = result$selected_loci$chr, F_X_block = NA_real_))
  }
  per_locus[, intersect(c("marker", "chr", "F_X_block", "theta", "se", "Q_min", "weight"), names(per_locus)), drop = FALSE]
}

validate_qtl2_scan1coef <- function(genoprobs,
                                    pheno,
                                    kinship = NULL,
                                    addcovar = NULL,
                                    contrasts = NULL,
                                    hsq = NULL,
                                    marker = NULL,
                                    sample_ids = NULL,
                                    standardize_pheno = FALSE,
                                    reml = TRUE,
                                    tolerance = 1e-6) {
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one finite non-negative numeric value", call. = FALSE)
  }
  qtl2mr_require_qtl2()
  comparison <- compare_founder_lmm_to_scan1coef(
    genoprobs = genoprobs,
    pheno = pheno,
    kinship = kinship,
    addcovar = addcovar,
    contrasts = contrasts,
    hsq = hsq,
    marker = marker,
    sample_ids = sample_ids,
    standardize_pheno = standardize_pheno,
    reml = reml
  )
  pass <- comparison$max_abs_beta_diff <= tolerance &&
    comparison$max_abs_se_diff <= tolerance
  data.frame(
    validation_status = if (isTRUE(pass)) "pass" else "fail",
    validation_available = TRUE,
    marker = comparison$direct$marker,
    lmm_backend = comparison$direct$lmm_backend,
    hsq = comparison$direct$hsq,
    max_beta_abs_diff = comparison$max_abs_beta_diff,
    max_se_abs_diff = comparison$max_abs_se_diff,
    tolerance = tolerance,
    qtl2_version = as.character(utils::packageVersion("qtl2")),
    stringsAsFactors = FALSE
  )
}

qtl2mr_write_figure_pair <- function(stem, plot_expr,
                                     width = 1800, height = 1300,
                                     pdf_width = 8.5, pdf_height = 6.5) {
  plot_expr <- substitute(plot_expr)
  qtl2mr_draw_explicit_device(paste0(stem, ".png"), eval(plot_expr, parent.frame()), width = width, height = height, pdf_width = pdf_width, pdf_height = pdf_height)
  qtl2mr_draw_explicit_device(paste0(stem, ".pdf"), eval(plot_expr, parent.frame()), pdf_width = pdf_width, pdf_height = pdf_height)
  invisible(c(paste0(stem, ".png"), paste0(stem, ".pdf")))
}

write_summary_md <- function(result, file) {
  lines <- summarise_qtl2mr_result(result)
  writeLines(c("# qtl2mr analysis summary", "", lines), file)
  invisible(file)
}

summarise_qtl2mr_result <- function(result) {
  s <- qtl2mr_add_q_fields(result$pilot$summary, fit = result$pilot$fit)
  warnings <- result$warnings
  main_warnings <- if (!is.null(warnings) && nrow(warnings) && "shown_in_main_report" %in% names(warnings)) {
    warnings[warnings$shown_in_main_report %in% TRUE, , drop = FALSE]
  } else warnings
  warning_ids <- if (!is.null(main_warnings) && nrow(main_warnings)) {
    paste(main_warnings$warning_id, collapse = ", ")
  } else "none"
  reported <- qtl2mr_summary_reported(s)
  so_theta <- qtl2mr_summary_scalar(s, "theta_so", qtl2mr_summary_scalar(s, "theta"))
  reliable <- qtl2mr_summary_scalar(s, "primary_reliable", NA)
  ci_type <- as.character(reported$ci_type %||% "")
  if (is.na(ci_type)) ci_type <- ""
  raw_differs <- !isTRUE(reliable) || !is.finite(reported$theta) || !is.finite(so_theta) ||
    !isTRUE(all.equal(as.numeric(reported$theta), as.numeric(so_theta)))
  lines <- c(
    "qtl2mr analysis complete", "",
    paste0("Exposure: ", s$exposure[[1L]]),
    paste0("Outcome: ", s$outcome[[1L]]),
    paste0("N samples: ", length(result$pilot$sample_ids)),
    paste0("N loci: ", s$n_loci[[1L]]),
    paste0("Reported estimator: ", reported$estimator),
    paste0("theta_reported: ", qtl2mr_fmt_num(reported$theta, 6)),
    paste0("reported 95% CI: [", qtl2mr_fmt_num(reported$ci_lower, 6), ", ", qtl2mr_fmt_num(reported$ci_upper, 6), "]", if (nzchar(ci_type)) paste0(" (", ci_type, ")") else ""),
    paste0("Numerically bounded SO fit: ", as.character(reliable), " (not a coverage or causal-validity guarantee)"),
    paste0("Effect units: ", qtl2mr_summary_scalar(s, "effect_units", "not recorded")),
    paste0("Sample policy: ", qtl2mr_summary_scalar(s, "sample_policy", "not recorded")),
    paste0("Null-model hsq: ", qtl2mr_hsq_line(result$pilot[["hsq"]])),
    paste0("Q statistic: ", qtl2mr_fmt_num(qtl2mr_summary_scalar(s, "Q"), 6),
      "; Q p-value: ", qtl2mr_q_reported_label(s, digits = 6),
      "; Q reference: ", qtl2mr_q_reference_label(s)),
    paste0("Actionable warning IDs: ", warning_ids),
    paste0("Outputs: ", if (isTRUE(result$metadata$make_html %||% TRUE)) "report.html; " else "",
      "summary.md; data/main_result.csv; data/warnings.csv; data/output_manifest.csv")
  )
  if (isTRUE(raw_differs)) {
    lines <- append(lines, paste0("Raw star-space diagnostic theta: ", qtl2mr_fmt_num(so_theta, 6), " (diagnostic only)"), after = 9L)
  }
  lines
}

qtl2mr_hsq_line <- function(hsq) {
  if (!is.data.frame(hsq) || !nrow(hsq) || !"hsq_source" %in% names(hsq)) return("not recorded")
  if (all(hsq$hsq_source == "no_kinship_ols")) return("not used (no kinship; ordinary least squares)")
  parts <- vapply(c("exposure", "outcome"), function(role) {
    d <- hsq[hsq$role == role, , drop = FALSE]
    if (!nrow(d)) return(NA_character_)
    r <- range(d$hsq, na.rm = TRUE)
    value <- if (isTRUE(all.equal(r[[1L]], r[[2L]]))) qtl2mr_fmt_num(r[[1L]], 3) else
      paste0(qtl2mr_fmt_num(r[[1L]], 3), "-", qtl2mr_fmt_num(r[[2L]], 3))
    n_boundary <- sum(d$hsq_status %in% c("boundary_zero", "boundary_one"))
    paste0(role, " ", value, if (n_boundary) paste0(" (", n_boundary, " on boundary)") else "")
  }, character(1))
  sources <- unique(hsq$hsq_source)
  paste0(paste(parts[!is.na(parts)], collapse = "; "), " [", paste(sources, collapse = ", "), "; data/hsq.csv]")
}

print.qtl2mr_analysis_result <- function(x, ...) {
  cat(paste(summarise_qtl2mr_result(x), collapse = "\n"), "\n")
  invisible(x)
}

# A one-row overview of the result; qtl2mr_primary_result() has every column.
summary.qtl2mr_analysis_result <- function(object, ...) {
  s <- qtl2mr_add_q_fields(object$pilot$summary, fit = object$pilot$fit)
  get <- function(name, default = NA) qtl2mr_summary_scalar(s, name, default)
  warnings <- object$warnings %||% data.frame()
  shown <- if (nrow(warnings) && "shown_in_main_report" %in% names(warnings)) {
    warnings$warning_id[warnings$shown_in_main_report %in% TRUE]
  } else character()
  data.frame(
    exposure = get("exposure", NA_character_), outcome = get("outcome", NA_character_),
    n_samples = length(object$pilot$sample_ids), n_loci = get("n_loci", NA_integer_),
    estimate = get("theta_reported", NA_real_), se = get("se_reported", NA_real_),
    ci_lower = get("reported_ci_lower", NA_real_), ci_upper = get("reported_ci_upper", NA_real_),
    effect_units = get("effect_units", NA_character_),
    reported_estimator = get("reported_estimator", NA_character_),
    primary_reliable = get("primary_reliable", NA),
    confidence_set_type = get("confidence_set_type", NA_character_),
    confidence_set_bounded = get("confidence_set_bounded", qtl2mr_confidence_set_bounded(object$pilot$fit)),
    Q = get("Q", NA_real_), Q_p = get("Q_p_reported", NA_real_), Q_p_method = get("Q_p_method", NA_character_),
    sample_policy = get("sample_policy", NA_character_),
    warnings = paste(unique(shown), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

qtl2mr_primary_result <- function(result) {
  if (!inherits(result, "qtl2mr_analysis_result")) stop("result must be a qtl2mr_analysis_result", call. = FALSE)
  qtl2mr_add_q_fields(result$primary %||% result$pilot$summary, fit = result$pilot$fit)
}

qtl2mr_warnings <- function(result) {
  if (!inherits(result, "qtl2mr_analysis_result")) stop("result must be a qtl2mr_analysis_result", call. = FALSE)
  result$warnings %||% qtl2mr_enrich_warnings(NULL)
}

qtl2mr_output_paths <- function(result, detail = FALSE) {
  if (!inherits(result, "qtl2mr_analysis_result")) stop("result must be a qtl2mr_analysis_result", call. = FALSE)
  manifest <- result$paths %||% result$output_manifest %||% data.frame()
  if (isTRUE(detail)) return(manifest)
  output_dir <- normalizePath(result$output_dir %||% ".", mustWork = FALSE)
  key_paths <- data.frame(
    item = c(
      "output_dir",
      "report_html",
      "summary_md",
      "main_result_csv",
      "warnings_csv",
      "output_manifest_csv",
      "output_manifest_html",
      "result_object_rds"
    ),
    path = c(
      output_dir,
      "report.html",
      "summary.md",
      "data/main_result.csv",
      "data/warnings.csv",
      "data/output_manifest.csv",
      "output_manifest.html",
      "result_object.rds"
    ),
    stringsAsFactors = FALSE
  )
  if (identical(result$metadata$make_html, FALSE)) {
    key_paths <- key_paths[key_paths$item != "report_html", , drop = FALSE]
  }
  if (nrow(manifest) && all(c("path_relative", "size_kb") %in% names(manifest))) {
    key_paths$size_kb <- NA_real_
    idx <- match(key_paths$path, manifest$path_relative)
    key_paths$size_kb[!is.na(idx)] <- manifest$size_kb[idx[!is.na(idx)]]
  }
  attr(key_paths, "note") <- "Use qtl2mr_output_paths(result, detail = TRUE) or data/output_manifest.csv for the full file manifest."
  key_paths
}

qtl2mr_method_comparison <- function(result) {
  if (!inherits(result, "qtl2mr_analysis_result")) stop("result must be a qtl2mr_analysis_result", call. = FALSE)
  result$tables$method_comparison %||% result$method_comparison %||% data.frame()
}

qtl2mr_add_stable_slots <- function(result) {
  s <- result$pilot$summary
  result$metadata <- list(
    package = "qtl2mr",
    version = tryCatch(as.character(utils::packageVersion("qtl2mr")), error = function(e) NA_character_),
    created_at = result$created_at %||% Sys.time(),
    input = result$input %||% NA_character_,
    exposure = result$exposure %||% s$exposure[[1L]],
    outcome = result$outcome %||% s$outcome[[1L]],
    make_html = result$metadata$make_html %||% TRUE
  )
  result$primary <- s
  result$diagnostics <- list(
    loo_mr = result$loo_mr %||% data.frame(),
    instrument_sensitivity = result$instrument_sensitivity %||% data.frame(),
    locus_diagnostics = result$pilot$locus_diagnostics,
    diagnostic_layer = result$pilot$diagnostic_layer
  )
  result$tables <- list(
    method_comparison = result$method_comparison %||% data.frame(),
    warnings = result$warnings %||% qtl2mr_enrich_warnings(NULL)
  )
  result$figures <- result$figure_manifest %||% data.frame()
  result$paths <- result$output_manifest %||% data.frame()
  result$raw_objects <- list(pilot = result$pilot)
  result
}

write_qtl2mr_html_report <- function(result, file) {
  result$pilot$summary <- qtl2mr_add_q_fields(result$pilot$summary, fit = result$pilot$fit)
  s <- result$pilot$summary
  reported <- qtl2mr_summary_reported(s)
  warnings <- result$warnings %||% data.frame()
  actionable <- if (nrow(warnings) && "shown_in_main_report" %in% names(warnings)) {
    warnings[warnings$shown_in_main_report %in% TRUE & warnings$severity %in% c("warning", "error"),
             intersect(c("warning_id", "severity", "message", "recommended_action"), names(warnings)), drop = FALSE]
  } else data.frame()
  methods <- result$method_comparison %||% data.frame()
  method_table <- data.frame()
  if (nrow(methods)) {
    method_table <- methods[, intersect(c("method", "estimate", "ci_lower", "ci_upper", "ci_type", "ci_level", "source"), names(methods)), drop = FALSE]
    level_label <- if ("ci_level" %in% names(method_table)) {
      ifelse(is.finite(method_table$ci_level), paste0(100 * method_table$ci_level, "% "), "")
    } else ""
    method_table$CI <- if (all(c("ci_lower", "ci_upper") %in% names(method_table))) {
      ifelse(
        is.finite(method_table$ci_lower) & is.finite(method_table$ci_upper),
        paste0(level_label, "[", qtl2mr_fmt_num(method_table$ci_lower, 5), ", ", qtl2mr_fmt_num(method_table$ci_upper, 5), "]"),
        "not available"
      )
    } else "not available"
    method_table$role <- "comparator"
    if ("source" %in% names(method_table)) {
      method_table$role[grepl("diagnostic", method_table$source, ignore.case = TRUE)] <- "diagnostic"
    }
    if (is.finite(reported$theta) && "method" %in% names(method_table)) {
      method_table$role[method_table$method == reported$estimator] <- "reported"
    }
    if ("ci_type" %in% names(method_table)) method_table$`CI method` <- method_table$ci_type
    method_table <- method_table[, intersect(c("method", "estimate", "CI", "CI method", "role"), names(method_table)), drop = FALSE]
  }
  loo <- result$loo_mr %||% data.frame()
  if (nrow(loo) && "status" %in% names(loo)) loo <- loo[loo$status == "ok", , drop = FALSE]
  loo <- if (nrow(loo)) utils::head(loo[order(loo$abs_delta_theta, decreasing = TRUE),
    intersect(c("marker", "theta_minus_m", "profile_ci_lower", "profile_ci_upper", "delta_theta"), names(loo)), drop = FALSE], 5L) else data.frame()
  sensitivity <- result$instrument_sensitivity %||% data.frame()
  if (nrow(sensitivity) && "status" %in% names(sensitivity)) sensitivity <- sensitivity[sensitivity$status == "ok", , drop = FALSE]
  sensitivity <- sensitivity[, intersect(c("set_name", "n_loci", "theta", "se", "profile_ci_lower", "profile_ci_upper"), names(sensitivity)), drop = FALSE]
  figs <- result$figure_manifest %||% data.frame()
  figure_html <- character()
  if (nrow(figs)) {
    pngs <- figs[figs$format == "png", , drop = FALSE]
    for (i in seq_len(nrow(pngs))) {
      path <- qtl2mr_relative_path(pngs$path[[i]], result$output_dir)
      figure_html <- c(figure_html, sprintf("<figure><img src=\"%s\" alt=\"%s\"><figcaption>%s</figcaption></figure>",
        qtl2mr_html_escape(path), qtl2mr_html_escape(pngs$label[[i]]), qtl2mr_html_escape(pngs$notes[[i]])))
    }
  }
  hsq_all <- result$pilot[["hsq"]]
  hsq_view <- if (is.data.frame(hsq_all) && nrow(hsq_all) && "hsq_source" %in% names(hsq_all)) do.call(rbind, lapply(split(hsq_all, paste(hsq_all$role, hsq_all$hsq_source)), function(d) data.frame(
    role = d$role[[1L]], trait = d$trait[[1L]], source = d$hsq_source[[1L]], n_chr = length(unique(d$chr)),
    hsq_min = min(d$hsq), hsq_median = stats::median(d$hsq), hsq_max = max(d$hsq),
    n_boundary = sum(d$hsq_status %in% c("boundary_zero", "boundary_one")), stringsAsFactors = FALSE))) else data.frame()
  ci_type <- as.character(reported$ci_type %||% "")
  if (is.na(ci_type)) ci_type <- ""
  reliable <- qtl2mr_summary_scalar(s, "primary_reliable", NA)
  so_theta <- qtl2mr_summary_scalar(s, "theta_so", qtl2mr_summary_scalar(s, "theta"))
  raw_note <- if (!isTRUE(reliable) || !is.finite(reported$theta) || !is.finite(so_theta) ||
                  !isTRUE(all.equal(as.numeric(reported$theta), as.numeric(so_theta)))) {
    paste0(" Raw star-space second-order results remain diagnostic only: theta = ", qtl2mr_fmt_num(so_theta, 6), ".")
  } else ""
  primary_table <- data.frame(
    item = c("Reported estimator", "Estimate", "95% CI", "Numerically bounded SO fit", "Effect units", "Confidence set", "Inference calibration", "Sample policy", "Q statistic", "Q p-value", "Q reference"),
    value = c(reported$estimator, qtl2mr_fmt_num(reported$theta, 6),
      paste0("[", qtl2mr_fmt_num(reported$ci_lower, 6), ", ", qtl2mr_fmt_num(reported$ci_upper, 6), "]", if (nzchar(ci_type)) paste0(" (", ci_type, ")") else ""),
      as.character(reliable), qtl2mr_summary_scalar(s, "effect_units", "not recorded"),
      qtl2mr_summary_scalar(s, "confidence_set_type", "not recorded"),
      qtl2mr_summary_scalar(s, "inference_status", "not established"),
      qtl2mr_summary_scalar(s, "sample_policy", "not recorded"),
      qtl2mr_fmt_num(qtl2mr_summary_scalar(s, "Q"), 6),
      qtl2mr_q_reported_label(s, digits = 5), qtl2mr_q_reference_label(s)),
    stringsAsFactors = FALSE
  )
  overview <- sprintf("<p><strong>Exposure:</strong> %s<br><strong>Outcome:</strong> %s<br><strong>Samples:</strong> %s &nbsp; <strong>Selected loci:</strong> %s</p>",
    qtl2mr_html_escape(s$exposure[[1L]]), qtl2mr_html_escape(s$outcome[[1L]]),
    length(result$pilot$sample_ids), s$n_loci[[1L]])
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\"><title>qtl2mr Analysis Report</title>",
    "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;color:#1d2528;line-height:1.5;}main{max-width:960px;margin:auto;}table{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0 22px;}th,td{border:1px solid #d9dfdd;padding:6px 8px;text-align:left;vertical-align:top;}th{background:#eef3f1;}figure{border:1px solid #ddd;padding:12px;margin:18px 0;}img{max-width:100%;height:auto;display:block;margin:auto;}code{background:#eef3f1;padding:1px 4px;}.primary{border:2px solid #00796b;background:#f2fbf8;padding:14px;margin:16px 0 22px;}.section-note{color:#445;max-width:900px;}</style>",
    "</head><body><main>",
    "<h1>qtl2mr analysis report</h1>",
    "<h2>Analysis overview and primary estimate</h2>", overview,
    "<div class=\"primary\">", qtl2mr_html_table(primary_table, digits = 6), "</div>",
    paste0("<p class=\"section-note\">Numerical success does not establish 95% coverage or causal validity. Weak instruments, same-sample instrument selection and dependence between loci can invalidate nominal uncertainty. Other methods are diagnostic comparators, not automatic replacements.", raw_note, "</p>"),
    "<h2>QC and actionable warnings</h2>",
    if (nrow(actionable)) c("<p class=\"section-note\">F&lt;10 is a heuristic for multi-founder block strength, not a calibrated Stock-Yogo cutoff.</p>", qtl2mr_html_table(actionable)) else "<p>No actionable warnings were recorded.</p>",
    if (nrow(hsq_view)) c("<h3>Polygenic background (null-model heritability)</h3>",
      "<p class=\"section-note\">hsq sets the kinship-based working covariance of each founder LMM. It is estimated by REML under a null model with the intercept and covariates (no locus), using the kinship for that chromosome, and is then treated as fixed. Per-chromosome values are in <code>data/hsq.csv</code>.</p>",
      qtl2mr_html_table(hsq_view, digits = 4)) else "",
    if (nrow(method_table) || nrow(loo) || nrow(sensitivity)) c(
      "<h2>Methods and available sensitivities</h2>",
      if (nrow(method_table)) qtl2mr_html_table(method_table, digits = 5) else "",
      if (nzchar(qtl2mr_diagnostic_ci_caption(methods))) paste0(
        "<p class=\"section-note\">", qtl2mr_html_escape(qtl2mr_diagnostic_ci_caption(methods, detailed = TRUE)), "</p>") else "",
      if (nrow(loo)) c("<h3>Leave-one-locus-out sensitivity</h3>", qtl2mr_html_table(loo, digits = 5)) else "",
      if (nrow(sensitivity)) c("<h3>Instrument-set sensitivity</h3>", qtl2mr_html_table(sensitivity, digits = 5)) else ""
    ) else "",
    if (length(figure_html)) c("<h2>Figures</h2>", figure_html) else "",
    "<h2>Outputs and provenance</h2>",
    "<p>Detailed machine-readable results are retained in <code>data/main_result.csv</code>, <code>data/method_comparison.csv</code>, <code>data/warnings.csv</code>, and <a href=\"data/output_manifest.csv\"><code>data/output_manifest.csv</code></a>. This report is a curated view; <a href=\"output_manifest.html\">output_manifest.html</a> lists all generated files.</p>",
    "<p>Session details: <code>data/session_info.txt</code>.</p>",
    "</main></body></html>"
  )
  writeLines(html[nzchar(html)], file)
  invisible(file)
}

write_qtl2mr_results <- function(result,
                                 output_dir = result$output_dir,
                                 overwrite = FALSE,
                                 make_html = TRUE,
                                 make_figures = TRUE,
                                 progress_callback = NULL,
                                 q_sim_n = NULL,
  q_sim_seed = 1L,
                                 .transaction = NULL) {
  transaction <- .transaction %||% qtl2mr_begin_output_transaction(output_dir, overwrite = overwrite)
  if (!identical(qtl2mr_assert_known_path_node(transaction$stage_dir, "transaction staging directory"), "directory")) {
    stop("Supplied qtl2mr transaction has no safe staging directory", call. = FALSE)
  }
  published <- FALSE
  on.exit({
    if (!published) qtl2mr_abort_output_transaction(transaction)
  }, add = TRUE)
  final_output_dir <- transaction$destination
  output_dir <- transaction$stage_dir
  result$output_dir <- final_output_dir
  result$metadata$make_html <- isTRUE(make_html)
  data_dir <- file.path(output_dir, "data")
  figure_dir <- file.path(output_dir, "figures")
  log_dir <- file.path(output_dir, "logs")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  qtl2mr_progress(progress_callback, 5L, "Write outputs", "prepare tables and output folders")

  fit_q_method <- result$pilot$fit$Q_p_method %||%
    if (length(result$pilot$locus_blocks %||% list()) == 1L) "chisq" else "simulation"
  if (is.null(result$pilot$fit$Q_p_method)) result$pilot$fit$Q_p_method <- fit_q_method
  if (is.null(result$pilot$fit$Q_chisq_df)) {
    bx_len <- sum(vapply(result$pilot$locus_blocks %||% list(), function(b) length(b$beta_x %||% numeric()), integer(1L)))
    result$pilot$fit$Q_chisq_df <- max(bx_len - 1L, 1L)
  }
  if (is.null(result$pilot$fit$Q_p_chisq) &&
      is.finite(result$pilot$fit$Q %||% NA_real_) &&
      is.finite(result$pilot$fit$Q_chisq_df %||% NA_integer_)) {
    result$pilot$fit$Q_p_chisq <- stats::pchisq(
      result$pilot$fit$Q,
      df = result$pilot$fit$Q_chisq_df,
      lower.tail = FALSE
    )
  }
  if ((is.null(result$pilot$fit$Q_p_reported) || !is.finite(result$pilot$fit$Q_p_reported)) &&
      identical(fit_q_method, "chisq")) {
    result$pilot$fit$Q_p_reported <- result$pilot$fit$Q_p_chisq %||% NA_real_
    result$pilot$fit$Q_p_note <- paste0("chi-square reference with df=", result$pilot$fit$Q_chisq_df %||% NA_integer_)
  }
  if (!is.null(q_sim_n) && identical(fit_q_method, "simulation") && isTRUE(as.integer(q_sim_n) > 0L) &&
      (is.null(result$pilot$fit$Q_p_sim) || !is.finite(result$pilot$fit$Q_p_sim))) {
    qtl2mr_progress(progress_callback, 5L, "Write outputs", sprintf("simulate heterogeneity Q calibration: nsim=%d", as.integer(q_sim_n)))
    q_cal <- calibrate_multilocus_q(
      locus_blocks = result$pilot$locus_blocks,
      fit = result$pilot$fit,
      nsim = as.integer(q_sim_n),
      seed = q_sim_seed,
      allow_cross_cov = TRUE,
      keep_simulated = FALSE,
      progress_callback = progress_callback
    )
    for (nm in setdiff(names(q_cal), "Q_sim_values")) {
      result$pilot$fit[[nm]] <- q_cal[[nm]]
    }
    result$pilot$fit$Q_p_reported <- result$pilot$fit$Q_p_sim
    result$pilot$fit$Q_p_note <- "conditional residual simulation under the fitted common-effect null"
  }
  if ((is.null(result$pilot$fit$Q_p_reported) || !is.finite(result$pilot$fit$Q_p_reported)) &&
      identical(fit_q_method, "simulation") &&
      is.finite(result$pilot$fit$Q_p_sim %||% NA_real_)) {
    result$pilot$fit$Q_p_reported <- result$pilot$fit$Q_p_sim
    result$pilot$fit$Q_p_note <- "conditional residual simulation under the fitted common-effect null"
  }
  if ((is.null(result$pilot$fit$Q_p_reported) || !is.finite(result$pilot$fit$Q_p_reported)) &&
      identical(fit_q_method, "none")) {
    result$pilot$fit$Q_p_reported <- NA_real_
    result$pilot$fit$Q_p_note <- "heterogeneity p-value not reported"
  }
  result$pilot$summary <- qtl2mr_add_q_fields(result$pilot$summary, fit = result$pilot$fit)
  method_comparison <- qtl2mr_method_comparison_table(result$pilot)
  block_first_stage <- qtl2mr_block_first_stage_table(result$pilot)
  psd_cols <- intersect(c("marker", "chr", "joint_psd_valid", "joint_psd_status", "joint_min_eigenvalue",
                          "joint_max_cross_singular_value", "star_space_psd_valid", "raw_space_psd_valid",
                          "raw_space_max_cross_singular_value", "star_space_max_cross_singular_value"),
                        names(result$pilot$locus_diagnostics))
  overlap_cols <- intersect(c("marker", "chr", "n_x", "n_y", "n_overlap",
                              "exposure_overlap_fraction", "outcome_overlap_fraction",
                              "jaccard_overlap_fraction", "complete_overlap",
                              "independent_samples", "rho_xy", "rho_xy_star",
                              "rho_xy_raw", "cross_cov_model"),
                            names(result$pilot$locus_diagnostics))
  warnings <- qtl2mr_warning_table(result$pilot)
  append_warning <- function(id, severity, message) {
    qtl2mr_enrich_warnings(rbind(
      warnings[, c("warning_id", "severity", "message"), drop = FALSE],
      data.frame(warning_id = id, severity = severity, message = message, stringsAsFactors = FALSE)
    ))
  }
  if (!is.null(result$loo_error)) {
    warnings <- append_warning("loo_mr_error", "warning", paste("LOO MR diagnostic failed:", result$loo_error))
  }
  if (!is.null(result$instrument_sensitivity_error)) {
    warnings <- append_warning("instrument_sensitivity_error", "warning", paste("Instrument-set sensitivity failed:", result$instrument_sensitivity_error))
  }
  result$warnings <- warnings
  result$method_comparison <- method_comparison
  result$first_stage_summary <- tryCatch(summarize_block_first_stage(block_first_stage), error = function(e) data.frame())

  data_paths <- c(
    main_result = file.path(data_dir, "main_result.csv"),
    method_comparison = file.path(data_dir, "method_comparison.csv"),
    per_locus_results = file.path(data_dir, "per_locus_results.csv"),
    block_first_stage = file.path(data_dir, "block_first_stage.csv"),
    psd_qc = file.path(data_dir, "psd_qc.csv"),
    overlap_qc = file.path(data_dir, "overlap_qc.csv"),
    loo_mr = file.path(data_dir, "loo_mr.csv"),
    instrument_sensitivity = file.path(data_dir, "instrument_sensitivity.csv"),
    q_calibration = file.path(data_dir, "q_calibration.csv"),
    first_stage_summary = file.path(data_dir, "first_stage_summary.csv"),
    founder_residuals = file.path(data_dir, "founder_residuals.csv"),
    founder_coordinate_ratios = file.path(data_dir, "founder_coordinate_ratios.csv"),
    warnings = file.path(data_dir, "warnings.csv")
  )
  data_paths <- c(data_paths,
    confidence_set = file.path(data_dir, "confidence_set.csv"),
    input_alignment = file.path(data_dir, "input_alignment.csv"),
    sample_alignment = file.path(data_dir, "sample_alignment.csv"),
    phenotype_scaling = file.path(data_dir, "phenotype_scaling.csv"),
    probability_normalization = file.path(data_dir, "probability_normalization.csv"),
    hsq = file.path(data_dir, "hsq.csv"))
  confidence_set <- result$pilot$fit$confidence_set %||% data.frame(lower = numeric(), upper = numeric())
  confidence_set$confidence_set_type <- rep(result$pilot$fit$confidence_set_type %||% "not_recorded", nrow(confidence_set))
  confidence_set$numerical_ok <- rep(isTRUE(result$pilot$fit$numerical_ok), nrow(confidence_set))
  confidence_set$inference_status <- rep(result$pilot$fit$inference_status %||% "not_calibrated", nrow(confidence_set))
  utils::write.csv(confidence_set, data_paths[["confidence_set"]], row.names = FALSE)
  utils::write.csv(result$pilot$input_alignment %||% data.frame(), data_paths[["input_alignment"]], row.names = FALSE)
  utils::write.csv(result$pilot$sample_alignment %||% data.frame(), data_paths[["sample_alignment"]], row.names = FALSE)
  utils::write.csv(result$pilot$phenotype_scaling %||% data.frame(), data_paths[["phenotype_scaling"]], row.names = FALSE)
  utils::write.csv(if (is.data.frame(result$pilot[["hsq"]])) result$pilot[["hsq"]] else data.frame(),
    data_paths[["hsq"]], row.names = FALSE)
  utils::write.csv(result$pilot$prepared_metadata$probability_normalization$adjustments %||% data.frame(),
    data_paths[["probability_normalization"]], row.names = FALSE)
  utils::write.csv(result$pilot$summary, data_paths[["main_result"]], row.names = FALSE)
  utils::write.csv(method_comparison, data_paths[["method_comparison"]], row.names = FALSE)
  utils::write.csv(result$pilot$diagnostic_layer$per_locus %||% data.frame(), data_paths[["per_locus_results"]], row.names = FALSE)
  utils::write.csv(block_first_stage, data_paths[["block_first_stage"]], row.names = FALSE)
  utils::write.csv(result$pilot$locus_diagnostics[, psd_cols, drop = FALSE], data_paths[["psd_qc"]], row.names = FALSE)
  overlap_qc <- if (!is.null(result$pilot$overlap_summary) && nrow(result$pilot$overlap_summary)) {
    result$pilot$overlap_summary
  } else {
    result$pilot$locus_diagnostics[, overlap_cols, drop = FALSE]
  }
  utils::write.csv(overlap_qc, data_paths[["overlap_qc"]], row.names = FALSE)
  utils::write.csv(result$loo_mr %||% data.frame(), data_paths[["loo_mr"]], row.names = FALSE)
  utils::write.csv(result$instrument_sensitivity %||% data.frame(), data_paths[["instrument_sensitivity"]], row.names = FALSE)
  q_calibration <- result$pilot$summary[, intersect(c(
    "Q", "Q_p_reported", "Q_p_method", "Q_p_note",
    "Q_p_chisq", "Q_chisq_df",
    "Q_p_sim", "Q_sim_n", "Q_sim_n_requested", "Q_sim_seed", "Q_sim_exceedances",
    "Q_sim_min_resolvable_p", "Q_sim_null_mean", "Q_sim_null_sd", "Q_sim_method"
  ), names(result$pilot$summary)), drop = FALSE]
  utils::write.csv(q_calibration, data_paths[["q_calibration"]], row.names = FALSE)
  utils::write.csv(result$first_stage_summary %||% data.frame(), data_paths[["first_stage_summary"]], row.names = FALSE)
  utils::write.csv(result$pilot$founder_residual_scores_global %||% data.frame(), data_paths[["founder_residuals"]], row.names = FALSE)
  utils::write.csv(result$pilot$founder_coordinate_ratios %||% data.frame(), data_paths[["founder_coordinate_ratios"]], row.names = FALSE)
  utils::write.csv(warnings, data_paths[["warnings"]], row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(data_dir, "session_info.txt"))
  qtl2mr_progress(progress_callback, 5L, "Write outputs", "data tables and session info")

  figure_rows <- list()
  add_fig <- function(label, stem, notes) {
    png <- paste0(stem, ".png")
    pdf <- paste0(stem, ".pdf")
    if (file.exists(png)) figure_rows[[length(figure_rows) + 1L]] <<- data.frame(label = label, path = png, format = "png", notes = notes, stringsAsFactors = FALSE)
    if (file.exists(pdf)) figure_rows[[length(figure_rows) + 1L]] <<- data.frame(label = paste(label, "PDF"), path = pdf, format = "pdf", notes = "PDF companion", stringsAsFactors = FALSE)
  }
  if (isTRUE(make_figures)) {
    qtl2mr_progress(progress_callback, 5L, "Write outputs", "figures")
    reported_for_fig <- qtl2mr_summary_reported(result$pilot$summary)
    qtl2mr_write_figure_pair(file.path(figure_dir, "method_comparison"), {
      plot_method_comparison(method_comparison, main = "qtl2mr method comparison", xlab = paste(result$pilot$summary$exposure[[1L]], "->", result$pilot$summary$outcome[[1L]], "effect"), show_value_labels = FALSE, primary_method = reported_for_fig)
    }, width = 2400, height = 1400, pdf_width = 10.5, pdf_height = 6.8)
    add_fig("Method comparison", file.path(figure_dir, "method_comparison"),
      trimws(paste("Reported qtl2mr estimate and comparator estimates.", qtl2mr_diagnostic_ci_caption(method_comparison))))
    if (!is.null(result$pilot$diagnostic_layer$per_locus)) {
      qtl2mr_write_figure_pair(file.path(figure_dir, "per_locus_forest"), {
        plot_per_locus_so_forest(result$pilot$diagnostic_layer$per_locus, estimates = result$pilot$diagnostic_layer$estimates,
                                 global_theta = reported_for_fig$theta, max_loci = 40L,
                                 order_by = "F",
                                 main = "Per-locus SO estimates", xlab = paste(result$pilot$summary$exposure[[1L]], "->", result$pilot$summary$outcome[[1L]], "effect"))
      }, width = 2600, height = 1800, pdf_width = 12, pdf_height = 9)
      add_fig("Per-locus diagnostic forest", file.path(figure_dir, "per_locus_forest"), "Top per-locus SO estimates with reported-estimate reference.")
    }
    if (!is.null(result$loo_mr) && nrow(result$loo_mr) && any(result$loo_mr$status == "ok", na.rm = TRUE)) {
      qtl2mr_write_figure_pair(file.path(figure_dir, "loo_mr_top25"), {
        plot_loo_multilocus_mr(result$loo_mr, top_n = 25L, show = "theta",
                               main = "Leave-one-locus-out MR sensitivity",
                               xlab = paste(result$pilot$summary$exposure[[1L]], "->", result$pilot$summary$outcome[[1L]], "effect"))
      }, width = 2400, height = 1800, pdf_width = 11, pdf_height = 8.5)
      add_fig("LOO MR sensitivity", file.path(figure_dir, "loo_mr_top25"), "Top leave-one-locus-out global refits by absolute theta change.")
    }
    if (!is.null(result$instrument_sensitivity) && nrow(result$instrument_sensitivity) && any(result$instrument_sensitivity$status == "ok", na.rm = TRUE)) {
      qtl2mr_write_figure_pair(file.path(figure_dir, "instrument_set_sensitivity"), {
        plot_instrument_set_sensitivity(result$instrument_sensitivity,
                                        global_theta = reported_for_fig$theta,
                                        main = "Instrument-set sensitivity",
                                        xlab = paste(result$pilot$summary$exposure[[1L]], "->", result$pilot$summary$outcome[[1L]], "effect"))
      }, width = 2200, height = 1500, pdf_width = 10.5, pdf_height = 7)
      add_fig("Instrument-set sensitivity", file.path(figure_dir, "instrument_set_sensitivity"), "Global refits across selected-locus subsets with reported-estimate reference.")
    }
    if (!is.null(result$pilot$founder_residual_scores_global) && nrow(result$pilot$founder_residual_scores_global)) {
      top_markers <- utils::head(result$pilot$selected_loci$marker, 30L)
      heat_df <- result$pilot$founder_residual_scores_global[result$pilot$founder_residual_scores_global$marker %in% top_markers, , drop = FALSE]
      qtl2mr_write_figure_pair(file.path(figure_dir, "founder_residual_heatmap"), {
        plot_founder_residual_heatmap(heat_df, value = "t_one_founder_vs_rest", row = "marker",
                                      main = "Founder residual scores, top exposure loci", label_threshold = 2)
      }, width = 2400, height = 1800, pdf_width = 11, pdf_height = 8)
      add_fig("Founder residual heatmap", file.path(figure_dir, "founder_residual_heatmap"), "Founder residual localization across selected exposure-locus ranks; descriptive only.")
    }
    qtl2mr_write_figure_pair(file.path(figure_dir, "weak_denominator_dashboard"), {
      plot_block_first_stage_dashboard(block_first_stage)
    }, width = 2200, height = 1100, pdf_width = 11, pdf_height = 5.5)
    add_fig("Weak-denominator dashboard", file.path(figure_dir, "weak_denominator_dashboard"), "Histogram and ranked view of block F_X heuristic diagnostics.")
  }
  figure_manifest <- if (length(figure_rows)) do.call(rbind, figure_rows) else data.frame(label = character(), path = character(), format = character(), notes = character())
  if (nrow(figure_manifest)) {
    # Files are created in staging, but saved objects and HTML describe the
    # published directory. Staging paths disappear when the transaction commits.
    relative_figures <- vapply(figure_manifest$path, qtl2mr_relative_path,
      character(1), root = output_dir)
    figure_manifest$path <- file.path(final_output_dir, relative_figures)
  }
  result$figure_manifest <- figure_manifest

  model_path <- file.path(output_dir, "result_object.rds")
  summary_path <- file.path(output_dir, "summary.md")
  html_path <- file.path(output_dir, "report.html")
  manifest_html_path <- file.path(output_dir, "output_manifest.html")
  session_path <- file.path(data_dir, "session_info.txt")
  run_log_path <- file.path(log_dir, "run.log")
  manifest_path <- file.path(data_dir, "output_manifest.csv")
  qtl2mr_progress(progress_callback, 5L, "Write outputs", "summary, HTML report, RDS object, and manifest")
  summary_value <- function(name, default = NA_real_) {
    if (!name %in% names(result$pilot$summary)) return(default)
    result$pilot$summary[[name]][[1L]]
  }
  reported_log <- qtl2mr_summary_reported(result$pilot$summary)
  writeLines(c(
    paste0("created_at\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("input\t", result$input),
    paste0("output_dir\t", result$output_dir),
    paste0("n_loci\t", result$pilot$summary$n_loci[[1L]]),
    paste0("reported_estimator\t", reported_log$estimator),
    paste0("theta_reported\t", reported_log$theta),
    paste0("se_reported\t", reported_log$se),
    paste0("reported_ci_lower\t", reported_log$ci_lower),
    paste0("reported_ci_upper\t", reported_log$ci_upper),
    paste0("reported_ci_type\t", reported_log$ci_type),
    paste0("reported_estimator_reason\t", reported_log$reason),
    paste0("theta_so\t", summary_value("theta_so", summary_value("theta"))),
    paste0("se_so\t", summary_value("se_so", summary_value("se"))),
    paste0("primary_reliable\t", summary_value("primary_reliable", NA)),
    paste0("primary_optimizer_status\t", summary_value("primary_optimizer_status", NA_character_)),
    paste0("primary_bracket_status\t", summary_value("primary_bracket_status", NA_character_)),
    paste0("Qmin\t", summary_value("Q")),
    paste0("Q_p_reported\t", summary_value("Q_p_reported")),
    paste0("Q_p_method\t", result$pilot$summary$Q_p_method[[1L]] %||% NA_character_),
    paste0("Q_p_note\t", result$pilot$summary$Q_p_note[[1L]] %||% NA_character_),
    paste0("Q_p_chisq\t", summary_value("Q_p_chisq")),
    paste0("Q_chisq_df\t", summary_value("Q_chisq_df")),
    paste0("Q_p_sim\t", summary_value("Q_p_sim")),
    paste0("Q_sim_n\t", summary_value("Q_sim_n")),
    paste0("Q_sim_n_requested\t", summary_value("Q_sim_n_requested")),
    paste0("Q_sim_exceedances\t", summary_value("Q_sim_exceedances")),
    paste0("Q_sim_seed\t", summary_value("Q_sim_seed")),
    paste0("Q_sim_min_resolvable_p\t", summary_value("Q_sim_min_resolvable_p")),
    paste0("Q_sim_method\t", result$pilot$summary$Q_sim_method[[1L]] %||% NA_character_),
    paste0("overlap_n_x_median\t", summary_value("n_x_median")),
    paste0("overlap_n_y_median\t", summary_value("n_y_median")),
    paste0("overlap_n_median\t", summary_value("n_overlap_median")),
    paste0("overlap_exposure_fraction_median\t", summary_value("exposure_overlap_fraction_median")),
    paste0("overlap_outcome_fraction_median\t", summary_value("outcome_overlap_fraction_median"))
  ), run_log_path)
  saveRDS(result, model_path)

  notes <- c(
    main_result = "Primary overlap-aware second-order Founder GLS-IVW estimate",
    method_comparison = "Primary and diagnostic sensitivity/comparator estimates",
    per_locus_results = "Per-locus second-order diagnostic table",
    block_first_stage = "Block first-stage heuristic diagnostics",
    psd_qc = "Star-space and raw-space PSD diagnostics; raw-space is diagnostic only",
    overlap_qc = "Per-locus overlap counts, overlap fractions, and residual-correlation diagnostics",
    loo_mr = "Leave-one-locus-out MR sensitivity table",
    instrument_sensitivity = "Instrument-set sensitivity refit table",
    q_calibration = "Reported global Q heterogeneity p-value with chi-square and simulation-calibration fields",
    first_stage_summary = "Summary of block F_X heuristic diagnostics",
    founder_residuals = "Descriptive founder residual localization diagnostics",
    founder_coordinate_ratios = "Diagnostic founder-coordinate Wald ratios",
    warnings = "Warnings with consequence and recommended action"
  )
  output_manifest <- qtl2mr_build_output_manifest(output_dir, notes = notes)
  result$output_manifest <- output_manifest
  result <- qtl2mr_add_stable_slots(result)
  utils::write.csv(output_manifest, manifest_path, row.names = FALSE)
  qtl2mr_write_output_manifest_html(output_manifest, manifest_html_path)
  write_summary_md(result, summary_path)
  if (isTRUE(make_html)) write_qtl2mr_html_report(result, html_path)

  output_manifest <- qtl2mr_build_output_manifest(output_dir, notes = notes)
  result$output_manifest <- output_manifest
  result$paths <- output_manifest
  utils::write.csv(output_manifest, manifest_path, row.names = FALSE)
  qtl2mr_write_output_manifest_html(output_manifest, manifest_html_path)
  if (isTRUE(make_html)) write_qtl2mr_html_report(result, html_path)
  saveRDS(result, model_path)
  qtl2mr_progress(progress_callback, 5L, "Write outputs", "publish staged result")
  qtl2mr_commit_output_transaction(transaction, keep_backup = TRUE)
  output_manifest <- qtl2mr_build_output_manifest(final_output_dir, notes = notes)
  result$output_manifest <- output_manifest
  result$paths <- output_manifest
  utils::write.csv(output_manifest, file.path(final_output_dir, "data", "output_manifest.csv"), row.names = FALSE)
  qtl2mr_write_output_manifest_html(output_manifest, file.path(final_output_dir, "output_manifest.html"))
  if (isTRUE(make_html)) write_qtl2mr_html_report(result, file.path(final_output_dir, "report.html"))
  saveRDS(result, file.path(final_output_dir, "result_object.rds"))
  qtl2mr_finalize_output_transaction(transaction)
  published <- TRUE
  result
}

run_qtl2mr_analysis <- function(input,
                                exposure = NULL,
                                outcome = NULL,
                                output_dir = NULL,
                                loci = NULL,
                                overwrite = FALSE,
                                make_html = TRUE,
                                make_figures = TRUE,
                                verbose = TRUE,
                                diagnostic_bootstrap_B = 500L,
                                diagnostic_bootstrap_seed = 1L,
                                reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
                                heterogeneity_method = c("auto", "simulation", "chisq", "none"),
                                heterogeneity_nsim = 1000L,
                                heterogeneity_seed = 1L,
                                run_loo_mr = TRUE,
                                run_instrument_sensitivity = TRUE,
                                progress = NULL,
                                ...) {
  dots <- list(...)
  allowed_dots <- c("standardize_phenotypes", "cross_cov_model", "ridge",
    "engine", "context_cache", "fit_cache")
  dot_names <- names(dots)
  if (length(dots) && (is.null(dot_names) || anyNA(dot_names) ||
      any(!nzchar(dot_names)) || anyDuplicated(dot_names) ||
      any(!dot_names %in% allowed_dots))) {
    bad <- if (is.null(dot_names)) "<unnamed>" else unique(ifelse(is.na(dot_names) | !nzchar(dot_names),
      "<unnamed>", dot_names)[!dot_names %in% allowed_dots | duplicated(dot_names)])
    stop("Unknown, unnamed or duplicate analysis argument in ...: ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed_dots, collapse = ", "), call. = FALSE)
  }
  if (!inherits(input, "qtl2mr_prepared_input")) {
    stop("input must be created by as_qtl2mr_prepared_input()", call. = FALSE)
  }
  reported_estimator_policy <- match.arg(reported_estimator_policy)
  heterogeneity_method <- match.arg(heterogeneity_method)
  exposure <- exposure %||% input$exposure
  outcome <- outcome %||% input$outcome
  if (is.null(exposure) || is.null(outcome)) {
    stop(
      "run_qtl2mr_analysis() requires exposure and outcome trait names. ",
      "Prepared inputs can store these names via as_qtl2mr_prepared_input(); ",
      "otherwise pass exposure and outcome explicitly.",
      call. = FALSE
    )
  }
  destination_dir <- qtl2mr_analysis_output_dir(exposure, outcome, output_dir)
  transaction <- qtl2mr_begin_output_transaction(destination_dir, overwrite = overwrite)
  on.exit({
    if (dir.exists(transaction$stage_dir) || (!is.null(transaction$backup_dir) && dir.exists(transaction$backup_dir))) {
      qtl2mr_abort_output_transaction(transaction)
    }
  }, add = TRUE)
  output_dir <- transaction$stage_dir
  data_dir <- file.path(output_dir, "data")
  figure_dir <- file.path(output_dir, "figures")
  log_dir <- file.path(output_dir, "logs")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  progress_callback <- qtl2mr_make_progress_reporter(
    progress = progress,
    verbose = verbose,
    log_file = file.path(log_dir, "progress.log")
  )
  qtl2mr_progress(
    progress_callback,
    1L,
    "Prepare input/output",
    sprintf("loci=%s, output=%s", qtl2mr_loci_label(loci), transaction$destination)
  )

  qtl2mr_progress(progress_callback, 2L, "Fit selected-locus founder LMM/MR blocks", "starting")
  cap <- qtl2mr_capture_messages(
    do.call(run_prepared_qtl2mr, c(list(
      input = input,
      exposure = exposure,
      outcome = outcome,
      loci = loci,
      output_dir = data_dir,
      write_outputs = TRUE,
      verbose = verbose,
      progress_callback = progress_callback,
      heterogeneity_method = heterogeneity_method,
      q_sim_n = heterogeneity_nsim,
      q_sim_seed = heterogeneity_seed,
      diagnostic_bootstrap_B = diagnostic_bootstrap_B,
      diagnostic_bootstrap_seed = diagnostic_bootstrap_seed,
      reported_estimator_policy = reported_estimator_policy,
      drop_nonestimable_loci = FALSE
    ), dots))
  )
  pilot <- cap$value
  writeLines(cap$messages, file.path(log_dir, "fit_messages.log"))

  loo_mr <- NULL
  loo_error <- NULL
  if (isTRUE(run_loo_mr)) {
    qtl2mr_progress(progress_callback, 3L, "Run leave-one-locus-out MR", "starting")
    loo_mr <- tryCatch(
      loo_multilocus_founder_mr(pilot$fit),
      error = function(e) {
        loo_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(loo_error)) writeLines(loo_error, file.path(log_dir, "loo_mr_error.log"))
    qtl2mr_progress(
      progress_callback,
      3L,
      "Run leave-one-locus-out MR",
      if (is.null(loo_error)) sprintf("complete: %d rows", nrow(loo_mr %||% data.frame())) else paste("failed:", loo_error)
    )
  } else {
    qtl2mr_progress(progress_callback, 3L, "Run leave-one-locus-out MR", "skipped")
  }

  instrument_sensitivity <- NULL
  instrument_sensitivity_error <- NULL
  if (isTRUE(run_instrument_sensitivity)) {
    qtl2mr_progress(progress_callback, 4L, "Run instrument-set sensitivity", "starting")
    instrument_sensitivity <- tryCatch(
      fit_instrument_set_sensitivity(
        structure(list(pilot = pilot), class = c("qtl2mr_analysis_result", "list")),
        loo = loo_mr
      ),
      error = function(e) {
        instrument_sensitivity_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(instrument_sensitivity_error)) writeLines(instrument_sensitivity_error, file.path(log_dir, "instrument_sensitivity_error.log"))
    qtl2mr_progress(
      progress_callback,
      4L,
      "Run instrument-set sensitivity",
      if (is.null(instrument_sensitivity_error)) sprintf("complete: %d sets", nrow(instrument_sensitivity %||% data.frame())) else paste("failed:", instrument_sensitivity_error)
    )
  } else {
    qtl2mr_progress(progress_callback, 4L, "Run instrument-set sensitivity", "skipped")
  }

  result <- structure(
    list(
      input = "qtl2mr_prepared_input",
      exposure = exposure,
      outcome = outcome,
      output_dir = transaction$destination,
      created_at = Sys.time(),
      pilot = pilot,
      loo_mr = loo_mr,
      instrument_sensitivity = instrument_sensitivity,
      loo_error = loo_error,
      instrument_sensitivity_error = instrument_sensitivity_error
    ),
    class = c("qtl2mr_analysis_result", "list")
  )
  result <- write_qtl2mr_results(
    result,
    output_dir = output_dir,
    overwrite = TRUE,
    make_html = make_html,
    make_figures = make_figures,
    progress_callback = progress_callback,
    .transaction = transaction
  )
  if (isTRUE(verbose)) print(result)
  invisible(result)
}

export_qtl2mr_lite_result <- function(result_dir,
                                      dest_zip = NULL,
                                      max_file_mb = 10,
                                      include_rds = FALSE) {
  if (inherits(result_dir, "qtl2mr_analysis_result")) result_dir <- result_dir$output_dir
  if (!is.character(result_dir) || length(result_dir) != 1L || !dir.exists(result_dir)) {
    stop("result_dir must be an existing qtl2mr result directory or qtl2mr_analysis_result", call. = FALSE)
  }
  result_state <- qtl2mr_assert_known_path_node(result_dir, "result_dir")
  if (!identical(result_state, "directory") || qtl2mr_tree_has_symlink(result_dir)) {
    stop("result_dir must not contain symlinks for lite export", call. = FALSE)
  }
  result_dir <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  dest_zip <- dest_zip %||% file.path(dirname(result_dir), paste0(basename(result_dir), "_lite.zip"))
  raw_dest_zip <- qtl2mr_trim_node_path(dest_zip)
  raw_dest_state <- qtl2mr_path_node_state(raw_dest_zip)
  if (identical(raw_dest_state, "link")) stop("dest_zip must not be a symlink", call. = FALSE)
  dest_parent <- qtl2mr_ensure_real_directory(dirname(raw_dest_zip), "dest_zip parent")
  dest_zip <- file.path(dest_parent, basename(raw_dest_zip))
  dest_state <- qtl2mr_assert_known_path_node(dest_zip, "dest_zip")
  if (identical(dest_state, "link")) stop("dest_zip must not be a symlink", call. = FALSE)
  if (!identical(dest_state, "absent")) stop("dest_zip already exists; choose a new destination", call. = FALSE)
  if (startsWith(dest_zip, paste0(sub("/+$", "", result_dir), "/"))) {
    stop("dest_zip must be outside result_dir", call. = FALSE)
  }
  zip_parent <- dirname(dest_zip)
  if (!dir.exists(zip_parent)) dir.create(zip_parent, recursive = TRUE, showWarnings = FALSE)
  if (!identical(qtl2mr_assert_known_path_node(zip_parent, "dest_zip parent"), "directory")) {
    stop("dest_zip parent must be a real writable directory", call. = FALSE)
  }
  lite_root <- qtl2mr_transaction_path(zip_parent, basename(result_dir), "lite")
  if (!dir.create(lite_root, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not create unique lite-export staging directory", call. = FALSE)
  }
  zip_lock <- qtl2mr_transaction_path(zip_parent, basename(dest_zip), "lite-lock")
  if (!file.create(zip_lock)) stop("Could not reserve qtl2mr lite ZIP destination", call. = FALSE)
  temp_zip <- tempfile(paste0(".", basename(dest_zip), ".qtl2mr-"), tmpdir = zip_parent, fileext = ".zip")
  published <- FALSE
  on.exit({
    if (identical(qtl2mr_path_node_state(lite_root), "directory")) unlink(lite_root, recursive = TRUE, force = TRUE)
    if (!published && identical(qtl2mr_path_node_state(temp_zip), "file")) unlink(temp_zip, force = TRUE)
    if (identical(qtl2mr_path_node_state(zip_lock), "file")) unlink(zip_lock, force = TRUE)
  }, add = TRUE)

  full_manifest <- qtl2mr_build_output_manifest(
    output_dir = result_dir,
    max_file_mb = max_file_mb,
    include_rds = include_rds
  )
  if (nrow(full_manifest) && !all(vapply(full_manifest$path_relative, qtl2mr_is_safe_relative_path, logical(1)))) {
    stop("Full result manifest contains an unsafe relative path", call. = FALSE)
  }
  full_manifest$lite_status <- ifelse(full_manifest$included_in_lite, "included", "excluded")
  full_manifest$exclusion_reason <- ""
  full_manifest$exclusion_reason[!full_manifest$included_in_lite & full_manifest$format == "rds"] <-
    "excluded from lite: complete R objects are large; retained in full local result"
  full_manifest$exclusion_reason[!full_manifest$included_in_lite & full_manifest$size_kb > max_file_mb * 1024] <-
    paste0("excluded from lite: file exceeds max_file_mb=", max_file_mb)
  full_manifest$exclusion_reason[startsWith(full_manifest$path_relative, "logs/")] <-
    "excluded from lite: local run and progress logs may contain absolute paths"
  full_manifest$exclusion_reason[basename(full_manifest$path_relative) == "output_manifest.csv"] <-
    "excluded from lite: regenerated as LITE_MANIFEST.csv"
  full_manifest$exclusion_reason[basename(full_manifest$path_relative) == "output_manifest.html"] <-
    "excluded from lite: regenerated from portable relative paths"

  include <- full_manifest[full_manifest$included_in_lite, , drop = FALSE]
  for (i in seq_len(nrow(include))) {
    src <- include$path_absolute[[i]]
    dst <- file.path(lite_root, include$path_relative[[i]])
    src_state <- qtl2mr_assert_known_path_node(src, "lite source member")
    if (!qtl2mr_is_safe_relative_path(include$path_relative[[i]]) || identical(src_state, "link")) {
      stop("Lite export encountered an unsafe source member", call. = FALSE)
    }
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dst, overwrite = FALSE)) stop("Could not copy lite-export member", call. = FALSE)
  }

  audit <- full_manifest[, setdiff(names(full_manifest), "path_absolute"), drop = FALSE]
  utils::write.csv(audit, file.path(lite_root, "LITE_EXPORT_AUDIT.csv"), row.names = FALSE)
  writeLines(c(
    "# qtl2mr Lite Result",
    "",
    "This is a portable lightweight qtl2mr result package for external discussion.",
    "",
    paste0("- max_file_mb: ", max_file_mb),
    paste0("- include_rds: ", include_rds),
    "",
    "Large complete R objects and large sample-level tables are excluded by default.",
    "See LITE_EXPORT_AUDIT.csv for inclusion/exclusion decisions.",
    "Use the full local result directory when complete R object reload is needed."
  ), file.path(lite_root, "README_LITE.md"))
  # Scrub copied source-derived text before generating the portable manifest so
  # manifest sizes and the HTML provenance view describe the shipped bytes.
  text_like <- list.files(lite_root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  text_like <- text_like[tolower(tools::file_ext(text_like)) %in% c("csv", "tsv", "txt", "log", "md", "html", "json", "yaml", "yml", "xml")]
  for (path in text_like) {
    writeLines(qtl2mr_scrub_portable_text(readLines(path, warn = FALSE, encoding = "UTF-8"), result_dir), path, useBytes = TRUE)
  }
  # Seed both regenerated provenance files before taking the portable inventory.
  # The final pass intentionally includes these self-referential members; their
  # membership is stable even though their recorded sizes reflect the pass just
  # before serialization.
  file.create(file.path(lite_root, "LITE_MANIFEST.csv"))
  file.create(file.path(lite_root, "output_manifest.html"))
  lite_manifest <- qtl2mr_portable_manifest(lite_root, max_file_mb = Inf, include_rds = TRUE)
  utils::write.csv(lite_manifest, file.path(lite_root, "LITE_MANIFEST.csv"), row.names = FALSE)
  qtl2mr_write_output_manifest_html(lite_manifest, file.path(lite_root, "output_manifest.html"))
  lite_manifest <- qtl2mr_portable_manifest(lite_root, max_file_mb = Inf, include_rds = TRUE)
  utils::write.csv(lite_manifest, file.path(lite_root, "LITE_MANIFEST.csv"), row.names = FALSE)
  qtl2mr_write_output_manifest_html(lite_manifest, file.path(lite_root, "output_manifest.html"))
  # Final pass is detect-only: never mutate files after their manifest rows and
  # HTML representation have been generated.
  text_like <- list.files(lite_root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  text_like <- text_like[tolower(tools::file_ext(text_like)) %in% c("csv", "tsv", "txt", "log", "md", "html", "json", "yaml", "yml", "xml")]
  for (path in text_like) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (qtl2mr_has_portable_path_leak(text, result_dir, lite_root)) {
      stop("Lite export portability scan found an unresolved absolute path", call. = FALSE)
    }
  }
  members <- list.files(lite_root, recursive = TRUE, full.names = FALSE, all.files = FALSE, no.. = TRUE)
  if (!length(members) || !all(vapply(members, qtl2mr_is_safe_relative_path, logical(1)))) {
    stop("Lite export has no safe relative members", call. = FALSE)
  }
  final_manifest <- utils::read.csv(file.path(lite_root, "LITE_MANIFEST.csv"), stringsAsFactors = FALSE)
  if (!identical(sort(members), sort(final_manifest$path_relative))) {
    stop("Lite export manifest does not cover every archive member", call. = FALSE)
  }

  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(lite_root)
  invisible(utils::capture.output(
    utils::zip(zipfile = temp_zip, files = members, flags = "-rq9X")
  ))
  final_dest_state <- qtl2mr_assert_known_path_node(dest_zip, "dest_zip")
  if (!identical(final_dest_state, "absent")) {
    stop("dest_zip appeared during lite export; existing destination was preserved", call. = FALSE)
  }
  # file.link is a same-filesystem, no-replace publication primitive: unlike
  # rename it fails if another process wins the destination race.
  if (!file.exists(temp_zip) || !file.link(temp_zip, dest_zip) || unlink(temp_zip, force = TRUE) != 0L) {
    stop("Could not publish qtl2mr lite ZIP", call. = FALSE)
  }
  published <- TRUE
  invisible(normalizePath(dest_zip, winslash = "/", mustWork = FALSE))
}

report_qtl2mr_locus <- function(result,
                                locus_id,
                                output_dir = NULL,
                                make_html = TRUE,
                                make_figures = TRUE) {
  if (is.character(result) && length(result) == 1L && file.exists(result)) result <- readRDS(result)
  if (!inherits(result, "qtl2mr_analysis_result")) stop("result must be a qtl2mr_analysis_result or path to result_object.rds", call. = FALSE)
  block <- result$pilot$locus_blocks[[locus_id]]
  if (is.null(block)) stop("locus_id not found: ", locus_id, call. = FALSE)
  output_dir <- output_dir %||% file.path(result$output_dir, "locus_reports", qtl2mr_slug(locus_id))
  data_dir <- file.path(output_dir, "data")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  theta <- result$pilot$summary$theta[[1L]]
  local_fit <- mr_block_so_diagnostics(block, theta = NULL, allow_cross_cov = TRUE)
  beta_table <- data.frame(
    coordinate = names(block$beta_x),
    beta_x = as.numeric(block$beta_x),
    beta_y = as.numeric(block$beta_y),
    se_x = sqrt(diag(block$Sigma_x)),
    se_y = sqrt(diag(block$Sigma_y)),
    stringsAsFactors = FALSE
  )
  cov_qc <- data.frame(
    marker = block$marker,
    chr = block$chr,
    theta_global = theta,
    theta_locus = local_fit$theta,
    Q_min = local_fit$Q_min,
    F_X_block = mr_block_first_stage_f(block),
    joint_psd_valid = block$joint_psd_valid %||% NA,
    star_space_psd_valid = block$star_space_psd_valid %||% NA,
    raw_space_psd_valid = block$raw_space_psd_valid %||% NA,
    stringsAsFactors = FALSE
  )
  residuals <- founder_residual_scores(block, theta = theta, allow_cross_cov = TRUE)
  ratios <- founder_coordinate_ratios(block, theta_ref = theta, allow_cross_cov = TRUE)
  locus_summary <- data.frame(
    marker = block$marker,
    chr = block$chr,
    theta_global = theta,
    theta_locus = local_fit$theta,
    se_locus = local_fit$se,
    profile_ci_lower_locus = local_fit$profile_ci_lower %||% NA_real_,
    profile_ci_upper_locus = local_fit$profile_ci_upper %||% NA_real_,
    Q_min_locus = local_fit$Q_min,
    F_X_block = mr_block_first_stage_f(block),
    n_overlap = block$n_overlap %||% NA_integer_,
    rho_xy = block$rho_xy %||% NA_real_,
    cross_cov_model = block$cross_cov_model %||% NA_character_,
    joint_psd_status = block$joint_psd_status %||% NA_character_,
    stringsAsFactors = FALSE
  )
  utils::write.csv(beta_table, file.path(data_dir, "locus_betas.csv"), row.names = FALSE)
  utils::write.csv(cov_qc, file.path(data_dir, "locus_covariance_qc.csv"), row.names = FALSE)
  utils::write.csv(residuals, file.path(data_dir, "founder_residuals.csv"), row.names = FALSE)
  utils::write.csv(ratios, file.path(data_dir, "founder_coordinate_ratios.csv"), row.names = FALSE)
  utils::write.csv(locus_summary, file.path(data_dir, "locus_summary.csv"), row.names = FALSE)
  if (isTRUE(make_figures)) {
    qtl2mr_write_figure_pair(file.path(figure_dir, "founder_beta_comparison"), {
      yy <- seq_len(nrow(beta_table))
      xr <- range(c(beta_table$beta_x, beta_table$beta_y, 0), finite = TRUE)
      pad <- diff(xr) * 0.12
      if (!is.finite(pad) || pad <= 0) pad <- 0.1
      graphics::plot(beta_table$beta_x, yy, xlim = xr + c(-pad, pad), yaxt = "n",
                     xlab = "Founder-coordinate effect", ylab = "",
                     main = paste0("Founder effects at ", locus_id), pch = 16, col = "#00796b")
      graphics::points(beta_table$beta_y, yy, pch = 17, col = "#4c63d2")
      graphics::abline(v = 0, lty = 3, col = "grey70")
      graphics::axis(2, at = yy, labels = beta_table$coordinate, las = 1, tick = FALSE)
      graphics::legend("topright", legend = c("Exposure", "Outcome"), pch = c(16, 17),
                       col = c("#00796b", "#4c63d2"), bty = "n")
    })
    qtl2mr_write_figure_pair(file.path(figure_dir, "founder_residual_scores"), {
      vals <- residuals$t_one_founder_vs_rest
      names(vals) <- residuals$founder
      cols <- ifelse(vals >= 0, "#00796b", "#b75a3c")
      graphics::barplot(vals, col = cols, border = NA, ylab = "Residual score",
                        main = paste0("Founder residual diagnostic at ", locus_id))
      graphics::abline(h = 0, col = "grey50")
    })
    qtl2mr_write_figure_pair(file.path(figure_dir, "founder_coordinate_ratio_forest"), {
      plot_founder_ratio_forest(ratios, main = paste0("Founder-coordinate ratios at ", locus_id))
    })
    qtl2mr_write_figure_pair(file.path(figure_dir, "covariance_psd_summary"), {
      vals <- c(
        star_psd_valid = as.numeric(isTRUE(block$star_space_psd_valid)),
        raw_psd_valid = as.numeric(isTRUE(block$raw_space_psd_valid)),
        joint_psd_valid = as.numeric(isTRUE(block$joint_psd_valid))
      )
      graphics::barplot(vals, ylim = c(0, 1), col = c("#00796b", "grey65", "#4c63d2"),
                        ylab = "Indicator", main = paste0("Covariance PSD summary at ", locus_id))
      graphics::mtext("Raw-space is diagnostic only.", side = 3, line = 0.2, cex = 0.8, col = "grey35")
    }, width = 1400, height = 1000, pdf_width = 7, pdf_height = 5)
  }
  if (isTRUE(make_html)) {
    html <- c(
      "<!doctype html><html><head><meta charset=\"utf-8\"><title>qtl2mr locus report</title>",
      "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:28px;line-height:1.45;}table{border-collapse:collapse;width:100%;}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left;}th{background:#eef3f1;}img{max-width:100%;}figure{border:1px solid #ddd;padding:12px;margin:16px 0;}code{background:#eef3f1;padding:1px 4px;}</style></head><body>",
      sprintf("<h1>Locus report: %s</h1>", qtl2mr_html_escape(locus_id)),
      "<h2>Locus Summary</h2>",
      qtl2mr_html_table(locus_summary, digits = 5),
      "<h2>Covariance QC</h2>",
      qtl2mr_html_table(cov_qc),
      "<h2>Founder-contrast effects</h2>",
      qtl2mr_html_table(beta_table, digits = 5),
      "<p>Founder-coordinate ratios and residual scores are diagnostic only. They are not separate causal estimands.</p>",
      if (isTRUE(make_figures)) paste0(
        "<h2>Figures</h2>",
        "<figure><img src=\"figures/founder_beta_comparison.png\"><figcaption>Founder-coordinate exposure and outcome effects.</figcaption></figure>",
        "<figure><img src=\"figures/founder_residual_scores.png\"><figcaption>Founder residual localization diagnostic.</figcaption></figure>",
        "<figure><img src=\"figures/founder_coordinate_ratio_forest.png\"><figcaption>Founder-coordinate Wald ratios are diagnostic only.</figcaption></figure>",
        "<figure><img src=\"figures/covariance_psd_summary.png\"><figcaption>Covariance PSD summary; raw-space is diagnostic only.</figcaption></figure>"
      ) else "",
      "</body></html>"
    )
    writeLines(html, file.path(output_dir, "report.html"))
  }
  invisible(normalizePath(output_dir, mustWork = FALSE))
}
