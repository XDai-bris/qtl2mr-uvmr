# Instrument subsets (0.6.1.9004). The instruments a pair is fitted on form the
# maximal set ("all"); each subset keeps the fitted loci of that pair that pass
# its rules and refits only the multilocus stage on their blocks, with the run's
# settings. The LMMs are not refitted, so a subset row equals the main row of a
# separate run whose instruments are that subset's loci.

qtl2mr_subset_rule_keys <- c("F_X_min", "steiger")

# Summary columns that describe the LMM stage of the pair rather than the
# subset's fit; subset rows carry their own n_loci instead.
qtl2mr_subset_dropped_summary <- c("exposure", "outcome", "n_loci", "n_loci_selected", "n_loci_dropped_nonestimable")

qtl2mr_subset_id_columns <- c("exposure", "outcome", "subset", "rules", "n_loci", "status", "same_loci_as",
  "error_message")

# Validates a subsets specification before anything is fitted. Returns NULL, or
# the named list with each subset's rules in the order F_X_min, steiger.
qtl2mr_check_subsets <- function(subsets) {
  if (is.null(subsets)) return(NULL)
  if (!is.list(subsets) || is.data.frame(subsets) || !length(subsets)) {
    stop("subsets must be NULL or a non-empty named list of rule lists, for example ",
      "list(F2 = list(F_X_min = 2), F1.5_SF = list(F_X_min = 1.5, steiger = TRUE))", call. = FALSE)
  }
  nms <- names(subsets)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) stop("every subset must have a non-empty name", call. = FALSE)
  bad <- nms[!grepl("^[A-Za-z0-9._]+$", nms)]
  if (length(bad)) {
    stop("subset names may use only letters, digits, '.' and '_': ", paste(unique(bad), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(nms)) stop("subset names must be unique: ", paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)
  if ("all" %in% nms) stop("the subset name 'all' is reserved for the maximal set", call. = FALSE)
  out <- lapply(nms, function(nm) {
    rules <- subsets[[nm]]
    if (!is.list(rules) || is.data.frame(rules)) {
      stop("subset ", nm, " must be a list of rules, for example list(F_X_min = 2)", call. = FALSE)
    }
    if (!length(rules)) stop("subset ", nm, " has no rules", call. = FALSE)
    keys <- names(rules)
    if (is.null(keys) || anyNA(keys) || any(!nzchar(keys))) stop("every rule of subset ", nm, " must be named", call. = FALSE)
    if (anyDuplicated(keys)) stop("subset ", nm, " repeats rule ", keys[duplicated(keys)][[1L]], call. = FALSE)
    unknown <- setdiff(keys, qtl2mr_subset_rule_keys)
    if (length(unknown)) {
      stop("subset ", nm, " has unknown rule(s) ", paste(unknown, collapse = ", "),
        "; the rules are F_X_min and steiger", call. = FALSE)
    }
    norm <- list()
    if ("F_X_min" %in% keys) {
      v <- rules[["F_X_min"]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
        stop("subset ", nm, ": F_X_min must be one finite number >= 0", call. = FALSE)
      }
      norm$F_X_min <- as.numeric(v)
    }
    if ("steiger" %in% keys) {
      v <- rules[["steiger"]]
      if (!is.logical(v) || length(v) != 1L || is.na(v)) stop("subset ", nm, ": steiger must be TRUE or FALSE", call. = FALSE)
      norm$steiger <- v
    }
    norm
  })
  names(out) <- nms
  out
}

qtl2mr_subset_rules_text <- function(rules) {
  parts <- c(if (!is.null(rules[["F_X_min"]])) paste0("F_X_min>=", as.character(rules[["F_X_min"]])),
    if (isTRUE(rules[["steiger"]])) "steiger")
  if (length(parts)) paste(parts, collapse = ";") else "none"
}

# Loci of a pair's first-stage table that pass a subset's rules. F_X_min keeps
# a finite F_X_block >= the value; steiger = TRUE keeps steiger_pass TRUE, so
# FALSE and NA are dropped.
qtl2mr_subset_members <- function(first_stage, rules) {
  keep <- rep(TRUE, nrow(first_stage))
  F_min <- rules[["F_X_min"]]
  if (!is.null(F_min)) keep <- keep & is.finite(first_stage$F_X_block) & first_stage$F_X_block >= F_min
  if (isTRUE(rules[["steiger"]])) keep <- keep & first_stage$steiger_pass %in% TRUE
  keep
}

# The multilocus stage and summary of a pair refitted on some of its blocks,
# returned as a pilot-like list (summary, fit, locus_blocks). Blocks keep the
# pair's order; phenotype scaling rows are per block, so the pair's rows for
# these blocks are the rows a separate fit would compute.
qtl2mr_subset_refit <- function(pilot, keep, settings, ridge = 1e-10, progress_callback = NULL) {
  blocks <- pilot$locus_blocks[keep]
  stage <- qtl2mr_multilocus_stage(blocks, ridge = ridge, heterogeneity_method = settings$heterogeneity_method,
    q_sim_n = settings$heterogeneity_nsim, q_sim_seed = settings$heterogeneity_seed,
    diagnostic_bootstrap_B = settings$diagnostic_bootstrap_B,
    diagnostic_bootstrap_seed = settings$diagnostic_bootstrap_seed,
    reported_estimator_policy = settings$reported_estimator_policy, progress_callback = progress_callback)
  scaling <- pilot$phenotype_scaling
  scaling <- scaling[scaling$marker %in% names(blocks), , drop = FALSE]
  rows <- qtl2mr_multilocus_summary(blocks, stage, scaling, qtl2mr_locus_first_stage(blocks, ridge = ridge),
    exposure = pilot$summary$exposure[[1L]], outcome = pilot$summary$outcome[[1L]],
    n_loci_selected = length(blocks), n_loci_dropped = 0L,
    standardize_phenotypes = settings$standardize_phenotypes, cross_cov_model = settings$cross_cov_model,
    diagnostic_bootstrap_B = settings$diagnostic_bootstrap_B,
    diagnostic_bootstrap_seed = settings$diagnostic_bootstrap_seed)
  list(summary = rows$summary, fit = stage$fit, locus_blocks = blocks)
}

# Subset rows of one fitted pair, in the order of `subsets`. Each row has the
# identification columns and the pair's summary columns (other than
# qtl2mr_subset_dropped_summary). A subset whose loci equal those of the maximal
# set or of an earlier subset copies that row (same_loci_as). `finish` maps a
# pilot-like list to its reported summary (identity in the batch; the single
# analysis completes the Q fields as for its main row). Returns rows, the
# membership of each fitted locus, and the errors of failed refits.
qtl2mr_pair_subsets <- function(pilot, subsets, settings, min_loci = 1L, ridge = 1e-10, progress_callback = NULL,
                                finish = function(p) p$summary) {
  blocks <- pilot$locus_blocks
  markers <- names(blocks)
  first <- pilot$locus_first_stage %||% qtl2mr_locus_first_stage(blocks, ridge = ridge)
  first <- first[match(markers, first$marker), , drop = FALSE]
  main <- finish(pilot)
  template <- main[NA_integer_, , drop = FALSE]
  keep_cols <- setdiff(names(main), qtl2mr_subset_dropped_summary)
  make_row <- function(nm, members, status, summary_row, same = NA_character_, message = "") {
    s <- summary_row[1L, names(template), drop = FALSE]
    rownames(s) <- NULL
    cbind(data.frame(exposure = main$exposure[[1L]], outcome = main$outcome[[1L]], subset = nm,
      rules = qtl2mr_subset_rules_text(subsets[[nm]]), n_loci = as.integer(sum(members)), status = status,
      same_loci_as = same, error_message = message, stringsAsFactors = FALSE), s[keep_cols])
  }
  done <- list(all = list(members = rep(TRUE, length(markers)), status = "ok", summary = main, message = "",
    error = NULL))
  rows <- list()
  membership <- list()
  errors <- list()
  for (nm in names(subsets)) {
    members <- tryCatch(qtl2mr_subset_members(first, subsets[[nm]]), error = function(e) e)
    if (inherits(members, "error")) {
      errors[[nm]] <- members
      membership[[nm]] <- rep(NA, length(markers))
      rows[[nm]] <- make_row(nm, logical(), "error", template, message = conditionMessage(members))
      rows[[nm]]$n_loci <- NA_integer_
      next
    }
    membership[[nm]] <- members
    same <- names(done)[vapply(done, function(d) identical(d$members, members), logical(1))]
    if (length(same)) {
      d <- done[[same[[1L]]]]
      if (!is.null(d$error)) errors[[nm]] <- d$error
      rows[[nm]] <- make_row(nm, members, d$status, d$summary, same = same[[1L]], message = d$message)
    } else if (sum(members) < min_loci) {
      message <- sprintf("too_few_loci: %d loci pass the subset rules, below min_loci = %d", sum(members), min_loci)
      done[[nm]] <- list(members = members, status = "too_few_loci", summary = template, message = message, error = NULL)
      rows[[nm]] <- make_row(nm, members, "too_few_loci", template, message = message)
    } else {
      refit <- tryCatch(finish(qtl2mr_subset_refit(pilot, members, settings, ridge = ridge,
        progress_callback = progress_callback)), error = function(e) e)
      if (inherits(refit, "error")) {
        errors[[nm]] <- refit
        done[[nm]] <- list(members = members, status = "error", summary = template,
          message = conditionMessage(refit), error = refit)
        rows[[nm]] <- make_row(nm, members, "error", template, message = conditionMessage(refit))
      } else {
        done[[nm]] <- list(members = members, status = "ok", summary = refit, message = "", error = NULL)
        rows[[nm]] <- make_row(nm, members, "ok", refit)
      }
    }
  }
  out <- do.call(rbind, unname(rows))
  rownames(out) <- NULL
  list(rows = out, markers = markers, membership = membership, errors = errors)
}

# Batch columns of subset_results: the identification columns, then the
# main_results summary columns that describe the fit.
qtl2mr_batch_subset_columns <- function() {
  c("pair_id", qtl2mr_subset_id_columns, setdiff(qtl2mr_batch_summary_columns, qtl2mr_subset_dropped_summary))
}

qtl2mr_empty_subset_results <- function() qtl2mr_empty_dataframe(qtl2mr_batch_subset_columns())

# Rows for every subset of a pair when the subset step itself fails (outside
# the guard of a single refit): status "error" with the message, no estimates.
# `template` is a one-row table whose columns follow the identification columns.
qtl2mr_subset_failure_rows <- function(subsets, exposure, outcome, message, template) {
  n <- length(subsets)
  ids <- data.frame(exposure = rep(exposure, n), outcome = rep(outcome, n), subset = names(subsets),
    rules = vapply(subsets, qtl2mr_subset_rules_text, character(1), USE.NAMES = FALSE),
    n_loci = rep(NA_integer_, n), status = rep("error", n), same_loci_as = rep(NA_character_, n),
    error_message = rep(message, n), stringsAsFactors = FALSE)
  rest <- template[rep(NA_integer_, n), setdiff(names(template), names(ids)), drop = FALSE]
  rownames(rest) <- NULL
  cbind(ids, rest)
}

qtl2mr_subset_error_rows <- function(pair, errors) {
  do.call(rbind, lapply(names(errors), function(nm) {
    row <- qtl2mr_error_row(pair, "subset", errors[[nm]])
    row$error_message <- paste0("subset ", nm, ": ", row$error_message)
    row
  }))
}

# A fitted batch pair: its subset rows, an in_subset_<name> column per subset on
# its locus rows (NA on rows that were not fitted), and an errors row with
# stage "subset" per failed subset. Everything is built first and attached only
# when complete; if any part of the step fails, every subset of the pair gets an
# "error" row and an errors row, and the pair's main row and locus rows are kept
# as they are.
qtl2mr_batch_add_subsets <- function(result, pair, pilot, subsets, settings) {
  parts <- tryCatch({
    sub <- qtl2mr_pair_subsets(pilot, subsets, settings, min_loci = settings$min_loci)
    rows <- sub$rows
    rows$pair_id <- rep(pair$pair_id, nrow(rows))
    rows <- rows[, qtl2mr_batch_subset_columns(), drop = FALSE]
    loci <- result$loci
    idx <- match(loci$marker, sub$markers)
    fitted <- loci$kept %in% TRUE & !is.na(idx)
    for (nm in names(subsets)) {
      value <- rep(NA, nrow(loci))
      value[fitted] <- sub$membership[[nm]][idx[fitted]]
      loci[[paste0("in_subset_", nm)]] <- value
    }
    list(rows = rows, loci = loci, errors = if (length(sub$errors)) qtl2mr_subset_error_rows(pair, sub$errors))
  }, error = function(e) {
    rows <- qtl2mr_subset_failure_rows(subsets, pair$exposure, pair$outcome, conditionMessage(e),
      qtl2mr_empty_subset_results())
    rows$pair_id <- rep(pair$pair_id, nrow(rows))
    list(rows = rows[, qtl2mr_batch_subset_columns(), drop = FALSE], loci = NULL,
      errors = qtl2mr_subset_error_rows(pair, stats::setNames(rep(list(e), length(subsets)), names(subsets))))
  })
  result$subsets <- parts$rows
  if (!is.null(parts$loci)) result$loci <- parts$loci
  if (!is.null(parts$errors)) result$subset_errors <- parts$errors
  result
}

# subset_results has one row per main_results pair and subset, in that order.
# Pairs without a fit (errors, including preflight and worker failures) get
# status "pair_not_fitted". locus_results ends with the in_subset_ columns.
qtl2mr_batch_finalize_subsets <- function(batch) {
  subsets <- batch$config$subsets
  have <- batch$subset_results %||% qtl2mr_empty_subset_results()
  main <- batch$main_results
  want_pair <- rep(main$pair_id, each = length(subsets))
  want_subset <- rep(names(subsets), times = nrow(main))
  idx <- match(paste(want_pair, want_subset, sep = "\r"), paste(have$pair_id, have$subset, sep = "\r"))
  out <- have[idx, , drop = FALSE]
  missing <- is.na(idx)
  if (any(missing)) {
    at <- match(want_pair[missing], main$pair_id)
    out$pair_id[missing] <- want_pair[missing]
    out$exposure[missing] <- main$exposure[at]
    out$outcome[missing] <- main$outcome[at]
    out$subset[missing] <- want_subset[missing]
    out$rules[missing] <- vapply(subsets[want_subset[missing]], qtl2mr_subset_rules_text, character(1), USE.NAMES = FALSE)
    out$status[missing] <- "pair_not_fitted"
    out$error_message[missing] <- ""
  }
  canonical <- qtl2mr_batch_subset_columns()
  for (nm in setdiff(canonical, names(out))) out[[nm]] <- rep(NA, nrow(out))
  out <- out[, c(canonical, setdiff(names(out), canonical)), drop = FALSE]
  rownames(out) <- NULL
  batch$subset_results <- out
  loci <- batch$locus_results
  cols <- paste0("in_subset_", names(subsets))
  for (nm in setdiff(cols, names(loci))) loci[[nm]] <- rep(NA, nrow(loci))
  batch$locus_results <- loci[, c(setdiff(names(loci), cols), cols), drop = FALSE]
  batch
}
