# By-trait planning for run_qtl2mr_batch(). Pairs are fitted exactly as before,
# one at a time and in order; planning only computes beforehand, once per trait
# and sample set, the LMM contexts and locus fits that the pairs will read, so
# that each pair finds them in the session registry. A pair never depends on
# whether its fits were planned: a fit that is missing is computed by the pair,
# and a planned fit is stored under the same exact key the pair looks up.
#
# Pairs are processed in waves (consecutive pairs, or whole exposure groups when
# groups run in parallel workers) sized so that one wave's fits and contexts take
# at most half the registry budget, measured as the registry measures them. After a wave, its entries that no later wave reads become
# the first to be released when the budget is reached; within the budget they
# stay available to later calls in the session.

# Per-pair demand: the two traits and the instruments. Each trait is fitted on
# its own complete samples, so a planned unit is keyed by the trait alone
# (samples = NULL and sig = ""; before 0.6.1.9006 a pair-common sample set could
# be planned here).
# With the outcome founder QC (min_founder_copies), the pair's QC result is kept
# (founder_qc) and the outcome needs only the loci that pass it.
qtl2mr_batch_plan <- function(resources, todo, instruments, min_founder_copies = NULL) {
  lapply(seq_len(nrow(todo)), function(i) {
    pair <- todo[i, , drop = FALSE]
    rows <- instruments[instruments$exposure == pair$exposure, , drop = FALSE]
    tryCatch({
      out <- list(traits = c(pair$exposure, pair$outcome), samples = NULL, sig = "",
        marker = as.character(rows$marker), chr = as.character(rows$chr))
      if (!is.null(min_founder_copies)) {
        out$founder_qc <- qtl2mr_pair_founder_qc(resources, rows, pair$exposure, pair$outcome,
          min_founder_copies)
      }
      out
    }, error = function(e) NULL)
  })
}

# The loci each trait of a planned pair is fitted at: every instrument for the
# exposure; for the outcome, those passing the outcome founder QC when it is on.
qtl2mr_plan_trait_loci <- function(p) {
  keep <- if (is.null(p$founder_qc)) rep(TRUE, length(p$marker)) else p$founder_qc$pass
  list(list(trait = p$traits[[1L]], marker = p$marker, chr = p$chr),
    list(trait = p$traits[[2L]], marker = p$marker[keep], chr = p$chr[keep]))
}

qtl2mr_batch_fit_items <- function(plan) {
  unique(unlist(lapply(plan, function(p) {
    if (is.null(p)) return(character())
    lapply(qtl2mr_plan_trait_loci(p), function(t) {
      if (!length(t$marker)) return(character())
      paste(t$trait, p$sig, t$marker, sep = "\r")
    })
  }), use.names = FALSE))
}

# Release keys of registry entries: trait and chromosome for contexts, plus the
# marker for fits (sample sets are not distinguished, so release is conservative).
qtl2mr_batch_release_items <- function(plan) {
  unique(unlist(lapply(plan, function(p) {
    if (is.null(p)) return(character())
    lapply(qtl2mr_plan_trait_loci(p), function(t) {
      if (!length(t$marker)) return(character())
      c(paste(t$trait, unique(t$chr), sep = "\r"), paste(t$trait, t$chr, t$marker, sep = "\r"))
    })
  }), use.names = FALSE))
}

qtl2mr_registry_item <- function(info) {
  if (is.null(info$trait) || is.null(info$chr)) return(NA_character_)
  if (is.null(info$marker) || is.na(info$marker)) paste(info$trait, info$chr, sep = "\r") else
    paste(info$trait, info$chr, info$marker, sep = "\r")
}

# Distinct contexts (trait, sample set, chromosome) a set of pairs reads.
qtl2mr_batch_context_items <- function(plan) {
  unique(unlist(lapply(plan, function(p) {
    if (is.null(p)) return(character())
    lapply(qtl2mr_plan_trait_loci(p), function(t) {
      if (!length(t$marker)) return(character())
      paste(t$trait, p$sig, unique(t$chr), sep = "\r")
    })
  }), use.names = FALSE))
}

# Waves of whole units (lists of pair indices), in order: a wave closes before
# a unit that would take its distinct fits and contexts, at the sizes the
# registry charges, past max_bytes, once it holds at least min_units units. A
# last wave with fewer than min_units units joins the one before, so that every
# wave of exposure groups runs them in parallel workers.
qtl2mr_batch_waves <- function(plan, units, max_bytes, sizes, min_units = 1L) {
  waves <- list()
  current <- list()
  fits <- contexts <- character()
  for (u in units) {
    u_fits <- qtl2mr_batch_fit_items(plan[u])
    u_contexts <- qtl2mr_batch_context_items(plan[u])
    merged_fits <- union(fits, u_fits)
    merged_contexts <- union(contexts, u_contexts)
    bytes <- length(merged_fits) * sizes$fit + length(merged_contexts) * sizes$context
    if (length(current) >= min_units && bytes > max_bytes) {
      waves[[length(waves) + 1L]] <- current
      current <- list()
      merged_fits <- u_fits
      merged_contexts <- u_contexts
    }
    current[[length(current) + 1L]] <- u
    fits <- merged_fits
    contexts <- merged_contexts
  }
  if (length(current)) {
    if (length(current) < min_units && length(waves)) {
      waves[[length(waves)]] <- c(waves[[length(waves)]], current)
    } else waves[[length(waves) + 1L]] <- current
  }
  lapply(waves, function(w) unlist(w, use.names = FALSE))
}

# object.size() of one locus fit and one stored (compact) context on all the
# samples, built with the package's own fitting code on a probe phenotype: the
# measure the registry charges, and an upper bound for the batch's entries,
# whose sample sets are subsets. Falls back to a conservative formula.
qtl2mr_batch_entry_sizes <- function(resources) {
  n <- max(length(resources$sample_ids), 1L)
  fallback <- list(fit = 8 * n * 100 + 50000, context = 8 * n * 60 + 50000)
  tryCatch(qtl2mr_quietly({
    ids <- resources$sample_ids
    probe <- stats::setNames(as.numeric(seq_along(ids) %% 7L) + 0.5, ids)
    ctx <- prepare_lmm_context(probe, kinship = NULL, addcovar = qtl2mr_prepared_covar(resources),
      sample_ids = ids, standardize_pheno = TRUE, trait_name = "size_probe", chr = "size_probe")
    gp <- resources$allele_probs[[1L]]
    marker <- dimnames(gp)[[3L]][[1L]]
    contrasts <- make_founder_contrasts(dimnames(gp)[[2L]])
    design <- qtl2mr_batch_marker_design(gp, marker, ctx$sample_ids, contrasts)
    fit <- qtl2mr_fit_founder_lmm_marker(design, marker, ctx, contrasts, dimnames(gp)[[2L]], NULL,
      paste0("FC", 1:7), 1e-12)
    list(fit = as.numeric(utils::object.size(fit)),
      context = as.numeric(utils::object.size(qtl2mr_compact_context(ctx))))
  }), error = function(e) fallback)
}

# Work units of a wave: one per (trait, sample set), with the markers it needs
# on each chromosome.
qtl2mr_batch_prefill_units <- function(plan) {
  units <- list()
  for (p in plan) {
    if (is.null(p)) next
    for (t in qtl2mr_plan_trait_loci(p)) {
      id <- paste(t$trait, p$sig, sep = "\r")
      u <- units[[id]] %||% list(trait = t$trait, samples = p$samples, markers = list())
      for (chr in unique(t$chr)) u$markers[[chr]] <- union(u$markers[[chr]], t$marker[t$chr == chr])
      units[[id]] <- u
    }
  }
  units
}

# Contexts and fits of the given units on one chromosome; nothing is signalled
# and errors are left for the pair that needs the fit to report.
qtl2mr_batch_prefill_chr <- function(input, unit, chr, settings, contrasts) {
  K <- if (is.list(input$kinship) && !is.matrix(input$kinship) && !is.data.frame(input$kinship)) {
    input$kinship[[chr]]
  } else input$kinship
  ids <- unit$samples %||% qtl2mr_prepared_role_filter(input, "exposure")
  tryCatch(qtl2mr_recording({
    ctx <- qtl2mr_prepare_lmm_context_cached(pheno = qtl2mr_prepared_trait(input, unit$trait), kinship = K,
      addcovar = qtl2mr_prepared_covar(input), hsq = NULL, sample_ids = ids,
      standardize_pheno = settings$standardize_phenotypes, reml = TRUE, trait_name = unit$trait, chr = chr, tol = 1e-12)
    fit_founder_lmm_blocks_batch(input$allele_probs[[chr]], unit$markers[[chr]], ctx, contrasts,
      fallback_genoprobs = input$allele_probs, allow_marker_errors = TRUE, role = "exposure")
  }, muffle = TRUE), error = function(e) NULL)
  invisible(NULL)
}

# Chromosome-major order, so traits with the same samples on a chromosome reuse
# its kinship decomposition.
qtl2mr_batch_prefill_run <- function(input, units, settings) {
  contrasts <- make_founder_contrasts(dimnames(input$allele_probs[[1L]])[[2L]])
  chrs <- unique(unlist(lapply(units, function(u) names(u$markers)), use.names = FALSE))
  for (chr in chrs) for (u in units) {
    if (length(u$markers[[chr]])) qtl2mr_batch_prefill_chr(input, u, chr, settings, contrasts)
  }
  invisible(NULL)
}

qtl2mr_batch_prefill <- function(resources, plan, settings, n_workers) {
  units <- qtl2mr_batch_prefill_units(plan)
  if (!length(units)) return(invisible(NULL))
  markers <- unique(unlist(lapply(units, function(u) unlist(u$markers, use.names = FALSE)), use.names = FALSE))
  input <- tryCatch(qtl2mr_recording({
    qtl2mr_prepared_probability_copy(qtl2mr_resources_input(resources, markers), markers)
  }, muffle = TRUE)$value, error = function(e) NULL)
  if (is.null(input)) return(invisible(NULL))
  # Units with the same missing phenotypes are adjacent, and so stay together
  # in one worker.
  pattern <- vapply(units, function(u) {
    v <- resources$pheno[[u$trait]]
    paste(c(u$samples %||% "", which(!is.finite(as.numeric(v)))), collapse = ",")
  }, character(1))
  units <- units[order(pattern, names(units))]
  # Entries are pinned as they are made, so that filling a wave never releases
  # an entry the same wave needs.
  reg <- qtl2mr_registry()
  reg$pin_new <- TRUE
  on.exit(reg$pin_new <- NULL, add = TRUE)
  if (!qtl2mr_can_fork(n_workers, length(units))) return(qtl2mr_batch_prefill_run(input, units, settings))
  cost <- vapply(units, function(u) 50 * length(u$markers) + length(unlist(u$markers)), numeric(1))
  n <- min(n_workers, length(units))
  chunk <- pmin(n, 1L + floor(n * (cumsum(cost) - cost) / sum(cost)))
  chunks <- split(seq_along(units), chunk)
  out <- parallel::mclapply(chunks, function(idx) {
    reg <- qtl2mr_registry()
    reg$journal <- new.env(parent = emptyenv())
    counters0 <- qtl2mr_cache_counters()
    qtl2mr_batch_prefill_run(input, units[idx], settings)
    ids <- ls(reg$journal, all.names = TRUE)
    entries <- mget(ids, envir = reg$entries, ifnotfound = list(NULL))
    list(entries = entries[!vapply(entries, is.null, logical(1))], computed = qtl2mr_cache_counters() - counters0,
      pid = Sys.getpid())
  }, mc.cores = length(chunks), mc.preschedule = FALSE, mc.set.seed = FALSE)
  for (o in out) {
    if (is.list(o) && !inherits(o, "try-error") && is.list(o$entries)) {
      qtl2mr_registry_import(o$entries)
      # A single chunk runs in this process, whose counters already include it.
      if (!identical(o$pid, Sys.getpid())) qtl2mr_cache_counters_add_computed(o$computed)
    }
  }
  invisible(NULL)
}

# Registry ids of entries read by the given plan.
qtl2mr_batch_registry_ids <- function(items) {
  qtl2mr_registry_ids_matching(function(info) qtl2mr_registry_item(info) %in% items)
}
