# By-trait planning for run_qtl2mr_batch(). Pairs are fitted exactly as before,
# one at a time and in order; planning only computes beforehand, once per trait
# and sample set, the LMM contexts and locus fits that the pairs will read, so
# that each pair finds them in the session registry. A pair never depends on
# whether its fits were planned: a fit that is missing is computed by the pair,
# and a planned fit is stored under the same exact key the pair looks up.
#
# Pairs are processed in waves (consecutive pairs, or whole exposure groups when
# groups run in parallel workers) sized so that one wave's fits stay well within
# the registry budget. After a wave, its entries that no later wave reads are
# released.

# Per-pair demand: the two traits, their sample set (NULL: each trait's own
# complete samples) and the instruments.
qtl2mr_batch_plan <- function(resources, todo, instruments) {
  lapply(seq_len(nrow(todo)), function(i) {
    pair <- todo[i, , drop = FALSE]
    rows <- instruments[instruments$exposure == pair$exposure, , drop = FALSE]
    list(traits = c(pair$exposure, pair$outcome), samples = NULL, sig = "",
      marker = as.character(rows$marker), chr = as.character(rows$chr))
  })
}

qtl2mr_batch_fit_items <- function(plan) {
  unique(unlist(lapply(plan, function(p) {
    if (is.null(p)) return(character())
    as.vector(outer(paste(p$traits, p$sig, sep = "\r"), p$marker, paste, sep = "\r"))
  }), use.names = FALSE))
}

# Release keys of registry entries: trait and chromosome for contexts, plus the
# marker for fits (sample sets are not distinguished, so release is conservative).
qtl2mr_batch_release_items <- function(plan) {
  unique(unlist(lapply(plan, function(p) {
    if (is.null(p)) return(character())
    c(as.vector(outer(p$traits, unique(p$chr), paste, sep = "\r")),
      as.vector(outer(p$traits, paste(p$chr, p$marker, sep = "\r"), paste, sep = "\r")))
  }), use.names = FALSE))
}

qtl2mr_registry_item <- function(info) {
  if (is.null(info$trait) || is.null(info$chr)) return(NA_character_)
  if (is.null(info$marker) || is.na(info$marker)) paste(info$trait, info$chr, sep = "\r") else
    paste(info$trait, info$chr, info$marker, sep = "\r")
}

# Waves of whole units (lists of pair indices), in order. A unit always fits in
# a wave of its own.
qtl2mr_batch_waves <- function(plan, units, max_fits) {
  waves <- list()
  current <- integer()
  items <- character()
  for (u in units) {
    merged <- union(items, qtl2mr_batch_fit_items(plan[u]))
    if (length(current) && length(merged) > max_fits) {
      waves[[length(waves) + 1L]] <- current
      current <- u
      items <- qtl2mr_batch_fit_items(plan[u])
    } else {
      current <- c(current, u)
      items <- merged
    }
  }
  if (length(current)) waves[[length(waves) + 1L]] <- current
  waves
}

# Fits per wave: half the registry budget, at the size of a fit for these samples.
qtl2mr_batch_max_fits <- function(resources) {
  n <- max(length(resources$sample_ids), 1L)
  max(1L, as.integer(0.5 * qtl2mr_cache_budget() / (8 * n * 45 + 20000)))
}

# Work units of a wave: one per (trait, sample set), with the markers it needs
# on each chromosome.
qtl2mr_batch_prefill_units <- function(plan) {
  units <- list()
  for (p in plan) {
    if (is.null(p)) next
    for (trait in p$traits) {
      id <- paste(trait, p$sig, sep = "\r")
      u <- units[[id]] %||% list(trait = trait, samples = p$samples, markers = list())
      for (chr in unique(p$chr)) u$markers[[chr]] <- union(u$markers[[chr]], p$marker[p$chr == chr])
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
  if (!qtl2mr_can_fork(n_workers, length(units))) return(qtl2mr_batch_prefill_run(input, units, settings))
  cost <- vapply(units, function(u) 50 * length(u$markers) + length(unlist(u$markers)), numeric(1))
  n <- min(n_workers, length(units))
  chunk <- pmin(n, 1L + floor(n * (cumsum(cost) - cost) / sum(cost)))
  chunks <- split(seq_along(units), chunk)
  out <- parallel::mclapply(chunks, function(idx) {
    reg <- qtl2mr_registry()
    reg$journal <- new.env(parent = emptyenv())
    qtl2mr_batch_prefill_run(input, units[idx], settings)
    ids <- ls(reg$journal, all.names = TRUE)
    entries <- mget(ids, envir = reg$entries, ifnotfound = list(NULL))
    entries[!vapply(entries, is.null, logical(1))]
  }, mc.cores = length(chunks), mc.preschedule = FALSE, mc.set.seed = FALSE)
  for (entries in out) if (is.list(entries) && !inherits(entries, "try-error")) qtl2mr_registry_import(entries)
  invisible(NULL)
}

# Registry ids of entries read by the given plan.
qtl2mr_batch_registry_ids <- function(items) {
  qtl2mr_registry_ids_matching(function(info) qtl2mr_registry_item(info) %in% items)
}
