qtl2mr_marker_map_from_genoprobs <- function(allele_probs) {
  rows <- lapply(names(allele_probs), function(chr) {
    markers <- dimnames(allele_probs[[chr]])[[3L]]
    data.frame(marker = markers, chr = as.character(chr), pos = seq_along(markers),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

qtl2mr_prepare_resources <- function(allele_probs,
                                     kinship,
                                     pheno,
                                     covar,
                                     marker_map = NULL,
                                     metadata = list(),
                                     probability_policy = c("normalize_rounding", "strict")) {
  probability_policy <- match.arg(probability_policy)
  if (!is.list(allele_probs) || !length(allele_probs)) {
    stop("allele_probs must be a non-empty qtl2-style list", call. = FALSE)
  }
  # Before 0.6.0 the fifth argument was the QTL scan.
  if (is.matrix(marker_map) && is.numeric(marker_map)) {
    stop("qtl_scan is no longer part of resources; pass the scan to qtl2mr_select_instruments()", call. = FALSE)
  }
  if (is.null(marker_map)) marker_map <- qtl2mr_marker_map_from_genoprobs(allele_probs)
  qtl2mr_validate_ids(rownames(pheno), "pheno row names")
  if (!is.null(covar)) qtl2mr_validate_ids(rownames(covar), "covar row names")
  pheno <- as.data.frame(pheno)
  covar <- if (is.null(covar)) NULL else as.data.frame(covar)
  marker_map <- qtl2mr_marker_map_table(marker_map)
  if (is.null(rownames(pheno)) || any(!nzchar(rownames(pheno)))) {
    stop("pheno must have sample IDs as row names", call. = FALSE)
  }
  if (!is.null(covar) && (is.null(rownames(covar)) || any(!nzchar(rownames(covar))))) {
    stop("covar must have sample IDs as row names", call. = FALSE)
  }
  if (!all(c("marker", "chr") %in% names(marker_map))) {
    stop("marker_map must contain marker and chr columns", call. = FALSE)
  }
  gp_ids <- dimnames(allele_probs[[1L]])[[1L]]
  if (is.null(gp_ids)) stop("allele_probs arrays must have sample row names", call. = FALSE)
  marker_map$chr <- as.character(marker_map$chr)
  marker_map$marker <- as.character(marker_map$marker)
  out <- list(
    allele_probs = allele_probs,
    kinship = kinship,
    pheno = pheno,
    covar = covar,
    marker_map = marker_map,
    sample_ids = gp_ids,
    metadata = metadata,
    probability_policy = probability_policy
  )
  class(out) <- c("qtl2mr_resources", "list")
  validate_qtl2mr_input(out)
  out <- qtl2mr_record_probability_policy(out, warn = TRUE)
  out$sample_alignment <- qtl2mr_input_alignment(out)
  qtl2mr_warn_alignment(out$sample_alignment)
  out
}

qtl2mr_validate_lod_thresh <- function(lod_thresh) {
  if (!is.numeric(lod_thresh) || is.logical(lod_thresh) || length(lod_thresh) != 1L ||
      !is.finite(lod_thresh) || lod_thresh <= 0) {
    stop("lod_thresh must be one finite numeric value > 0", call. = FALSE)
  }
  as.numeric(lod_thresh)
}

qtl2mr_check_count <- function(x, name, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  as.integer(x)
}

# A scan identifies its traits through list names or column names; a plain LOD
# vector or an unnamed column does not.
qtl2mr_scan_has_traits <- function(qtl_scan) {
  if (is.list(qtl_scan) && !is.data.frame(qtl_scan)) return(TRUE)
  traits <- colnames(qtl_scan)
  !is.null(traits) && !anyNA(traits) && all(nzchar(traits))
}

qtl2mr_scan_traits <- function(qtl_scan) {
  traits <- if (is.list(qtl_scan) && !is.data.frame(qtl_scan)) names(qtl_scan) else colnames(qtl_scan)
  if (is.null(traits) || !length(traits) || anyNA(traits) || any(!nzchar(traits))) {
    stop("qtl_scan has no trait names; pass exposures explicitly", call. = FALSE)
  }
  traits
}

qtl2mr_lod_vector_from_scan <- function(qtl_scan, trait) {
  if (length(trait) != 1L || is.na(trait) || !nzchar(trait)) {
    stop("trait must be one non-empty trait name", call. = FALSE)
  }
  scan <- qtl_scan
  trait_validated_by_list <- FALSE
  if (is.list(qtl_scan) && !is.data.frame(qtl_scan) && !is.matrix(qtl_scan)) {
    scan_names <- names(qtl_scan)
    if (is.null(scan_names) || any(is.na(scan_names) | !nzchar(scan_names))) {
      stop("named-list qtl_scan must use non-empty trait names", call. = FALSE)
    }
    idx <- which(scan_names == trait)
    if (!length(idx)) stop("qtl_scan is missing trait: ", trait, call. = FALSE)
    if (length(idx) > 1L) stop("qtl_scan has duplicated trait: ", trait, call. = FALSE)
    scan <- qtl_scan[[idx[[1L]]]]
    trait_validated_by_list <- TRUE
  }
  if (is.data.frame(scan) || is.matrix(scan)) {
    scan <- as.matrix(scan)
    if (!is.numeric(scan)) stop("qtl_scan LOD values must be numeric", call. = FALSE)
    marker_ids <- rownames(scan)
    if (is.null(marker_ids) || any(is.na(marker_ids) | !nzchar(marker_ids)) || anyDuplicated(marker_ids)) {
      stop("qtl_scan must have non-empty, unique marker row names", call. = FALSE)
    }
    n_traits <- ncol(scan)
    trait_names <- colnames(scan)
    if (n_traits > 1L) {
      if (is.null(trait_names) || any(is.na(trait_names) | !nzchar(trait_names))) {
        stop("multi-column qtl_scan must have non-empty trait column names", call. = FALSE)
      }
      if (anyDuplicated(trait_names)) {
        stop("multi-column qtl_scan has duplicated trait column names", call. = FALSE)
      }
      idx <- match(trait, trait_names)
      if (is.na(idx)) stop("qtl_scan is missing trait: ", trait, call. = FALSE)
      lod <- as.numeric(scan[, idx])
    } else {
      if (!is.null(trait_names) && length(trait_names) && !is.na(trait_names[[1L]]) && nzchar(trait_names[[1L]])) {
        if (!isTRUE(trait_validated_by_list) && !identical(trait_names[[1L]], trait)) {
          stop("single-column qtl_scan trait name does not match the requested trait: ", trait, call. = FALSE)
        }
      } else {
        warning(
          "single-column qtl_scan has no trait column name; the trait could not be validated",
          call. = FALSE
        )
      }
      lod <- as.numeric(scan[, 1L])
    }
    names(lod) <- marker_ids
  } else {
    if (!is.numeric(scan)) stop("qtl_scan LOD values must be numeric", call. = FALSE)
    lod <- as.numeric(scan)
    names(lod) <- names(scan)
  }
  marker_ids <- names(lod)
  if (is.null(marker_ids) || any(is.na(marker_ids) | !nzchar(marker_ids)) || anyDuplicated(marker_ids)) {
    stop("qtl_scan LOD vector must have non-empty, unique marker names", call. = FALSE)
  }
  lod
}

# LOD peaks of one trait: the maximum of each contiguous above-threshold region,
# defined in map order (positions, or map rows when positions are missing) within
# a chromosome, never across a chromosome boundary.
qtl2mr_scan_candidates <- function(qtl_scan, marker_map, trait, lod_thresh, exclude_chr) {
  lod <- qtl2mr_lod_vector_from_scan(qtl_scan, trait)
  mm <- marker_map
  map_idx <- match(names(lod), mm$marker)
  if (anyNA(map_idx)) stop("marker_map is missing qtl_scan markers", call. = FALSE)
  peak_idx <- unlist(lapply(unique(mm$chr[map_idx]), function(chr) {
    idx <- which(mm$chr[map_idx] == chr)
    pos <- mm$pos[map_idx[idx]]
    idx <- if (all(is.finite(pos))) idx[order(pos, map_idx[idx])] else idx[order(map_idx[idx])]
    idx[qtl2mr_get_tophits(lod[idx], lod_thresh = lod_thresh)]
  }), use.names = FALSE)
  idx <- match(names(lod)[peak_idx], mm$marker)
  out <- data.frame(exposure = rep(trait, length(peak_idx)), marker = names(lod)[peak_idx],
    chr = as.character(mm$chr[idx]), pos = mm$pos[idx], lod = as.numeric(lod[peak_idx]),
    stringsAsFactors = FALSE)
  out <- out[!out$chr %in% exclude_chr, , drop = FALSE]
  out <- out[order(out$lod, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

qtl2mr_is_autosome <- function(chr) {
  !toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE)) %in% c("X", "Y", "M", "MT")
}

# Prepared input built from resources, shared by instrument selection, batch
# pairs and qtl2mr_pair_input() so that every route fits the same samples.
qtl2mr_resources_input <- function(resources, selected_loci, exposure = NULL, outcome = NULL,
                                   metadata = resources$metadata) {
  as_qtl2mr_prepared_input(
    allele_probs = resources$allele_probs, pheno = resources$pheno, covar = resources$covar,
    kinship = resources$kinship, marker_map = resources$marker_map, selected_loci = selected_loci,
    exposure = exposure, outcome = outcome, metadata = metadata %||% list(),
    probability_policy = resources$probability_policy %||% "strict"
  )
}

# Exposure founder LMM at each candidate locus, fitted exactly as in the MR stage
# (same samples, covariates, LOCO kinship and null-model hsq), so F_X_block here
# equals the block F_X of any pair using this exposure; with
# standardize_phenotypes = FALSE in the MR stage they agree up to rounding.
qtl2mr_exposure_first_stage <- function(input, exposure, candidates, ridge = 1e-10) {
  out <- qtl2mr_empty_first_stage(nrow(candidates))
  if (!nrow(candidates)) return(out)
  x <- qtl2mr_as_named_numeric(qtl2mr_prepared_trait(input, exposure), "exposure")
  covar <- qtl2mr_prepared_covar(input)
  ids <- qtl2mr_prepared_role_filter(input, "exposure")
  contrasts <- make_founder_contrasts(dimnames(input$allele_probs[[1L]])[[2L]])
  kinship <- input$kinship
  for (chr in unique(candidates$chr)) {
    rows <- which(candidates$chr == chr)
    K <- if (is.list(kinship) && !is.matrix(kinship) && !is.data.frame(kinship)) kinship[[chr]] else kinship
    fits <- tryCatch({
      ctx <- qtl2mr_prepare_lmm_context_cached(pheno = x, kinship = K, addcovar = covar,
        hsq = qtl2mr_prepared_hsq(input, chr, "exposure", exposure), sample_ids = ids,
        standardize_pheno = TRUE, reml = TRUE, trait_name = exposure, chr = chr)
      fit_founder_lmm_blocks_batch(input$allele_probs[[chr]], candidates$marker[rows], ctx, contrasts,
        fallback_genoprobs = input$allele_probs, allow_marker_errors = TRUE, role = "exposure")
    }, error = function(e) e)
    for (j in seq_along(rows)) {
      fit <- if (inherits(fits, "error")) {
        qtl2mr_lmm_marker_error(candidates$marker[[rows[[j]]]], chr, "exposure", fits)
      } else fits[[j]]
      if (qtl2mr_is_lmm_marker_error(fit)) {
        out$fit_error[[rows[[j]]]] <- if (identical(qtl2mr_marker_drop_reason(fit), "nonestimable_founder_contrast")) {
          "nonestimable_x"
        } else "x_fit_error"
        next
      }
      out$F_X_block[[rows[[j]]]] <- mr_block_first_stage_f(list(beta_x = fit$beta_fc,
        Sigma_x = qtl2mr_sym(fit$vcov_fc)), ridge = ridge)
      out$n_x[[rows[[j]]]] <- as.integer(fit$n)
      out$hsq_x[[rows[[j]]]] <- as.numeric(fit$hsq)
    }
  }
  out
}

qtl2mr_empty_first_stage <- function(n) {
  data.frame(F_X_block = rep(NA_real_, n), n_x = rep(NA_integer_, n), hsq_x = rep(NA_real_, n),
    fit_error = rep("", n), stringsAsFactors = FALSE)
}

qtl2mr_empty_instruments <- function() {
  data.frame(exposure = character(), marker = character(), chr = character(), pos = numeric(),
    lod = numeric(), F_X_block = numeric(), n_x = integer(), hsq_x = numeric(), selected = logical(),
    drop_reason = character(), locus_rank = integer(), selected_by = character(),
    stringsAsFactors = FALSE)
}

# Candidates arrive sorted by LOD. Non-autosomal and non-estimable loci can never
# be used; the F rule and then max_loci decide among the rest.
qtl2mr_instrument_rows <- function(candidates, first_stage, F_X_min, max_loci) {
  if (!nrow(candidates)) return(qtl2mr_empty_instruments())
  reason <- rep("", nrow(candidates))
  reason[!qtl2mr_is_autosome(candidates$chr)] <- "non_autosomal"
  open <- !nzchar(reason) & nzchar(first_stage$fit_error)
  reason[open] <- first_stage$fit_error[open]
  if (!is.null(F_X_min)) {
    weak <- !nzchar(reason) & !(is.finite(first_stage$F_X_block) & first_stage$F_X_block >= F_X_min)
    reason[weak] <- "F"
  }
  passing <- which(!nzchar(reason))
  if (!is.null(max_loci) && length(passing) > max_loci) reason[passing[-seq_len(max_loci)]] <- "max_loci"
  selected <- !nzchar(reason)
  data.frame(candidates, F_X_block = first_stage$F_X_block, n_x = first_stage$n_x,
    hsq_x = first_stage$hsq_x, selected = selected, drop_reason = reason,
    locus_rank = ifelse(selected, cumsum(selected), NA_integer_),
    selected_by = rep("contiguous_lod_island_max", nrow(candidates)), stringsAsFactors = FALSE)
}

qtl2mr_select_instruments <- function(resources,
                                      qtl_scan,
                                      exposures = NULL,
                                      lod_thresh = 4,
                                      exclude_chr = "X",
                                      max_loci = NULL,
                                      F_X_min = NULL,
                                      n_workers = 1L,
                                      progress = TRUE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  lod_thresh <- qtl2mr_validate_lod_thresh(lod_thresh)
  max_loci <- qtl2mr_check_count(max_loci, "max_loci", allow_null = TRUE)
  if (!is.null(F_X_min) && (!is.numeric(F_X_min) || length(F_X_min) != 1L || !is.finite(F_X_min) || F_X_min < 0)) {
    stop("F_X_min must be NULL or one non-negative number", call. = FALSE)
  }
  n_workers <- qtl2mr_check_count(n_workers, "n_workers")
  progress <- qtl2mr_check_flag(progress, "progress")
  exclude_chr <- as.character(exclude_chr %||% character())
  exposures <- if (is.null(exposures)) qtl2mr_scan_traits(qtl_scan) else as.character(exposures)
  if (!length(exposures) || anyNA(exposures) || any(!nzchar(exposures)) || anyDuplicated(exposures)) {
    stop("exposures must be unique, non-empty trait names", call. = FALSE)
  }
  missing <- setdiff(exposures, names(resources$pheno))
  if (length(missing)) stop("pheno is missing exposure(s): ", paste(missing, collapse = ", "), call. = FALSE)
  if (length(exposures) > 1L && !qtl2mr_scan_has_traits(qtl_scan)) {
    stop("qtl_scan has no trait names, so it cannot be matched to several exposures; ",
      "use a scan with trait column names or one exposure per call", call. = FALSE)
  }
  candidates <- lapply(exposures, function(trait) {
    qtl2mr_scan_candidates(qtl_scan, resources$marker_map, trait, lod_thresh, exclude_chr)
  })
  markers <- unique(unlist(lapply(candidates, function(d) d$marker[qtl2mr_is_autosome(d$chr)]), use.names = FALSE))
  input <- NULL
  if (length(markers)) {
    input <- qtl2mr_resources_input(resources, markers)
    input <- qtl2mr_prepared_probability_copy(input, markers)
  }
  select_one <- function(i) {
    cand <- candidates[[i]]
    if (progress) message(sprintf("[qtl2mr select] %d/%d %s: %d candidate loci", i, length(exposures),
      exposures[[i]], nrow(cand)))
    first_stage <- qtl2mr_empty_first_stage(nrow(cand))
    auto <- qtl2mr_is_autosome(cand$chr)
    if (any(auto)) {
      first_stage[auto, ] <- qtl2mr_exposure_first_stage(input, exposures[[i]], cand[auto, , drop = FALSE])
    }
    qtl2mr_instrument_rows(cand, first_stage, F_X_min, max_loci)
  }
  index <- seq_along(exposures)
  if (n_workers > 1L && length(exposures) > 1L && identical(.Platform$OS.type, "unix")) {
    chunks <- split(index, (index - 1L) %% min(n_workers, length(exposures)))
    out <- parallel::mclapply(chunks, function(idx) lapply(idx, select_one),
      mc.cores = length(chunks), mc.preschedule = FALSE)
    failed <- !vapply(out, is.list, logical(1))
    if (any(failed)) {
      stop("qtl2mr_select_instruments() worker failed: ", as.character(out[failed][[1L]]), call. = FALSE)
    }
    rows <- unlist(out, recursive = FALSE)[order(unlist(chunks, use.names = FALSE))]
  } else {
    rows <- lapply(index, select_one)
  }
  out <- do.call(rbind, c(list(qtl2mr_empty_instruments()), rows))
  out$lod_thresh <- rep(lod_thresh, nrow(out))
  out$F_X_min <- rep(F_X_min %||% NA_real_, nrow(out))
  out$max_loci <- rep(max_loci %||% NA_integer_, nrow(out))
  out$exclude_chr <- rep(paste(exclude_chr, collapse = ","), nrow(out))
  rownames(out) <- NULL
  out
}

# Instrument tables are keyed by exposure. Rows with selected = FALSE are kept
# only as a record; markers must exist in marker_map and allele_probs.
qtl2mr_normalize_instruments <- function(instruments, resources) {
  if (is.character(instruments) && length(instruments) == 1L && !is.na(instruments)) {
    if (!file.exists(instruments) || dir.exists(instruments)) {
      stop("instruments file does not exist: ", instruments, call. = FALSE)
    }
    instruments <- utils::read.csv(instruments, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(instruments)) {
    stop("instruments must be a data frame or CSV path with exposure and marker columns, for example from qtl2mr_select_instruments()", call. = FALSE)
  }
  missing <- setdiff(c("exposure", "marker"), names(instruments))
  if (length(missing)) stop("instruments must contain column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  if ("selected" %in% names(instruments)) {
    keep <- instruments$selected
    if (!is.logical(keep) || anyNA(keep)) stop("instruments$selected must be TRUE or FALSE", call. = FALSE)
    instruments <- instruments[keep, , drop = FALSE]
  }
  exposure <- as.character(instruments$exposure)
  marker <- as.character(instruments$marker)
  if (anyNA(exposure) || anyNA(marker) || any(!nzchar(exposure)) || any(!nzchar(marker))) {
    stop("instruments contains missing or blank exposure or marker values", call. = FALSE)
  }
  dup <- duplicated(paste(exposure, marker, sep = "\r"))
  if (any(dup)) {
    stop("instruments lists marker ", marker[dup][[1L]], " more than once for exposure ", exposure[dup][[1L]], call. = FALSE)
  }
  map <- resources$marker_map
  idx <- match(marker, map$marker)
  if (anyNA(idx)) {
    stop("instruments markers are not in marker_map: ", paste(utils::head(unique(marker[is.na(idx)]), 5L), collapse = ", "), call. = FALSE)
  }
  chr <- as.character(map$chr[idx])
  if ("chr" %in% names(instruments)) {
    given <- as.character(instruments$chr)
    bad <- !is.na(given) & given != chr
    if (any(bad)) {
      stop("instruments chr does not match marker_map for marker(s): ", paste(utils::head(unique(marker[bad]), 5L), collapse = ", "), call. = FALSE)
    }
  }
  in_probs <- mapply(function(m, c) {
    a <- resources$allele_probs[[c]]
    !is.null(a) && m %in% dimnames(a)[[3L]]
  }, marker, chr)
  if (!all(in_probs)) {
    stop("instruments markers are not in allele_probs: ", paste(utils::head(unique(marker[!in_probs]), 5L), collapse = ", "), call. = FALSE)
  }
  data.frame(exposure = exposure, marker = marker, chr = chr, pos = map$pos[idx], stringsAsFactors = FALSE)
}

qtl2mr_exposure_instruments <- function(instruments, exposure) {
  rows <- instruments[instruments$exposure == exposure, , drop = FALSE]
  if (!nrow(rows)) stop("no_instruments: instruments has no selected loci for exposure ", exposure, call. = FALSE)
  rows
}

qtl2mr_pair_prepared <- function(resources, rows, exposure, outcome, metadata = resources$metadata) {
  qtl2mr_resources_input(resources, data.frame(marker = rows$marker, chr = rows$chr, stringsAsFactors = FALSE),
    exposure = exposure, outcome = outcome, metadata = metadata)
}

qtl2mr_pair_input <- function(resources, instruments, exposure, outcome) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  for (name in c("exposure", "outcome")) {
    if (!qtl2mr_nonempty_string(get(name))) stop(name, " must be one non-empty trait name", call. = FALSE)
  }
  rows <- qtl2mr_exposure_instruments(qtl2mr_normalize_instruments(instruments, resources), exposure)
  qtl2mr_pair_prepared(resources, rows, exposure, outcome)
}

qtl2mr_resources_kinship <- function(kinship, chr) {
  if (is.null(kinship)) return(NULL)
  if (is.list(kinship) && !is.matrix(kinship) && !is.data.frame(kinship)) {
    K <- kinship[[as.character(chr)]]
    if (is.null(K)) stop("kinship list is missing selected chromosome: ", chr, call. = FALSE)
    return(as.matrix(K))
  }
  as.matrix(kinship)
}

# Samples with finite values for every listed trait, complete covariates, and
# genotype probabilities and kinship rows on every selected chromosome.
qtl2mr_complete_samples <- function(resources, traits, selected_chr = NULL) {
  pheno <- resources$pheno
  covar <- resources$covar
  ids <- Reduce(intersect, list(resources$sample_ids, rownames(pheno),
    if (is.null(covar)) resources$sample_ids else rownames(covar)))
  for (chr in unique(as.character(selected_chr))) {
    gp <- resources$allele_probs[[chr]]
    if (is.null(gp)) stop("allele_probs is missing selected chromosome: ", chr, call. = FALSE)
    ids <- ids[ids %in% dimnames(gp)[[1L]]]
    K <- qtl2mr_resources_kinship(resources$kinship, chr)
    if (!is.null(K)) ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
  }
  keep <- qtl2mr_complete_rows(pheno[ids, traits, drop = FALSE])
  if (!is.null(covar) && ncol(covar)) keep <- keep & qtl2mr_complete_rows(covar[ids, , drop = FALSE])
  ids[keep]
}

# One row per trait and exclusion reason. Each trait is fitted on its own
# complete samples; reasons can overlap.
qtl2mr_batch_sample_exclusions <- function(resources, pairs, instruments) {
  if (!nrow(pairs)) return(qtl2mr_empty_batch_table("sample_exclusions"))
  pheno <- resources$pheno
  covar <- resources$covar
  genotyped <- resources$sample_ids
  ids <- union(genotyped, rownames(pheno))
  in_pheno <- ids %in% rownames(pheno)
  missing_covar <- if (is.null(covar)) rep(FALSE, length(ids)) else !ids %in% rownames(covar)
  incomplete_covar <- rep(FALSE, length(ids))
  if (!is.null(covar) && ncol(covar)) {
    matched <- which(!missing_covar)
    incomplete_covar[matched] <- !qtl2mr_complete_rows(covar[ids[matched], , drop = FALSE])
  }
  traits <- unique(c(pairs$exposure, pairs$outcome))
  do.call(rbind, lapply(traits, function(trait) {
    value <- rep(NA_real_, length(ids))
    value[in_pheno] <- as.numeric(pheno[ids[in_pheno], trait])
    exposures <- pairs$exposure[pairs$exposure == trait | pairs$outcome == trait]
    missing_kinship <- rep(FALSE, length(ids))
    if (!is.null(resources$kinship)) {
      for (chr in unique(instruments$chr[instruments$exposure %in% exposures])) {
        K <- qtl2mr_resources_kinship(resources$kinship, chr)
        missing_kinship <- missing_kinship | !(ids %in% rownames(K) & ids %in% colnames(K))
      }
    }
    reasons <- list(
      missing_genoprob_id = !ids %in% genotyped,
      missing_phenotype_id = !in_pheno,
      nonfinite_value = in_pheno & !is.finite(value),
      missing_covar_id = missing_covar,
      incomplete_covar = incomplete_covar,
      missing_kinship_id = missing_kinship
    )
    do.call(rbind, lapply(names(reasons), function(reason) data.frame(
      trait = trait, reason = reason, n_input_samples = length(ids),
      n_excluded = sum(reasons[[reason]]),
      excluded_ids = paste(ids[reasons[[reason]]], collapse = ";"),
      stringsAsFactors = FALSE
    )))
  }))
}

qtl2mr_pair_id <- function(exposure, outcome) {
  paste0(make.names(exposure), "__to__", make.names(outcome))
}

qtl2mr_normalize_pairs <- function(pairs) {
  if (is.character(pairs) && length(pairs) == 1L) {
    if (is.na(pairs) || !file.exists(pairs) || dir.exists(pairs)) stop("pairs file does not exist: ", pairs, call. = FALSE)
    pairs <- utils::read.csv(pairs, stringsAsFactors = FALSE)
  }
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("exposure", "outcome") %in% names(pairs))) {
    stop("pairs must contain exposure and outcome columns", call. = FALSE)
  }
  retired <- intersect(c("instrument_trait", "lod_thresh", "exclude_chr", "max_loci", "selected_loci_file"), names(pairs))
  if (length(retired)) {
    stop("pairs column(s) ", paste(retired, collapse = ", "), " were removed in qtl2mr 0.6.0. ",
      "Select instruments with qtl2mr_select_instruments() and pass the table as `instruments`.", call. = FALSE)
  }
  if (!"pair_id" %in% names(pairs)) {
    pairs$pair_id <- mapply(qtl2mr_pair_id, pairs$exposure, pairs$outcome, USE.NAMES = FALSE)
  }
  if (!nrow(pairs)) stop("pairs must contain at least one pair", call. = FALSE)
  for (nm in c("pair_id", "exposure", "outcome")) {
    pairs[[nm]] <- as.character(pairs[[nm]])
    if (anyNA(pairs[[nm]]) || any(!nzchar(pairs[[nm]]))) stop("pairs contains missing or blank ", nm, call. = FALSE)
  }
  if (anyDuplicated(pairs$pair_id)) stop("pair_id must be unique; supply explicit IDs for distinct requests", call. = FALSE)
  pairs
}

qtl2mr_nonempty_string <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

qtl2mr_batch_probability_table <- function(pair, pilot) {
  adjustments <- pilot$prepared_metadata$probability_normalization$adjustments
  if (is.null(adjustments) || !nrow(adjustments)) return(qtl2mr_empty_batch_table("probability_normalization"))
  cbind(pair_id = pair$pair_id, adjustments, stringsAsFactors = FALSE)
}

qtl2mr_batch_preflight <- function(resources, pairs, instruments, min_loci, min_samples, fail_fast) {
  errors <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    tryCatch({
      rows <- qtl2mr_exposure_instruments(instruments, pair$exposure)
      if (nrow(rows) < min_loci) {
        stop(sprintf("too_few_loci: %d instruments for exposure %s, below min_loci = %d",
          nrow(rows), pair$exposure, min_loci), call. = FALSE)
      }
      probe <- resources
      probe$selected_loci <- rows$marker
      validate_qtl2mr_input(probe, exposure = pair$exposure, outcome = pair$outcome, require_selected_loci = TRUE)
      n <- vapply(c(pair$exposure, pair$outcome), function(trait) {
        length(qtl2mr_complete_samples(resources, trait, rows$chr))
      }, integer(1))
      if (any(n < min_samples)) {
        stop(sprintf("too_few_samples: %d complete samples for the exposure and %d for the outcome, below min_samples = %d",
          n[[1L]], n[[2L]], min_samples), call. = FALSE)
      }
      rows$marker
    }, error = function(e) {
      if (isTRUE(fail_fast)) stop(e)
      qtl2mr_error_row(pair, "preflight", e)
    })
  })
  invalid <- vapply(errors, is.data.frame, logical(1))
  markers <- unique(unlist(errors[!invalid], use.names = FALSE))
  list(pairs = pairs, valid = pairs[!invalid, , drop = FALSE], errors = do.call(rbind, errors[invalid]),
    markers = markers)
}

qtl2mr_append_preflight_errors <- function(batch, preflight) {
  errors <- preflight$errors
  if (!is.null(errors) && nrow(errors)) {
    for (nm in c("main_results", "errors")) {
      if ("pair_id" %in% names(batch[[nm]])) batch[[nm]] <- batch[[nm]][!batch[[nm]]$pair_id %in% errors$pair_id, , drop = FALSE]
    }
    pairs <- preflight$pairs[match(errors$pair_id, preflight$pairs$pair_id), , drop = FALSE]
    batch <- qtl2mr_batch_append(batch, "main_results", qtl2mr_batch_error_main(pairs, errors$error_message))
    batch <- qtl2mr_batch_append(batch, "errors", errors)
  }
  batch$pairs <- preflight$pairs
  batch
}

qtl2mr_fingerprint_object <- function(object) {
  path <- tempfile("qtl2mr-fingerprint-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, compress = FALSE, version = 2L)
  unname(tools::md5sum(path))
}

qtl2mr_batch_fingerprint <- function(resources, pairs, instruments, config, markers = NULL) {
  config$created_at <- NULL
  if (!is.null(markers)) {
    resources$allele_probs <- lapply(resources$allele_probs, function(a) {
      a[, , dimnames(a)[[3L]] %in% markers, drop = FALSE]
    })
  }
  instruments <- instruments[instruments$exposure %in% pairs$exposure, , drop = FALSE]
  rownames(instruments) <- NULL
  env <- environment(qtl2mr_batch_fingerprint)
  fn <- sort(ls(env, all.names = TRUE))
  fn <- fn[vapply(fn, function(nm) is.function(get(nm, env, inherits = FALSE)), logical(1))]
  # Deparsed text, not language objects: byte compilation and source references
  # can change a closure's serialized form without changing its code.
  implementation <- lapply(fn, function(nm) {
    c(nm, deparse(get(nm, env, inherits = FALSE), control = c("keepInteger", "keepNA")))
  })
  list(schema = "qtl2mr-batch-v060-1", data = qtl2mr_fingerprint_object(resources),
       request = qtl2mr_fingerprint_object(list(pairs, instruments, config)),
       implementation = qtl2mr_fingerprint_object(implementation))
}

qtl2mr_resume_abort <- function(message) {
  stop(structure(list(message = message, call = NULL),
    class = c("qtl2mr_resume_error", "error", "condition")))
}

# Paths written by run_qtl2mr_batch(). Names used only by earlier versions are
# kept so overwrite = TRUE can clear old folders.
qtl2mr_batch_owned_paths <- function(output_dir) {
  if (is.null(output_dir)) return(character())
  tables <- c("main_results", "method_comparison", "locus_results", "overlap_qc", "warnings", "errors",
    "run_manifest", "sample_exclusions", "probability_normalization", "input_alignment",
    "selected_loci", "first_stage_summary", "sample_qc", "first_stage_filter_results",
    "mr_block_cache_manifest", "exact_group_plan", "sample_groups", "selected_loci_summary", "workload_summary")
  paths <- file.path(output_dir, c(paste0("batch_", tables, ".csv"),
    "batch_config.json", "batch_resume_state.rds", "batch_pair_records", "mr_block_cache", "by_exposure"))
  links <- Sys.readlink(paths)
  paths[file.exists(paths) | (!is.na(links) & nzchar(links))]
}

# Batch writes follow links, so a linked table or checkpoint could redirect a
# write outside output_dir. Such folders are refused before anything is written.
qtl2mr_batch_linked_paths <- function(output_dir) {
  paths <- qtl2mr_batch_owned_paths(output_dir)
  records <- qtl2mr_batch_records_dir(output_dir)
  if (dir.exists(records)) paths <- c(paths, list.files(records, all.files = TRUE, no.. = TRUE, full.names = TRUE))
  links <- Sys.readlink(paths)
  paths[!is.na(links) & nzchar(links)]
}

qtl2mr_check_batch_destination <- function(output_dir, resume = FALSE, overwrite = FALSE) {
  if (is.null(output_dir)) return(invisible(TRUE))
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be NULL or one directory path", call. = FALSE)
  }
  if (file.exists(output_dir) && !dir.exists(output_dir)) stop("output_dir exists and is not a directory", call. = FALSE)
  linked <- qtl2mr_batch_linked_paths(output_dir)
  if (length(linked)) {
    stop("Batch output contains symbolic links (", paste(basename(linked), collapse = ", "),
      "); qtl2mr does not write through links. Remove them or use a new output_dir.", call. = FALSE)
  }
  occupied <- length(qtl2mr_batch_owned_paths(output_dir)) > 0L
  if (occupied && !isTRUE(resume) && !isTRUE(overwrite)) {
    stop("Batch output already exists; use a new output_dir, resume = TRUE for the identical run, or overwrite = TRUE to replace it", call. = FALSE)
  }
  invisible(TRUE)
}

qtl2mr_empty_dataframe <- function(columns) {
  out <- stats::setNames(lapply(columns, function(x) logical()), columns)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# Main-result columns shared with the single-pair main_result table: same names,
# same values for the same pair, instruments and settings.
qtl2mr_batch_summary_columns <- c(
  "exposure", "outcome", "lmm_backend", "n_loci", "n_loci_selected", "n_loci_dropped_nonestimable",
  "n_loci_steiger_fail", "n_x_median", "n_y_median", "n_overlap_median",
  "theta_reported", "se_reported", "reported_ci_lower", "reported_ci_upper", "reported_ci_type",
  "reported_estimator", "reported_estimator_source", "reported_estimator_policy", "reported_estimator_reason",
  "theta", "se", "profile_ci_lower", "profile_ci_upper", "profile_p",
  "confidence_set_type", "confidence_set_bounded", "primary_optimizer_status", "primary_bracket_status",
  "primary_bracket_n_expand", "primary_profile_ci_finite", "primary_reliable", "numerical_ok",
  "inference_status", "effect_units", "exposure_sd", "outcome_sd", "raw_conversion_available",
  "theta_raw_units", "se_raw_units", "first_order_theta", "first_order_se",
  "theta_independence", "theta_overlap_delta", "theta_two_step_ivw", "theta_weighted_median",
  "theta_weighted_mode", "Q", "Q_chisq_df", "Q_p_chisq", "Q_p_reported", "Q_p_method",
  "cross_cov_model", "standardize_phenotypes"
)

qtl2mr_empty_batch_table <- function(name) {
  columns <- switch(
    name,
    main_results = c("pair_id", "exposure", "outcome", "status", "error_message",
      setdiff(qtl2mr_batch_summary_columns, c("exposure", "outcome")), "locus_fit_status", "runtime_sec"),
    method_comparison = c(
      "pair_id", "exposure", "outcome", "method", "estimate", "se", "ci_lower", "ci_upper", "ci_type",
      "source", "note", "ci_level", "bootstrap_B", "bootstrap_seed", "bootstrap_complete_B",
      "bootstrap_failed_B"
    ),
    locus_results = c(
      "pair_id", "exposure", "outcome", "marker", "chr", "pos", "kept", "drop_reason",
      "F_X_block", "F_Y_block", "partial_R2_X_block", "partial_R2_Y_block", "steiger_pass",
      "n_x", "n_y", "hsq_x", "hsq_y"
    ),
    overlap_qc = c(
      "pair_id", "locus_index", "marker", "chr", "n_x", "n_y", "n_overlap",
      "exposure_overlap_fraction", "outcome_overlap_fraction", "jaccard_overlap_fraction",
      "complete_overlap", "independent_samples"
    ),
    warnings = c(
      "pair_id", "warning_id", "severity", "message", "consequence",
      "recommended_action", "shown_in_main_report"
    ),
    errors = c(
      "pair_id", "exposure", "outcome", "stage", "error_class",
      "error_message", "traceback_available", "runtime_sec"
    ),
    run_manifest = c(
      "run_id", "created_at", "n_pairs_requested", "n_pairs_ok", "n_pairs_error",
      "n_pairs_resumed", "nonestimable_locus_policy", "n_loci_dropped_nonestimable",
      "parallel_backend", "n_workers", "context_cache_hits", "context_cache_misses",
      "lmm_fit_cache_hits", "lmm_fit_cache_misses", "wall_sec", "mean_wall_sec_per_pair", "lmm_backend"
    ),
    sample_exclusions = c("trait", "reason", "n_input_samples", "n_excluded", "excluded_ids"),
    probability_normalization = c("pair_id", "chr", "marker", "n_normalized", "max_abs_sum_deviation"),
    character()
  )
  qtl2mr_empty_dataframe(columns)
}

qtl2mr_batch_tables <- c("main_results", "method_comparison", "locus_results", "overlap_qc", "warnings",
  "errors", "run_manifest", "sample_exclusions", "probability_normalization")

qtl2mr_empty_batch <- function(resources, pairs, instruments, config) {
  config$probability_policy <- resources$probability_policy %||% "strict"
  out <- stats::setNames(lapply(qtl2mr_batch_tables, qtl2mr_empty_batch_table), qtl2mr_batch_tables)
  out$sample_exclusions <- qtl2mr_batch_sample_exclusions(resources, pairs, instruments)
  out$input_alignment <- resources$sample_alignment %||% qtl2mr_input_alignment(resources)
  out$config <- config
  out$resources_metadata <- resources$metadata
  out$pairs <- pairs
  class(out) <- c("qtl2mr_batch_result", "list")
  out
}

qtl2mr_batch_append <- function(batch, name, rows) {
  if (is.null(rows) || !nrow(rows)) return(batch)
  # Point single-analysis file references in warnings to the batch tables.
  if (identical(name, "warnings") && "warning_id" %in% names(rows)) {
    normalized <- rows$warning_id == "probability_rounding_normalized"
    covariates <- rows$warning_id == "incomplete_covariates"
    for (column in intersect(c("message", "recommended_action"), names(rows))) {
      rows[[column]][normalized] <- gsub("data/probability_normalization.csv",
        "batch_probability_normalization.csv", rows[[column]][normalized], fixed = TRUE)
      rows[[column]][covariates] <- gsub("Per-locus, per-role counts and sample IDs are recorded in data/sample_alignment.csv; counts across loci are not additive.",
        "Per-trait counts and sample IDs are recorded in batch_sample_exclusions.csv.", rows[[column]][covariates], fixed = TRUE)
    }
  }
  batch[[name]] <- if (nrow(batch[[name]])) {
    common <- union(names(batch[[name]]), names(rows))
    for (nm in setdiff(common, names(batch[[name]]))) batch[[name]][[nm]] <- NA
    for (nm in setdiff(common, names(rows))) rows[[nm]] <- NA
    rbind(batch[[name]][common], rows[common])
  } else {
    rows
  }
  batch
}

# Every instrument of the pair with its first-stage statistics; loci dropped as
# non-estimable keep their row with kept = FALSE and the failing role.
qtl2mr_batch_locus_results <- function(pair, pilot) {
  original <- pilot$selected_loci_original %||% data.frame()
  if (!nrow(original)) return(qtl2mr_empty_batch_table("locus_results"))
  first <- pilot$locus_first_stage %||% qtl2mr_locus_first_stage(pilot$locus_blocks)
  dropped <- pilot$dropped_loci %||% data.frame()
  idx <- match(original$marker, first$marker)
  reason <- vapply(original$marker, function(marker) {
    d <- if (nrow(dropped)) dropped[dropped$marker == marker, , drop = FALSE] else dropped
    if (!nrow(d)) return("")
    kind <- ifelse(d$drop_reason == "nonestimable_founder_contrast", "nonestimable_", "fit_error_")
    paste(unique(paste0(kind, d$role)), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
  data.frame(
    pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
    marker = original$marker, chr = as.character(original$chr),
    pos = if ("pos" %in% names(original)) original$pos else NA_real_,
    kept = !is.na(idx), drop_reason = reason,
    F_X_block = first$F_X_block[idx], F_Y_block = first$F_Y_block[idx],
    partial_R2_X_block = first$partial_R2_X_block[idx], partial_R2_Y_block = first$partial_R2_Y_block[idx],
    steiger_pass = first$steiger_pass[idx], n_x = first$n_x[idx], n_y = first$n_y[idx],
    hsq_x = first$hsq_x[idx], hsq_y = first$hsq_y[idx],
    stringsAsFactors = FALSE
  )
}

qtl2mr_batch_main_row <- function(pair, pilot, runtime_sec) {
  s <- pilot$summary
  missing <- setdiff(qtl2mr_batch_summary_columns, names(s))
  if (length(missing)) stop("Internal error: summary is missing ", paste(missing, collapse = ", "), call. = FALSE)
  row <- s[1L, qtl2mr_batch_summary_columns, drop = FALSE]
  rownames(row) <- NULL
  cbind(
    data.frame(pair_id = pair$pair_id, stringsAsFactors = FALSE),
    row,
    data.frame(
      status = "ok", error_message = "",
      locus_fit_status = if (isTRUE((s$n_loci_dropped_nonestimable %||% 0L) > 0L)) "ok_with_dropped_nonestimable_loci" else "ok",
      runtime_sec = runtime_sec, stringsAsFactors = FALSE
    )
  )
}

qtl2mr_batch_error_main <- function(pairs, message, runtime_sec = NA_real_) {
  data.frame(pair_id = pairs$pair_id, exposure = pairs$exposure, outcome = pairs$outcome,
    status = "error", error_message = message, runtime_sec = runtime_sec, stringsAsFactors = FALSE)
}

qtl2mr_error_row <- function(pair, stage, error, runtime_sec = NA_real_) {
  data.frame(
    pair_id = pair$pair_id %||% NA_character_,
    exposure = pair$exposure %||% NA_character_,
    outcome = pair$outcome %||% NA_character_,
    stage = stage,
    error_class = class(error)[[1L]],
    error_message = conditionMessage(error),
    traceback_available = FALSE,
    runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
}

qtl2mr_json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  gsub("\"", "\\\\\"", x)
}

qtl2mr_json_value <- function(x) {
  if (is.null(x)) return("null")
  if ((is.atomic(x) || is.factor(x)) && length(x) > 1L) {
    return(paste0("[", paste(vapply(as.list(x), qtl2mr_json_value, character(1)), collapse = ","), "]"))
  }
  if (is.logical(x)) return(ifelse(is.na(x), "null", ifelse(x, "true", "false")))
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), "null", as.character(x)))
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      return(paste0("[", paste(vapply(x, qtl2mr_json_value, character(1)), collapse = ","), "]"))
    }
    fields <- sprintf("\"%s\":%s", qtl2mr_json_escape(nms), vapply(x, qtl2mr_json_value, character(1)))
    return(paste0("{", paste(fields, collapse = ","), "}"))
  }
  paste0("\"", qtl2mr_json_escape(x), "\"")
}

write_qtl2mr_batch_outputs <- function(batch, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- c(qtl2mr_batch_tables, if (!is.null(batch$input_alignment)) "input_alignment")
  for (name in tables) {
    utils::write.csv(batch[[name]], file.path(output_dir, paste0("batch_", name, ".csv")), row.names = FALSE)
  }
  writeLines(qtl2mr_json_value(batch$config), file.path(output_dir, "batch_config.json"))
  invisible(c(paste0("batch_", tables, ".csv"), "batch_config.json"))
}

qtl2mr_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(name, " must be TRUE or FALSE", call. = FALSE)
  x
}

qtl2mr_batch_records_dir <- function(output_dir) file.path(output_dir, "batch_pair_records")

qtl2mr_atomic_save_rds <- function(object, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 2L)
  if (!file.rename(tmp, path)) stop("Cannot write checkpoint file: ", path, call. = FALSE)
  invisible(path)
}

qtl2mr_write_resume_state <- function(output_dir, fingerprint, table_hashes = NULL) {
  qtl2mr_atomic_save_rds(list(schema = "qtl2mr-batch-v060", fingerprint = fingerprint,
    tables = table_hashes), file.path(output_dir, "batch_resume_state.rds"))
}

qtl2mr_pair_record_path <- function(output_dir, pair_index) {
  file.path(qtl2mr_batch_records_dir(output_dir), sprintf("pair_%06d.rds", as.integer(pair_index)))
}

qtl2mr_write_pair_record <- function(output_dir, pair_index, pair_id, fingerprint, payload) {
  qtl2mr_atomic_save_rds(list(
    schema = "qtl2mr-batch-pair-v060", pair_index = as.integer(pair_index), pair_id = pair_id,
    fingerprint = fingerprint, payload = payload, payload_md5 = qtl2mr_fingerprint_object(payload)
  ), qtl2mr_pair_record_path(output_dir, pair_index))
}

# Validate a checkpoint before anything in output_dir is modified. Completed pairs
# are stored as one record per pair, so a run interrupted at any point can resume
# from the records that were fully written. Paths are relative to output_dir and
# the worker count is not part of the fingerprint.
qtl2mr_load_resume_records <- function(output_dir, fingerprint, pairs) {
  state_path <- file.path(output_dir, "batch_resume_state.rds")
  if (!file.exists(state_path)) {
    qtl2mr_resume_abort("Cannot resume: legacy or incomplete output has no compatible resume fingerprint")
  }
  state <- tryCatch(readRDS(state_path), error = function(e) {
    qtl2mr_resume_abort(paste("Cannot resume: unreadable checkpoint:", conditionMessage(e)))
  })
  if (!is.list(state) || !identical(state$schema, "qtl2mr-batch-v060") || !identical(state$fingerprint, fingerprint)) {
    qtl2mr_resume_abort("Cannot resume: input, parameters, instruments or implementation fingerprint changed")
  }
  if (length(state$tables)) {
    files <- file.path(output_dir, names(state$tables))
    if (!all(file.exists(files)) || !identical(unname(tools::md5sum(files)), unname(state$tables))) {
      qtl2mr_resume_abort("Cannot resume: checkpoint output tables changed or are incomplete")
    }
  }
  files <- list.files(qtl2mr_batch_records_dir(output_dir), pattern = "^pair_[0-9]+\\.rds$", full.names = TRUE)
  records <- lapply(files, function(path) {
    record <- tryCatch(readRDS(path), error = function(e) {
      qtl2mr_resume_abort(paste0("Cannot resume: unreadable pair record ", basename(path),
        "; remove it or use overwrite = TRUE"))
    })
    index <- record$pair_index %||% NA_integer_
    valid <- is.list(record) && identical(record$schema, "qtl2mr-batch-pair-v060") &&
      identical(record$fingerprint, fingerprint) && length(index) == 1L && !is.na(index) &&
      index >= 1L && index <= nrow(pairs) && identical(record$pair_id, pairs$pair_id[[index]]) &&
      identical(qtl2mr_fingerprint_object(record$payload), record$payload_md5)
    if (!isTRUE(valid)) {
      qtl2mr_resume_abort(paste0("Cannot resume: pair record ", basename(path),
        " changed or does not belong to this run"))
    }
    record
  })
  names(records) <- vapply(records, function(r) r$pair_id, character(1))
  records
}

qtl2mr_read_pair_record <- function(output_dir, pair_index, pair_id, fingerprint) {
  if (is.null(output_dir)) return(NULL)
  path <- qtl2mr_pair_record_path(output_dir, pair_index)
  if (!file.exists(path)) return(NULL)
  record <- tryCatch(readRDS(path), error = function(e) NULL)
  ok <- is.list(record) && identical(record$fingerprint, fingerprint) && identical(record$pair_id, pair_id) &&
    identical(qtl2mr_fingerprint_object(record$payload), record$payload_md5)
  if (isTRUE(ok)) record$payload else NULL
}

qtl2mr_clear_batch_outputs <- function(output_dir) {
  paths <- qtl2mr_batch_owned_paths(output_dir)
  unlink(paths, recursive = TRUE)
  invisible(paths)
}

qtl2mr_batch_append_pair <- function(batch, result) {
  batch <- qtl2mr_batch_append(batch, "main_results", result$main)
  batch <- qtl2mr_batch_append(batch, "method_comparison", result$methods)
  batch <- qtl2mr_batch_append(batch, "locus_results", result$loci)
  batch <- qtl2mr_batch_append(batch, "overlap_qc", result$overlap)
  batch <- qtl2mr_batch_append(batch, "warnings", result$warnings)
  batch <- qtl2mr_batch_append(batch, "errors", result$error)
  qtl2mr_batch_append(batch, "probability_normalization", result$probability)
}

# Fit one exposure-outcome pair exactly as run_qtl2mr_analysis() would fit
# qtl2mr_pair_input(resources, instruments, exposure, outcome). Errors are
# returned as error rows; the caller applies fail_fast.
qtl2mr_batch_fit_pair <- function(resources, pair, instruments, settings, context_cache, fit_cache) {
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    rows <- qtl2mr_exposure_instruments(instruments, pair$exposure)
    prepared <- qtl2mr_pair_prepared(resources, rows, pair$exposure, pair$outcome,
      metadata = c(resources$metadata, list(pair_id = pair$pair_id)))
    pilot <- qtl2mr_fit_prepared(prepared, pair$exposure, pair$outcome, settings, verbose = FALSE,
      context_cache = context_cache, fit_cache = fit_cache)
    if (pilot$summary$n_loci[[1L]] < settings$min_loci) {
      stop(sprintf("too_few_loci: %d loci remain after dropping non-estimable loci, below min_loci = %d",
        pilot$summary$n_loci[[1L]], settings$min_loci), call. = FALSE)
    }
    runtime <- proc.time()[["elapsed"]] - t0
    warn <- qtl2mr_warning_table(pilot)
    if (nrow(warn)) warn <- cbind(pair_id = pair$pair_id, warn)
    overlap <- pilot$overlap_summary
    list(
      main = qtl2mr_batch_main_row(pair, pilot, runtime),
      methods = cbind(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
        qtl2mr_method_comparison_table(pilot)),
      loci = qtl2mr_batch_locus_results(pair, pilot),
      overlap = if (nrow(overlap)) cbind(pair_id = pair$pair_id, overlap) else NULL,
      warnings = warn,
      probability = qtl2mr_batch_probability_table(pair, pilot), error = NULL
    )
  }, error = function(e) {
    runtime <- proc.time()[["elapsed"]] - t0
    list(
      main = qtl2mr_batch_error_main(pair, conditionMessage(e), runtime),
      error = qtl2mr_error_row(pair, "pair", e, runtime_sec = runtime)
    )
  })
}

qtl2mr_batch_worker_failure_rows <- function(pairs, message) {
  list(
    main = qtl2mr_batch_error_main(pairs, message),
    error = data.frame(pair_id = pairs$pair_id, exposure = pairs$exposure, outcome = pairs$outcome,
      stage = "parallel_worker", error_class = "qtl2mr_worker_failure", error_message = message,
      traceback_available = FALSE, runtime_sec = NA_real_, stringsAsFactors = FALSE)
  )
}

# Tables are rebuilt in pair order from resumed records and new results.
qtl2mr_batch_assemble <- function(batch, records, results, todo, pairs, preflight) {
  by_id <- c(lapply(records, function(r) r$payload), stats::setNames(results, todo$pair_id))
  for (id in pairs$pair_id) if (!is.null(by_id[[id]])) batch <- qtl2mr_batch_append_pair(batch, by_id[[id]])
  batch <- qtl2mr_append_preflight_errors(batch, preflight)
  qtl2mr_batch_finalize(batch, pairs)
}

# On a fail_fast stop, completed pairs are published so the user can inspect them;
# their records already allow resume = TRUE to continue.
qtl2mr_batch_publish_partial <- function(output_dir, records, results, todo, batch, pairs, preflight) {
  if (is.null(output_dir)) return(invisible(NULL))
  done <- !vapply(results, is.null, logical(1))
  batch <- qtl2mr_batch_assemble(batch, records, results[done], todo[done, , drop = FALSE], pairs, preflight)
  batch$run_manifest <- data.frame(created_at = as.character(Sys.time()), n_pairs_requested = nrow(pairs),
    n_pairs_ok = sum(batch$main_results$status == "ok", na.rm = TRUE),
    n_pairs_error = sum(batch$main_results$status == "error", na.rm = TRUE), status = "stopped_fail_fast",
    stringsAsFactors = FALSE)
  write_qtl2mr_batch_outputs(batch, output_dir)
  invisible(output_dir)
}

# Every table gets its documented columns first, in a fixed order, whether or
# not it has rows; main results follow the order of pairs.
qtl2mr_batch_finalize <- function(batch, pairs) {
  for (name in setdiff(qtl2mr_batch_tables, "run_manifest")) {
    table <- batch[[name]]
    canonical <- names(qtl2mr_empty_batch_table(name))
    for (nm in setdiff(canonical, names(table))) table[[nm]] <- rep(NA, nrow(table))
    table <- table[, c(canonical, setdiff(names(table), canonical)), drop = FALSE]
    rownames(table) <- NULL
    batch[[name]] <- table
  }
  main <- batch$main_results
  if (nrow(main)) main <- main[order(match(main$pair_id, pairs$pair_id), na.last = TRUE), , drop = FALSE]
  rownames(main) <- NULL
  batch$main_results <- main
  batch
}

run_qtl2mr_batch <- function(resources,
                             pairs,
                             instruments,
                             on_nonestimable_locus = c("fail", "drop"),
                             reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
                             heterogeneity_method = c("chisq", "simulation", "auto", "none"),
                             heterogeneity_nsim = 1000L,
                             heterogeneity_seed = 1L,
                             diagnostic_bootstrap_B = 0L,
                             diagnostic_bootstrap_seed = 1L,
                             standardize_phenotypes = TRUE,
                             cross_cov_model = c("star", "raw", "zero"),
                             min_samples = 50L,
                             min_loci = 1L,
                             n_workers = 1L,
                             fail_fast = FALSE,
                             resume = FALSE,
                             overwrite = FALSE,
                             output_dir = NULL,
                             progress = TRUE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  if (missing(instruments)) {
    stop("instruments is required; select them with qtl2mr_select_instruments() or supply a table with exposure and marker columns", call. = FALSE)
  }
  fail_fast <- qtl2mr_check_flag(fail_fast, "fail_fast")
  resume <- qtl2mr_check_flag(resume, "resume")
  overwrite <- qtl2mr_check_flag(overwrite, "overwrite")
  progress <- qtl2mr_check_flag(progress, "progress")
  n_workers <- qtl2mr_check_count(n_workers, "n_workers")
  settings <- qtl2mr_fit_settings(
    on_nonestimable_locus = match.arg(on_nonestimable_locus),
    reported_estimator_policy = match.arg(reported_estimator_policy),
    heterogeneity_method = match.arg(heterogeneity_method),
    heterogeneity_nsim = heterogeneity_nsim, heterogeneity_seed = heterogeneity_seed,
    diagnostic_bootstrap_B = diagnostic_bootstrap_B, diagnostic_bootstrap_seed = diagnostic_bootstrap_seed,
    standardize_phenotypes = standardize_phenotypes, cross_cov_model = match.arg(cross_cov_model)
  )
  settings$min_samples <- qtl2mr_check_count(min_samples, "min_samples")
  settings$min_loci <- qtl2mr_check_count(min_loci, "min_loci")
  pairs <- qtl2mr_normalize_pairs(pairs)
  instruments <- qtl2mr_normalize_instruments(instruments, resources)
  config <- settings

  # Destination, pair preflight and the checkpoint are all checked before
  # anything in output_dir is written or removed.
  qtl2mr_check_batch_destination(output_dir, resume, overwrite)
  preflight <- qtl2mr_batch_preflight(resources, pairs, instruments, settings$min_loci,
    settings$min_samples, fail_fast)
  valid <- preflight$valid
  batch <- qtl2mr_empty_batch(resources, valid, instruments, c(config, list(
    batch_traits = unique(c(valid$exposure, valid$outcome)), created_at = as.character(Sys.time()))))
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  t_run <- proc.time()[["elapsed"]]

  fingerprint <- if (is.null(output_dir)) NULL else
    qtl2mr_batch_fingerprint(resources, pairs, instruments, config, markers = preflight$markers)
  occupied <- !is.null(output_dir) && length(qtl2mr_batch_owned_paths(output_dir)) > 0L
  records <- list()
  if (occupied && resume && !overwrite) records <- qtl2mr_load_resume_records(output_dir, fingerprint, pairs)
  if (occupied && overwrite) qtl2mr_clear_batch_outputs(output_dir)
  if (!is.null(output_dir)) {
    dir.create(qtl2mr_batch_records_dir(output_dir), recursive = TRUE, showWarnings = FALSE)
    qtl2mr_write_resume_state(output_dir, fingerprint)
  }
  todo <- valid[!valid$pair_id %in% names(records), , drop = FALSE]
  pair_index <- match(todo$pair_id, pairs$pair_id)
  fit_one <- function(i, context_cache, fit_cache) {
    pair <- todo[i, , drop = FALSE]
    result <- qtl2mr_batch_fit_pair(resources, pair, instruments, settings, context_cache, fit_cache)
    if (!is.null(output_dir) && is.null(result$error)) {
      qtl2mr_write_pair_record(output_dir, pair_index[[i]], pair$pair_id, fingerprint, result)
    }
    result
  }
  results <- vector("list", nrow(todo))
  cache_stats <- list()
  first_failure <- NULL
  backend <- "sequential"
  groups <- split(seq_len(nrow(todo)), factor(todo$exposure, levels = unique(todo$exposure)))
  use_fork <- n_workers > 1L && length(groups) > 1L && identical(.Platform$OS.type, "unix")
  if (use_fork) {
    backend <- "parallel::mclapply"
    # Each exposure group runs in one worker. With fail_fast, a worker stops at
    # its first failing pair; groups already running elsewhere finish first.
    run_group <- function(idx) {
      context_cache <- qtl2mr_context_cache_new()
      fit_cache <- qtl2mr_context_cache_new()
      if (progress) message(sprintf("[qtl2mr batch] exposure %s: %d pairs", todo$exposure[[idx[[1L]]]], length(idx)))
      out <- vector("list", length(idx))
      for (k in seq_along(idx)) {
        out[[k]] <- fit_one(idx[[k]], context_cache, fit_cache)
        if (fail_fast && !is.null(out[[k]]$error)) break
      }
      list(results = out, stats = list(qtl2mr_context_cache_stats(context_cache), qtl2mr_context_cache_stats(fit_cache)))
    }
    group_out <- parallel::mclapply(groups, function(idx) tryCatch(run_group(idx), error = function(e) e),
      mc.cores = min(n_workers, length(groups)), mc.preschedule = FALSE)
    for (g in seq_along(groups)) {
      value <- group_out[[g]]
      idx <- groups[[g]]
      if (is.list(value) && !inherits(value, "condition") && is.list(value$results)) {
        results[idx] <- value$results
        cache_stats[[length(cache_stats) + 1L]] <- value$stats
        next
      }
      # A worker that died or errored outside a pair: keep its fully written
      # records and report the remaining pairs of the group as errors.
      message_text <- if (inherits(value, "condition")) conditionMessage(value) else
        if (inherits(value, "try-error")) as.character(value) else "parallel worker returned no result"
      failure <- qtl2mr_batch_worker_failure_rows(todo[idx, , drop = FALSE],
        paste("parallel worker failed:", message_text))
      for (k in seq_along(idx)) {
        saved <- qtl2mr_read_pair_record(output_dir, pair_index[[idx[[k]]]], todo$pair_id[[idx[[k]]]], fingerprint)
        results[[idx[[k]]]] <- saved %||% list(main = failure$main[k, , drop = FALSE], error = failure$error[k, , drop = FALSE])
      }
    }
    if (fail_fast) {
      failed <- which(vapply(results, function(r) !is.null(r$error), logical(1)))
      if (length(failed)) first_failure <- results[[failed[[1L]]]]$error
    }
  } else {
    context_cache <- qtl2mr_context_cache_new()
    fit_cache <- qtl2mr_context_cache_new()
    for (i in seq_len(nrow(todo))) {
      if (progress) message(sprintf("[qtl2mr batch] %d/%d %s: %s -> %s", i, nrow(todo), todo$pair_id[[i]],
        todo$exposure[[i]], todo$outcome[[i]]))
      results[[i]] <- fit_one(i, context_cache, fit_cache)
      if (fail_fast && !is.null(results[[i]]$error)) {
        first_failure <- results[[i]]$error
        break
      }
    }
    cache_stats[[1L]] <- list(qtl2mr_context_cache_stats(context_cache), qtl2mr_context_cache_stats(fit_cache))
  }
  if (!is.null(first_failure)) {
    qtl2mr_batch_publish_partial(output_dir, records, results, todo, batch, pairs, preflight)
    stop(sprintf("Pair %s failed (fail_fast = TRUE): %s", first_failure$pair_id[[1L]],
      first_failure$error_message[[1L]]), call. = FALSE)
  }

  batch <- qtl2mr_batch_assemble(batch, records, results, todo, pairs, preflight)
  stat_sum <- function(which, field) sum(vapply(cache_stats, function(s) as.numeric(s[[which]][[field]] %||% 0), numeric(1)))
  wall_sec <- proc.time()[["elapsed"]] - t_run
  batch$run_manifest <- data.frame(
    run_id = run_id, created_at = as.character(Sys.time()),
    n_pairs_requested = nrow(pairs),
    n_pairs_ok = sum(batch$main_results$status == "ok", na.rm = TRUE),
    n_pairs_error = sum(batch$main_results$status == "error", na.rm = TRUE),
    n_pairs_resumed = length(records),
    nonestimable_locus_policy = if (identical(settings$on_nonestimable_locus, "drop")) "drop_audited" else "fail_closed",
    n_loci_dropped_nonestimable = sum(batch$main_results$n_loci_dropped_nonestimable %||% 0L, na.rm = TRUE),
    parallel_backend = backend, n_workers = n_workers,
    context_cache_hits = stat_sum(1L, "hits"), context_cache_misses = stat_sum(1L, "misses"),
    lmm_fit_cache_hits = stat_sum(2L, "hits"), lmm_fit_cache_misses = stat_sum(2L, "misses"),
    wall_sec = wall_sec, mean_wall_sec_per_pair = wall_sec / max(nrow(todo), 1L),
    lmm_backend = qtl2mr_backend_label(batch$main_results$lmm_backend[batch$main_results$status == "ok"]),
    stringsAsFactors = FALSE
  )
  if (!is.null(output_dir)) {
    files <- write_qtl2mr_batch_outputs(batch, output_dir)
    hashes <- tools::md5sum(file.path(output_dir, files))
    names(hashes) <- files
    qtl2mr_write_resume_state(output_dir, fingerprint, hashes)
  }
  batch
}
