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
                                     qtl_scan = NULL,
                                     marker_map = NULL,
                                     metadata = list(),
                                     probability_policy = c("normalize_rounding", "strict")) {
  probability_policy <- match.arg(probability_policy)
  if (!is.list(allele_probs) || !length(allele_probs)) {
    stop("allele_probs must be a non-empty qtl2-style list", call. = FALSE)
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
  if (!"pos" %in% names(marker_map)) marker_map$pos <- seq_len(nrow(marker_map))
  gp_ids <- dimnames(allele_probs[[1L]])[[1L]]
  if (is.null(gp_ids)) stop("allele_probs arrays must have sample row names", call. = FALSE)
  marker_map$chr <- as.character(marker_map$chr)
  marker_map$marker <- as.character(marker_map$marker)
  out <- list(
    allele_probs = allele_probs,
    kinship = kinship,
    pheno = pheno,
    covar = covar,
    qtl_scan = qtl_scan,
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

qtl2mr_pair_lod_thresh <- function(pair, default) {
  if (!"lod_thresh" %in% names(pair)) return(qtl2mr_validate_lod_thresh(default))
  value <- pair$lod_thresh
  # An absent per-pair value inherits the validated batch default. Any supplied
  # non-missing value is validated rather than silently coerced or ignored.
  if (length(value) == 1L && is.na(value) && (!is.numeric(value) || !is.nan(value))) {
    return(qtl2mr_validate_lod_thresh(default))
  }
  qtl2mr_validate_lod_thresh(value)
}

qtl2mr_pair_max_loci <- function(pair, default) {
  value <- if ("max_loci" %in% names(pair)) pair$max_loci else NULL
  if (is.null(value) || (length(value) == 1L && is.na(value) && !is.nan(value))) return(default)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 1 || value != as.integer(value)) stop("max_loci must be a positive integer", call. = FALSE)
  as.integer(value)
}

qtl2mr_lod_vector_from_scan <- function(qtl_scan, instrument_trait) {
  if (length(instrument_trait) != 1L || is.na(instrument_trait) || !nzchar(instrument_trait)) {
    stop("instrument_trait must be one non-empty trait name", call. = FALSE)
  }
  scan <- qtl_scan
  trait_validated_by_list <- FALSE
  if (is.list(qtl_scan) && !is.data.frame(qtl_scan) && !is.matrix(qtl_scan)) {
    scan_names <- names(qtl_scan)
    if (is.null(scan_names) || any(is.na(scan_names) | !nzchar(scan_names))) {
      stop("named-list qtl_scan must use non-empty trait names", call. = FALSE)
    }
    idx <- which(scan_names == instrument_trait)
    if (!length(idx)) stop("qtl_scan is missing instrument trait: ", instrument_trait, call. = FALSE)
    if (length(idx) > 1L) stop("qtl_scan has duplicated instrument trait: ", instrument_trait, call. = FALSE)
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
      idx <- match(instrument_trait, trait_names)
      if (is.na(idx)) stop("qtl_scan is missing instrument trait: ", instrument_trait, call. = FALSE)
      lod <- as.numeric(scan[, idx])
    } else {
      if (!is.null(trait_names) && length(trait_names) && !is.na(trait_names[[1L]]) && nzchar(trait_names[[1L]])) {
        if (!isTRUE(trait_validated_by_list) && !identical(trait_names[[1L]], instrument_trait)) {
          stop("single-column qtl_scan trait name does not match instrument_trait: ", instrument_trait, call. = FALSE)
        }
      } else {
        warning(
          "single-column qtl_scan has no trait column name; instrument_trait could not be validated",
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

qtl2mr_select_loci_from_scan <- function(qtl_scan,
                                         marker_map,
                                         instrument_trait,
                                         lod_thresh = 4,
                                         exclude_chr = "X",
                                         max_loci = NULL) {
  lod_thresh <- qtl2mr_validate_lod_thresh(lod_thresh)
  lod <- qtl2mr_lod_vector_from_scan(qtl_scan, instrument_trait)
  mm <- qtl2mr_marker_map_table(marker_map)
  map_idx <- match(names(lod), mm$marker)
  if (anyNA(map_idx)) stop("marker_map is missing qtl_scan markers", call. = FALSE)
  if (!is.null(max_loci) && (length(max_loci) != 1L || !is.numeric(max_loci) || !is.finite(max_loci) || max_loci < 1 || max_loci != as.integer(max_loci))) stop("max_loci must be a positive integer", call. = FALSE)
  # Islands are defined in map order within a chromosome, never across a
  # boundary in the concatenated scan matrix.
  peak_idx <- unlist(lapply(unique(mm$chr[map_idx]), function(chr) {
    idx <- which(mm$chr[map_idx] == chr)
    pos <- mm$pos[map_idx[idx]]
    if (all(is.finite(pos))) idx <- idx[order(pos, map_idx[idx])]
    idx[qtl2mr_get_tophits(lod[idx], lod_thresh = lod_thresh)]
  }), use.names = FALSE)
  if (!length(peak_idx)) {
    return(data.frame(instrument_trait = character(), marker = character(), chr = character(), pos = numeric(),
                      exposure_lod = numeric(), locus_rank = integer(),
                      selected_by = character(), lod_thresh = numeric(),
                      excluded_chr_rule = character(), stringsAsFactors = FALSE))
  }
  idx <- match(names(lod)[peak_idx], mm$marker)
  selected <- data.frame(
    marker = names(lod)[peak_idx],
    chr = mm$chr[idx],
    pos = if ("pos" %in% names(mm)) mm$pos[idx] else peak_idx,
    exposure_lod = as.numeric(lod[peak_idx]),
    stringsAsFactors = FALSE
  )
  selected <- selected[!is.na(selected$chr), , drop = FALSE]
  if (!is.null(exclude_chr) && length(exclude_chr)) {
    selected <- selected[!selected$chr %in% as.character(exclude_chr), , drop = FALSE]
  }
  selected <- selected[order(selected$exposure_lod, decreasing = TRUE), , drop = FALSE]
  if (!is.null(max_loci)) selected <- utils::head(selected, as.integer(max_loci))
  selected$locus_rank <- seq_len(nrow(selected))
  selected$instrument_trait <- instrument_trait
  selected$selected_by <- "contiguous_lod_island_max"
  selected$lod_thresh <- lod_thresh
  selected$excluded_chr_rule <- paste(exclude_chr %||% character(), collapse = ",")
  rownames(selected) <- NULL
  selected[, c("instrument_trait", setdiff(names(selected), "instrument_trait")), drop = FALSE]
}

qtl2mr_select_loci_for_traits <- function(qtl_scan,
                                          marker_map,
                                          instrument_traits,
                                          lod_thresh = 4,
                                          exclude_chr = "X",
                                          max_loci = NULL) {
  if (is.null(instrument_traits) || !length(instrument_traits)) {
    stop("instrument_traits must contain at least one trait name", call. = FALSE)
  }
  rows <- lapply(as.character(instrument_traits), function(trait) {
    selected <- qtl2mr_select_loci_from_scan(
      qtl_scan = qtl_scan,
      marker_map = marker_map,
      instrument_trait = trait,
      lod_thresh = lod_thresh,
      exclude_chr = exclude_chr,
      max_loci = max_loci
    )
    if (!nrow(selected)) {
      selected$instrument_trait <- character()
      return(selected)
    }
    if (!"instrument_trait" %in% names(selected)) selected$instrument_trait <- trait
    selected[, c("instrument_trait", setdiff(names(selected), "instrument_trait")), drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

qtl2mr_make_common_complete_samples <- function(resources,
                                                exposure,
                                                outcome,
                                                selected_chr = NULL,
                                                sample_rule = "common_complete") {
  mode <- if (is.list(sample_rule)) sample_rule$mode %||% "common_complete" else sample_rule
  if (identical(mode, "pair_common_complete")) mode <- "common_complete"
  if (identical(mode, "batch_common_complete")) {
    stop("Use qtl2mr_make_batch_common_complete_samples() for batch_common_complete sample sets", call. = FALSE)
  }
  if (!identical(mode, "common_complete")) {
    stop("Supported batch sample modes are 'common_complete', 'pair_common_complete', and 'batch_common_complete'", call. = FALSE)
  }
  requested <- if (is.list(sample_rule)) sample_rule$sample_ids else NULL
  if (!is.null(requested)) {
    qtl2mr_validate_ids(requested, "sample_rule sample_ids")
    if (length(setdiff(requested, resources$sample_ids))) stop("sample_rule sample_ids are absent from allele_probs", call. = FALSE)
  }
  pheno <- resources$pheno
  covar <- resources$covar
  ids <- Reduce(intersect, list(resources$sample_ids, rownames(pheno), if (is.null(covar)) resources$sample_ids else rownames(covar)))
  if (!is.null(requested)) ids <- ids[ids %in% requested]
  if (!exposure %in% names(pheno)) stop("Exposure trait not found in pheno: ", exposure, call. = FALSE)
  if (!outcome %in% names(pheno)) stop("Outcome trait not found in pheno: ", outcome, call. = FALSE)
  if (!is.numeric(pheno[[exposure]]) || !is.numeric(pheno[[outcome]])) stop("Exposure and outcome phenotypes must be numeric", call. = FALSE)
  if (!is.null(resources$kinship)) {
    if (is.list(resources$kinship) && !is.matrix(resources$kinship) && !is.data.frame(resources$kinship)) {
      for (chr in unique(as.character(selected_chr))) {
        K <- resources$kinship[[chr]]
        if (is.null(K)) stop("kinship list is missing selected chromosome: ", chr, call. = FALSE)
        ids <- ids[ids %in% rownames(as.matrix(K)) & ids %in% colnames(as.matrix(K))]
      }
    } else {
      K <- as.matrix(resources$kinship)
      ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
    }
  }
  keep <- is.finite(pheno[ids, exposure]) & is.finite(pheno[ids, outcome])
  if (!is.null(covar) && ncol(covar)) keep <- keep & qtl2mr_complete_rows(covar[ids, , drop = FALSE])
  ids[keep]
}

qtl2mr_make_batch_common_complete_samples <- function(resources,
                                                      traits,
                                                      selected_chr = NULL) {
  pheno <- resources$pheno
  covar <- resources$covar
  traits <- unique(as.character(traits))
  missing <- setdiff(traits, colnames(pheno))
  if (length(missing)) stop("pheno is missing batch traits: ", paste(missing, collapse = ", "), call. = FALSE)
  if (any(!vapply(pheno[traits], is.numeric, logical(1)))) stop("Batch phenotypes must be numeric", call. = FALSE)
  ids <- Reduce(intersect, list(resources$sample_ids, rownames(pheno), if (is.null(covar)) resources$sample_ids else rownames(covar)))
  if (length(selected_chr)) {
    selected_chr <- unique(as.character(selected_chr))
    for (chr in selected_chr) {
      gp <- resources$allele_probs[[chr]]
      if (is.null(gp)) stop("allele_probs is missing selected chromosome: ", chr, call. = FALSE)
      ids <- ids[ids %in% dimnames(gp)[[1L]]]
    }
    if (is.list(resources$kinship) && !is.matrix(resources$kinship) && !is.data.frame(resources$kinship)) {
      for (chr in selected_chr) {
        K <- resources$kinship[[chr]]
        if (is.null(K)) stop("kinship list is missing selected chromosome: ", chr, call. = FALSE)
        ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
      }
    } else if (!is.null(resources$kinship)) {
      K <- as.matrix(resources$kinship)
      ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
    }
  }
  keep <- qtl2mr_complete_rows(pheno[ids, traits, drop = FALSE])
  if (!is.null(covar) && ncol(covar)) keep <- keep & qtl2mr_complete_rows(covar[ids, , drop = FALSE])
  ids[keep]
}

qtl2mr_make_batch_eligible_samples <- function(resources, selected_chr = NULL) {
  covar <- resources$covar
  ids <- if (is.null(covar)) resources$sample_ids else intersect(resources$sample_ids, rownames(covar))
  if (length(selected_chr)) {
    selected_chr <- unique(as.character(selected_chr))
    for (chr in selected_chr) {
      gp <- resources$allele_probs[[chr]]
      if (is.null(gp)) stop("allele_probs is missing selected chromosome: ", chr, call. = FALSE)
      ids <- ids[ids %in% dimnames(gp)[[1L]]]
    }
    if (is.list(resources$kinship) && !is.matrix(resources$kinship) && !is.data.frame(resources$kinship)) {
      for (chr in selected_chr) {
        K <- resources$kinship[[chr]]
        if (is.null(K)) stop("kinship list is missing selected chromosome: ", chr, call. = FALSE)
        ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
      }
    } else if (!is.null(resources$kinship)) {
      K <- as.matrix(resources$kinship)
      ids <- ids[ids %in% rownames(K) & ids %in% colnames(K)]
    }
  }
  if (!is.null(covar) && ncol(covar)) ids <- ids[qtl2mr_complete_rows(covar[ids, , drop = FALSE])]
  ids
}

qtl2mr_selected_loci_for_pairs <- function(resources, pairs, instrument_rule) {
  out <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    lod_thresh <- qtl2mr_pair_lod_thresh(pair, instrument_rule$lod_thresh)
    exclude_chr <- if ("exclude_chr" %in% names(pair) && qtl2mr_nonempty_string(pair$exclude_chr)) strsplit(pair$exclude_chr, ",")[[1L]] else instrument_rule$exclude_chr
    max_loci <- qtl2mr_pair_max_loci(pair, instrument_rule$max_loci)
    if ("selected_loci_file" %in% names(pair) &&
        qtl2mr_nonempty_string(pair$selected_loci_file)) {
      selected <- qtl2mr_read_selected_loci_file(pair$selected_loci_file)
    } else {
      selected <- qtl2mr_select_loci_from_scan(
        resources$qtl_scan,
        resources$marker_map,
        instrument_trait = pair$instrument_trait,
        lod_thresh = lod_thresh,
        exclude_chr = exclude_chr,
        max_loci = max_loci
      )
    }
    if (!nrow(selected)) return(data.frame())
    cbind(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
          instrument_trait = pair$instrument_trait, selected)
  })
  columns <- unique(unlist(lapply(out, names), use.names = FALSE))
  out <- lapply(out, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- rep(NA, nrow(x))
    x[columns]
  })
  out <- do.call(rbind, out)
  if (is.null(out)) out <- data.frame()
  rownames(out) <- NULL
  out
}

qtl2mr_resolve_instrument_request <- function(pair, instrument_rule) {
  lod_thresh <- qtl2mr_pair_lod_thresh(pair, instrument_rule$lod_thresh)
  exclude_chr <- if ("exclude_chr" %in% names(pair) && qtl2mr_nonempty_string(pair$exclude_chr)) {
    strsplit(pair$exclude_chr, ",")[[1L]]
  } else {
    instrument_rule$exclude_chr
  }
  max_loci <- qtl2mr_pair_max_loci(pair, instrument_rule$max_loci)
  selected_loci_file <- if ("selected_loci_file" %in% names(pair) &&
                            qtl2mr_nonempty_string(pair$selected_loci_file)) {
    as.character(pair$selected_loci_file)
  } else {
    NA_character_
  }
  list(
    instrument_trait = as.character(pair$instrument_trait),
    lod_thresh = lod_thresh,
    exclude_chr = as.character(exclude_chr %||% character()),
    max_loci = max_loci,
    selected_loci_file = selected_loci_file
  )
}

qtl2mr_select_loci_for_request <- function(resources, request) {
  if (!is.na(request$selected_loci_file)) {
    qtl2mr_read_selected_loci_file(request$selected_loci_file)
  } else {
    qtl2mr_select_loci_from_scan(
      resources$qtl_scan,
      resources$marker_map,
      instrument_trait = request$instrument_trait,
      lod_thresh = request$lod_thresh,
      exclude_chr = request$exclude_chr,
      max_loci = request$max_loci
    )
  }
}

qtl2mr_pair_id <- function(exposure, outcome) {
  paste0(make.names(exposure), "__to__", make.names(outcome))
}

qtl2mr_normalize_pairs <- function(pairs) {
  if (is.character(pairs) && length(pairs) == 1L && file.exists(pairs)) {
    pairs <- utils::read.csv(pairs, stringsAsFactors = FALSE)
  }
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("exposure", "outcome") %in% names(pairs))) {
    stop("pairs must contain exposure and outcome columns", call. = FALSE)
  }
  if (!"pair_id" %in% names(pairs)) {
    pairs$pair_id <- mapply(qtl2mr_pair_id, pairs$exposure, pairs$outcome, USE.NAMES = FALSE)
  }
  if (!"instrument_trait" %in% names(pairs)) pairs$instrument_trait <- pairs$exposure
  if (!nrow(pairs)) stop("pairs must contain at least one pair", call. = FALSE)
  for (nm in c("pair_id", "exposure", "outcome", "instrument_trait")) {
    pairs[[nm]] <- as.character(pairs[[nm]])
    if (anyNA(pairs[[nm]]) || any(!nzchar(pairs[[nm]]))) stop("pairs contains missing or blank ", nm, call. = FALSE)
  }
  if (anyDuplicated(pairs$pair_id)) stop("pair_id must be unique; supply explicit IDs for distinct requests", call. = FALSE)
  pairs
}

qtl2mr_nonempty_string <- function(x) {
  length(x) == 1L && !is.na(x) && nzchar(x)
}

qtl2mr_read_selected_loci_file <- function(path) {
  if (!qtl2mr_nonempty_string(path) || !file.exists(path) || dir.exists(path)) stop("selected_loci_file does not exist or is not a file: ", path, call. = FALSE)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("marker", "chr") %in% names(d))) {
    stop("selected_loci_file must contain marker and chr columns", call. = FALSE)
  }
  if (!"locus_rank" %in% names(d)) d$locus_rank <- seq_len(nrow(d))
  d
}

qtl2mr_batch_defaults <- function(x, defaults) {
  if (is.null(x)) return(defaults)
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))) stop("Options must be a uniquely named list", call. = FALSE)
  unknown <- setdiff(names(x), names(defaults))
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if ("source" %in% names(x) && !identical(x$source, "exposure_qtl_scan")) stop("Unsupported instrument_rule source", call. = FALSE)
  if ("peak_method" %in% names(x) && !identical(x$peak_method, "contiguous_lod_island_max")) stop("Unsupported peak_method", call. = FALSE)
  utils::modifyList(defaults, x, keep.null = TRUE)
}

qtl2mr_batch_fit_options <- function(options) {
  qtl2mr_batch_defaults(options, list(
    standardize_phenotypes = TRUE, cross_cov_model = "star",
    heterogeneity_method = "none", q_sim_n = 0L, q_sim_seed = 1L,
    diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L,
    reported_estimator_policy = "diagnostic_only", min_samples = 50L,
    min_loci = 1L, engine = "batch",
    cache_contexts = TRUE, cache_lmm_fits = TRUE, verbose = FALSE
  ))
}

qtl2mr_batch_sample_policy <- function(policy = NULL) {
  if (identical(policy, "batch_common_complete")) "batch_common_complete" else "pair_common_complete"
}

qtl2mr_batch_sample_exclusions <- function(resources, pairs, sample_policy = NULL,
                                           batch_traits = NULL, sample_ids = NULL) {
  if (!nrow(pairs)) return(qtl2mr_empty_batch_table("sample_exclusions"))
  ids <- resources$sample_ids
  pheno <- resources$pheno
  covar <- resources$covar
  policy <- qtl2mr_batch_sample_policy(sample_policy)
  missing_pheno <- !ids %in% rownames(pheno)
  missing_covar <- if (is.null(covar)) rep(FALSE, length(ids)) else !ids %in% rownames(covar)
  incomplete_covar <- rep(FALSE, length(ids))
  if (!is.null(covar) && ncol(covar)) {
    matched <- which(!missing_covar)
    incomplete_covar[matched] <- !qtl2mr_complete_rows(covar[ids[matched], , drop = FALSE])
  }
  nonfinite_trait <- function(trait) {
    bad <- rep(FALSE, length(ids))
    if (trait %in% names(pheno) && is.numeric(pheno[[trait]])) {
      matched <- which(!missing_pheno)
      bad[matched] <- !is.finite(pheno[ids[matched], trait])
    }
    bad
  }
  batch_traits <- batch_traits %||% unique(c(pairs$exposure, pairs$outcome))
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    reasons <- list(
      missing_phenotype_id = missing_pheno,
      nonfinite_exposure = nonfinite_trait(pair$exposure),
      nonfinite_outcome = nonfinite_trait(pair$outcome),
      missing_covar_id = missing_covar,
      incomplete_covar = incomplete_covar
    )
    if (!is.null(sample_ids)) reasons$not_in_sample_ids <- !ids %in% as.character(sample_ids)
    if (identical(policy, "batch_common_complete")) {
      other_traits <- setdiff(batch_traits, c(pair$exposure, pair$outcome))
      reasons$nonfinite_other_batch_trait <- Reduce(`|`, lapply(other_traits, nonfinite_trait),
        init = rep(FALSE, length(ids)))
    }
    do.call(rbind, lapply(names(reasons), function(reason) data.frame(
      pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
      sample_policy = policy, reason = reason, n_input_samples = length(ids),
      n_excluded = sum(reasons[[reason]]),
      excluded_ids = paste(ids[reasons[[reason]]], collapse = ";"),
      counts_overlap = TRUE, stringsAsFactors = FALSE
    )))
  })
  do.call(rbind, rows)
}

qtl2mr_batch_covariate_warnings <- function(exclusions) {
  rows <- exclusions[exclusions$reason %in% c("missing_covar_id", "incomplete_covar") &
    exclusions$n_excluded > 0L, , drop = FALSE]
  if (!nrow(rows)) return(qtl2mr_empty_batch_table("warnings"))
  do.call(rbind, lapply(unique(rows$pair_id), function(id) {
    d <- rows[rows$pair_id == id, , drop = FALSE]
    data.frame(pair_id = id, warning_id = "batch_covariate_sample_exclusions",
      severity = "warning",
      message = sprintf("%d input samples are ineligible for this pair because covariate IDs or finite covariate values are missing.", sum(d$n_excluded)),
      consequence = "The fitted cohort excludes these samples; exclusion reasons can overlap phenotype exclusions.",
      recommended_action = "Review batch_sample_exclusions.csv and fitted sample counts; no imputation is performed.",
      shown_in_main_report = TRUE, stringsAsFactors = FALSE)
  }))
}

qtl2mr_batch_probability_table <- function(pair, pilot) {
  adjustments <- pilot$prepared_metadata$probability_normalization$adjustments
  if (is.null(adjustments) || !nrow(adjustments)) return(qtl2mr_empty_batch_table("probability_normalization"))
  cbind(pair_id = pair$pair_id, adjustments, stringsAsFactors = FALSE)
}

qtl2mr_batch_preflight <- function(resources, pairs, instrument_rule, fit_options, fail_fast,
                                   sample_ids = NULL) {
  rule <- qtl2mr_batch_defaults(instrument_rule, list(source = "exposure_qtl_scan", lod_thresh = 4,
    peak_method = "contiguous_lod_island_max", exclude_chr = "X", max_loci = NULL))
  errors <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    tryCatch({
      selected <- qtl2mr_select_loci_for_request(resources, qtl2mr_resolve_instrument_request(pair, rule))
      if (nrow(selected) < fit_options$min_loci) stop("no_selected_loci: selected loci below min_loci", call. = FALSE)
      probe <- resources
      probe$selected_loci <- selected
      validate_qtl2mr_input(probe, exposure = pair$exposure, outcome = pair$outcome, require_selected_loci = TRUE)
      samples <- qtl2mr_make_common_complete_samples(resources, pair$exposure, pair$outcome, selected$chr,
        sample_rule = list(mode = "common_complete", sample_ids = sample_ids))
      if (length(samples) < fit_options$min_samples) stop("too_few_samples: complete samples below min_samples", call. = FALSE)
      as.character(selected$marker)
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
    main <- data.frame(pair_id = pairs$pair_id, exposure = pairs$exposure, outcome = pairs$outcome,
      instrument_trait = pairs$instrument_trait, status = "error", error_message = errors$error_message,
      runtime_sec = NA_real_, sample_policy = qtl2mr_batch_sample_policy(batch$config$sample_policy),
      stringsAsFactors = FALSE)
    batch <- qtl2mr_batch_append(batch, "main_results", main)
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

qtl2mr_batch_fingerprint <- function(resources, pairs, config, markers = NULL) {
  config$created_at <- NULL
  if (!is.null(markers)) {
    resources$allele_probs <- lapply(resources$allele_probs, function(a) {
      a[, , dimnames(a)[[3L]] %in% markers, drop = FALSE]
    })
  }
  env <- environment(qtl2mr_batch_fingerprint)
  fn <- sort(ls(env, all.names = TRUE))
  fn <- fn[vapply(fn, function(nm) is.function(get(nm, env, inherits = FALSE)), logical(1))]
  # Deparsed text, not language objects: byte compilation and source references
  # can change a closure's serialized form without changing its code.
  implementation <- lapply(fn, function(nm) {
    c(nm, deparse(get(nm, env, inherits = FALSE), control = c("keepInteger", "keepNA")))
  })
  selected_files <- if ("selected_loci_file" %in% names(pairs)) lapply(pairs$selected_loci_file, function(path) {
    if (!qtl2mr_nonempty_string(path) || !file.exists(path)) return(list(path = path, missing = TRUE))
    list(path = path, hash = unname(tools::md5sum(path)))
  }) else NULL
  list(schema = "qtl2mr-batch-v050-1", data = qtl2mr_fingerprint_object(resources),
       request = qtl2mr_fingerprint_object(list(pairs, config, selected_files)),
       implementation = qtl2mr_fingerprint_object(implementation))
}

qtl2mr_resume_abort <- function(message) {
  stop(structure(list(message = message, call = NULL),
    class = c("qtl2mr_resume_error", "error", "condition")))
}

qtl2mr_batch_owned_paths <- function(output_dir) {
  if (is.null(output_dir)) return(character())
  tables <- c("main_results", "method_comparison", "selected_loci", "overlap_qc",
    "first_stage_summary", "first_stage_filter_results", "mr_block_cache_manifest",
    "warnings", "errors", "run_manifest", "sample_qc", "sample_exclusions", "probability_normalization", "input_alignment",
    "exact_group_plan", "sample_groups", "selected_loci_summary", "workload_summary")
  paths <- file.path(output_dir, c(paste0("batch_", tables, ".csv"),
    "batch_config.json", "batch_resume_state.rds", "batch_pair_records", "mr_block_cache", "by_exposure"))
  links <- Sys.readlink(paths)
  paths[file.exists(paths) | (!is.na(links) & nzchar(links))]
}

qtl2mr_check_batch_destination <- function(output_dir, resume = FALSE, force = FALSE) {
  if (is.null(output_dir)) return(invisible(TRUE))
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be NULL or one directory path", call. = FALSE)
  }
  if (file.exists(output_dir) && !dir.exists(output_dir)) stop("output_dir exists and is not a directory", call. = FALSE)
  occupied <- length(qtl2mr_batch_owned_paths(output_dir)) > 0L
  if (occupied && !isTRUE(resume) && !isTRUE(force)) {
    stop("Batch output already exists; use a new output_dir, resume = TRUE for the identical run, or force = TRUE to explicitly overwrite", call. = FALSE)
  }
  invisible(TRUE)
}

qtl2mr_empty_dataframe <- function(columns) {
  out <- stats::setNames(lapply(columns, function(x) logical()), columns)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

qtl2mr_empty_batch_table <- function(name) {
  columns <- switch(
    name,
    main_results = c(
      "pair_id", "exposure", "outcome", "instrument_trait", "status", "error_message",
      "lmm_backend", "locus_fit_status",
      "n_loci", "n_loci_selected", "n_loci_dropped_nonestimable",
      "n_x_median", "n_y_median", "n_overlap_median",
      "theta_so", "se_so", "ci_lower_so", "ci_upper_so",
      "theta_reported", "se_reported", "reported_ci_lower", "reported_ci_upper",
      "reported_ci_type", "reported_estimator", "reported_estimator_source",
      "reported_estimator_policy", "reported_estimator_reason",
      "theta", "se", "profile_ci_lower", "profile_ci_upper", "profile_p",
      "primary_optimizer_status", "primary_bracket_status",
      "primary_bracket_n_expand", "primary_profile_ci_finite", "primary_reliable",
      "numerical_ok", "confidence_set_type", "inference_status", "effect_units",
      "exposure_sd", "outcome_sd", "raw_conversion_available", "theta_raw_units", "se_raw_units",
      "first_order_theta", "first_order_se", "theta_independence",
      "theta_overlap_delta", "theta_two_step_ivw", "theta_weighted_median",
      "theta_weighted_mode", "Q", "Q_chisq_df", "Q_p_chisq", "Q_p_reported",
      "Q_p_method", "cross_cov_model", "standardize_phenotypes", "lod_thresh",
      "exclude_chr", "runtime_sec", "sample_policy"
    ),
    method_comparison = c(
      "pair_id", "exposure", "outcome", "method", "estimate", "se",
      "ci_lower", "ci_upper", "ci_type", "source", "note"
    ),
    selected_loci = c(
      "pair_id", "exposure", "outcome", "instrument_trait", "marker", "chr",
      "pos", "exposure_lod", "selected_by", "lod_thresh", "excluded_chr_rule",
      "locus_rank", "locus_status", "drop_roles", "drop_reason",
      "error_class", "error_message"
    ),
    overlap_qc = c(
      "pair_id", "marker", "chr", "n_x", "n_y", "n_overlap", "overlap_frac_x",
      "overlap_frac_y", "overlap_jaccard", "rho_xy", "rho_xy_star"
    ),
    first_stage_summary = c(
      "pair_id", "marker", "chr", "F_X_block", "F_Y_block",
      "partial_R2_X_block", "partial_R2_Y_block",
      "n_x", "n_y", "df_x", "df_y", "theta", "se", "Q_min", "weight",
      "hsq_x", "hsq_y", "hsq_status_x", "hsq_status_y"
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
      "n_pairs_resumed", "sample_policy", "batch_common_n", "sample_retention_ratio",
      "sample_retention_min", "sample_retention_gate_status", "nonestimable_locus_policy",
      "n_loci_dropped_nonestimable", "parallel_backend", "n_workers",
      "context_cache_enabled", "context_cache_hits", "context_cache_misses",
      "lmm_fit_cache_enabled", "lmm_fit_cache_hits", "lmm_fit_cache_misses",
      "wall_sec",
      "mean_wall_sec_per_pair", "lmm_backend"
    ),
    sample_exclusions = c("pair_id", "exposure", "outcome", "sample_policy", "reason",
      "n_input_samples", "n_excluded", "excluded_ids", "counts_overlap"),
    probability_normalization = c("pair_id", "chr", "marker", "n_normalized", "max_abs_sum_deviation"),
    sample_qc = c(
      "n_total_eligible", "n_batch_common", "sample_retention_ratio",
      "sample_retention_min", "sample_retention_gate_status",
      "n_traits", "n_pairs", "n_selected_rows", "n_selected_chr"
    ),
    character()
  )
  qtl2mr_empty_dataframe(columns)
}

qtl2mr_empty_batch <- function(resources, pairs, config) {
  config$sample_policy <- qtl2mr_batch_sample_policy(config$sample_policy)
  config$probability_policy <- resources$probability_policy %||% "strict"
  exclusions <- qtl2mr_batch_sample_exclusions(resources, pairs, config$sample_policy,
    batch_traits = config$batch_traits, sample_ids = config$sample_ids)
  out <- list(
    main_results = qtl2mr_empty_batch_table("main_results"),
    method_comparison = qtl2mr_empty_batch_table("method_comparison"),
    selected_loci = qtl2mr_empty_batch_table("selected_loci"),
    overlap_qc = qtl2mr_empty_batch_table("overlap_qc"),
    first_stage_summary = qtl2mr_empty_batch_table("first_stage_summary"),
    warnings = qtl2mr_batch_covariate_warnings(exclusions),
    errors = qtl2mr_empty_batch_table("errors"),
    run_manifest = qtl2mr_empty_batch_table("run_manifest"),
    sample_qc = qtl2mr_empty_batch_table("sample_qc"),
    sample_exclusions = exclusions,
    probability_normalization = qtl2mr_empty_batch_table("probability_normalization"),
    input_alignment = resources$sample_alignment %||% qtl2mr_input_alignment(resources),
    config = config,
    resources_metadata = resources$metadata,
    pairs = pairs
  )
  class(out) <- c("qtl2mr_batch_result", "list")
  out
}

qtl2mr_block_wald_f <- function(beta, Sigma, ridge = 1e-10) {
  beta <- matrix(as.numeric(beta), ncol = 1L)
  k <- nrow(beta)
  if (!k) return(NA_real_)
  out <- tryCatch({
    Sinv <- qtl2mr_mat_inv(as.matrix(Sigma), ridge = ridge)
    as.numeric(t(beta) %*% Sinv %*% beta) / k
  }, error = function(e) NA_real_)
  if (is.finite(out)) out else NA_real_
}

qtl2mr_partial_r2_from_block_f <- function(F_value, k, df) {
  if (!is.finite(F_value) || !is.finite(k) || !is.finite(df) || k <= 0 || df <= 0) {
    return(NA_real_)
  }
  stat <- F_value * k
  stat / (stat + df)
}

qtl2mr_strip_locus_block <- function(block,
                                     per_locus_row = NULL,
                                     selected_row = NULL,
                                     ridge = 1e-10) {
  if (is.null(block$beta_x) || is.null(block$beta_y) ||
      is.null(block$Sigma_x) || is.null(block$Sigma_y) || is.null(block$Sigma_xy)) {
    stop("Cannot strip locus block: beta/covariance fields are missing", call. = FALSE)
  }
  k <- length(block$beta_x)
  df_x <- block$fit_x$df %||% NA_real_
  df_y <- block$fit_y$df %||% NA_real_
  F_x <- per_locus_row$F_X_block %||% NA_real_
  if (!is.finite(F_x)) F_x <- qtl2mr_block_wald_f(block$beta_x, block$Sigma_x, ridge = ridge)
  F_y <- qtl2mr_block_wald_f(block$beta_y, block$Sigma_y, ridge = ridge)
  marker <- block$marker %||% selected_row$marker %||% NA_character_
  chr <- block$chr %||% selected_row$chr %||% NA_character_
  out <- list(
    marker = marker,
    chr = chr,
    pos = selected_row$pos %||% NA_real_,
    exposure_lod = selected_row$exposure_lod %||% NA_real_,
    locus_rank = selected_row$locus_rank %||% NA_integer_,
    beta_x = block$beta_x,
    beta_y = block$beta_y,
    Sigma_x = block$Sigma_x,
    Sigma_y = block$Sigma_y,
    Sigma_xy = block$Sigma_xy,
    n_overlap = block$n_overlap %||% NA_integer_,
    O_xy_summary = block$O_xy_summary %||% list(),
    cross_cov_model = block$cross_cov_model %||% NA_character_,
    rho_xy = block$rho_xy %||% NA_real_,
    rho_xy_raw = block$rho_xy_raw %||% NA_real_,
    rho_xy_star = block$rho_xy_star %||% NA_real_,
    joint_psd_valid = block$joint_psd_valid %||% NA,
    joint_psd_status = block$joint_psd_status %||% NA_character_,
    joint_max_cross_singular_value = block$joint_max_cross_singular_value %||% NA_real_,
    raw_space_psd_valid = block$raw_space_psd_valid %||% NA,
    raw_space_psd_status = block$raw_space_psd_status %||% NA_character_,
    raw_space_max_cross_singular_value = block$raw_space_max_cross_singular_value %||% NA_real_,
    star_space_psd_valid = block$star_space_psd_valid %||% NA,
    star_space_psd_status = block$star_space_psd_status %||% NA_character_,
    star_space_max_cross_singular_value = block$star_space_max_cross_singular_value %||% NA_real_,
    founders = block$founders %||% rownames(block$founder_H) %||% NULL,
    contrasts = block$contrasts %||% NULL,
    founder_H = block$founder_H %||% NULL,
    hsq_x = block$hsq_x %||% block$fit_x$hsq %||% NA_real_,
    hsq_y = block$hsq_y %||% block$fit_y$hsq %||% NA_real_,
    hsq_status_x = block$fit_x$hsq_status %||% NA_character_,
    hsq_status_y = block$fit_y$hsq_status %||% NA_character_,
    fit_x = list(
      n = block$fit_x$n %||% block$O_xy_summary$n_x %||% NA_integer_,
      df = df_x,
      sigma2 = block$fit_x$sigma2 %||% NA_real_,
      phenotype_mean = block$fit_x$phenotype_mean %||% NA_real_,
      phenotype_sd = block$fit_x$phenotype_sd %||% NA_real_
    ),
    fit_y = list(
      n = block$fit_y$n %||% block$O_xy_summary$n_y %||% NA_integer_,
      df = df_y,
      sigma2 = block$fit_y$sigma2 %||% NA_real_,
      phenotype_mean = block$fit_y$phenotype_mean %||% NA_real_,
      phenotype_sd = block$fit_y$phenotype_sd %||% NA_real_
    ),
    F_X_block = F_x,
    F_Y_block = F_y,
    partial_R2_X_block = qtl2mr_partial_r2_from_block_f(F_x, k, df_x),
    partial_R2_Y_block = qtl2mr_partial_r2_from_block_f(F_y, k, df_y),
    local_theta = per_locus_row$theta %||% NA_real_,
    local_se = per_locus_row$se %||% NA_real_,
    local_Q_min = per_locus_row$Q_min %||% NA_real_,
    local_weight = per_locus_row$weight %||% NA_real_
  )
  class(out) <- c("qtl2mr_stripped_locus_block", "list")
  out
}

# Per-locus first-stage summaries (block F, partial R2, hsq) for the batch tables.
qtl2mr_make_mr_block_cache <- function(pair, pilot, ridge = 1e-10) {
  blocks <- pilot$locus_blocks
  if (!is.list(blocks) || !length(blocks)) stop("pilot does not contain locus_blocks", call. = FALSE)
  per_locus <- pilot$diagnostic_layer$per_locus %||% data.frame()
  selected <- pilot$selected_loci %||% data.frame()
  stripped <- lapply(seq_along(blocks), function(i) {
    block <- blocks[[i]]
    marker <- block$marker %||% names(blocks)[[i]] %||% NA_character_
    per_row <- if (nrow(per_locus) && "marker" %in% names(per_locus)) {
      per_locus[match(marker, per_locus$marker), , drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(per_row) != 1L || is.na(per_row$marker[[1L]])) per_row <- data.frame()
    selected_row <- if (nrow(selected) && "marker" %in% names(selected)) {
      selected[match(marker, selected$marker), , drop = FALSE]
    } else {
      data.frame()
    }
    if (nrow(selected_row) != 1L || is.na(selected_row$marker[[1L]])) selected_row <- data.frame()
    qtl2mr_strip_locus_block(block, per_locus_row = per_row, selected_row = selected_row, ridge = ridge)
  })
  names(stripped) <- vapply(stripped, function(x) as.character(x$marker %||% NA_character_), character(1))
  if (anyNA(names(stripped)) || any(!nzchar(names(stripped)))) names(stripped) <- paste0("locus_", seq_along(stripped))
  first_stage <- qtl2mr_mr_block_cache_first_stage_table(stripped)
  out <- list(
    pair_id = pair$pair_id,
    exposure = pair$exposure,
    outcome = pair$outcome,
    instrument_trait = pair$instrument_trait,
    cache_level = "stripped",
    created_at = as.character(Sys.time()),
    summary = pilot$summary,
    selected_loci = selected,
    overlap_summary = pilot$overlap_summary %||% data.frame(),
    first_stage = first_stage,
    blocks = stripped
  )
  class(out) <- c("qtl2mr_mr_block_cache", "list")
  out
}

qtl2mr_mr_block_cache_first_stage_table <- function(cache_or_blocks) {
  blocks <- if (inherits(cache_or_blocks, "qtl2mr_mr_block_cache")) cache_or_blocks$blocks else cache_or_blocks
  if (!is.list(blocks) || !length(blocks)) return(qtl2mr_empty_batch_table("first_stage_summary"))
  rows <- lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    data.frame(
      marker = b$marker %||% names(blocks)[[i]] %||% paste0("locus_", i),
      chr = b$chr %||% NA_character_,
      F_X_block = b$F_X_block %||% NA_real_,
      F_Y_block = b$F_Y_block %||% NA_real_,
      partial_R2_X_block = b$partial_R2_X_block %||% NA_real_,
      partial_R2_Y_block = b$partial_R2_Y_block %||% NA_real_,
      n_x = b$O_xy_summary$n_x %||% b$fit_x$n %||% NA_integer_,
      n_y = b$O_xy_summary$n_y %||% b$fit_y$n %||% NA_integer_,
      df_x = b$fit_x$df %||% NA_real_,
      df_y = b$fit_y$df %||% NA_real_,
      theta = b$local_theta %||% NA_real_,
      se = b$local_se %||% NA_real_,
      Q_min = b$local_Q_min %||% NA_real_,
      weight = b$local_weight %||% NA_real_,
      hsq_x = as.numeric(b$hsq_x %||% NA_real_),
      hsq_y = as.numeric(b$hsq_y %||% NA_real_),
      hsq_status_x = b$hsq_status_x %||% NA_character_,
      hsq_status_y = b$hsq_status_y %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

qtl2mr_batch_append <- function(batch, name, rows) {
  if (is.null(rows) || !nrow(rows)) return(batch)
  if (identical(name, "warnings") && "warning_id" %in% names(rows)) {
    normalized <- rows$warning_id == "probability_rounding_normalized"
    for (column in intersect(c("message", "recommended_action"), names(rows))) {
      rows[[column]][normalized] <- gsub("data/probability_normalization.csv",
        "batch_probability_normalization.csv", rows[[column]][normalized], fixed = TRUE)
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

qtl2mr_batch_locus_audit <- function(pair, pilot) {
  original <- pilot$selected_loci_original %||% pilot$selected_loci %||% data.frame()
  if (!nrow(original)) return(qtl2mr_empty_batch_table("selected_loci"))
  retained <- as.character((pilot$selected_loci %||% data.frame())$marker %||% character())
  dropped <- pilot$dropped_loci %||% data.frame()
  collapse_for_marker <- function(marker, column) {
    if (!nrow(dropped) || !column %in% names(dropped)) return("")
    values <- unique(as.character(dropped[dropped$marker == marker, column, drop = TRUE]))
    values <- values[!is.na(values) & nzchar(values)]
    paste(values, collapse = ";")
  }
  original$locus_status <- ifelse(
    original$marker %in% retained,
    "retained",
    "dropped_nonestimable"
  )
  original$drop_roles <- vapply(original$marker, collapse_for_marker, character(1), column = "role")
  original$drop_reason <- vapply(original$marker, collapse_for_marker, character(1), column = "drop_reason")
  original$error_class <- vapply(original$marker, collapse_for_marker, character(1), column = "error_class")
  original$error_message <- vapply(original$marker, collapse_for_marker, character(1), column = "error_message")
  cbind(
    pair_id = pair$pair_id,
    exposure = pair$exposure,
    outcome = pair$outcome,
    instrument_trait = pair$instrument_trait,
    original,
    stringsAsFactors = FALSE
  )
}

qtl2mr_compact_main_row <- function(pair, pilot, status, error_message, runtime_sec,
                                    instrument_rule, sample_policy) {
  s <- pilot$summary
  data.frame(
    pair_id = pair$pair_id,
    exposure = pair$exposure,
    outcome = pair$outcome,
    instrument_trait = pair$instrument_trait,
    status = status,
    error_message = error_message %||% "",
    lmm_backend = s$lmm_backend %||% "not_available",
    locus_fit_status = if (isTRUE((s$n_loci_dropped_nonestimable %||% 0L) > 0L)) {
      "ok_with_dropped_nonestimable_loci"
    } else {
      "ok"
    },
    n_loci = s$n_loci[[1L]],
    n_loci_selected = s$n_loci_selected %||% s$n_loci[[1L]],
    n_loci_dropped_nonestimable = s$n_loci_dropped_nonestimable %||% 0L,
    n_x_median = s$n_x_median %||% NA_real_,
    n_y_median = s$n_y_median %||% NA_real_,
    n_overlap_median = s$n_overlap_median %||% NA_real_,
    theta_so = s$theta_so %||% s$theta[[1L]],
    se_so = s$se_so %||% s$se[[1L]],
    ci_lower_so = s$ci_lower_so %||% s$profile_ci_lower[[1L]],
    ci_upper_so = s$ci_upper_so %||% s$profile_ci_upper[[1L]],
    theta_reported = s$theta_reported %||% s$theta[[1L]],
    se_reported = s$se_reported %||% s$se[[1L]],
    reported_ci_lower = s$reported_ci_lower %||% s$profile_ci_lower[[1L]],
    reported_ci_upper = s$reported_ci_upper %||% s$profile_ci_upper[[1L]],
    reported_ci_type = s$reported_ci_type %||% "profile",
    reported_estimator = s$reported_estimator %||% "Star-space qtl2MR",
    reported_estimator_source = s$reported_estimator_source %||% "primary",
    reported_estimator_policy = s$reported_estimator_policy %||% "diagnostic_only",
    reported_estimator_reason = s$reported_estimator_reason %||% "",
    theta = s$theta[[1L]],
    se = s$se[[1L]],
    profile_ci_lower = s$profile_ci_lower[[1L]],
    profile_ci_upper = s$profile_ci_upper[[1L]],
    profile_p = s$profile_p[[1L]],
    primary_optimizer_status = s$primary_optimizer_status %||% NA_character_,
    primary_bracket_status = s$primary_bracket_status %||% NA_character_,
    primary_bracket_n_expand = s$primary_bracket_n_expand %||% NA_integer_,
    primary_profile_ci_finite = s$primary_profile_ci_finite %||% NA,
    primary_reliable = s$primary_reliable %||% NA,
    numerical_ok = s$numerical_ok %||% NA,
    confidence_set_type = s$confidence_set_type %||% "unresolved",
    inference_status = s$inference_status %||% "not_calibrated",
    effect_units = s$effect_units %||% NA_character_,
    exposure_sd = s$exposure_sd %||% NA_real_,
    outcome_sd = s$outcome_sd %||% NA_real_,
    raw_conversion_available = s$raw_conversion_available %||% FALSE,
    theta_raw_units = s$theta_raw_units %||% NA_real_,
    se_raw_units = s$se_raw_units %||% NA_real_,
    first_order_theta = s$first_order_theta %||% NA_real_,
    first_order_se = s$first_order_se %||% NA_real_,
    theta_independence = s$theta_independence[[1L]],
    theta_overlap_delta = s$theta_overlap_delta[[1L]],
    theta_two_step_ivw = s$theta_two_step_ivw[[1L]],
    theta_weighted_median = s$theta_weighted_median[[1L]],
    theta_weighted_mode = s$theta_weighted_mode[[1L]],
    Q = s$Q[[1L]],
    Q_chisq_df = s$Q_chisq_df[[1L]],
    Q_p_chisq = s$Q_p_chisq[[1L]],
    Q_p_reported = s$Q_p_reported[[1L]],
    Q_p_method = s$Q_p_method[[1L]],
    cross_cov_model = s$cross_cov_model[[1L]],
    standardize_phenotypes = s$standardize_phenotypes[[1L]],
    lod_thresh = instrument_rule$lod_thresh,
    exclude_chr = paste(instrument_rule$exclude_chr %||% character(), collapse = ","),
    runtime_sec = runtime_sec,
    sample_policy = qtl2mr_batch_sample_policy(sample_policy),
    stringsAsFactors = FALSE
  )
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
  utils::write.csv(batch$main_results, file.path(output_dir, "batch_main_results.csv"), row.names = FALSE)
  utils::write.csv(batch$method_comparison, file.path(output_dir, "batch_method_comparison.csv"), row.names = FALSE)
  utils::write.csv(batch$selected_loci, file.path(output_dir, "batch_selected_loci.csv"), row.names = FALSE)
  utils::write.csv(batch$overlap_qc, file.path(output_dir, "batch_overlap_qc.csv"), row.names = FALSE)
  utils::write.csv(batch$first_stage_summary, file.path(output_dir, "batch_first_stage_summary.csv"), row.names = FALSE)
  utils::write.csv(batch$warnings, file.path(output_dir, "batch_warnings.csv"), row.names = FALSE)
  utils::write.csv(batch$errors, file.path(output_dir, "batch_errors.csv"), row.names = FALSE)
  utils::write.csv(batch$run_manifest, file.path(output_dir, "batch_run_manifest.csv"), row.names = FALSE)
  utils::write.csv(batch$sample_qc, file.path(output_dir, "batch_sample_qc.csv"), row.names = FALSE)
  utils::write.csv(batch$sample_exclusions %||% qtl2mr_empty_batch_table("sample_exclusions"),
    file.path(output_dir, "batch_sample_exclusions.csv"), row.names = FALSE)
  utils::write.csv(batch$probability_normalization %||% qtl2mr_empty_batch_table("probability_normalization"),
    file.path(output_dir, "batch_probability_normalization.csv"), row.names = FALSE)
  if (!is.null(batch$input_alignment)) utils::write.csv(batch$input_alignment, file.path(output_dir, "batch_input_alignment.csv"), row.names = FALSE)
  writeLines(qtl2mr_json_value(batch$config), file.path(output_dir, "batch_config.json"))
  tables <- c("main_results", "method_comparison", "selected_loci", "overlap_qc", "first_stage_summary",
    "warnings", "errors", "run_manifest",
    "sample_qc", "sample_exclusions", "probability_normalization",
    if (!is.null(batch$input_alignment)) "input_alignment")
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
  qtl2mr_atomic_save_rds(list(schema = "qtl2mr-batch-v050", fingerprint = fingerprint,
    tables = table_hashes), file.path(output_dir, "batch_resume_state.rds"))
}

qtl2mr_pair_record_path <- function(output_dir, pair_index) {
  file.path(qtl2mr_batch_records_dir(output_dir), sprintf("pair_%06d.rds", as.integer(pair_index)))
}

qtl2mr_write_pair_record <- function(output_dir, pair_index, pair_id, fingerprint, payload) {
  qtl2mr_atomic_save_rds(list(
    schema = "qtl2mr-batch-pair-v050", pair_index = as.integer(pair_index), pair_id = pair_id,
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
  if (!is.list(state) || !identical(state$schema, "qtl2mr-batch-v050") || !identical(state$fingerprint, fingerprint)) {
    qtl2mr_resume_abort("Cannot resume: input, parameters, selected loci or implementation fingerprint changed")
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
        "; remove it or use force = TRUE"))
    })
    index <- record$pair_index %||% NA_integer_
    valid <- is.list(record) && identical(record$schema, "qtl2mr-batch-pair-v050") &&
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
  batch <- qtl2mr_batch_append(batch, "selected_loci", result$selected)
  batch <- qtl2mr_batch_append(batch, "overlap_qc", result$overlap)
  batch <- qtl2mr_batch_append(batch, "first_stage_summary", result$first)
  batch <- qtl2mr_batch_append(batch, "warnings", result$warnings)
  batch <- qtl2mr_batch_append(batch, "errors", result$error)
  qtl2mr_batch_append(batch, "probability_normalization", result$probability)
}

qtl2mr_batch_select_for_pair <- function(resources, pair, instrument_rule) {
  qtl2mr_select_loci_for_request(resources, qtl2mr_resolve_instrument_request(pair, instrument_rule))
}

# Fit one exposure-outcome pair. Errors are returned as error rows; the caller
# applies fail_fast.
qtl2mr_batch_fit_pair <- function(resources, pair, instrument_rule, sample_policy,
                                  on_nonestimable_locus, fit_options, batch_samples,
                                  sample_ids, context_cache, fit_cache) {
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    selected <- qtl2mr_batch_select_for_pair(resources, pair, instrument_rule)
    if (nrow(selected) < fit_options$min_loci) stop("no_selected_loci: selected loci below min_loci", call. = FALSE)
    samples <- if (identical(sample_policy, "batch_common_complete")) batch_samples else
      qtl2mr_make_common_complete_samples(resources, exposure = pair$exposure, outcome = pair$outcome,
        selected_chr = selected$chr, sample_rule = list(mode = "common_complete", sample_ids = sample_ids))
    if (length(samples) < fit_options$min_samples) stop("too_few_samples: common complete samples below min_samples", call. = FALSE)
    prepared <- as_qtl2mr_prepared_input(
      allele_probs = resources$allele_probs, pheno = resources$pheno, covar = resources$covar,
      kinship = resources$kinship, marker_map = resources$marker_map, selected_loci = selected,
      exposure = pair$exposure, outcome = pair$outcome,
      exposure_samples = samples, outcome_samples = samples,
      probability_policy = resources$probability_policy %||% "strict",
      metadata = c(resources$metadata, list(pair_id = pair$pair_id))
    )
    pilot <- run_prepared_qtl2mr(
      prepared, exposure = pair$exposure, outcome = pair$outcome, loci = NULL, write_outputs = FALSE,
      verbose = isTRUE(fit_options$verbose),
      standardize_phenotypes = isTRUE(fit_options$standardize_phenotypes),
      cross_cov_model = fit_options$cross_cov_model,
      heterogeneity_method = fit_options$heterogeneity_method,
      q_sim_n = fit_options$q_sim_n, q_sim_seed = fit_options$q_sim_seed,
      diagnostic_bootstrap_B = fit_options$diagnostic_bootstrap_B,
      diagnostic_bootstrap_seed = fit_options$diagnostic_bootstrap_seed,
      reported_estimator_policy = fit_options$reported_estimator_policy,
      engine = fit_options$engine, context_cache = context_cache, fit_cache = fit_cache,
      drop_nonestimable_loci = identical(on_nonestimable_locus, "drop")
    )
    runtime <- proc.time()[["elapsed"]] - t0
    main <- qtl2mr_compact_main_row(pair, pilot, "ok", "", runtime, instrument_rule, sample_policy)
    methods <- cbind(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
      qtl2mr_method_comparison_table(pilot))
    cache <- qtl2mr_make_mr_block_cache(pair = pair, pilot = pilot)
    warn <- qtl2mr_warning_table(pilot, sample_policy = sample_policy)
    if (nrow(warn)) warn <- cbind(pair_id = pair$pair_id, warn)
    list(
      main = main, methods = methods, selected = qtl2mr_batch_locus_audit(pair, pilot),
      overlap = cbind(pair_id = pair$pair_id, pilot$overlap_summary),
      first = cbind(pair_id = pair$pair_id, qtl2mr_mr_block_cache_first_stage_table(cache)),
      warnings = warn,
      probability = qtl2mr_batch_probability_table(pair, pilot), error = NULL
    )
  }, error = function(e) {
    runtime <- proc.time()[["elapsed"]] - t0
    list(
      main = data.frame(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
        instrument_trait = pair$instrument_trait, status = "error", error_message = conditionMessage(e),
        runtime_sec = runtime, sample_policy = sample_policy, stringsAsFactors = FALSE),
      error = qtl2mr_error_row(pair, "pair", e, runtime_sec = runtime)
    )
  })
}

qtl2mr_batch_worker_failure_rows <- function(pairs, message, sample_policy) {
  list(
    main = data.frame(pair_id = pairs$pair_id, exposure = pairs$exposure, outcome = pairs$outcome,
      instrument_trait = pairs$instrument_trait, status = "error", error_message = message,
      runtime_sec = NA_real_, sample_policy = sample_policy, stringsAsFactors = FALSE),
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
  qtl2mr_batch_finalize_main(batch, pairs)
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

qtl2mr_batch_finalize_main <- function(batch, pairs) {
  main <- batch$main_results
  canonical <- names(qtl2mr_empty_batch_table("main_results"))
  for (nm in setdiff(canonical, names(main))) main[[nm]] <- rep(NA, nrow(main))
  main <- main[, c(canonical, setdiff(names(main), canonical)), drop = FALSE]
  if (nrow(main)) main <- main[order(match(main$pair_id, pairs$pair_id), na.last = TRUE), , drop = FALSE]
  rownames(main) <- NULL
  batch$main_results <- main
  batch
}

run_qtl2mr_batch <- function(resources,
                             pairs,
                             instrument_rule = list(),
                             sample_policy = c("pair_common_complete", "batch_common_complete"),
                             on_nonestimable_locus = c("fail", "drop"),
                             min_sample_retention = 0.85,
                             fit_options = list(),
                             n_workers = 1L,
                             fail_fast = FALSE,
                             resume = FALSE,
                             force = FALSE,
                             output_dir = NULL,
                             progress = TRUE,
                             ...) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  legacy <- list(...)
  if (length(legacy) && (is.null(names(legacy)) || any(is.na(names(legacy)) | !nzchar(names(legacy))))) {
    stop("Arguments passed through ... must be named; only the deprecated sample_rule and output_mode are accepted",
      call. = FALSE)
  }
  unknown <- setdiff(names(legacy), c("sample_rule", "output_mode"))
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if (!is.null(legacy$output_mode)) {
    qtl2mr_deprecated("run_qtl2mr_batch(output_mode)", "it had no effect on results and is ignored")
  }
  sample_ids <- NULL
  if (!is.null(legacy$sample_rule)) {
    if (!missing(sample_policy)) stop("Supply sample_policy or the deprecated sample_rule, not both", call. = FALSE)
    qtl2mr_deprecated("run_qtl2mr_batch(sample_rule)", "use sample_policy = \"pair_common_complete\" or \"batch_common_complete\"")
    # 0.4.x run_qtl2mr_batch() had no retention gate; keep that unless one is requested.
    if (missing(min_sample_retention)) min_sample_retention <- 0
    rule <- legacy$sample_rule
    mode <- if (is.list(rule)) rule$mode %||% "common_complete" else rule
    if (is.list(rule)) {
      bad <- setdiff(names(rule), c("mode", "sample_ids"))
      if (length(bad)) stop("Unknown sample_rule option(s): ", paste(bad, collapse = ", "), call. = FALSE)
      sample_ids <- rule$sample_ids
    }
    if (mode %in% c("common_complete", "pair_common_complete")) sample_policy <- "pair_common_complete"
    else if (identical(mode, "batch_common_complete")) sample_policy <- "batch_common_complete"
    else stop("Unsupported sample_rule mode", call. = FALSE)
  }
  sample_policy <- match.arg(sample_policy)
  on_nonestimable_locus <- match.arg(on_nonestimable_locus)
  fail_fast <- qtl2mr_check_flag(fail_fast, "fail_fast")
  resume <- qtl2mr_check_flag(resume, "resume")
  force <- qtl2mr_check_flag(force, "force")
  progress <- qtl2mr_check_flag(progress, "progress")
  if (!is.numeric(min_sample_retention) || length(min_sample_retention) != 1L || !is.finite(min_sample_retention) ||
      min_sample_retention < 0 || min_sample_retention > 1) {
    stop("min_sample_retention must be one number between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(n_workers) || length(n_workers) != 1L || !is.finite(n_workers) || n_workers < 1 ||
      n_workers != as.integer(n_workers)) {
    stop("n_workers must be a positive integer", call. = FALSE)
  }
  n_workers <- as.integer(n_workers)
  if (!is.null(sample_ids)) qtl2mr_validate_ids(as.character(sample_ids), "sample_rule sample_ids")
  if ("drop_nonestimable_loci" %in% names(fit_options)) {
    qtl2mr_deprecated("fit_options$drop_nonestimable_loci", "use on_nonestimable_locus = \"drop\" or \"fail\"")
    if (isTRUE(fit_options$drop_nonestimable_loci)) on_nonestimable_locus <- "drop"
    fit_options$drop_nonestimable_loci <- NULL
  }
  pairs <- qtl2mr_normalize_pairs(pairs)
  instrument_rule <- qtl2mr_batch_defaults(instrument_rule, list(
    source = "exposure_qtl_scan", lod_thresh = 4, peak_method = "contiguous_lod_island_max",
    exclude_chr = "X", max_loci = NULL))
  fit_options <- qtl2mr_batch_fit_options(fit_options)
  config <- list(
    instrument_rule = instrument_rule,
    sample_policy = sample_policy,
    on_nonestimable_locus = on_nonestimable_locus,
    min_sample_retention = if (identical(sample_policy, "batch_common_complete")) min_sample_retention else NA_real_,
    sample_ids = sample_ids,
    fit_options = fit_options
  )

  # Destination, pair preflight, the sample-retention gate and the checkpoint are
  # all checked before anything in output_dir is written or removed.
  qtl2mr_check_batch_destination(output_dir, resume, force)
  preflight <- qtl2mr_batch_preflight(resources, pairs, instrument_rule, fit_options, fail_fast,
    sample_ids = sample_ids)
  valid <- preflight$valid
  batch <- qtl2mr_empty_batch(resources, pairs, c(config, list(
    batch_traits = unique(c(valid$exposure, valid$outcome)), created_at = as.character(Sys.time()))))
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  t_run <- proc.time()[["elapsed"]]

  batch_samples <- NULL
  if (identical(sample_policy, "batch_common_complete") && nrow(valid)) {
    selected <- qtl2mr_selected_loci_for_pairs(resources, valid, instrument_rule)
    selected_chr <- unique(as.character(selected$chr))
    traits <- unique(c(valid$exposure, valid$outcome))
    eligible <- qtl2mr_make_batch_eligible_samples(resources, selected_chr = selected_chr)
    # A requested sample_ids restriction is the user's choice, not sample loss.
    if (!is.null(sample_ids)) eligible <- intersect(eligible, as.character(sample_ids))
    batch_samples <- qtl2mr_make_batch_common_complete_samples(resources, traits = traits, selected_chr = selected_chr)
    if (!is.null(sample_ids)) {
      if (length(setdiff(sample_ids, batch_samples))) stop("sample_rule sample_ids must be complete for every batch trait and input", call. = FALSE)
      batch_samples <- as.character(sample_ids)
    }
    ratio <- if (length(eligible)) length(batch_samples) / length(eligible) else NA_real_
    gate <- if (is.finite(ratio) && ratio >= min_sample_retention) "pass" else "fail"
    batch$sample_qc <- data.frame(
      n_total_eligible = length(eligible), n_batch_common = length(batch_samples),
      sample_retention_ratio = ratio, sample_retention_min = min_sample_retention,
      sample_retention_gate_status = gate, n_traits = length(traits), n_pairs = nrow(valid),
      n_selected_rows = nrow(selected), n_selected_chr = length(selected_chr), stringsAsFactors = FALSE)
    if (identical(gate, "fail")) {
      stop(sprintf(paste0("batch_common_complete keeps %d of %d eligible samples (retention %.3f), below ",
        "min_sample_retention = %.3f. Use sample_policy = \"pair_common_complete\", or lower ",
        "min_sample_retention to accept the sample loss."), length(batch_samples), length(eligible),
        ratio, min_sample_retention), call. = FALSE)
    }
    if (length(batch_samples) < fit_options$min_samples) {
      stop("too_few_samples: batch common complete samples below min_samples", call. = FALSE)
    }
  }

  fingerprint <- if (is.null(output_dir)) NULL else
    qtl2mr_batch_fingerprint(resources, pairs, config, markers = preflight$markers)
  occupied <- !is.null(output_dir) && length(qtl2mr_batch_owned_paths(output_dir)) > 0L
  records <- list()
  if (occupied && resume && !force) records <- qtl2mr_load_resume_records(output_dir, fingerprint, pairs)
  if (occupied && force) qtl2mr_clear_batch_outputs(output_dir)
  if (!is.null(output_dir)) {
    dir.create(qtl2mr_batch_records_dir(output_dir), recursive = TRUE, showWarnings = FALSE)
    qtl2mr_write_resume_state(output_dir, fingerprint)
  }
  todo <- valid[!valid$pair_id %in% names(records), , drop = FALSE]
  pair_index <- match(todo$pair_id, pairs$pair_id)
  fit_one <- function(i, context_cache, fit_cache) {
    pair <- todo[i, , drop = FALSE]
    result <- qtl2mr_batch_fit_pair(resources, pair, instrument_rule, sample_policy, on_nonestimable_locus,
      fit_options, batch_samples, sample_ids, context_cache, fit_cache)
    if (!is.null(output_dir) && is.null(result$error)) {
      qtl2mr_write_pair_record(output_dir, pair_index[[i]], pair$pair_id, fingerprint, result)
    }
    result
  }
  new_cache <- function(enabled) if (isTRUE(enabled)) qtl2mr_context_cache_new() else NULL
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
      context_cache <- new_cache(fit_options$cache_contexts)
      fit_cache <- new_cache(fit_options$cache_lmm_fits)
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
        paste("parallel worker failed:", message_text), sample_policy)
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
    context_cache <- new_cache(fit_options$cache_contexts)
    fit_cache <- new_cache(fit_options$cache_lmm_fits)
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
    sample_policy = sample_policy,
    batch_common_n = length(batch_samples %||% character()),
    sample_retention_ratio = batch$sample_qc$sample_retention_ratio[1L] %||% NA_real_,
    sample_retention_min = config$min_sample_retention,
    sample_retention_gate_status = if (identical(sample_policy, "batch_common_complete")) batch$sample_qc$sample_retention_gate_status[1L] %||% NA_character_ else "not_applicable",
    nonestimable_locus_policy = if (identical(on_nonestimable_locus, "drop")) "drop_audited" else "fail_closed",
    n_loci_dropped_nonestimable = sum(batch$main_results$n_loci_dropped_nonestimable %||% 0L, na.rm = TRUE),
    parallel_backend = backend, n_workers = n_workers,
    context_cache_enabled = isTRUE(fit_options$cache_contexts),
    context_cache_hits = stat_sum(1L, "hits"), context_cache_misses = stat_sum(1L, "misses"),
    lmm_fit_cache_enabled = isTRUE(fit_options$cache_lmm_fits),
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

qtl2mr_deprecated <- function(what, instead) {
  warning(what, " is deprecated in qtl2mr 0.5.0 and will be removed in 0.6.0; ", instead, ".", call. = FALSE)
}

run_qtl2mr_batch_exact_grouped <- function(resources, pairs, plan = NULL, instrument_rule = list(),
                                           fit_options = list(), output_mode = "compact", fail_fast = FALSE,
                                           n_workers = 1L, output_dir = NULL, progress = TRUE, force = FALSE) {
  qtl2mr_deprecated("run_qtl2mr_batch_exact_grouped()", "use run_qtl2mr_batch(), which reuses fits for shared sample sets automatically")
  if (!is.null(plan)) warning("plan is ignored; selection and sample settings now come from instrument_rule and fit_options (min_samples, min_loci).", call. = FALSE)
  run_qtl2mr_batch(resources, pairs, instrument_rule = instrument_rule, sample_policy = "pair_common_complete",
    on_nonestimable_locus = "fail", fit_options = fit_options, n_workers = n_workers, fail_fast = fail_fast,
    force = force, output_dir = output_dir, progress = progress)
}

run_qtl2mr_batch_fast <- function(resources, pairs, instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
                                  sample_retention_min = 0.85, allow_sample_loss = FALSE, fit_options = list(),
                                  output_dir = NULL, progress = TRUE, fail_fast = FALSE, resume = FALSE,
                                  force = FALSE) {
  qtl2mr_deprecated("run_qtl2mr_batch_fast()",
    "use run_qtl2mr_batch(sample_policy = \"batch_common_complete\", on_nonestimable_locus = \"drop\")")
  drop <- if (is.null(fit_options$drop_nonestimable_loci)) TRUE else isTRUE(fit_options$drop_nonestimable_loci)
  fit_options$drop_nonestimable_loci <- NULL
  run_qtl2mr_batch(resources, pairs, instrument_rule = instrument_rule, sample_policy = "batch_common_complete",
    on_nonestimable_locus = if (drop) "drop" else "fail",
    min_sample_retention = if (isTRUE(allow_sample_loss)) 0 else sample_retention_min,
    fit_options = fit_options, fail_fast = fail_fast, resume = resume, force = force,
    output_dir = output_dir, progress = progress)
}

run_qtl2mr_batch_parallel_by_exposure <- function(resources, pairs, fast = TRUE, n_workers = 2L,
                                                  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
                                                  sample_retention_min = 0.85, allow_sample_loss = FALSE,
                                                  fit_options = list(), output_dir = NULL, progress = TRUE,
                                                  fail_fast = FALSE, resume = FALSE, force = FALSE) {
  qtl2mr_deprecated("run_qtl2mr_batch_parallel_by_exposure()", "use run_qtl2mr_batch(n_workers = ...)")
  drop <- if (is.null(fit_options$drop_nonestimable_loci)) isTRUE(fast) else isTRUE(fit_options$drop_nonestimable_loci)
  fit_options$drop_nonestimable_loci <- NULL
  run_qtl2mr_batch(resources, pairs, instrument_rule = instrument_rule,
    sample_policy = if (isTRUE(fast)) "batch_common_complete" else "pair_common_complete",
    on_nonestimable_locus = if (drop) "drop" else "fail",
    min_sample_retention = if (isTRUE(allow_sample_loss)) 0 else sample_retention_min,
    fit_options = fit_options, n_workers = n_workers, fail_fast = fail_fast, resume = resume,
    force = force, output_dir = output_dir, progress = progress)
}
