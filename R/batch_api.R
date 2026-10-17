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

qtl2mr_subset_allele_probs_to_loci <- function(allele_probs,
                                               selected_loci,
                                               marker_col = "marker",
                                               chr_col = "chr",
                                               drop_empty_chr = TRUE) {
  if (!is.list(allele_probs) || !length(allele_probs)) {
    stop("allele_probs must be a non-empty qtl2-style list", call. = FALSE)
  }
  if (!is.data.frame(selected_loci)) selected_loci <- as.data.frame(selected_loci, stringsAsFactors = FALSE)
  if (!all(c(marker_col, chr_col) %in% names(selected_loci))) {
    stop("selected_loci must contain marker and chr columns", call. = FALSE)
  }
  loci <- data.frame(
    marker = as.character(selected_loci[[marker_col]]),
    chr = as.character(selected_loci[[chr_col]]),
    stringsAsFactors = FALSE
  )
  keep <- !is.na(loci$marker) & !is.na(loci$chr) & nzchar(loci$marker) & nzchar(loci$chr)
  loci <- loci[keep, , drop = FALSE]
  loci <- unique(loci)
  if (!nrow(loci)) stop("selected_loci contains no non-empty marker/chr rows", call. = FALSE)
  missing_chr <- setdiff(unique(loci$chr), names(allele_probs))
  if (length(missing_chr)) {
    stop("allele_probs is missing selected chromosomes: ", paste(missing_chr, collapse = ", "), call. = FALSE)
  }
  out <- allele_probs
  missing_marker <- character()
  kept <- 0L
  for (chr in names(allele_probs)) {
    arr <- allele_probs[[chr]]
    markers <- dimnames(arr)[[3L]]
    wanted <- unique(loci$marker[loci$chr == as.character(chr)])
    if (!length(wanted)) {
      if (isTRUE(drop_empty_chr)) {
        out[[chr]] <- NULL
      } else {
        out[[chr]] <- arr[, , integer(), drop = FALSE]
      }
      next
    }
    idx <- match(wanted, markers)
    if (anyNA(idx)) missing_marker <- c(missing_marker, wanted[is.na(idx)])
    idx <- idx[!is.na(idx)]
    out[[chr]] <- arr[, , idx, drop = FALSE]
    kept <- kept + length(idx)
  }
  if (length(missing_marker)) {
    shown <- utils::head(unique(missing_marker), 10L)
    stop(
      "allele_probs is missing selected markers: ",
      paste(shown, collapse = ", "),
      if (length(unique(missing_marker)) > length(shown)) " ...",
      call. = FALSE
    )
  }
  if (!length(out) || kept == 0L) stop("No selected markers were retained in allele_probs subset", call. = FALSE)
  class(out) <- class(allele_probs)
  attr(out, "qtl2mr_subset") <- list(
    n_selected_rows = nrow(selected_loci),
    n_unique_loci = nrow(loci),
    n_retained_loci = kept,
    n_retained_chr = length(out),
    created_at = as.character(Sys.time())
  )
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

qtl2mr_batch_sample_qc <- function(resources,
                                  pairs,
                                  instrument_rule,
                                  sample_retention_min = 0.85,
                                  allow_sample_loss = FALSE) {
  pairs <- qtl2mr_normalize_pairs(pairs)
  request_cols <- c("instrument_trait", "lod_thresh", "exclude_chr", "max_loci", "selected_loci_file")
  for (nm in setdiff(request_cols, names(pairs))) pairs[[nm]] <- NA
  request_key <- do.call(paste, c(pairs[request_cols], sep = "\r"))
  request_pairs <- pairs[!duplicated(request_key), , drop = FALSE]
  selected <- qtl2mr_selected_loci_for_pairs(resources, request_pairs, instrument_rule)
  selected_chr <- unique(as.character(selected$chr %||% names(resources$allele_probs)))
  traits <- unique(c(as.character(pairs$exposure), as.character(pairs$outcome)))
  eligible <- qtl2mr_make_batch_eligible_samples(resources, selected_chr = selected_chr)
  batch_common <- qtl2mr_make_batch_common_complete_samples(resources, traits = traits, selected_chr = selected_chr)
  n_total <- length(eligible)
  n_common <- length(batch_common)
  ratio <- if (n_total > 0L) n_common / n_total else NA_real_
  gate_pass <- is.finite(ratio) && ratio >= sample_retention_min
  gate_status <- if (gate_pass) "pass" else if (isTRUE(allow_sample_loss)) "override" else "fail"
  data.frame(
    n_total_eligible = n_total,
    n_batch_common = n_common,
    sample_retention_ratio = ratio,
    sample_retention_min = sample_retention_min,
    sample_retention_gate_status = gate_status,
    allow_sample_loss = isTRUE(allow_sample_loss),
    n_traits = length(traits),
    n_pairs = nrow(pairs),
    n_selected_rows = nrow(selected),
    n_selected_chr = length(selected_chr),
    sample_approximation = TRUE,
    stringsAsFactors = FALSE
  )
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

qtl2mr_instrument_request_key <- function(request) {
  paste(
    request$instrument_trait,
    signif(as.numeric(request$lod_thresh), 12),
    paste(request$exclude_chr %||% character(), collapse = ","),
    if (is.null(request$max_loci)) "NA" else as.character(request$max_loci),
    request$selected_loci_file %||% NA_character_,
    sep = "\r"
  )
}

qtl2mr_sample_group_key <- function(exposure, instrument_request_id, sample_ids) {
  paste(exposure, instrument_request_id, paste(sort(as.character(sample_ids)), collapse = "\r"), sep = "\f")
}

qtl2mr_sample_signature_label <- function(sample_ids) {
  ids <- sort(as.character(sample_ids))
  n <- length(ids)
  if (!n) return("n0")
  raw <- paste(ids, collapse = "\r")
  codes <- utf8ToInt(raw)
  checksum <- if (length(codes)) sum((seq_along(codes) %% 997L + 1L) * codes) %% 1000000007 else 0
  paste0("n", n, "_", ids[[1L]], "_", ids[[n]], "_", format(checksum, scientific = FALSE))
}

qtl2mr_collapse_unique <- function(x, max_items = 12L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  shown <- utils::head(x, max_items)
  paste0(paste(shown, collapse = ";"), if (length(x) > length(shown)) ";..." else "")
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

plan_qtl2mr_batch_exact_groups <- function(resources,
                                          pairs,
                                          instrument_rule = list(),
                                          min_samples = 50L,
                                          min_loci = 1L) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  pairs <- qtl2mr_normalize_pairs(pairs)
  instrument_rule <- qtl2mr_batch_defaults(instrument_rule, list(
    source = "exposure_qtl_scan",
    lod_thresh = 4,
    peak_method = "contiguous_lod_island_max",
    exclude_chr = "X",
    max_loci = NULL
  ))
  request_fingerprint <- qtl2mr_batch_fingerprint(resources, pairs,
    list(instrument_rule = instrument_rule, min_samples = min_samples, min_loci = min_loci))
  requests <- lapply(seq_len(nrow(pairs)), function(i) qtl2mr_resolve_instrument_request(pairs[i, , drop = FALSE], instrument_rule))
  request_keys <- vapply(requests, qtl2mr_instrument_request_key, character(1))
  request_ids <- paste0("inst", sprintf("%04d", match(request_keys, unique(request_keys))))
  pairs$instrument_request_id <- request_ids

  unique_idx <- match(unique(request_keys), request_keys)
  instrument_requests <- do.call(rbind, lapply(seq_along(unique_idx), function(i) {
    req <- requests[[unique_idx[[i]]]]
    data.frame(
      instrument_request_id = paste0("inst", sprintf("%04d", i)),
      instrument_trait = req$instrument_trait,
      lod_thresh = req$lod_thresh,
      exclude_chr = paste(req$exclude_chr %||% character(), collapse = ","),
      max_loci = if (is.null(req$max_loci)) NA_integer_ else as.integer(req$max_loci),
      selected_loci_file = req$selected_loci_file,
      request_key = unique(request_keys)[[i]],
      stringsAsFactors = FALSE
    )
  }))

  selected_by_request <- vector("list", nrow(instrument_requests))
  names(selected_by_request) <- instrument_requests$instrument_request_id
  selected_rows <- lapply(seq_len(nrow(instrument_requests)), function(i) {
    req <- list(
      instrument_trait = instrument_requests$instrument_trait[[i]],
      lod_thresh = instrument_requests$lod_thresh[[i]],
      exclude_chr = strsplit(instrument_requests$exclude_chr[[i]], ",", fixed = TRUE)[[1L]],
      max_loci = if (is.na(instrument_requests$max_loci[[i]])) NULL else instrument_requests$max_loci[[i]],
      selected_loci_file = instrument_requests$selected_loci_file[[i]]
    )
    selected <- qtl2mr_select_loci_for_request(resources, req)
    if (nrow(selected) < min_loci) {
      stop("no_selected_loci: selected loci below min_loci for instrument_trait ", req$instrument_trait, call. = FALSE)
    }
    selected_by_request[[instrument_requests$instrument_request_id[[i]]]] <<- selected
    cbind(instrument_request_id = instrument_requests$instrument_request_id[[i]],
          instrument_trait = req$instrument_trait,
          selected)
  })
  selected_columns <- unique(unlist(lapply(selected_rows, names), use.names = FALSE))
  selected_rows <- lapply(selected_rows, function(x) {
    for (nm in setdiff(selected_columns, names(x))) x[[nm]] <- rep(NA, nrow(x))
    x[selected_columns]
  })
  selected_loci <- do.call(rbind, selected_rows)
  rownames(selected_loci) <- NULL

  selected_loci_summary <- do.call(rbind, lapply(names(selected_by_request), function(id) {
    selected <- selected_by_request[[id]]
    req <- instrument_requests[instrument_requests$instrument_request_id == id, , drop = FALSE]
    data.frame(
      instrument_request_id = id,
      instrument_trait = req$instrument_trait,
      n_loci = nrow(selected),
      n_chr = length(unique(as.character(selected$chr))),
      selected_chr = paste(unique(as.character(selected$chr)), collapse = ","),
      top_marker = if (nrow(selected)) selected$marker[[1L]] else NA_character_,
      top_lod = if ("exposure_lod" %in% names(selected) && nrow(selected)) selected$exposure_lod[[1L]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  sample_sets <- list()
  sample_group_map <- new.env(parent = emptyenv())
  pair_plan_rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    selected <- selected_by_request[[pair$instrument_request_id]]
    samples <- qtl2mr_make_common_complete_samples(
      resources,
      exposure = pair$exposure,
      outcome = pair$outcome,
      selected_chr = selected$chr,
      sample_rule = list(mode = "common_complete")
    )
    status <- if (length(samples) < min_samples) "too_few_samples" else "planned"
    group_key <- qtl2mr_sample_group_key(pair$exposure, pair$instrument_request_id, samples)
    if (!exists(group_key, envir = sample_group_map, inherits = FALSE)) {
      group_id <- paste0("grp", sprintf("%05d", length(ls(sample_group_map, all.names = TRUE)) + 1L))
      assign(group_key, group_id, envir = sample_group_map)
      sample_sets[[group_id]] <<- samples
    } else {
      group_id <- get(group_key, envir = sample_group_map, inherits = FALSE)
    }
    data.frame(
      pair_index = i,
      pair_id = pair$pair_id,
      exposure = pair$exposure,
      outcome = pair$outcome,
      instrument_trait = pair$instrument_trait,
      instrument_request_id = pair$instrument_request_id,
      sample_group_id = group_id,
      n_samples = length(samples),
      n_loci = nrow(selected),
      sample_signature = qtl2mr_sample_signature_label(samples),
      status = status,
      stringsAsFactors = FALSE
    )
  })
  pair_plan <- do.call(rbind, pair_plan_rows)
  rownames(pair_plan) <- NULL

  planned_pairs <- pair_plan[pair_plan$status == "planned", , drop = FALSE]
  sample_groups <- do.call(rbind, lapply(unique(planned_pairs$sample_group_id), function(group_id) {
    rows <- planned_pairs[planned_pairs$sample_group_id == group_id, , drop = FALSE]
    data.frame(
      sample_group_id = group_id,
      exposure = rows$exposure[[1L]],
      instrument_trait = rows$instrument_trait[[1L]],
      instrument_request_id = rows$instrument_request_id[[1L]],
      n_pairs = nrow(rows),
      n_samples = rows$n_samples[[1L]],
      n_loci = rows$n_loci[[1L]],
      outcomes = qtl2mr_collapse_unique(rows$outcome),
      pair_ids = qtl2mr_collapse_unique(rows$pair_id),
      sample_signature = rows$sample_signature[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(sample_groups)) sample_groups <- data.frame()
  rownames(sample_groups) <- NULL
  sample_groups <- sample_groups[order(sample_groups$exposure, sample_groups$instrument_request_id, sample_groups$sample_group_id), , drop = FALSE]

  exposure_marker_fits_naive <- sum(planned_pairs$n_loci)
  exposure_marker_fits_grouped <- if (nrow(sample_groups)) sum(sample_groups$n_loci) else 0L
  outcome_fit_keys <- character()
  for (i in seq_len(nrow(planned_pairs))) {
    row <- planned_pairs[i, , drop = FALSE]
    selected <- selected_by_request[[row$instrument_request_id]]
    outcome_fit_keys <- c(outcome_fit_keys, paste(row$outcome, row$sample_group_id, selected$marker, sep = "\r"))
  }
  outcome_marker_fits_unique <- length(unique(outcome_fit_keys))
  workload_summary <- data.frame(
    n_pairs = nrow(pairs),
    n_planned_pairs = nrow(planned_pairs),
    n_error_pairs = sum(pair_plan$status != "planned"),
    n_instrument_requests = nrow(instrument_requests),
    n_exact_sample_groups = nrow(sample_groups),
    n_selected_loci_rows = nrow(selected_loci),
    n_union_loci = nrow(unique(selected_loci[, c("marker", "chr")])),
    exposure_marker_fits_naive = exposure_marker_fits_naive,
    exposure_marker_fits_grouped = exposure_marker_fits_grouped,
    exposure_marker_fit_reduction = if (exposure_marker_fits_naive > 0) 1 - exposure_marker_fits_grouped / exposure_marker_fits_naive else NA_real_,
    outcome_marker_fits_unique = outcome_marker_fits_unique,
    total_marker_fits_estimated = exposure_marker_fits_grouped + outcome_marker_fits_unique,
    sample_approximation = FALSE,
    stringsAsFactors = FALSE
  )

  out <- list(
    pairs = pairs,
    instrument_rule = instrument_rule,
    instrument_requests = instrument_requests,
    selected_loci = selected_loci,
    selected_loci_by_request = selected_by_request,
    selected_loci_summary = selected_loci_summary,
    pair_plan = pair_plan,
    sample_groups = sample_groups,
    sample_sets = sample_sets,
    workload_summary = workload_summary,
    min_samples = min_samples,
    min_loci = min_loci,
    request_fingerprint = request_fingerprint,
    created_at = as.character(Sys.time())
  )
  class(out) <- c("qtl2mr_exact_group_plan", "list")
  out
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
  allowed <- names(defaults)
  if ("standardize_phenotypes" %in% allowed) allowed <- union(allowed, c(
    "save_mr_block_cache", "mr_block_cache_level", "mr_block_cache_compress",
    "run_first_stage_filter_grid", "first_stage_filter_grid", "first_stage_filter_min_loci", "steiger_filter"))
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  if ("source" %in% names(x) && !identical(x$source, "exposure_qtl_scan")) stop("Unsupported instrument_rule source", call. = FALSE)
  if ("peak_method" %in% names(x) && !identical(x$peak_method, "contiguous_lod_island_max")) stop("Unsupported peak_method", call. = FALSE)
  utils::modifyList(defaults, x, keep.null = TRUE)
}

qtl2mr_batch_fit_options <- function(options, drop_nonestimable_loci = FALSE) {
  qtl2mr_batch_defaults(options, list(
    standardize_phenotypes = TRUE, cross_cov_model = "star",
    heterogeneity_method = "none", q_sim_n = 0L, q_sim_seed = 1L,
    diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L,
    reported_estimator_policy = "diagnostic_only", min_samples = 50L,
    min_loci = 1L, engine = "batch", drop_nonestimable_loci = drop_nonestimable_loci,
    cache_contexts = TRUE, cache_lmm_fits = TRUE, save_mr_block_cache = FALSE,
    mr_block_cache_level = "stripped", mr_block_cache_compress = "xz",
    run_first_stage_filter_grid = FALSE,
    first_stage_filter_grid = c(NA_real_, 1, 2, 3, 5, 10),
    first_stage_filter_min_loci = 1L, steiger_filter = FALSE, verbose = FALSE
  ))
}

qtl2mr_batch_sample_policy <- function(sample_rule = NULL) {
  mode <- if (is.list(sample_rule)) sample_rule$mode else sample_rule
  if (identical(mode, "batch_common_complete")) "batch_common_complete" else "pair_common_complete"
}

qtl2mr_batch_sample_exclusions <- function(resources, pairs, sample_rule = NULL) {
  if (!nrow(pairs)) return(qtl2mr_empty_batch_table("sample_exclusions"))
  ids <- resources$sample_ids
  pheno <- resources$pheno
  covar <- resources$covar
  policy <- qtl2mr_batch_sample_policy(sample_rule)
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
  batch_traits <- unique(c(pairs$exposure, pairs$outcome))
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    pair <- pairs[i, , drop = FALSE]
    reasons <- list(
      missing_phenotype_id = missing_pheno,
      nonfinite_exposure = nonfinite_trait(pair$exposure),
      nonfinite_outcome = nonfinite_trait(pair$outcome),
      missing_covar_id = missing_covar,
      incomplete_covar = incomplete_covar
    )
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

qtl2mr_batch_preflight <- function(resources, pairs, instrument_rule, fit_options, fail_fast) {
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
      samples <- qtl2mr_make_common_complete_samples(resources, pair$exposure, pair$outcome, selected$chr)
      if (length(samples) < fit_options$min_samples) stop("too_few_samples: complete samples below min_samples", call. = FALSE)
      NULL
    }, error = function(e) {
      if (isTRUE(fail_fast)) stop(e)
      qtl2mr_error_row(pair, "preflight", e)
    })
  })
  invalid <- !vapply(errors, is.null, logical(1))
  list(pairs = pairs, valid = pairs[!invalid, , drop = FALSE], errors = do.call(rbind, errors))
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
      sample_policy = qtl2mr_batch_sample_policy(batch$config$sample_rule),
      runtime_sec = NA_real_, stringsAsFactors = FALSE)
    batch <- qtl2mr_batch_append(batch, "main_results", main)
    batch <- qtl2mr_batch_append(batch, "errors", errors)
  }
  batch$pairs <- preflight$pairs
  if (!nrow(batch$run_manifest)) batch$run_manifest <- data.frame(created_at = as.character(Sys.time()))
  batch$run_manifest$n_pairs_requested <- nrow(preflight$pairs)
  batch$run_manifest$n_pairs_ok <- sum(batch$main_results$status == "ok")
  batch$run_manifest$n_pairs_error <- sum(batch$main_results$status == "error")
  batch$run_manifest$sample_policy <- qtl2mr_batch_sample_policy(batch$config$sample_rule)
  qtl2mr_order_batch_main(batch, preflight$pairs)
}

qtl2mr_fingerprint_object <- function(object) {
  path <- tempfile("qtl2mr-fingerprint-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, compress = FALSE, version = 2L)
  unname(tools::md5sum(path))
}

qtl2mr_batch_fingerprint <- function(resources, pairs, config) {
  config$created_at <- NULL
  env <- environment(qtl2mr_batch_fingerprint)
  fn <- sort(ls(env, all.names = TRUE))
  fn <- fn[vapply(fn, function(nm) is.function(get(nm, env, inherits = FALSE)), logical(1))]
  implementation <- lapply(fn, function(nm) {
    f <- get(nm, env, inherits = FALSE)
    list(name = nm, formals = formals(f), body = body(f))
  })
  selected_files <- if ("selected_loci_file" %in% names(pairs)) lapply(pairs$selected_loci_file, function(path) {
    if (!qtl2mr_nonempty_string(path) || !file.exists(path)) return(list(path = path, missing = TRUE))
    list(path = path, hash = unname(tools::md5sum(path)))
  }) else NULL
  list(schema = "qtl2mr-batch-v040-1", data = qtl2mr_fingerprint_object(resources),
       request = qtl2mr_fingerprint_object(list(pairs, config, selected_files)),
       implementation = qtl2mr_fingerprint_object(implementation))
}

qtl2mr_resume_abort <- function(message) {
  stop(structure(list(message = message, call = NULL),
    class = c("qtl2mr_resume_error", "error", "condition")))
}

qtl2mr_validate_resume <- function(output_dir, fingerprint) {
  if (!length(qtl2mr_batch_owned_paths(output_dir))) return(invisible(FALSE))
  if (!file.exists(file.path(output_dir, "batch_main_results.csv"))) qtl2mr_resume_abort("Cannot resume: incomplete checkpoint is missing batch_main_results.csv")
  path <- file.path(output_dir, "batch_resume_state.rds")
  if (!file.exists(path)) qtl2mr_resume_abort("Cannot resume: legacy or incomplete output has no compatible resume fingerprint")
  state <- tryCatch(readRDS(path), error = function(e) {
    qtl2mr_resume_abort(paste("Cannot resume: unreadable checkpoint:", conditionMessage(e)))
  })
  if (!is.list(state) || !identical(state$fingerprint, fingerprint)) qtl2mr_resume_abort("Cannot resume: input, parameters, selected loci or implementation fingerprint changed")
  files <- file.path(output_dir, names(state$tables))
  if (!length(files) || !all(file.exists(files)) || !identical(unname(tools::md5sum(files)), unname(state$tables))) qtl2mr_resume_abort("Cannot resume: checkpoint output tables changed or are incomplete")
  if (length(state$artifacts)) {
    if (!all(file.exists(names(state$artifacts))) ||
        !identical(unname(tools::md5sum(names(state$artifacts))), unname(state$artifacts))) {
      qtl2mr_resume_abort("Cannot resume: referenced block-cache artifacts changed or are incomplete (including worker checkpoints)")
    }
  }
  invisible(TRUE)
}

qtl2mr_batch_owned_paths <- function(output_dir) {
  if (is.null(output_dir)) return(character())
  tables <- c("main_results", "method_comparison", "selected_loci", "overlap_qc",
    "first_stage_summary", "first_stage_filter_results", "mr_block_cache_manifest",
    "warnings", "errors", "run_manifest", "sample_qc", "sample_exclusions", "probability_normalization", "input_alignment",
    "exact_group_plan", "sample_groups", "selected_loci_summary", "workload_summary")
  paths <- file.path(output_dir, c(paste0("batch_", tables, ".csv"),
    "batch_config.json", "batch_resume_state.rds", "mr_block_cache", "by_exposure"))
  links <- Sys.readlink(paths)
  paths[file.exists(paths) | (!is.na(links) & nzchar(links))]
}

qtl2mr_check_batch_destination <- function(output_dir, resume = FALSE, force = FALSE) {
  occupied <- length(qtl2mr_batch_owned_paths(output_dir)) > 0L
  if (occupied && !isTRUE(resume) && !isTRUE(force)) {
    stop("Batch output already exists; use a new output_dir, resume = TRUE for the identical run, or force = TRUE to explicitly overwrite", call. = FALSE)
  }
  if (occupied && isTRUE(resume) && !isTRUE(force) &&
      !file.exists(file.path(output_dir, "batch_main_results.csv"))) {
    qtl2mr_resume_abort("Cannot resume: incomplete checkpoint is missing batch_main_results.csv")
  }
  invisible(TRUE)
}

qtl2mr_batch_existing <- function(output_dir) {
  if (is.null(output_dir)) return(NULL)
  path <- file.path(output_dir, "batch_main_results.csv")
  if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else NULL
}

qtl2mr_load_existing_batch_tables <- function(output_dir) {
  if (is.null(output_dir) || !dir.exists(output_dir)) return(NULL)
  read_one <- function(file, table_name) {
    path <- file.path(output_dir, file)
    if (!file.exists(path) || file.info(path)$size == 0) return(qtl2mr_empty_batch_table(table_name))
    utils::read.csv(path, stringsAsFactors = FALSE)
  }
  list(
    main_results = read_one("batch_main_results.csv", "main_results"),
    method_comparison = read_one("batch_method_comparison.csv", "method_comparison"),
    selected_loci = read_one("batch_selected_loci.csv", "selected_loci"),
    overlap_qc = read_one("batch_overlap_qc.csv", "overlap_qc"),
    first_stage_summary = read_one("batch_first_stage_summary.csv", "first_stage_summary"),
    first_stage_filter_results = read_one("batch_first_stage_filter_results.csv", "first_stage_filter_results"),
    mr_block_cache_manifest = read_one("batch_mr_block_cache_manifest.csv", "mr_block_cache_manifest"),
    warnings = read_one("batch_warnings.csv", "warnings"),
    errors = read_one("batch_errors.csv", "errors"),
    probability_normalization = read_one("batch_probability_normalization.csv", "probability_normalization"),
    run_manifest = read_one("batch_run_manifest.csv", "run_manifest"),
    sample_qc = read_one("batch_sample_qc.csv", "sample_qc")
  )
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
      "lmm_backend",
      "locus_fit_status",
      "n_loci", "n_loci_selected", "n_loci_dropped_nonestimable",
      "n_x_median", "n_y_median", "n_overlap_median",
      "theta_so", "se_so", "ci_lower_so", "ci_upper_so",
      "theta_reported", "se_reported", "reported_ci_lower", "reported_ci_upper",
      "reported_ci_type", "reported_estimator", "reported_estimator_source",
      "reported_estimator_policy", "reported_estimator_reason",
      "theta", "se",
      "profile_ci_lower", "profile_ci_upper", "profile_p",
      "primary_optimizer_status", "primary_bracket_status",
      "primary_bracket_n_expand", "primary_profile_ci_finite",
      "primary_reliable", "first_order_theta", "first_order_se",
      "theta_independence",
      "theta_overlap_delta", "theta_two_step_ivw", "theta_weighted_median",
      "theta_weighted_mode", "Q", "Q_chisq_df", "Q_p_chisq", "Q_p_reported",
      "Q_p_method", "cross_cov_model", "standardize_phenotypes", "lod_thresh",
      "exclude_chr", "sample_rule", "sample_policy", "runtime_sec"
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
      "n_x", "n_y", "df_x", "df_y", "theta", "se", "Q_min", "weight"
    ),
    first_stage_filter_results = c(
      "pair_id", "exposure", "outcome", "instrument_trait", "filter_id",
      "F_X_min", "steiger_filter", "n_loci_input", "n_loci_kept",
      "n_loci_dropped_F", "n_loci_dropped_steiger", "theta", "se",
      "profile_ci_lower", "profile_ci_upper", "Q", "Q_p_reported",
      "Q_p_method", "theta_independence", "theta_overlap_delta",
      "F_X_block_min", "F_X_block_median", "F_X_block_max",
      "partial_R2_X_block_median", "partial_R2_Y_block_median",
      "status", "notes", "cache_file"
    ),
    mr_block_cache_manifest = c(
      "pair_id", "exposure", "outcome", "instrument_trait", "cache_file",
      "cache_level", "n_loci", "n_loci_with_finite_F_X",
      "n_loci_F_X_ge_2", "n_loci_steiger_pass", "file_size_bytes",
      "created_at"
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
      "run_id", "created_at", "n_pairs_requested", "n_pairs_ok",
      "n_pairs_error", "output_mode", "diagnostic_level", "sample_mode", "sample_policy",
      "batch_common_n", "context_cache_enabled", "context_cache_hits",
      "context_cache_misses", "context_cache_size", "lmm_fit_cache_enabled",
      "lmm_fit_cache_hits", "lmm_fit_cache_misses", "lmm_fit_cache_size",
      "sample_approximation", "sample_retention_ratio", "sample_retention_min",
      "sample_retention_gate_status", "parallel_strategy", "parallel_backend",
      "n_workers", "n_groups", "wall_sec", "mean_wall_sec_per_pair",
      "nonestimable_locus_policy", "n_loci_dropped_nonestimable", "lmm_backend"
    ),
    sample_exclusions = c("pair_id", "exposure", "outcome", "sample_policy", "reason",
      "n_input_samples", "n_excluded", "excluded_ids", "counts_overlap"),
    probability_normalization = c("pair_id", "chr", "marker", "n_normalized", "max_abs_sum_deviation"),
    sample_qc = c(
      "n_total_eligible", "n_batch_common", "sample_retention_ratio",
      "sample_retention_min", "sample_retention_gate_status",
      "allow_sample_loss", "n_traits", "n_pairs", "n_selected_rows",
      "n_selected_chr", "sample_approximation"
    ),
    character()
  )
  qtl2mr_empty_dataframe(columns)
}

qtl2mr_empty_batch <- function(resources, pairs, config) {
  config$sample_policy <- qtl2mr_batch_sample_policy(config$sample_rule)
  config$probability_policy <- resources$probability_policy %||% "strict"
  exclusions <- qtl2mr_batch_sample_exclusions(resources, pairs, config$sample_rule)
  out <- list(
    main_results = qtl2mr_empty_batch_table("main_results"),
    method_comparison = qtl2mr_empty_batch_table("method_comparison"),
    selected_loci = qtl2mr_empty_batch_table("selected_loci"),
    overlap_qc = qtl2mr_empty_batch_table("overlap_qc"),
    first_stage_summary = qtl2mr_empty_batch_table("first_stage_summary"),
    first_stage_filter_results = qtl2mr_empty_batch_table("first_stage_filter_results"),
    mr_block_cache_manifest = qtl2mr_empty_batch_table("mr_block_cache_manifest"),
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

qtl2mr_make_mr_block_cache <- function(pair,
                                       pilot,
                                       cache_level = c("stripped"),
                                       ridge = 1e-10) {
  cache_level <- match.arg(cache_level)
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
    cache_level = cache_level,
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
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

qtl2mr_write_mr_block_cache <- function(cache,
                                        output_dir,
                                        pair_id = cache$pair_id,
                                        compress = "xz") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(qtl2mr_slug(pair_id), ".rds"))
  saveRDS(cache, path, compress = compress)
  path
}

qtl2mr_cache_manifest_row <- function(cache, path = NA_character_) {
  fs <- qtl2mr_mr_block_cache_first_stage_table(cache)
  file_size <- if (!is.na(path) && file.exists(path)) file.info(path)$size else NA_real_
  data.frame(
    pair_id = cache$pair_id,
    exposure = cache$exposure,
    outcome = cache$outcome,
    instrument_trait = cache$instrument_trait,
    cache_file = path,
    cache_level = cache$cache_level,
    n_loci = length(cache$blocks),
    n_loci_with_finite_F_X = sum(is.finite(fs$F_X_block)),
    n_loci_F_X_ge_2 = sum(is.finite(fs$F_X_block) & fs$F_X_block >= 2),
    n_loci_steiger_pass = sum(is.finite(fs$partial_R2_X_block) & is.finite(fs$partial_R2_Y_block) &
                                fs$partial_R2_X_block >= fs$partial_R2_Y_block),
    file_size_bytes = file_size,
    created_at = cache$created_at %||% as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
}

qtl2mr_existing_file_count <- function(path) {
  if (is.null(path) || !length(path)) return(0L)
  path <- as.character(path)
  path <- path[!is.na(path) & nzchar(path)]
  if (!length(path)) return(0L)
  sum(file.exists(path))
}

qtl2mr_file_size_sum <- function(size) {
  if (is.null(size) || !length(size)) return(0)
  size <- suppressWarnings(as.numeric(size))
  sum(size[is.finite(size)], na.rm = TRUE)
}

qtl2mr_read_mr_block_cache <- function(cache) {
  if (inherits(cache, "qtl2mr_mr_block_cache")) return(cache)
  if (is.character(cache) && length(cache) == 1L && file.exists(cache)) return(readRDS(cache))
  stop("cache must be a qtl2mr_mr_block_cache object or an existing .rds path", call. = FALSE)
}

refit_qtl2mr_mr_block_cache <- function(cache,
                                        F_X_min = NA_real_,
                                        steiger_filter = FALSE,
                                        min_loci = 1L,
                                        allow_cross_cov = TRUE,
                                        ridge = 1e-10,
                                        heterogeneity_method = "none",
                                        q_sim_n = 0L,
                                        q_sim_seed = 1L,
                                        cache_file = NA_character_) {
  cache <- qtl2mr_read_mr_block_cache(cache)
  fs <- qtl2mr_mr_block_cache_first_stage_table(cache)
  keep <- rep(TRUE, nrow(fs))
  dropped_F <- 0L
  dropped_steiger <- 0L
  filter_label <- if (is.finite(F_X_min)) paste0("F_X_ge_", F_X_min) else "all_loci"
  if (is.finite(F_X_min)) {
    by_F <- is.finite(fs$F_X_block) & fs$F_X_block >= F_X_min
    dropped_F <- sum(!by_F)
    keep <- keep & by_F
  }
  if (isTRUE(steiger_filter)) {
    by_steiger <- is.finite(fs$partial_R2_X_block) & is.finite(fs$partial_R2_Y_block) &
      fs$partial_R2_X_block >= fs$partial_R2_Y_block
    dropped_steiger <- sum(keep & !by_steiger)
    keep <- keep & by_steiger
    filter_label <- paste0(filter_label, "_steiger")
  }
  n_kept <- sum(keep)
  base <- data.frame(
    pair_id = cache$pair_id,
    exposure = cache$exposure,
    outcome = cache$outcome,
    instrument_trait = cache$instrument_trait,
    filter_id = filter_label,
    F_X_min = if (is.finite(F_X_min)) F_X_min else NA_real_,
    steiger_filter = isTRUE(steiger_filter),
    n_loci_input = nrow(fs),
    n_loci_kept = n_kept,
    n_loci_dropped_F = dropped_F,
    n_loci_dropped_steiger = dropped_steiger,
    stringsAsFactors = FALSE
  )
  if (n_kept < min_loci) {
    return(cbind(base, data.frame(
      theta = NA_real_, se = NA_real_, profile_ci_lower = NA_real_,
      profile_ci_upper = NA_real_, Q = NA_real_, Q_p_reported = NA_real_,
      Q_p_method = heterogeneity_method, theta_independence = NA_real_,
      theta_overlap_delta = NA_real_, F_X_block_min = NA_real_,
      F_X_block_median = NA_real_, F_X_block_max = NA_real_,
      partial_R2_X_block_median = NA_real_, partial_R2_Y_block_median = NA_real_,
      status = "too_few_loci", notes = "Filter left fewer loci than min_loci.",
      cache_file = cache_file, stringsAsFactors = FALSE
    )))
  }
  blocks <- cache$blocks[keep]
  fit <- tryCatch(
    fit_multilocus_founder_mr(
      blocks,
      allow_cross_cov = allow_cross_cov,
      ridge = ridge,
      compute_independence = TRUE,
      heterogeneity_method = heterogeneity_method,
      q_sim_n = q_sim_n,
      q_sim_seed = q_sim_seed
    ),
    error = function(e) e
  )
  f_keep <- fs$F_X_block[keep]
  r2x <- fs$partial_R2_X_block[keep]
  r2y <- fs$partial_R2_Y_block[keep]
  if (inherits(fit, "error")) {
    return(cbind(base, data.frame(
      theta = NA_real_, se = NA_real_, profile_ci_lower = NA_real_,
      profile_ci_upper = NA_real_, Q = NA_real_, Q_p_reported = NA_real_,
      Q_p_method = heterogeneity_method, theta_independence = NA_real_,
      theta_overlap_delta = NA_real_, F_X_block_min = min(f_keep, na.rm = TRUE),
      F_X_block_median = stats::median(f_keep, na.rm = TRUE),
      F_X_block_max = max(f_keep, na.rm = TRUE),
      partial_R2_X_block_median = stats::median(r2x, na.rm = TRUE),
      partial_R2_Y_block_median = stats::median(r2y, na.rm = TRUE),
      status = "failed", notes = conditionMessage(fit),
      cache_file = cache_file, stringsAsFactors = FALSE
    )))
  }
  cbind(base, data.frame(
    theta = fit$theta,
    se = fit$se,
    profile_ci_lower = fit$profile_ci_lower,
    profile_ci_upper = fit$profile_ci_upper,
    Q = fit$Q,
    Q_p_reported = fit$Q_p_reported %||% NA_real_,
    Q_p_method = fit$Q_p_method %||% heterogeneity_method,
    theta_independence = fit$independence$theta %||% NA_real_,
    theta_overlap_delta = fit$theta_overlap_delta %||% NA_real_,
    F_X_block_min = min(f_keep, na.rm = TRUE),
    F_X_block_median = stats::median(f_keep, na.rm = TRUE),
    F_X_block_max = max(f_keep, na.rm = TRUE),
    partial_R2_X_block_median = stats::median(r2x, na.rm = TRUE),
    partial_R2_Y_block_median = stats::median(r2y, na.rm = TRUE),
    status = "ok",
    notes = "",
    cache_file = cache_file,
    stringsAsFactors = FALSE
  ))
}

refit_qtl2mr_mr_block_cache_grid <- function(cache,
                                             F_X_min_grid = c(NA_real_, 1, 2, 3, 5, 10),
                                             steiger_filter = FALSE,
                                             min_loci = 1L,
                                             cache_file = NULL,
                                             ...) {
  cache_file <- cache_file %||% if (is.character(cache) && length(cache) == 1L) cache else NA_character_
  cache <- qtl2mr_read_mr_block_cache(cache)
  rows <- lapply(F_X_min_grid, function(fmin) {
    refit_qtl2mr_mr_block_cache(
      cache,
      F_X_min = as.numeric(fmin),
      steiger_filter = steiger_filter,
      min_loci = min_loci,
      cache_file = cache_file,
      ...
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
                                    instrument_rule, sample_rule) {
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
    sample_rule = if (is.list(sample_rule)) sample_rule$mode %||% "common_complete" else sample_rule,
    sample_policy = qtl2mr_batch_sample_policy(sample_rule),
    runtime_sec = runtime_sec,
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
  utils::write.csv(batch$first_stage_filter_results, file.path(output_dir, "batch_first_stage_filter_results.csv"), row.names = FALSE)
  utils::write.csv(batch$mr_block_cache_manifest, file.path(output_dir, "batch_mr_block_cache_manifest.csv"), row.names = FALSE)
  utils::write.csv(batch$warnings, file.path(output_dir, "batch_warnings.csv"), row.names = FALSE)
  utils::write.csv(batch$errors, file.path(output_dir, "batch_errors.csv"), row.names = FALSE)
  utils::write.csv(batch$run_manifest, file.path(output_dir, "batch_run_manifest.csv"), row.names = FALSE)
  utils::write.csv(batch$sample_qc, file.path(output_dir, "batch_sample_qc.csv"), row.names = FALSE)
  utils::write.csv(batch$sample_exclusions %||% qtl2mr_empty_batch_table("sample_exclusions"),
    file.path(output_dir, "batch_sample_exclusions.csv"), row.names = FALSE)
  utils::write.csv(batch$probability_normalization %||% qtl2mr_empty_batch_table("probability_normalization"),
    file.path(output_dir, "batch_probability_normalization.csv"), row.names = FALSE)
  if (!is.null(batch$input_alignment)) utils::write.csv(batch$input_alignment, file.path(output_dir, "batch_input_alignment.csv"), row.names = FALSE)
  if (!is.null(batch$exact_group_plan)) {
    utils::write.csv(batch$exact_group_plan, file.path(output_dir, "batch_exact_group_plan.csv"), row.names = FALSE)
  }
  if (!is.null(batch$sample_groups)) {
    utils::write.csv(batch$sample_groups, file.path(output_dir, "batch_sample_groups.csv"), row.names = FALSE)
  }
  if (!is.null(batch$selected_loci_summary)) {
    utils::write.csv(batch$selected_loci_summary, file.path(output_dir, "batch_selected_loci_summary.csv"), row.names = FALSE)
  }
  if (!is.null(batch$workload_summary)) {
    utils::write.csv(batch$workload_summary, file.path(output_dir, "batch_workload_summary.csv"), row.names = FALSE)
  }
  writeLines(qtl2mr_json_value(batch$config), file.path(output_dir, "batch_config.json"))
  if (!is.null(batch$resume_fingerprint)) {
    files <- c(list.files(output_dir, pattern = "^batch_.*\\.csv$", full.names = TRUE),
      file.path(output_dir, "batch_config.json"))
    hashes <- tools::md5sum(files)
    names(hashes) <- basename(files)
    artifacts <- c(batch$mr_block_cache_manifest$cache_file %||% character(),
      batch$resume_artifacts %||% character())
    artifacts <- unique(artifacts[!is.na(artifacts) & nzchar(artifacts)])
    if (length(artifacts) && !all(file.exists(artifacts))) stop("Cannot checkpoint: referenced block-cache artifact is missing", call. = FALSE)
    artifact_hashes <- if (length(artifacts)) tools::md5sum(normalizePath(artifacts, mustWork = TRUE)) else character()
    saveRDS(list(fingerprint = batch$resume_fingerprint, tables = hashes, artifacts = artifact_hashes),
            file.path(output_dir, "batch_resume_state.rds"), version = 2L)
  }
  invisible(output_dir)
}

run_qtl2mr_batch_exact_grouped <- function(resources,
                                           pairs,
                                           plan = NULL,
                                           instrument_rule = list(),
                                           fit_options = list(),
                                           output_mode = c("compact", "qc"),
                                           fail_fast = FALSE,
                                           n_workers = 1L,
                                           output_dir = NULL,
                                           progress = TRUE,
                                           force = FALSE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  qtl2mr_check_batch_destination(output_dir, force = force)
  output_mode <- match.arg(output_mode)
  if (as.integer(n_workers) != 1L) {
    stop("run_qtl2mr_batch_exact_grouped() currently runs sequential exact groups; use an external exposure-block runner for parallel execution.", call. = FALSE)
  }
  pairs <- qtl2mr_normalize_pairs(pairs)
  supplied_fit_options <- names(fit_options)
  fit_options <- qtl2mr_batch_fit_options(fit_options)
  if (isTRUE(fit_options$drop_nonestimable_loci)) {
    stop(
      "run_qtl2mr_batch_exact_grouped() is confirmatory and requires drop_nonestimable_loci = FALSE; use run_qtl2mr_batch_fast() for audited screening drops.",
      call. = FALSE
    )
  }
  if (!is.null(plan)) {
    if (!inherits(plan, "qtl2mr_exact_group_plan")) stop("plan must come from plan_qtl2mr_batch_exact_groups()", call. = FALSE)
    if (!length(instrument_rule)) instrument_rule <- plan$instrument_rule
    for (nm in c("min_samples", "min_loci")) if (!nm %in% supplied_fit_options) fit_options[[nm]] <- plan[[nm]]
    expected <- qtl2mr_batch_fingerprint(resources, pairs, list(
      instrument_rule = qtl2mr_batch_defaults(instrument_rule, list(source = "exposure_qtl_scan", lod_thresh = 4,
        peak_method = "contiguous_lod_island_max", exclude_chr = "X", max_loci = NULL)),
      min_samples = fit_options$min_samples, min_loci = fit_options$min_loci))
    if (!identical(plan$request_fingerprint, expected)) stop("Exact-group plan does not match current inputs or requests; rebuild the plan", call. = FALSE)
  }
  preflight <- qtl2mr_batch_preflight(resources, pairs, instrument_rule, fit_options, fail_fast)
  if (!nrow(preflight$valid)) {
    batch <- qtl2mr_append_preflight_errors(qtl2mr_empty_batch(resources, pairs, list()), preflight)
    if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
    return(batch)
  }
  pairs <- preflight$valid
  if (is.null(plan)) {
    plan <- plan_qtl2mr_batch_exact_groups(
      resources = resources,
      pairs = pairs,
      instrument_rule = instrument_rule,
      min_samples = fit_options$min_samples,
      min_loci = fit_options$min_loci
    )
  }
  if (!inherits(plan, "qtl2mr_exact_group_plan")) {
    stop("plan must be NULL or an object returned by plan_qtl2mr_batch_exact_groups()", call. = FALSE)
  }
  config <- list(
    output_mode = output_mode,
    instrument_rule = plan$instrument_rule,
    sample_rule = list(mode = "common_complete"),
    fit_options = fit_options,
    execution_strategy = "exact_sample_grouped",
    created_at = as.character(Sys.time())
  )
  batch <- qtl2mr_empty_batch(resources, plan$pairs, config)
  batch$exact_group_plan <- plan$pair_plan
  batch$sample_groups <- plan$sample_groups
  batch$selected_loci_summary <- plan$selected_loci_summary
  batch$workload_summary <- plan$workload_summary

  invalid <- plan$pair_plan[plan$pair_plan$status != "planned", , drop = FALSE]
  if (nrow(invalid) && isTRUE(fail_fast)) {
    stop("Exact grouped plan contains non-planned pairs: ", paste(invalid$pair_id, collapse = ", "), call. = FALSE)
  }
  if (nrow(invalid)) {
    for (i in seq_len(nrow(invalid))) {
      row <- invalid[i, , drop = FALSE]
      main <- data.frame(
        pair_id = row$pair_id, exposure = row$exposure, outcome = row$outcome,
        instrument_trait = row$instrument_trait, status = "error",
        error_message = row$status, runtime_sec = NA_real_,
        stringsAsFactors = FALSE
      )
      err <- data.frame(
        pair_id = row$pair_id, exposure = row$exposure, outcome = row$outcome,
        stage = "plan", error_class = row$status, error_message = row$status,
        traceback_available = FALSE, runtime_sec = NA_real_,
        stringsAsFactors = FALSE
      )
      batch <- qtl2mr_batch_append(batch, "main_results", main)
      batch <- qtl2mr_batch_append(batch, "errors", err)
    }
  }

  context_cache <- if (isTRUE(fit_options$cache_contexts)) qtl2mr_context_cache_new() else NULL
  fit_cache <- if (isTRUE(fit_options$cache_lmm_fits)) qtl2mr_context_cache_new() else NULL
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  t_run <- proc.time()[["elapsed"]]
  planned_groups <- plan$sample_groups
  if (nrow(planned_groups)) {
    for (g in seq_len(nrow(planned_groups))) {
      group <- planned_groups[g, , drop = FALSE]
      group_pairs <- plan$pair_plan[
        plan$pair_plan$sample_group_id == group$sample_group_id & plan$pair_plan$status == "planned",
        ,
        drop = FALSE
      ]
      selected <- plan$selected_loci_by_request[[group$instrument_request_id]]
      common_samples <- plan$sample_sets[[group$sample_group_id]]
      if (isTRUE(progress)) {
        message(sprintf(
          "[qtl2MR grouped] group %d/%d %s exposure=%s pairs=%d samples=%d loci=%d",
          g, nrow(planned_groups), group$sample_group_id, group$exposure,
          nrow(group_pairs), length(common_samples), nrow(selected)
        ))
      }
      for (j in seq_len(nrow(group_pairs))) {
        row <- group_pairs[j, , drop = FALSE]
        pair <- plan$pairs[plan$pairs$pair_id == row$pair_id, , drop = FALSE]
        t0 <- proc.time()[["elapsed"]]
        pair_result <- tryCatch({
          prepared <- as_qtl2mr_prepared_input(
            allele_probs = resources$allele_probs,
            pheno = resources$pheno,
            covar = resources$covar,
            kinship = resources$kinship,
            marker_map = resources$marker_map,
            selected_loci = selected,
            exposure = pair$exposure,
            outcome = pair$outcome,
            exposure_samples = common_samples,
            outcome_samples = common_samples,
            probability_policy = resources$probability_policy %||% "strict",
            metadata = c(resources$metadata, list(pair_id = pair$pair_id, sample_group_id = group$sample_group_id))
          )
          pilot <- run_prepared_qtl2mr(
            prepared,
            exposure = pair$exposure,
            outcome = pair$outcome,
            loci = NULL,
            write_outputs = FALSE,
            verbose = isTRUE(fit_options$verbose),
            standardize_phenotypes = isTRUE(fit_options$standardize_phenotypes),
            cross_cov_model = fit_options$cross_cov_model,
            heterogeneity_method = fit_options$heterogeneity_method,
            q_sim_n = fit_options$q_sim_n,
            q_sim_seed = fit_options$q_sim_seed,
            diagnostic_bootstrap_B = fit_options$diagnostic_bootstrap_B,
            diagnostic_bootstrap_seed = fit_options$diagnostic_bootstrap_seed,
            reported_estimator_policy = fit_options$reported_estimator_policy,
            engine = fit_options$engine,
            context_cache = context_cache,
            fit_cache = fit_cache,
            drop_nonestimable_loci = FALSE
          )
          runtime <- proc.time()[["elapsed"]] - t0
          main <- qtl2mr_compact_main_row(pair, pilot, "ok", "", runtime, plan$instrument_rule, list(mode = "common_complete"))
          methods <- qtl2mr_method_comparison_table(pilot)
          methods <- cbind(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome, methods)
          selected_out <- qtl2mr_batch_locus_audit(pair, pilot)
          overlap <- cbind(pair_id = pair$pair_id, pilot$overlap_summary)
          cache <- qtl2mr_make_mr_block_cache(
            pair = pair,
            pilot = pilot,
            cache_level = fit_options$mr_block_cache_level
          )
          cache_file <- NA_character_
          cache_manifest <- data.frame()
          if (isTRUE(fit_options$save_mr_block_cache) && !is.null(output_dir)) {
            cache_file <- qtl2mr_write_mr_block_cache(
              cache,
              output_dir = file.path(output_dir, "mr_block_cache"),
              pair_id = pair$pair_id,
              compress = fit_options$mr_block_cache_compress
            )
            cache_manifest <- qtl2mr_cache_manifest_row(cache, cache_file)
          } else if (isTRUE(fit_options$save_mr_block_cache) && is.null(output_dir)) {
            cache_manifest <- qtl2mr_cache_manifest_row(cache, NA_character_)
          }
          first <- cbind(pair_id = pair$pair_id, qtl2mr_mr_block_cache_first_stage_table(cache))
          filter <- data.frame()
          if (isTRUE(fit_options$run_first_stage_filter_grid)) {
            filter <- refit_qtl2mr_mr_block_cache_grid(
              cache,
              F_X_min_grid = fit_options$first_stage_filter_grid,
              steiger_filter = isTRUE(fit_options$steiger_filter),
              min_loci = fit_options$first_stage_filter_min_loci,
              allow_cross_cov = TRUE,
              heterogeneity_method = "none",
              q_sim_n = 0L,
              cache_file = cache_file
            )
          }
          warn <- qtl2mr_warning_table(pilot, adapter = "prepared")
          warn <- if (nrow(warn)) cbind(pair_id = pair$pair_id, warn) else warn
          list(
            main = main, methods = methods, selected = selected_out,
            overlap = overlap, first = first, filter = filter,
            cache_manifest = cache_manifest, warnings = warn, error = NULL
          )
        }, error = function(e) {
          runtime <- proc.time()[["elapsed"]] - t0
          if (isTRUE(fail_fast)) stop(e)
          err <- qtl2mr_error_row(pair, "pair", e, runtime_sec = runtime)
          main <- data.frame(
            pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
            instrument_trait = pair$instrument_trait, status = "error",
            error_message = conditionMessage(e), runtime_sec = runtime,
            stringsAsFactors = FALSE
          )
          list(main = main, methods = data.frame(), selected = data.frame(), overlap = data.frame(),
               first = data.frame(), filter = data.frame(), cache_manifest = data.frame(),
               warnings = data.frame(), error = err)
        })
        batch <- qtl2mr_batch_append(batch, "main_results", pair_result$main)
        batch <- qtl2mr_batch_append(batch, "method_comparison", pair_result$methods)
        batch <- qtl2mr_batch_append(batch, "selected_loci", pair_result$selected)
        batch <- qtl2mr_batch_append(batch, "overlap_qc", pair_result$overlap)
        batch <- qtl2mr_batch_append(batch, "first_stage_summary", pair_result$first)
        batch <- qtl2mr_batch_append(batch, "first_stage_filter_results", pair_result$filter)
        batch <- qtl2mr_batch_append(batch, "mr_block_cache_manifest", pair_result$cache_manifest)
        batch <- qtl2mr_batch_append(batch, "warnings", pair_result$warnings)
        batch <- qtl2mr_batch_append(batch, "errors", pair_result$error)
        if (is.null(pair_result$error)) batch <- qtl2mr_batch_append(batch,
          "probability_normalization", qtl2mr_batch_probability_table(pair, pilot))
      }
    }
  }
  wall_sec <- proc.time()[["elapsed"]] - t_run
  context_stats <- qtl2mr_context_cache_stats(context_cache)
  fit_stats <- qtl2mr_context_cache_stats(fit_cache)
  batch <- qtl2mr_order_batch_main(batch, plan$pairs)
  batch$run_manifest <- data.frame(
    run_id = run_id,
    created_at = as.character(Sys.time()),
    n_pairs_requested = nrow(plan$pairs),
    n_pairs_ok = sum(batch$main_results$status == "ok", na.rm = TRUE),
    n_pairs_error = sum(batch$main_results$status == "error", na.rm = TRUE),
    output_mode = output_mode,
    diagnostic_level = output_mode,
    sample_mode = "common_complete",
    batch_common_n = 0L,
    context_cache_enabled = isTRUE(fit_options$cache_contexts),
    context_cache_hits = context_stats$hits,
    context_cache_misses = context_stats$misses,
    context_cache_size = context_stats$size,
    lmm_fit_cache_enabled = isTRUE(fit_options$cache_lmm_fits),
    lmm_fit_cache_hits = fit_stats$hits,
    lmm_fit_cache_misses = fit_stats$misses,
    lmm_fit_cache_size = fit_stats$size,
    sample_approximation = FALSE,
    sample_retention_ratio = NA_real_,
    sample_retention_min = NA_real_,
    sample_retention_gate_status = "not_applicable",
    execution_strategy = "exact_sample_grouped",
    n_instrument_requests = plan$workload_summary$n_instrument_requests[[1L]],
    n_exact_sample_groups = plan$workload_summary$n_exact_sample_groups[[1L]],
    expected_marker_fits = plan$workload_summary$total_marker_fits_estimated[[1L]],
    observed_lmm_fit_misses = fit_stats$misses,
    save_mr_block_cache = isTRUE(fit_options$save_mr_block_cache),
    n_mr_block_cache_files = qtl2mr_existing_file_count(batch$mr_block_cache_manifest$cache_file),
    mr_block_cache_bytes = qtl2mr_file_size_sum(batch$mr_block_cache_manifest$file_size_bytes),
    run_first_stage_filter_grid = isTRUE(fit_options$run_first_stage_filter_grid),
    first_stage_filter_grid = paste(fit_options$first_stage_filter_grid, collapse = ","),
    steiger_filter = isTRUE(fit_options$steiger_filter),
    wall_sec = wall_sec,
    mean_wall_sec_per_pair = wall_sec / max(nrow(plan$pairs), 1L),
    nonestimable_locus_policy = "fail_closed",
    n_loci_dropped_nonestimable = 0L,
    lmm_backend = qtl2mr_backend_label(batch$main_results$lmm_backend[batch$main_results$status == "ok"]),
    stringsAsFactors = FALSE
  )
  batch <- qtl2mr_append_preflight_errors(batch, preflight)
  if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
  batch
}

run_qtl2mr_batch <- function(resources,
                             pairs,
                             instrument_rule = list(),
                             sample_rule = list(mode = "common_complete"),
                             fit_options = list(),
                             output_mode = c("compact", "qc"),
                             fail_fast = FALSE,
                             resume = FALSE,
                             force = FALSE,
                             output_dir = NULL,
                             progress = TRUE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  output_mode <- match.arg(output_mode)
  pairs <- qtl2mr_normalize_pairs(pairs)
  qtl2mr_check_batch_destination(output_dir, resume, force)
  if (is.list(sample_rule)) {
    unknown <- setdiff(names(sample_rule), c("mode", "sample_ids"))
    if (length(unknown)) stop("Unknown sample_rule option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  sample_mode <- if (is.list(sample_rule)) sample_rule$mode %||% "common_complete" else sample_rule
  if (identical(sample_mode, "pair_common_complete")) sample_mode <- "common_complete"
  if (!sample_mode %in% c("common_complete", "batch_common_complete")) stop("Unsupported sample_rule mode", call. = FALSE)
  instrument_rule <- qtl2mr_batch_defaults(instrument_rule, list(
    source = "exposure_qtl_scan",
    lod_thresh = 4,
    peak_method = "contiguous_lod_island_max",
    exclude_chr = "X",
    max_loci = NULL
  ))
  fit_options <- qtl2mr_batch_fit_options(fit_options)
  if (isTRUE(fit_options$drop_nonestimable_loci) && !identical(sample_mode, "batch_common_complete")) {
    stop(
      "drop_nonestimable_loci = TRUE is restricted to explicit batch_common_complete screening; exact pair-specific batches fail closed.",
      call. = FALSE
    )
  }
  config <- list(
    output_mode = output_mode,
    instrument_rule = instrument_rule,
    sample_rule = sample_rule,
    fit_options = fit_options,
    created_at = as.character(Sys.time())
  )
  batch <- qtl2mr_empty_batch(resources, pairs, config)
  batch$resume_fingerprint <- if (!is.null(output_dir)) qtl2mr_batch_fingerprint(resources, pairs, config) else NULL
  if (isTRUE(resume) && !isTRUE(force)) qtl2mr_validate_resume(output_dir, batch$resume_fingerprint)
  existing <- if (isTRUE(resume) && !isTRUE(force)) qtl2mr_load_existing_batch_tables(output_dir) else NULL
  if (!is.null(existing) && anyDuplicated(existing$main_results$pair_id)) stop("Cannot resume: duplicate pair_id in checkpoint", call. = FALSE)
  done_ids <- if (!is.null(existing) && nrow(existing$main_results)) existing$main_results$pair_id[existing$main_results$status == "ok"] else character()
  if (length(done_ids)) {
    for (nm in intersect(names(existing), names(batch))) {
      tab <- existing[[nm]]
      if ("pair_id" %in% names(tab)) tab <- tab[tab$pair_id %in% done_ids, , drop = FALSE]
      batch[[nm]] <- tab
    }
    pending_covariate_warnings <- qtl2mr_batch_covariate_warnings(batch$sample_exclusions)
    pending_covariate_warnings <- pending_covariate_warnings[
      !pending_covariate_warnings$pair_id %in% done_ids, , drop = FALSE]
    batch <- qtl2mr_batch_append(batch, "warnings", pending_covariate_warnings)
  }

  batch_common_samples <- NULL
  if (identical(sample_mode, "batch_common_complete")) {
    batch_common_samples <- if (is.list(sample_rule) && !is.null(sample_rule$sample_ids)) {
      requested <- as.character(sample_rule$sample_ids)
      qtl2mr_validate_ids(requested, "sample_rule sample_ids")
      complete <- qtl2mr_make_batch_common_complete_samples(resources,
        traits = unique(c(pairs$exposure, pairs$outcome)), selected_chr = names(resources$allele_probs))
      if (length(setdiff(requested, complete))) stop("sample_rule sample_ids must be complete for every batch trait and input", call. = FALSE)
      requested
    } else {
      qtl2mr_make_batch_common_complete_samples(
        resources,
        traits = unique(c(pairs$exposure, pairs$outcome)),
        selected_chr = names(resources$allele_probs)
      )
    }
    if (length(batch_common_samples) < fit_options$min_samples) {
      stop("too_few_samples: batch common complete samples below min_samples", call. = FALSE)
    }
  }
  context_cache <- if (isTRUE(fit_options$cache_contexts)) qtl2mr_context_cache_new() else NULL
  fit_cache <- if (isTRUE(fit_options$cache_lmm_fits)) qtl2mr_context_cache_new() else NULL
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, , drop = FALSE]
    if (pair$pair_id %in% done_ids) next
    t0 <- proc.time()[["elapsed"]]
    if (isTRUE(progress)) message(sprintf("[qtl2MR batch] %d/%d %s: %s -> %s", i, nrow(pairs), pair$pair_id, pair$exposure, pair$outcome))
    pair_result <- tryCatch({
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
      if (nrow(selected) < fit_options$min_loci) {
        stop("no_selected_loci: selected loci below min_loci", call. = FALSE)
      }
      common_samples <- if (identical(sample_mode, "batch_common_complete")) {
        batch_common_samples
      } else {
        qtl2mr_make_common_complete_samples(
          resources,
          exposure = pair$exposure,
          outcome = pair$outcome,
          selected_chr = selected$chr,
          sample_rule = sample_rule
        )
      }
      if (length(common_samples) < fit_options$min_samples) {
        stop("too_few_samples: common complete samples below min_samples", call. = FALSE)
      }
      prepared <- as_qtl2mr_prepared_input(
        allele_probs = resources$allele_probs,
        pheno = resources$pheno,
        covar = resources$covar,
        kinship = resources$kinship,
        marker_map = resources$marker_map,
        selected_loci = selected,
        exposure = pair$exposure,
        outcome = pair$outcome,
        exposure_samples = common_samples,
        outcome_samples = common_samples,
        probability_policy = resources$probability_policy %||% "strict",
        metadata = c(resources$metadata, list(pair_id = pair$pair_id))
      )
      pilot <- run_prepared_qtl2mr(
        prepared,
        exposure = pair$exposure,
        outcome = pair$outcome,
        loci = NULL,
        write_outputs = FALSE,
        verbose = isTRUE(fit_options$verbose),
        standardize_phenotypes = isTRUE(fit_options$standardize_phenotypes),
        cross_cov_model = fit_options$cross_cov_model,
        heterogeneity_method = fit_options$heterogeneity_method,
        q_sim_n = fit_options$q_sim_n,
        q_sim_seed = fit_options$q_sim_seed,
        diagnostic_bootstrap_B = fit_options$diagnostic_bootstrap_B,
        diagnostic_bootstrap_seed = fit_options$diagnostic_bootstrap_seed,
        reported_estimator_policy = fit_options$reported_estimator_policy,
        engine = fit_options$engine,
        context_cache = context_cache,
        fit_cache = fit_cache,
        drop_nonestimable_loci = isTRUE(fit_options$drop_nonestimable_loci)
      )
      runtime <- proc.time()[["elapsed"]] - t0
      main <- qtl2mr_compact_main_row(pair, pilot, "ok", "", runtime, instrument_rule, sample_rule)
      methods <- qtl2mr_method_comparison_table(pilot)
      methods <- cbind(pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome, methods)
      selected_out <- qtl2mr_batch_locus_audit(pair, pilot)
      overlap <- cbind(pair_id = pair$pair_id, pilot$overlap_summary)
      cache <- qtl2mr_make_mr_block_cache(
        pair = pair,
        pilot = pilot,
        cache_level = fit_options$mr_block_cache_level
      )
      cache_file <- NA_character_
      cache_manifest <- data.frame()
      if (isTRUE(fit_options$save_mr_block_cache) && !is.null(output_dir)) {
        cache_file <- qtl2mr_write_mr_block_cache(
          cache,
          output_dir = file.path(output_dir, "mr_block_cache"),
          pair_id = pair$pair_id,
          compress = fit_options$mr_block_cache_compress
        )
        cache_manifest <- qtl2mr_cache_manifest_row(cache, cache_file)
      } else if (isTRUE(fit_options$save_mr_block_cache) && is.null(output_dir)) {
        cache_manifest <- qtl2mr_cache_manifest_row(cache, NA_character_)
      }
      first <- cbind(pair_id = pair$pair_id, qtl2mr_mr_block_cache_first_stage_table(cache))
      filter <- data.frame()
      if (isTRUE(fit_options$run_first_stage_filter_grid)) {
        filter <- refit_qtl2mr_mr_block_cache_grid(
          cache,
          F_X_min_grid = fit_options$first_stage_filter_grid,
          steiger_filter = isTRUE(fit_options$steiger_filter),
          min_loci = fit_options$first_stage_filter_min_loci,
          allow_cross_cov = TRUE,
          heterogeneity_method = "none",
          q_sim_n = 0L,
          cache_file = cache_file
        )
      }
      warn <- qtl2mr_warning_table(pilot, adapter = "prepared")
      warn <- if (nrow(warn)) cbind(pair_id = pair$pair_id, warn) else warn
      list(
        main = main, methods = methods, selected = selected_out,
        overlap = overlap, first = first, filter = filter,
        cache_manifest = cache_manifest, warnings = warn, error = NULL
      )
    }, error = function(e) {
      runtime <- proc.time()[["elapsed"]] - t0
      if (isTRUE(fail_fast)) stop(e)
      err <- qtl2mr_error_row(pair, "pair", e, runtime_sec = runtime)
      main <- data.frame(
        pair_id = pair$pair_id, exposure = pair$exposure, outcome = pair$outcome,
        instrument_trait = pair$instrument_trait, status = "error",
        sample_policy = qtl2mr_batch_sample_policy(sample_rule),
        error_message = conditionMessage(e), runtime_sec = runtime,
        stringsAsFactors = FALSE
      )
      list(main = main, methods = data.frame(), selected = data.frame(), overlap = data.frame(),
           first = data.frame(), filter = data.frame(), cache_manifest = data.frame(),
           warnings = data.frame(), error = err)
    })
    batch <- qtl2mr_batch_append(batch, "main_results", pair_result$main)
    batch <- qtl2mr_batch_append(batch, "method_comparison", pair_result$methods)
    batch <- qtl2mr_batch_append(batch, "selected_loci", pair_result$selected)
    batch <- qtl2mr_batch_append(batch, "overlap_qc", pair_result$overlap)
    batch <- qtl2mr_batch_append(batch, "first_stage_summary", pair_result$first)
    batch <- qtl2mr_batch_append(batch, "first_stage_filter_results", pair_result$filter)
    batch <- qtl2mr_batch_append(batch, "mr_block_cache_manifest", pair_result$cache_manifest)
    batch <- qtl2mr_batch_append(batch, "warnings", pair_result$warnings)
    batch <- qtl2mr_batch_append(batch, "errors", pair_result$error)
    if (is.null(pair_result$error)) batch <- qtl2mr_batch_append(batch,
      "probability_normalization", qtl2mr_batch_probability_table(pair, pilot))
    context_stats <- qtl2mr_context_cache_stats(context_cache)
    fit_stats <- qtl2mr_context_cache_stats(fit_cache)
    batch$run_manifest <- data.frame(
      run_id = run_id,
      created_at = as.character(Sys.time()),
      n_pairs_requested = nrow(pairs),
      n_pairs_ok = sum(batch$main_results$status == "ok", na.rm = TRUE),
      n_pairs_error = sum(batch$main_results$status == "error", na.rm = TRUE),
      output_mode = output_mode,
      diagnostic_level = output_mode,
      sample_mode = sample_mode,
      sample_policy = qtl2mr_batch_sample_policy(sample_rule),
      batch_common_n = length(batch_common_samples %||% character()),
      context_cache_enabled = isTRUE(fit_options$cache_contexts),
      context_cache_hits = context_stats$hits,
      context_cache_misses = context_stats$misses,
      context_cache_size = context_stats$size,
      lmm_fit_cache_enabled = isTRUE(fit_options$cache_lmm_fits),
      lmm_fit_cache_hits = fit_stats$hits,
      lmm_fit_cache_misses = fit_stats$misses,
      lmm_fit_cache_size = fit_stats$size,
      sample_approximation = identical(sample_mode, "batch_common_complete"),
      sample_retention_ratio = NA_real_,
      sample_retention_min = NA_real_,
      sample_retention_gate_status = if (identical(sample_mode, "batch_common_complete")) "not_checked" else "not_applicable",
      save_mr_block_cache = isTRUE(fit_options$save_mr_block_cache),
      n_mr_block_cache_files = qtl2mr_existing_file_count(batch$mr_block_cache_manifest$cache_file),
      mr_block_cache_bytes = qtl2mr_file_size_sum(batch$mr_block_cache_manifest$file_size_bytes),
      run_first_stage_filter_grid = isTRUE(fit_options$run_first_stage_filter_grid),
      first_stage_filter_grid = paste(fit_options$first_stage_filter_grid, collapse = ","),
      steiger_filter = isTRUE(fit_options$steiger_filter),
      nonestimable_locus_policy = if (isTRUE(fit_options$drop_nonestimable_loci)) "drop_audited_screening" else "fail_closed",
      n_loci_dropped_nonestimable = sum(batch$main_results$n_loci_dropped_nonestimable %||% 0L, na.rm = TRUE),
      lmm_backend = qtl2mr_backend_label(batch$main_results$lmm_backend[batch$main_results$status == "ok"]),
      stringsAsFactors = FALSE
    )
    if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
  }
  batch
}

run_qtl2mr_batch_fast <- function(resources,
                                  pairs,
                                  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
                                  sample_retention_min = 0.85,
                                  allow_sample_loss = FALSE,
                                  fit_options = list(),
                                  output_dir = NULL,
                                  progress = TRUE,
                                  fail_fast = FALSE,
                                  resume = FALSE,
                                  force = FALSE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  qtl2mr_check_batch_destination(output_dir, resume, force)
  pairs <- qtl2mr_normalize_pairs(pairs)
  instrument_rule <- qtl2mr_batch_defaults(instrument_rule, list(
    source = "exposure_qtl_scan",
    lod_thresh = 4,
    peak_method = "contiguous_lod_island_max",
    exclude_chr = "X",
    max_loci = NULL
  ))
  fit_options <- qtl2mr_batch_fit_options(fit_options, drop_nonestimable_loci = TRUE)
  preflight <- qtl2mr_batch_preflight(resources, pairs, instrument_rule, fit_options, fail_fast)
  resumed_checkpoint <- FALSE
  if (isTRUE(resume) && !isTRUE(force) && !is.null(output_dir)) {
    resume_config <- list(output_mode = "compact", instrument_rule = instrument_rule,
      sample_rule = list(mode = "batch_common_complete"), fit_options = fit_options)
    resumed_checkpoint <- qtl2mr_validate_resume(output_dir,
      qtl2mr_batch_fingerprint(resources, preflight$valid, resume_config))
  }
  if (!nrow(preflight$valid)) {
    batch <- qtl2mr_append_preflight_errors(qtl2mr_empty_batch(resources, pairs,
      list(sample_rule = list(mode = "batch_common_complete"))), preflight)
    if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
    return(batch)
  }
  pairs <- preflight$valid
  sample_qc <- qtl2mr_batch_sample_qc(
    resources = resources,
    pairs = pairs,
    instrument_rule = instrument_rule,
    sample_retention_min = sample_retention_min,
    allow_sample_loss = allow_sample_loss
  )
  gate_failed <- identical(sample_qc$sample_retention_gate_status[[1L]], "fail")
  if (gate_failed) {
    if (isTRUE(resumed_checkpoint)) qtl2mr_resume_abort(
      "Cannot resume: current sample-retention gate fails; use a new output_dir for the changed screening request")
    config <- list(
      output_mode = "compact",
      instrument_rule = instrument_rule,
      sample_rule = list(mode = "batch_common_complete"),
      fit_options = fit_options,
      created_at = as.character(Sys.time())
    )
    batch <- qtl2mr_empty_batch(resources, pairs, config)
    batch$sample_qc <- sample_qc
    batch$run_manifest <- data.frame(
      run_id = format(Sys.time(), "%Y%m%d_%H%M%S"),
      created_at = as.character(Sys.time()),
      n_pairs_requested = nrow(pairs),
      n_pairs_ok = 0L,
      n_pairs_error = 0L,
      output_mode = "compact",
      diagnostic_level = "compact",
      sample_mode = "batch_common_complete",
      sample_policy = "batch_common_complete",
      batch_common_n = sample_qc$n_batch_common[[1L]],
      context_cache_enabled = TRUE,
      context_cache_hits = 0L,
      context_cache_misses = 0L,
      context_cache_size = 0L,
      lmm_fit_cache_enabled = TRUE,
      lmm_fit_cache_hits = 0L,
      lmm_fit_cache_misses = 0L,
      lmm_fit_cache_size = 0L,
      sample_approximation = TRUE,
      sample_retention_ratio = sample_qc$sample_retention_ratio[[1L]],
      sample_retention_min = sample_retention_min,
      sample_retention_gate_status = "fail",
      nonestimable_locus_policy = if (isTRUE(fit_options$drop_nonestimable_loci)) "drop_audited_screening" else "fail_closed",
      n_loci_dropped_nonestimable = 0L,
      lmm_backend = "not_available",
      stringsAsFactors = FALSE
    )
    if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
    stop(
      sprintf(
        "batch_common_complete sample retention %.3f is below sample_retention_min %.3f; use run_qtl2mr_batch() for pair-specific complete samples or set allow_sample_loss = TRUE to override.",
        sample_qc$sample_retention_ratio[[1L]], sample_retention_min
      ),
      call. = FALSE
    )
  }
  batch <- run_qtl2mr_batch(
    resources = resources,
    pairs = pairs,
    instrument_rule = instrument_rule,
    sample_rule = list(mode = "batch_common_complete"),
    fit_options = fit_options,
    output_mode = "compact",
    fail_fast = fail_fast,
    resume = resume,
    force = force,
    output_dir = output_dir,
    progress = progress
  )
  batch$sample_qc <- sample_qc
  batch$run_manifest$sample_approximation <- TRUE
  batch$run_manifest$sample_retention_ratio <- sample_qc$sample_retention_ratio[[1L]]
  batch$run_manifest$sample_retention_min <- sample_retention_min
  batch$run_manifest$sample_retention_gate_status <- sample_qc$sample_retention_gate_status[[1L]]
  batch$run_manifest$nonestimable_locus_policy <- if (isTRUE(fit_options$drop_nonestimable_loci)) "drop_audited_screening" else "fail_closed"
  batch$run_manifest$n_loci_dropped_nonestimable <- sum(batch$main_results$n_loci_dropped_nonestimable %||% 0L, na.rm = TRUE)
  if (identical(sample_qc$sample_retention_gate_status[[1L]], "override")) {
    warn <- data.frame(
      pair_id = NA_character_,
      warning_id = "batch_common_sample_retention_low",
      severity = "warning",
      message = sprintf(
        "batch_common_complete sample retention %.3f is below %.3f but allow_sample_loss = TRUE.",
        sample_qc$sample_retention_ratio[[1L]], sample_retention_min
      ),
      consequence = "Fast batch estimates may differ materially from pair-specific complete-case estimates.",
      recommended_action = "Use run_qtl2mr_batch() for confirmatory analyses or review batch_sample_qc.csv.",
      shown_in_main_report = TRUE,
      stringsAsFactors = FALSE
    )
    batch <- qtl2mr_batch_append(batch, "warnings", warn)
  }
  batch <- qtl2mr_append_preflight_errors(batch, preflight)
  if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
  batch
}

qtl2mr_order_batch_main <- function(batch, pairs) {
  if (nrow(batch$main_results) && "pair_id" %in% names(batch$main_results)) {
    ord <- match(batch$main_results$pair_id, pairs$pair_id)
    batch$main_results <- batch$main_results[order(ord, na.last = TRUE), , drop = FALSE]
    rownames(batch$main_results) <- NULL
  }
  batch
}

qtl2mr_merge_batch_group_results <- function(resources,
                                             pairs,
                                             group_results,
                                             config,
                                             run_manifest,
                                             sample_qc = NULL) {
  batch <- qtl2mr_empty_batch(resources, pairs, config)
  for (res in group_results) {
    batch <- qtl2mr_batch_append(batch, "main_results", res$main_results)
    batch <- qtl2mr_batch_append(batch, "method_comparison", res$method_comparison)
    batch <- qtl2mr_batch_append(batch, "selected_loci", res$selected_loci)
    batch <- qtl2mr_batch_append(batch, "overlap_qc", res$overlap_qc)
    batch <- qtl2mr_batch_append(batch, "first_stage_summary", res$first_stage_summary)
    warnings <- res$warnings
    if (nrow(warnings)) warnings <- warnings[
      warnings$warning_id != "batch_covariate_sample_exclusions", , drop = FALSE]
    batch <- qtl2mr_batch_append(batch, "warnings", warnings)
    batch <- qtl2mr_batch_append(batch, "errors", res$errors)
    batch <- qtl2mr_batch_append(batch, "first_stage_filter_results", res$first_stage_filter_results)
    batch <- qtl2mr_batch_append(batch, "mr_block_cache_manifest", res$mr_block_cache_manifest)
    batch <- qtl2mr_batch_append(batch, "probability_normalization", res$probability_normalization)
  }
  batch <- qtl2mr_order_batch_main(batch, pairs)
  batch$run_manifest <- run_manifest
  batch$run_manifest$sample_policy <- rep(qtl2mr_batch_sample_policy(config$sample_rule), nrow(run_manifest))
  if (!is.null(sample_qc)) batch$sample_qc <- sample_qc
  batch
}

qtl2mr_parallel_group_failure <- function(group, error) {
  structure(
    list(
      exposure_group = as.character(group$key),
      n_pairs = nrow(group$pairs),
      resume_integrity_failure = inherits(error, "qtl2mr_resume_error"),
      error_class = if (inherits(error, "try-error")) "try-error" else class(error)[[1L]],
      error_message = if (inherits(error, "try-error")) as.character(error) else conditionMessage(error)
    ),
    class = c("qtl2mr_parallel_group_failure", "list")
  )
}

qtl2mr_collect_parallel_group_results <- function(group_results, groups, fail_fast = FALSE) {
  failures <- vector("list", 0L)
  successes <- vector("list", 0L)
  for (i in seq_along(group_results)) {
    value <- group_results[[i]]
    if (inherits(value, "qtl2mr_parallel_group_failure")) {
      failures[[length(failures) + 1L]] <- value
    } else if (inherits(value, "try-error")) {
      failures[[length(failures) + 1L]] <- qtl2mr_parallel_group_failure(groups[[i]], value)
    } else {
      successes[[length(successes) + 1L]] <- value
    }
  }
  if (length(failures)) {
    failure_table <- do.call(rbind, lapply(failures, function(x) data.frame(
      exposure_group = x$exposure_group,
      n_pairs = as.integer(x$n_pairs),
      error_class = x$error_class,
      error_message = x$error_message,
      stringsAsFactors = FALSE
    )))
    integrity <- vapply(failures, function(x) isTRUE(x$resume_integrity_failure), logical(1))
    if (any(integrity) || isTRUE(fail_fast)) {
      first <- failure_table[if (any(integrity)) which(integrity)[[1L]] else 1L, , drop = FALSE]
      stop(
        sprintf("parallel exposure group '%s' failed: %s", first$exposure_group[[1L]], first$error_message[[1L]]),
        call. = FALSE
      )
    }
  } else {
    failure_table <- data.frame(
      exposure_group = character(), n_pairs = integer(), error_class = character(),
      error_message = character(), stringsAsFactors = FALSE
    )
  }
  list(successes = successes, failures = failure_table)
}

run_qtl2mr_batch_parallel_by_exposure <- function(resources,
                                                  pairs,
                                                  fast = TRUE,
                                                  n_workers = 2L,
                                                  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
                                                  sample_retention_min = 0.85,
                                                  allow_sample_loss = FALSE,
                                                  fit_options = list(),
                                                  output_dir = NULL,
                                                  progress = TRUE,
                                                  fail_fast = FALSE,
                                                  resume = FALSE,
                                                  force = FALSE) {
  if (!inherits(resources, "qtl2mr_resources")) stop("resources must come from qtl2mr_prepare_resources()", call. = FALSE)
  qtl2mr_check_batch_destination(output_dir, resume, force)
  pairs <- qtl2mr_normalize_pairs(pairs)
  instrument_rule <- qtl2mr_batch_defaults(instrument_rule, list(
    source = "exposure_qtl_scan",
    lod_thresh = 4,
    peak_method = "contiguous_lod_island_max",
    exclude_chr = "X",
    max_loci = NULL
  ))
  fit_options <- qtl2mr_batch_fit_options(fit_options, drop_nonestimable_loci = isTRUE(fast))
  parent_config <- list(output_mode = "compact", instrument_rule = instrument_rule,
    sample_rule = list(mode = if (isTRUE(fast)) "batch_common_complete" else "common_complete"),
    fit_options = fit_options, parallel_strategy = "by_exposure", fast = isTRUE(fast),
    n_workers = as.integer(n_workers), sample_retention_min = sample_retention_min,
    allow_sample_loss = allow_sample_loss)
  parent_fingerprint <- if (is.null(output_dir)) NULL else
    qtl2mr_batch_fingerprint(resources, pairs, parent_config)
  # Reject a changed or damaged parent before preflight/gate branches can write.
  if (isTRUE(resume) && !isTRUE(force)) qtl2mr_validate_resume(output_dir, parent_fingerprint)
  preflight <- qtl2mr_batch_preflight(resources, pairs, instrument_rule, fit_options, fail_fast)
  if (!nrow(preflight$valid)) {
    batch <- qtl2mr_append_preflight_errors(qtl2mr_empty_batch(resources, pairs, parent_config), preflight)
    batch$resume_fingerprint <- parent_fingerprint
    if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
    return(batch)
  }
  pairs <- preflight$valid
  sample_qc <- NULL
  batch_common_samples <- NULL
  if (isTRUE(fast)) {
    sample_qc <- qtl2mr_batch_sample_qc(
      resources = resources,
      pairs = pairs,
      instrument_rule = instrument_rule,
      sample_retention_min = sample_retention_min,
      allow_sample_loss = allow_sample_loss
    )
    if (identical(sample_qc$sample_retention_gate_status[[1L]], "fail")) {
      config <- list(
        output_mode = "compact",
        instrument_rule = instrument_rule,
        sample_rule = list(mode = "batch_common_complete"),
        fit_options = fit_options,
        parallel_strategy = "by_exposure",
        created_at = as.character(Sys.time())
      )
      manifest <- data.frame(
        run_id = format(Sys.time(), "%Y%m%d_%H%M%S"),
        created_at = as.character(Sys.time()),
        n_pairs_requested = nrow(pairs),
        n_pairs_ok = 0L,
        n_pairs_error = 0L,
        output_mode = "compact",
        diagnostic_level = "compact",
        sample_mode = "batch_common_complete",
        batch_common_n = sample_qc$n_batch_common[[1L]],
        context_cache_enabled = TRUE,
        context_cache_hits = 0L,
        context_cache_misses = 0L,
        context_cache_size = 0L,
        lmm_fit_cache_enabled = TRUE,
        lmm_fit_cache_hits = 0L,
        lmm_fit_cache_misses = 0L,
        lmm_fit_cache_size = 0L,
        sample_approximation = TRUE,
        sample_retention_ratio = sample_qc$sample_retention_ratio[[1L]],
        sample_retention_min = sample_retention_min,
        sample_retention_gate_status = "fail",
        parallel_strategy = "by_exposure",
        n_workers = as.integer(n_workers),
        n_groups = 0L,
        nonestimable_locus_policy = if (isTRUE(fit_options$drop_nonestimable_loci)) "drop_audited_screening" else "fail_closed",
        n_loci_dropped_nonestimable = 0L,
        lmm_backend = "not_available",
        stringsAsFactors = FALSE
      )
      batch <- qtl2mr_merge_batch_group_results(resources, pairs, list(), config, manifest, sample_qc)
      batch$resume_fingerprint <- parent_fingerprint
      if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
      stop(
        sprintf(
          "batch_common_complete sample retention %.3f is below sample_retention_min %.3f; use run_qtl2mr_batch() or set allow_sample_loss = TRUE.",
          sample_qc$sample_retention_ratio[[1L]], sample_retention_min
        ),
        call. = FALSE
      )
    }
    selected <- qtl2mr_selected_loci_for_pairs(resources, pairs, instrument_rule)
    batch_common_samples <- qtl2mr_make_batch_common_complete_samples(
      resources,
      traits = unique(c(pairs$exposure, pairs$outcome)),
      selected_chr = unique(as.character(selected$chr))
    )
  }

  split_key <- as.character(pairs$exposure)
  group_index <- split(seq_len(nrow(pairs)), split_key)
  groups <- lapply(names(group_index), function(key) {
    list(key = key, pairs = pairs[group_index[[key]], , drop = FALSE],
      output_dir = if (is.null(output_dir)) NULL else file.path(output_dir, "by_exposure", make.names(key)))
  })
  names(groups) <- names(group_index)
  if (!is.null(output_dir) && anyDuplicated(vapply(groups, function(group) make.names(group$key), character(1)))) {
    stop("Exposure names map to the same output directory; use distinct syntactic trait names", call. = FALSE)
  }
  worker_sample_rule <- if (isTRUE(fast)) {
    list(mode = "batch_common_complete", sample_ids = batch_common_samples)
  } else list(mode = "common_complete")
  # Validate every existing child before any child can update its checkpoint.
  if (isTRUE(resume) && !isTRUE(force) && !is.null(output_dir)) {
    worker_config <- list(output_mode = "compact", instrument_rule = instrument_rule,
      sample_rule = worker_sample_rule, fit_options = fit_options)
    for (group in groups) {
      fingerprint <- qtl2mr_batch_fingerprint(resources, group$pairs, worker_config)
      qtl2mr_validate_resume(group$output_dir, fingerprint)
    }
  }
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  run_one <- function(group) {
    group_dir <- group$output_dir
    if (isTRUE(progress)) message(sprintf("[qtl2MR parallel] exposure=%s pairs=%d", group$key, nrow(group$pairs)))
    if (isTRUE(fast)) {
      out <- run_qtl2mr_batch(
        resources = resources,
        pairs = group$pairs,
        instrument_rule = instrument_rule,
        sample_rule = list(mode = "batch_common_complete", sample_ids = batch_common_samples),
        fit_options = fit_options,
        output_mode = "compact",
        fail_fast = fail_fast,
        resume = resume,
        force = force,
        output_dir = group_dir,
        progress = FALSE
      )
      out$sample_qc <- sample_qc
      out$run_manifest$sample_retention_ratio <- sample_qc$sample_retention_ratio[[1L]]
      out$run_manifest$sample_retention_min <- sample_retention_min
      out$run_manifest$sample_retention_gate_status <- sample_qc$sample_retention_gate_status[[1L]]
      out$run_manifest$sample_approximation <- TRUE
      if (!is.null(group_dir)) write_qtl2mr_batch_outputs(out, group_dir)
      out
    } else {
      run_qtl2mr_batch(
        resources = resources,
        pairs = group$pairs,
        instrument_rule = instrument_rule,
        sample_rule = list(mode = "common_complete"),
        fit_options = fit_options,
        output_mode = "compact",
        fail_fast = fail_fast,
        resume = resume,
        force = force,
        output_dir = group_dir,
        progress = FALSE
      )
    }
  }
  n_workers <- max(1L, min(as.integer(n_workers), length(groups)))
  backend <- "sequential"
  t0 <- proc.time()[["elapsed"]]
  safe_run_one <- function(group) {
    tryCatch(run_one(group), error = function(e) qtl2mr_parallel_group_failure(group, e))
  }
  if (n_workers > 1L && identical(.Platform$OS.type, "unix")) {
    backend <- "parallel::mclapply"
    group_results <- parallel::mclapply(groups, safe_run_one, mc.cores = n_workers)
  } else {
    group_results <- lapply(groups, safe_run_one)
  }
  collected <- qtl2mr_collect_parallel_group_results(group_results, groups, fail_fast = fail_fast)
  group_results <- collected$successes
  group_failures <- collected$failures
  wall_sec <- proc.time()[["elapsed"]] - t0
  context_hits <- sum(vapply(group_results, function(x) x$run_manifest$context_cache_hits[[1L]] %||% 0L, numeric(1)))
  context_misses <- sum(vapply(group_results, function(x) x$run_manifest$context_cache_misses[[1L]] %||% 0L, numeric(1)))
  context_size <- sum(vapply(group_results, function(x) x$run_manifest$context_cache_size[[1L]] %||% 0L, numeric(1)))
  fit_hits <- sum(vapply(group_results, function(x) x$run_manifest$lmm_fit_cache_hits[[1L]] %||% 0L, numeric(1)))
  fit_misses <- sum(vapply(group_results, function(x) x$run_manifest$lmm_fit_cache_misses[[1L]] %||% 0L, numeric(1)))
  fit_size <- sum(vapply(group_results, function(x) x$run_manifest$lmm_fit_cache_size[[1L]] %||% 0L, numeric(1)))
  n_ok <- sum(vapply(group_results, function(x) sum(x$main_results$status == "ok", na.rm = TRUE), integer(1)))
  n_err <- sum(vapply(group_results, function(x) sum(x$main_results$status == "error", na.rm = TRUE), integer(1))) +
    sum(group_failures$n_pairs)
  config <- parent_config
  config$created_at <- as.character(Sys.time())
  manifest <- data.frame(
    run_id = format(Sys.time(), "%Y%m%d_%H%M%S"),
    created_at = as.character(Sys.time()),
    n_pairs_requested = nrow(pairs),
    n_pairs_ok = n_ok,
    n_pairs_error = n_err,
    output_mode = "compact",
    diagnostic_level = "compact",
    sample_mode = if (isTRUE(fast)) "batch_common_complete" else "common_complete",
    batch_common_n = if (isTRUE(fast)) length(batch_common_samples) else 0L,
    context_cache_enabled = isTRUE(fit_options$cache_contexts),
    context_cache_hits = context_hits,
    context_cache_misses = context_misses,
    context_cache_size = context_size,
    lmm_fit_cache_enabled = isTRUE(fit_options$cache_lmm_fits),
    lmm_fit_cache_hits = fit_hits,
    lmm_fit_cache_misses = fit_misses,
    lmm_fit_cache_size = fit_size,
    sample_approximation = isTRUE(fast),
    sample_retention_ratio = if (isTRUE(fast)) sample_qc$sample_retention_ratio[[1L]] else NA_real_,
    sample_retention_min = if (isTRUE(fast)) sample_retention_min else NA_real_,
    sample_retention_gate_status = if (isTRUE(fast)) sample_qc$sample_retention_gate_status[[1L]] else "not_applicable",
    parallel_strategy = "by_exposure",
    parallel_backend = backend,
    n_workers = n_workers,
    n_groups = length(groups),
    wall_sec = wall_sec,
    mean_wall_sec_per_pair = wall_sec / max(nrow(pairs), 1L),
    nonestimable_locus_policy = if (isTRUE(fit_options$drop_nonestimable_loci)) "drop_audited_screening" else "fail_closed",
    n_loci_dropped_nonestimable = sum(vapply(
      group_results,
      function(x) sum(x$main_results$n_loci_dropped_nonestimable %||% 0L, na.rm = TRUE),
      numeric(1)
    )),
    lmm_backend = qtl2mr_backend_label(unlist(lapply(
      group_results,
      function(x) x$main_results$lmm_backend[x$main_results$status == "ok"]
    ), use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  batch <- qtl2mr_merge_batch_group_results(resources, pairs, group_results, config, manifest, sample_qc)
  batch$resume_fingerprint <- parent_fingerprint
  if (!is.null(output_dir)) {
    batch$resume_artifacts <- unique(unlist(lapply(groups, function(group) {
      paths <- qtl2mr_batch_owned_paths(group$output_dir)
      unlist(lapply(paths, function(path) {
        if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else path
      }), use.names = FALSE)
    }), use.names = FALSE))
  }
  if (nrow(group_failures)) {
    batch$parallel_group_failures <- group_failures
    affected <- pairs[pairs$exposure %in% group_failures$exposure_group, , drop = FALSE]
    failures <- group_failures[match(affected$exposure, group_failures$exposure_group), , drop = FALSE]
    error_rows <- data.frame(
      pair_id = affected$pair_id,
      exposure = affected$exposure,
      outcome = affected$outcome,
      stage = "parallel_exposure_group",
      error_class = failures$error_class,
      error_message = failures$error_message,
      traceback_available = FALSE,
      runtime_sec = NA_real_,
      exposure_group = failures$exposure_group,
      n_pairs_affected = failures$n_pairs,
      stringsAsFactors = FALSE
    )
    batch <- qtl2mr_batch_append(batch, "errors", error_rows)
    main <- data.frame(pair_id = affected$pair_id, exposure = affected$exposure,
      outcome = affected$outcome, instrument_trait = affected$instrument_trait,
      status = "error", error_message = failures$error_message,
      sample_policy = qtl2mr_batch_sample_policy(config$sample_rule),
      runtime_sec = NA_real_, stringsAsFactors = FALSE)
    batch <- qtl2mr_batch_append(batch, "main_results", main)
  }
  batch <- qtl2mr_append_preflight_errors(batch, preflight)
  if (!is.null(output_dir)) write_qtl2mr_batch_outputs(batch, output_dir)
  batch
}
