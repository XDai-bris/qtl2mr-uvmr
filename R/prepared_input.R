as_qtl2mr_prepared_input <- function(allele_probs,
                                     pheno,
                                     covar = NULL,
                                     kinship = NULL,
                                     marker_map,
                                     selected_loci = NULL,
                                     exposure = NULL,
                                     outcome = NULL,
                                     exposure_samples = NULL,
                                     outcome_samples = NULL,
                                     hsq = NULL,
                                     sample_ids = NULL,
                                     metadata = list(),
                                     use_all_loci = FALSE,
                                     probability_policy = c("normalize_rounding", "strict")) {
  probability_policy <- match.arg(probability_policy)
  if (!is.null(metadata) && !is.list(metadata)) stop("metadata must be a list", call. = FALSE)
  if (!is.logical(use_all_loci) || length(use_all_loci) != 1L || is.na(use_all_loci)) {
    stop("use_all_loci must be TRUE or FALSE", call. = FALSE)
  }
  input <- list(
    allele_probs = allele_probs,
    pheno = pheno,
    covar = covar,
    kinship = kinship,
    marker_map = marker_map,
    selected_loci = selected_loci,
    use_all_loci = use_all_loci,
    exposure = exposure,
    outcome = outcome,
    exposure_samples = exposure_samples,
    outcome_samples = outcome_samples,
    hsq = hsq,
    sample_ids = sample_ids,
    metadata = metadata %||% list(),
    probability_policy = probability_policy
  )
  class(input) <- c("qtl2mr_prepared_input", "list")
  validate_qtl2mr_input(input)
  input <- qtl2mr_record_probability_policy(input, warn = TRUE)
  input$sample_alignment <- qtl2mr_input_alignment(input)
  qtl2mr_warn_alignment(input$sample_alignment)
  input
}

qtl2mr_validate_ids <- function(ids, label) {
  if (is.null(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(label, " must contain non-empty, unique sample/marker IDs", call. = FALSE)
  }
  invisible(TRUE)
}

qtl2mr_input_alignment <- function(input) {
  ids <- qtl2mr_prepared_sample_ids(input)
  sources <- list(pheno = rownames(input$pheno))
  if (!is.null(input$covar)) sources$covar <- rownames(input$covar)
  do.call(rbind, lapply(names(sources), function(source) {
    source_ids <- sources[[source]]
    data.frame(source = source, n_requested = length(ids), n_available = length(source_ids),
               n_matched = sum(ids %in% source_ids),
               n_missing = sum(!ids %in% source_ids),
               missing_ids = paste(setdiff(ids, source_ids), collapse = ";"),
               stringsAsFactors = FALSE)
  }))
}

qtl2mr_warn_alignment <- function(alignment) {
  dropped <- alignment[alignment$n_missing > 0L, , drop = FALSE]
  if (nrow(dropped)) warning(
    "Sample ID mismatch: ", paste(paste0(dropped$source, " missing ", dropped$n_missing,
                                         " requested IDs"), collapse = "; "),
    ". Unmatched samples cannot be fitted; inspect sample_alignment.", call. = FALSE)
  invisible(alignment)
}

qtl2mr_complete_rows <- function(x) {
  x <- as.data.frame(x)
  keep <- rep(TRUE, nrow(x))
  for (column in x) keep <- keep & if (is.numeric(column)) is.finite(column) else !is.na(column)
  keep
}

qtl2mr_marker_map_table <- function(marker_map) {
  if (is.null(marker_map)) stop("marker_map is required", call. = FALSE)
  if (is.data.frame(marker_map)) {
    out <- marker_map
    if (!"marker" %in% names(out)) {
      if (!is.null(rownames(out)) && all(nzchar(rownames(out)))) {
        out$marker <- rownames(out)
      } else {
        stop("marker_map data frame must contain a marker column or marker row names", call. = FALSE)
      }
    }
    if (!"chr" %in% names(out)) {
      chr_col <- intersect(c("chromosome", "chr_id", "chrom"), names(out))
      if (length(chr_col)) out$chr <- out[[chr_col[[1L]]]] else stop("marker_map must contain a chr column", call. = FALSE)
    }
    out$marker <- as.character(out$marker)
    out$chr <- as.character(out$chr)
    if (!"pos" %in% names(out)) out$pos <- NA_real_
    if (!is.numeric(out$pos)) stop("marker_map pos must be numeric", call. = FALSE)
    qtl2mr_validate_ids(out$marker, "marker_map markers")
    if (anyNA(out$chr) || any(!nzchar(out$chr))) stop("marker_map contains missing chromosomes", call. = FALSE)
    return(out[, unique(c("marker", "chr", "pos", names(out))), drop = FALSE])
  }
  if (is.list(marker_map)) {
    if (is.null(names(marker_map)) || anyNA(names(marker_map)) || any(!nzchar(names(marker_map))) || anyDuplicated(names(marker_map))) {
      stop("qtl2-style marker_map must be named by chromosome", call. = FALSE)
    }
    rows <- list()
    for (chr in names(marker_map) %||% as.character(seq_along(marker_map))) {
      mm <- marker_map[[chr]]
      if (is.null(names(mm))) {
        marker <- as.character(mm)
        pos <- NA_real_
      } else {
        marker <- names(mm)
        pos <- as.numeric(mm)
      }
      rows[[length(rows) + 1L]] <- data.frame(marker = marker, chr = as.character(chr), pos = pos, stringsAsFactors = FALSE)
    }
    return(qtl2mr_marker_map_table(do.call(rbind, rows)))
  }
  if (!is.null(names(marker_map))) {
    return(qtl2mr_marker_map_table(data.frame(marker = names(marker_map), chr = as.character(marker_map), pos = NA_real_, stringsAsFactors = FALSE)))
  }
  stop("marker_map must be a data frame, qtl2-style list, or named vector", call. = FALSE)
}

qtl2mr_allele_prob_markers <- function(allele_probs, markers_to_validate = NULL,
                                     probability_policy = c("strict", "normalize_rounding")) {
  probability_policy <- match.arg(probability_policy)
  if (!is.list(allele_probs) || !length(allele_probs)) stop("allele_probs must be a non-empty qtl2-style list of 3D arrays", call. = FALSE)
  qtl2mr_validate_ids(names(allele_probs), "allele_probs chromosome names")
  chr_names <- names(allele_probs) %||% as.character(seq_along(allele_probs))
  rows <- list()
  adjustments <- list()
  for (i in seq_along(allele_probs)) {
    arr <- allele_probs[[i]]
    if (length(dim(arr)) != 3L) stop("Each allele_probs element must be a sample x founder x marker array", call. = FALSE)
    dn <- dimnames(arr)
    for (d in seq_len(3L)) qtl2mr_validate_ids(dn[[d]], paste0("allele_probs dimension ", d))
    if (!is.numeric(arr)) stop("allele_probs must be numeric", call. = FALSE)
    # Bound temporary allocations for genome-wide inputs.
    checked <- if (is.null(markers_to_validate)) seq_len(dim(arr)[3L]) else which(dn[[3L]] %in% markers_to_validate)
    for (first in seq.int(1L, max(length(checked), 1L), by = 128L)) {
      if (!length(checked)) break
      chunk <- arr[, , checked[seq.int(first, min(first + 127L, length(checked)))], drop = FALSE]
      if (any(!is.finite(chunk))) stop("allele_probs must contain finite founder allele probabilities", call. = FALSE)
      if (any(chunk < -1e-8 | chunk > 1 + 1e-8)) {
        stop(sprintf("allele_probs must lie in [0, 1] up to floating-point error (1e-8); observed range [%.10g, %.10g]",
          min(chunk), max(chunk)), call. = FALSE)
      }
      deviations <- abs(apply(chunk, c(1L, 3L), sum) - 1)
      sum_tolerance <- if (probability_policy == "strict") 1e-6 else dim(arr)[2L] * 5e-5
      if (any(deviations > sum_tolerance + 16 * .Machine$double.eps)) {
        stop(sprintf(paste0("allele_probs must sum to 1 over founders: maximum observed |sum - 1|=%.8g; ",
          "allowed tolerance=%.8g under probability_policy='%s'. Check the probability export, rounding precision, and normalization; this exceeds the selected policy's rounding allowance."),
          max(deviations), sum_tolerance, probability_policy), call. = FALSE)
      }
      adjustments[[length(adjustments) + 1L]] <- data.frame(chr = chr_names[[i]],
        marker = dimnames(chunk)[[3L]], n_normalized = colSums(deviations > 1e-6 + 16 * .Machine$double.eps),
        max_abs_sum_deviation = apply(deviations, 2L, max), stringsAsFactors = FALSE)
    }
    if (i > 1L && !identical(dn[[2L]], dimnames(allele_probs[[1L]])[[2L]])) {
      stop("allele_probs founder names and order must match across chromosomes", call. = FALSE)
    }
    markers <- dn[[3L]]
    rows[[length(rows) + 1L]] <- data.frame(marker = markers, chr = chr_names[[i]], stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  qtl2mr_validate_ids(out$marker, "allele_probs marker names")
  attr(out, "probability_normalization") <- do.call(rbind, adjustments)
  out
}

qtl2mr_probability_policy_disclosed <- function(previous, current) {
  if (!is.list(previous) || !identical(previous$policy, current$policy) ||
      !identical(previous$strict_sum_tolerance, current$strict_sum_tolerance) ||
      !identical(previous$rounding_sum_tolerance, current$rounding_sum_tolerance)) return(FALSE)
  old <- previous$adjustments
  now <- current$adjustments
  if (!is.data.frame(old) || !all(names(now) %in% names(old))) return(FALSE)
  selected <- match(now$marker, old$marker)
  if (anyNA(selected)) return(FALSE)
  old <- old[selected, names(now), drop = FALSE]
  rownames(old) <- rownames(now) <- NULL
  identical(old, now)
}

qtl2mr_record_probability_policy <- function(input, markers = NULL, warn = FALSE) {
  previous <- input$metadata$probability_normalization
  policy <- input$probability_policy %||% "strict"
  if (is.null(markers)) markers <- if (is.data.frame(input$selected_loci)) input$selected_loci$marker else input$selected_loci
  marker_info <- qtl2mr_allele_prob_markers(input$allele_probs, markers, probability_policy = policy)
  adjustment <- attr(marker_info, "probability_normalization") %||% data.frame(
    chr = character(), marker = character(), n_normalized = integer(),
    max_abs_sum_deviation = numeric(), stringsAsFactors = FALSE)
  adjustment <- adjustment[adjustment$n_normalized > 0L, , drop = FALSE]
  input$metadata$probability_normalization <- list(policy = policy,
    strict_sum_tolerance = 1e-6, rounding_sum_tolerance = dim(input$allele_probs[[1L]])[2L] * 5e-5,
    original_arrays_preserved = TRUE, adjustments = adjustment)
  if (isTRUE(warn) && nrow(adjustment) &&
      !qtl2mr_probability_policy_disclosed(previous, input$metadata$probability_normalization)) warning(
    "Founder probability rounding: ", sum(adjustment$n_normalized),
    " sample-marker vectors will be normalized for analysis under probability_policy='normalize_rounding'; original arrays are retained. Inspect metadata$probability_normalization.",
    call. = FALSE)
  input
}

qtl2mr_prepared_probability_copy <- function(input, markers) {
  previous <- input$metadata$probability_normalization
  input <- qtl2mr_record_probability_policy(input, markers)
  record <- input$metadata$probability_normalization
  if (nrow(record$adjustments) && !qtl2mr_probability_policy_disclosed(previous, record)) {
    warning("Founder probability rounding normalized for analysis; inspect prepared_metadata$probability_normalization. Original input arrays are unchanged.", call. = FALSE)
  }
  for (chr in unique(record$adjustments$chr)) {
    arr <- input$allele_probs[[chr]]
    for (marker in record$adjustments$marker[record$adjustments$chr == chr]) {
      values <- arr[, , marker, drop = FALSE]
      sums <- apply(values, 1L, sum)
      adjust <- abs(sums - 1) > record$strict_sum_tolerance + 16 * .Machine$double.eps
      arr[adjust, , marker] <- sweep(matrix(values[adjust, , , drop = FALSE], nrow = sum(adjust)), 1L, sums[adjust], "/")
    }
    input$allele_probs[[chr]] <- arr
  }
  input
}

qtl2mr_prepared_sample_ids <- function(input) {
  gp_ids <- unique(unlist(lapply(input$allele_probs, function(x) dimnames(x)[[1L]]), use.names = FALSE))
  if (is.null(gp_ids) || !length(gp_ids) || any(!nzchar(gp_ids))) {
    stop("allele_probs arrays must have sample row names", call. = FALSE)
  }
  ids <- input$sample_ids %||% gp_ids
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(ids))) stop("sample_ids must not contain missing or blank IDs", call. = FALSE)
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) stop("sample_ids contains duplicate IDs: ", paste(dup, collapse = ", "), call. = FALSE)
  if (length(setdiff(ids, gp_ids))) stop("sample_ids contains IDs absent from allele_probs", call. = FALSE)
  ids
}

validate_qtl2mr_input <- function(input,
                                  founder_dim = 8L,
                                  require_selected_loci = FALSE,
                                  exposure = NULL,
                                  outcome = NULL) {
  if (!inherits(input, "qtl2mr_prepared_input") && !is.list(input)) {
    stop("input must be a qtl2mr_prepared_input or compatible list", call. = FALSE)
  }
  exposure <- exposure %||% input$exposure
  outcome <- outcome %||% input$outcome
  for (role in c("exposure", "outcome")) {
    trait <- if (role == "exposure") exposure else outcome
    if (!is.null(trait) && (!is.character(trait) || length(trait) != 1L || is.na(trait) || !nzchar(trait))) {
      stop(role, " must be one non-empty trait name", call. = FALSE)
    }
  }
  if (!is.null(exposure) && !is.null(outcome) && identical(exposure, outcome)) {
    stop("exposure and outcome must be different traits; both are '", exposure, "'", call. = FALSE)
  }
  needed <- c("allele_probs", "pheno", "marker_map")
  miss <- setdiff(needed, names(input))
  if (length(miss)) stop("prepared input is missing: ", paste(miss, collapse = ", "), call. = FALSE)
  marker_tab <- qtl2mr_marker_map_table(input$marker_map)
  checked_markers <- if (is.data.frame(input$selected_loci)) input$selected_loci$marker else input$selected_loci
  gp_markers <- qtl2mr_allele_prob_markers(input$allele_probs, markers_to_validate = checked_markers,
    probability_policy = input$probability_policy %||% "strict")
  missing_map <- setdiff(gp_markers$marker, marker_tab$marker)
  if (length(missing_map)) stop("marker_map is missing allele probability markers: ", paste(utils::head(missing_map, 5L), collapse = ", "), call. = FALSE)
  founder_n <- unique(vapply(input$allele_probs, function(x) dim(x)[2L], integer(1)))
  if (length(founder_n) != 1L || founder_n[[1L]] != founder_dim) {
    stop("founder dimension must be ", founder_dim, " for the autosomal additive founder model", call. = FALSE)
  }
  ids <- qtl2mr_prepared_sample_ids(input)
  pheno <- as.data.frame(input$pheno)
  qtl2mr_validate_prepared_hsq(input, unique(marker_tab$chr), names(pheno))
  qtl2mr_validate_ids(rownames(input$pheno), "pheno row names")
  if (is.null(rownames(pheno)) || any(!nzchar(rownames(pheno)))) {
    stop("pheno must have row names matching sample IDs", call. = FALSE)
  }
  if (!is.null(exposure) && !exposure %in% names(pheno)) stop("Exposure trait not found in prepared pheno: ", exposure, call. = FALSE)
  if (!is.null(outcome) && !outcome %in% names(pheno)) stop("Outcome trait not found in prepared pheno: ", outcome, call. = FALSE)
  traits <- unique(c(exposure, outcome))
  for (trait in traits) {
    if (!is.numeric(pheno[[trait]]) || is.factor(pheno[[trait]])) {
      stop("Phenotype ", trait, " must be numeric; factor/character values are not silently recoded", call. = FALSE)
    }
  }
  if (!length(intersect(ids, rownames(pheno)))) stop("pheno and allele_probs have no overlapping sample IDs", call. = FALSE)
  for (role in c("exposure", "outcome")) {
    samples <- input[[paste0(role, "_samples")]]
    if (is.null(samples)) next
    samples <- as.character(samples)
    if (anyNA(samples) || any(!nzchar(samples))) stop(role, "_samples must not contain missing or blank IDs", call. = FALSE)
    dup <- unique(samples[duplicated(samples)])
    if (length(dup)) stop(role, "_samples contains duplicate IDs: ", paste(dup, collapse = ", "), call. = FALSE)
    missing_gp <- setdiff(samples, ids)
    if (length(missing_gp)) stop(role, "_samples are missing from allele_probs/sample_ids: ", paste(utils::head(missing_gp, 5L), collapse = ", "), call. = FALSE)
    missing_pheno <- setdiff(samples, rownames(pheno))
    if (length(missing_pheno)) stop(role, "_samples are missing from pheno: ", paste(utils::head(missing_pheno, 5L), collapse = ", "), call. = FALSE)
  }
  if (!is.null(input$covar)) {
    covar <- as.data.frame(input$covar)
    qtl2mr_validate_ids(rownames(input$covar), "covar row names")
    if (!length(Reduce(intersect, list(ids, rownames(pheno), rownames(covar))))) stop("covar, pheno and allele_probs have no overlapping sample IDs", call. = FALSE)
    if (is.null(rownames(covar)) || any(!nzchar(rownames(covar)))) {
      stop("covar must have row names matching sample IDs", call. = FALSE)
    }
  }
  selected <- input$selected_loci
  selected_markers <- NULL
  if (is.null(selected)) {
    if (isTRUE(require_selected_loci) && !isTRUE(input$use_all_loci)) stop("selected_loci is required; choose loci explicitly or set use_all_loci = TRUE", call. = FALSE)
  } else {
    if (is.data.frame(selected) && !"marker" %in% names(selected)) stop("selected_loci must contain a marker column", call. = FALSE)
    selected_markers <- if (is.data.frame(selected)) as.character(selected$marker) else as.character(selected)
    if (!length(selected_markers)) {
      stop("selected_loci contains no loci; supply at least one marker, or use selected_loci = NULL with use_all_loci = TRUE", call. = FALSE)
    }
    if (anyNA(selected_markers) || any(!nzchar(selected_markers))) stop("selected_loci contains missing or blank markers", call. = FALSE)
    dup_m <- unique(selected_markers[duplicated(selected_markers)])
    if (length(dup_m)) stop("selected_loci contains duplicate markers: ", paste(dup_m, collapse = ", "), call. = FALSE)
    missing_sel <- setdiff(selected_markers, gp_markers$marker)
    if (length(missing_sel)) stop("selected_loci markers are not present in allele_probs: ", paste(missing_sel, collapse = ", "), call. = FALSE)
    selected_chr <- marker_tab$chr[match(selected_markers, marker_tab$marker)]
    if (any(toupper(sub("^chr", "", selected_chr, ignore.case = TRUE)) %in% c("X", "Y", "M", "MT"))) {
      stop("Only autosomal selected_loci are supported; sex-chromosome and mitochondrial models are not implemented", call. = FALSE)
    }
    if (is.data.frame(selected) && "chr" %in% names(selected) && !identical(as.character(selected$chr), selected_chr)) {
      stop("selected_loci chromosome assignments disagree with marker_map", call. = FALSE)
    }
  }
  if (any(gp_markers$chr != marker_tab$chr[match(gp_markers$marker, marker_tab$marker)])) stop("marker_map chromosome assignments disagree with allele_probs", call. = FALSE)
  if (!is.null(input$kinship)) {
    kinships <- qtl2mr_prepared_kinship_list(input$kinship)
    selected_chr <- unique(qtl2mr_marker_map_table(input$marker_map)$chr)
    if (!is.null(selected_markers)) {
      map0 <- qtl2mr_marker_map_table(input$marker_map)
      selected_chr <- unique(map0$chr[match(selected_markers, map0$marker)])
    }
    if (is.list(input$kinship) && !is.matrix(input$kinship) && !is.data.frame(input$kinship)) {
      missing_chr <- setdiff(as.character(selected_chr), names(kinships))
      if (length(missing_chr)) {
        stop("kinship list is missing selected chromosome(s): ", paste(missing_chr, collapse = ", "), call. = FALSE)
      }
    }
    role_ids <- unique(c(
      qtl2mr_prepared_candidate_ids(input, exposure, "exposure"),
      qtl2mr_prepared_candidate_ids(input, outcome, "outcome")
    ))
    if (!length(role_ids)) role_ids <- intersect(ids, rownames(pheno))
    for (nm in names(kinships) %||% as.character(seq_along(kinships))) {
      K <- as.matrix(kinships[[nm]])
      if (!is.numeric(K) || nrow(K) != ncol(K) || any(!is.finite(K))) stop("kinship must be a finite numeric square matrix", call. = FALSE)
      qtl2mr_validate_ids(rownames(K), "kinship row names")
      qtl2mr_validate_ids(colnames(K), "kinship column names")
      if (!identical(rownames(K), colnames(K))) stop("kinship row and column names must have the same order", call. = FALSE)
      if (!isTRUE(all.equal(K, t(K), tolerance = 1e-8))) stop("kinship must be symmetric", call. = FALSE)
      if (is.null(rownames(K)) || is.null(colnames(K))) stop("kinship matrices must have row and column names", call. = FALSE)
      missing_k <- setdiff(role_ids, intersect(rownames(K), colnames(K)))
      if (length(missing_k)) stop("kinship ", nm, " is missing sample IDs: ", paste(utils::head(missing_k, 5L), collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

qtl2mr_prepared_role_filter <- function(input, role) {
  role_samples <- input[[paste0(role, "_samples")]]
  if (!is.null(role_samples)) return(as.character(role_samples))
  input$sample_ids %||% qtl2mr_prepared_sample_ids(input)
}

qtl2mr_prepared_candidate_ids <- function(input, trait = NULL, role = c("exposure", "outcome")) {
  role <- match.arg(role)
  ids <- qtl2mr_prepared_role_filter(input, role)
  pheno <- as.data.frame(input$pheno)
  ids <- ids[ids %in% rownames(pheno)]
  if (!is.null(trait) && trait %in% names(pheno)) {
    z <- as.numeric(pheno[[trait]])
    names(z) <- rownames(pheno)
    ids <- ids[is.finite(z[ids])]
  }
  if (!is.null(input$covar)) {
    covar <- as.data.frame(input$covar)
    ids <- ids[ids %in% rownames(covar)]
  }
  ids
}

qtl2mr_prepared_selected_loci <- function(input, loci = NULL) {
  if (!is.null(loci)) {
    valid_count <- is.numeric(loci) && length(loci) == 1L && is.finite(loci) && loci >= 1 && loci == as.integer(loci)
    valid_names <- is.character(loci) && length(loci) > 0L && !anyNA(loci) && all(nzchar(loci)) && !anyDuplicated(loci)
    if (!valid_count && !valid_names) stop("loci must be a positive integer count or unique, non-empty marker names", call. = FALSE)
  }
  marker_tab <- qtl2mr_marker_map_table(input$marker_map)
  if (is.null(input$selected_loci)) {
    explicit_markers <- is.character(loci) && length(loci) > 0L
    if (!isTRUE(input$use_all_loci) && !explicit_markers) stop("selected_loci is required; provide marker names in loci or explicitly set use_all_loci = TRUE", call. = FALSE)
    selected <- marker_tab
  } else if (is.data.frame(input$selected_loci)) {
    selected <- merge(input$selected_loci, marker_tab[, c("marker", "chr", "pos"), drop = FALSE],
                      by = "marker", all.x = TRUE, suffixes = c("", ".map"), sort = FALSE)
    if ("chr.map" %in% names(selected)) {
      selected$chr <- selected$chr %||% selected$chr.map
      selected$chr[is.na(selected$chr)] <- selected$chr.map[is.na(selected$chr)]
      selected$chr.map <- NULL
    }
  } else {
    selected <- marker_tab[match(as.character(input$selected_loci), marker_tab$marker), , drop = FALSE]
  }
  if (!"locus_rank" %in% names(selected)) selected$locus_rank <- seq_len(nrow(selected))
  if (!is.null(loci)) {
    if (is.numeric(loci) && length(loci) == 1L) selected <- utils::head(selected, as.integer(loci))
    else selected <- selected[match(as.character(loci), selected$marker), , drop = FALSE]
  }
  if (!nrow(selected) || anyNA(selected$marker)) stop("No selected loci remain after loci filtering", call. = FALSE)
  if (any(toupper(sub("^chr", "", selected$chr, ignore.case = TRUE)) %in% c("X", "Y", "M", "MT"))) stop("Only autosomal selected_loci are supported", call. = FALSE)
  selected
}

qtl2mr_prepared_trait <- function(input, trait) {
  pheno <- as.data.frame(input$pheno)
  if (!trait %in% names(pheno)) stop("Trait not found in prepared pheno: ", trait, call. = FALSE)
  if (!is.numeric(pheno[[trait]])) stop("Phenotype ", trait, " must be numeric", call. = FALSE)
  out <- as.numeric(pheno[[trait]])
  names(out) <- rownames(pheno)
  out
}

qtl2mr_prepared_covar <- function(input) {
  if (is.null(input$covar)) return(NULL)
  as.data.frame(input$covar)
}

qtl2mr_prepared_kinship_list <- function(kinship) {
  if (is.null(kinship)) return(NULL)
  if (is.matrix(kinship) || is.data.frame(kinship)) return(list(`1` = as.matrix(kinship)))
  if (!is.list(kinship)) stop("kinship must be a matrix or chromosome-named list of matrices", call. = FALSE)
  if (any(names(kinship) %in% c("exposure", "outcome"))) {
    stop("prepared inputs use one kinship source for exposure and outcome; separate exposure/outcome kinships are not supported", call. = FALSE)
  }
  if (is.null(names(kinship)) || any(!nzchar(names(kinship)))) {
    stop("kinship lists must be named by chromosome (for example, qtl2::calc_kinship(type = \"loco\"))", call. = FALSE)
  }
  if (!all(vapply(kinship, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    stop("each kinship list element must be a matrix", call. = FALSE)
  }
  lapply(kinship, as.matrix)
}

qtl2mr_prepared_kinship <- function(input, chr) {
  if (is.null(input$kinship)) return(NULL)
  if (is.matrix(input$kinship) || is.data.frame(input$kinship)) return(as.matrix(input$kinship))
  kinships <- qtl2mr_prepared_kinship_list(input$kinship)
  K <- kinships[[as.character(chr)]]
  if (is.null(K)) stop("kinship list is missing selected chromosome: ", chr, call. = FALSE)
  K
}

qtl2mr_prepared_hsq <- function(input, chr, trait_role, trait_name) {
  h <- input$hsq
  if (is.null(h)) return(NULL)
  if (is.numeric(h) && is.null(names(h)) && length(h) == 1L) return(qtl2mr_validate_hsq(h))
  chr <- as.character(chr)
  get1 <- function(x, key) {
    if (is.null(x) || is.null(key) || is.null(names(x)) || !key %in% names(x)) return(NULL)
    x[[key]]
  }
  candidates <- list(
    get1(get1(h, chr), trait_role),
    get1(get1(h, chr), trait_name),
    get1(get1(h, trait_role), chr),
    get1(get1(h, trait_name), chr),
    get1(h, chr),
    get1(h, trait_role),
    get1(h, trait_name)
  )
  candidates <- Filter(function(z) !is.null(z) && !is.list(z) && is.null(names(z)), candidates)
  if (!length(candidates)) return(NULL)
  values <- vapply(candidates, qtl2mr_validate_hsq, numeric(1))
  if (length(unique(values)) > 1L) stop("Conflicting hsq values for chromosome ", chr, " and ", trait_role, "/", trait_name, call. = FALSE)
  values[[1L]]
}

qtl2mr_validate_prepared_hsq <- function(input, chromosomes, traits) {
  h <- input$hsq
  if (is.null(h)) return(invisible(TRUE))
  if (is.numeric(h) && is.null(names(h)) && length(h) == 1L) {
    qtl2mr_validate_hsq(h)
    if (is.null(input$kinship) && h != 0) stop("Nonzero hsq requires kinship", call. = FALSE)
    return(invisible(TRUE))
  }
  roles <- unique(c("exposure", "outcome", traits))
  keys <- function(x, allowed) {
    if ((!is.list(x) && !is.numeric(x)) || is.object(x) || is.null(names(x)) ||
        !length(x) || anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)) ||
        any(!names(x) %in% allowed)) {
      stop("hsq must be a scalar or a uniquely named numeric vector/list keyed by chromosome, trait, or exposure/outcome role", call. = FALSE)
    }
  }
  leaf <- function(x) {
    if (is.null(x)) return(invisible(NULL))
    value <- qtl2mr_validate_hsq(x)
    if (is.null(input$kinship) && value != 0) stop("Nonzero hsq requires kinship", call. = FALSE)
  }
  keys(h, unique(c(chromosomes, roles)))
  for (key in names(h)) {
    value <- h[[key]]
    if (is.null(value)) next
    if (is.list(value) || !is.null(names(value))) {
      if (key %in% chromosomes && key %in% roles) stop("Ambiguous hsq key: ", key, call. = FALSE)
      keys(value, if (key %in% chromosomes) roles else chromosomes)
      for (entry in value) leaf(entry)
    } else leaf(value)
  }
  for (chr in chromosomes) for (role in c("exposure", "outcome")) {
    trait <- input[[role]]
    if (!is.null(trait)) qtl2mr_prepared_hsq(input, chr, role, trait)
  }
  invisible(TRUE)
}

qtl2mr_locus_overlap_table <- function(blocks) {
  if (!length(blocks)) return(data.frame())
  do.call(rbind, lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    s <- b$O_xy_summary %||% list()
    n_x <- as.integer(s$n_x %||% b$fit_x$n %||% NA_integer_)
    n_y <- as.integer(s$n_y %||% b$fit_y$n %||% NA_integer_)
    n_overlap <- as.integer(s$n_overlap %||% b$n_overlap %||% NA_integer_)
    data.frame(
      locus_index = i,
      marker = b$marker %||% names(blocks)[[i]],
      chr = b$chr %||% NA_character_,
      n_x = n_x,
      n_y = n_y,
      n_overlap = n_overlap,
      exposure_overlap_fraction = as.numeric(s$exposure_overlap_fraction %||% if (is.finite(n_x) && n_x > 0) n_overlap / n_x else NA_real_),
      outcome_overlap_fraction = as.numeric(s$outcome_overlap_fraction %||% if (is.finite(n_y) && n_y > 0) n_overlap / n_y else NA_real_),
      jaccard_overlap_fraction = as.numeric(s$jaccard_overlap_fraction %||% {
        denom <- n_x + n_y - n_overlap
        if (is.finite(denom) && denom > 0) n_overlap / denom else NA_real_
      }),
      complete_overlap = isTRUE(s$complete_overlap),
      independent_samples = isTRUE(n_overlap == 0L),
      stringsAsFactors = FALSE
    )
  }))
}

qtl2mr_reported_estimate <- function(fit,
                                     diagnostic_layer,
                                     policy = c("diagnostic_only", "strict_so", "auto_fallback")) {
  policy <- match.arg(policy)
  so <- list(
    theta = fit$theta %||% NA_real_,
    se = fit$se %||% NA_real_,
    ci_lower = fit$profile_ci_lower %||% NA_real_,
    ci_upper = fit$profile_ci_upper %||% NA_real_,
    ci_type = "profile",
    estimator = "Star-space qtl2mr",
    reason = "The numerical search returned a finite estimate and bounded confidence interval; statistical calibration and causal assumptions require separate assessment.",
    source = "primary"
  )
  so_reliable <- isTRUE(fit$primary_reliable)
  diagnostic_estimates <- diagnostic_layer$estimates %||% data.frame()
  two_step <- diagnostic_estimates[diagnostic_estimates$estimator == "Two-step IVW", , drop = FALSE]
  two_step_value <- function(col, default = NA_real_) {
    if (!nrow(two_step) || !col %in% names(two_step)) return(default)
    two_step[[col]][[1L]]
  }
  has_two_step <- nrow(two_step) &&
    is.finite(two_step$theta[[1L]]) &&
    is.finite(two_step$se[[1L]])

  if (identical(policy, "strict_so") || isTRUE(so_reliable)) {
    so$reason <- if (isTRUE(so_reliable)) {
      "The numerical search returned a finite estimate and bounded confidence interval; statistical calibration and causal assumptions require separate assessment."
    } else {
      "Strict SO policy requested; the primary estimate is reported despite numerical or confidence-set warnings."
    }
    so$policy <- policy
    return(so)
  }

  if (identical(policy, "diagnostic_only")) {
    return(list(
      theta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      ci_type = NA_character_,
      estimator = "not_reported_unreliable_so",
      reason = "The primary estimate lacks an acceptable bounded numerical result; diagnostic_only does not substitute a fallback causal estimate.",
      source = "none",
      policy = policy
    ))
  }

  if (isTRUE(has_two_step)) {
    return(list(
      theta = two_step$theta[[1L]],
      se = two_step$se[[1L]],
      ci_lower = two_step_value("ci_lower", NA_real_),
      ci_upper = two_step_value("ci_upper", NA_real_),
      ci_type = two_step_value("ci_type", "Wald"),
      estimator = "Two-step IVW",
      reason = "Primary star-space second-order estimator is unreliable; auto_fallback reports Two-step IVW diagnostic estimate.",
      source = "diagnostic_fallback",
      policy = policy
    ))
  }

  list(
    theta = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    ci_type = NA_character_,
    estimator = "not_reported_no_fallback",
    reason = "Primary star-space second-order estimator is unreliable and no finite Two-step IVW fallback was available.",
    source = "none",
    policy = policy
  )
}

# TRUE when every confidence-set component is finite (a bounded interval or
# several finite pieces), FALSE when a component reaches -Inf/Inf, NA when the
# set could not be resolved.
qtl2mr_confidence_set_bounded <- function(fit) {
  type <- fit$confidence_set_type %||% "unresolved"
  set <- fit$confidence_set
  if (identical(type, "unresolved") || !is.data.frame(set) || !nrow(set)) return(NA)
  all(is.finite(c(set$lower, set$upper)))
}

qtl2mr_hsq_table <- function(blocks, exposure, outcome) {
  rows <- lapply(blocks, function(block) {
    do.call(rbind, lapply(c("exposure", "outcome"), function(role) {
      fit <- if (role == "exposure") block$fit_x else block$fit_y
      status <- fit$hsq_status %||% NA_character_
      data.frame(
        role = role,
        trait = if (role == "exposure") exposure else outcome,
        chr = as.character(block$chr %||% NA_character_),
        hsq = as.numeric(fit$hsq %||% NA_real_),
        hsq_status = if (identical(status, "fixed")) "user_supplied" else status,
        hsq_source = if (identical(status, "not_applicable")) "no_kinship_ols" else
          if (identical(status, "fixed")) "user_supplied" else "reml_null_model",
        n_samples = as.integer(fit$n %||% NA_integer_),
        marker = block$marker %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }))
  })
  d <- do.call(rbind, rows)
  if (is.null(d) || !nrow(d)) {
    return(data.frame(role = character(), trait = character(), chr = character(), hsq = numeric(),
      hsq_status = character(), hsq_source = character(), n_samples = integer(),
      n_loci = integer(), markers = character(), stringsAsFactors = FALSE))
  }
  key <- paste(d$role, d$chr, format(d$hsq, digits = 15), d$hsq_status, d$n_samples, sep = "\r")
  out <- d[!duplicated(key), setdiff(names(d), "marker"), drop = FALSE]
  out_key <- key[!duplicated(key)]
  out$n_loci <- as.integer(vapply(out_key, function(k) sum(key == k), numeric(1)))
  out$markers <- vapply(out_key, function(k) paste(d$marker[key == k], collapse = ";"), character(1))
  rownames(out) <- NULL
  out
}

run_prepared_qtl2mr <- function(input,
                                exposure,
                                outcome,
                                loci = NULL,
                                output_dir = NULL,
                                write_outputs = FALSE,
                                verbose = TRUE,
                                standardize_phenotypes = TRUE,
                                cross_cov_model = "star",
                                ridge = 1e-10,
                                progress_callback = NULL,
                                heterogeneity_method = c("auto", "simulation", "chisq", "none"),
                                q_sim_n = 0L,
                                q_sim_seed = 1L,
                                diagnostic_bootstrap_B = 500L,
                                diagnostic_bootstrap_seed = 1L,
                                reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
                                engine = c("batch", "locus"),
                                context_cache = NULL,
                                fit_cache = NULL,
                                drop_nonestimable_loci = FALSE,
                                ...) {
  if (length(list(...))) stop("Unused arguments in run_prepared_qtl2mr: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  heterogeneity_method <- match.arg(heterogeneity_method)
  reported_estimator_policy <- match.arg(reported_estimator_policy)
  engine <- match.arg(engine)
  exposure <- exposure %||% input$exposure
  outcome <- outcome %||% input$outcome
  if (is.null(exposure) || is.null(outcome)) {
    stop("prepared analyses require exposure and outcome trait names", call. = FALSE)
  }
  validation_input <- input
  if (is.null(input$selected_loci) && is.character(loci) && length(loci)) validation_input$selected_loci <- loci
  validate_qtl2mr_input(validation_input, exposure = exposure, outcome = outcome, require_selected_loci = TRUE)
  sample_alignment <- qtl2mr_input_alignment(input)
  if (!identical(sample_alignment, input$sample_alignment)) qtl2mr_warn_alignment(sample_alignment)
  selected <- qtl2mr_prepared_selected_loci(input, loci = loci)
  input <- qtl2mr_prepared_probability_copy(input, selected$marker)
  x <- qtl2mr_prepared_trait(input, exposure)
  y <- qtl2mr_prepared_trait(input, outcome)
  covar <- qtl2mr_prepared_covar(input)
  founders <- dimnames(input$allele_probs[[1L]])[[2L]]
  contrasts <- make_founder_contrasts(founders)
  sample_ids_x_filter <- qtl2mr_prepared_role_filter(input, "exposure")
  sample_ids_y_filter <- qtl2mr_prepared_role_filter(input, "outcome")
  if (identical(engine, "batch")) {
    chrs <- unique(as.character(selected$chr))
    hsq_x <- stats::setNames(lapply(chrs, function(chr) {
      qtl2mr_prepared_hsq(input, chr, "exposure", exposure)
    }), chrs)
    hsq_y <- stats::setNames(lapply(chrs, function(chr) {
      qtl2mr_prepared_hsq(input, chr, "outcome", outcome)
    }), chrs)
    blocks <- fit_mr_locus_blocks_batch(
      genoprobs = input$allele_probs,
      selected_loci = selected,
      exposure = x,
      outcome = y,
      kinship = input$kinship,
      addcovar = covar,
      contrasts = contrasts,
      hsq_x = hsq_x,
      hsq_y = hsq_y,
      sample_ids_x = sample_ids_x_filter,
      sample_ids_y = sample_ids_y_filter,
      standardize_phenotypes = standardize_phenotypes,
      cross_cov_model = cross_cov_model,
      ridge = ridge,
      progress_callback = progress_callback,
      context_cache = context_cache,
	      fit_cache = fit_cache,
	      exposure_name = exposure,
	      outcome_name = outcome,
	      drop_nonestimable_loci = drop_nonestimable_loci
	    )
	  } else {
    blocks <- list()
    for (i in seq_len(nrow(selected))) {
      marker <- selected$marker[[i]]
      chr <- selected$chr[[i]]
      qtl2mr_progress(
        progress_callback,
        2L,
        "Fit selected-locus founder LMM/MR blocks",
        sprintf("%d/%d %s chr%s", i, nrow(selected), marker, chr)
      )
      if (isTRUE(verbose)) message(sprintf("[%d/%d] %s chr%s", i, nrow(selected), marker, chr))
      K <- qtl2mr_prepared_kinship(input, chr)
      blocks[[marker]] <- fit_mr_locus_block_overlap(
        genoprobs = input$allele_probs,
        exposure = x,
        outcome = y,
        kinship_x = K,
        kinship_y = K,
        addcovar_x = covar,
        addcovar_y = covar,
        contrasts = contrasts,
        hsq_x = qtl2mr_prepared_hsq(input, chr, "exposure", exposure),
        hsq_y = qtl2mr_prepared_hsq(input, chr, "outcome", outcome),
        marker = marker,
        sample_ids_x = sample_ids_x_filter,
        sample_ids_y = sample_ids_y_filter,
        standardize_phenotypes = standardize_phenotypes,
        cross_cov_model = cross_cov_model,
        ridge = ridge
      )
      os <- blocks[[marker]]$O_xy_summary
      qtl2mr_progress(
        progress_callback,
        2L,
        "Fit selected-locus founder LMM/MR blocks",
        sprintf(
          "%d/%d %s overlap=%d/%d X, %d Y",
          i, nrow(selected), marker,
          os$n_overlap %||% 0L,
          os$n_x %||% NA_integer_,
          os$n_y %||% NA_integer_
        )
	      )
	    }
	  }
  dropped_loci <- attr(blocks, "dropped_loci") %||% data.frame(
    marker = character(),
    chr = character(),
    role = character(),
    drop_reason = character(),
    error_class = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
  n_loci_selected <- attr(blocks, "n_loci_selected") %||% nrow(selected)
  n_loci_dropped <- attr(blocks, "n_loci_dropped") %||% length(unique(dropped_loci$marker))
  selected_retained <- selected[selected$marker %in% names(blocks), , drop = FALSE]
  if (nrow(dropped_loci)) {
    qtl2mr_progress(
      progress_callback,
      2L,
      "Fit selected-locus founder LMM/MR blocks",
      sprintf("dropped %d non-estimable loci; retained %d/%d",
              n_loci_dropped, length(blocks), n_loci_selected)
    )
  }
	  fit <- fit_multilocus_founder_mr(
	    blocks,
    allow_cross_cov = TRUE,
    ridge = ridge,
    heterogeneity_method = heterogeneity_method,
    q_sim_n = q_sim_n,
    q_sim_seed = q_sim_seed,
    progress_callback = progress_callback
  )
  diagnostic_layer <- fit_multilocus_diagnostic_layer(
    blocks,
    global_fit = fit,
    ridge = ridge,
    bootstrap_B = diagnostic_bootstrap_B,
    bootstrap_seed = diagnostic_bootstrap_seed
  )
  reported <- qtl2mr_reported_estimate(
    fit = fit,
    diagnostic_layer = diagnostic_layer,
    policy = reported_estimator_policy
  )
  get_est <- function(name, col) {
    rows <- diagnostic_layer$estimates$estimator == name
    if (!any(rows)) return(NA_real_)
    diagnostic_layer$estimates[[col]][which(rows)[1L]]
  }
  get_est_chr <- function(name, col) {
    rows <- diagnostic_layer$estimates$estimator == name
    if (!any(rows) || !col %in% names(diagnostic_layer$estimates)) return(NA_character_)
    as.character(diagnostic_layer$estimates[[col]][which(rows)[1L]])
  }
  backend_provenance <- qtl2mr_lmm_backend_provenance(blocks)
  phenotype_scaling <- do.call(rbind, lapply(blocks, function(block) {
    do.call(rbind, lapply(c("exposure", "outcome"), function(role) {
      fitted <- if (role == "exposure") block$fit_x else block$fit_y
      values <- if (role == "exposure") x[fitted$sample_ids] else y[fitted$sample_ids]
      data.frame(marker = block$marker, role = role,
                 trait = if (role == "exposure") exposure else outcome,
                 n = length(values), mean = mean(values), sd = stats::sd(values),
                 standardized = isTRUE(standardize_phenotypes), stringsAsFactors = FALSE)
    }))
  }))
  sample_alignment <- do.call(rbind, lapply(blocks, function(block) {
    do.call(rbind, lapply(c("exposure", "outcome"), function(role) {
      requested <- qtl2mr_prepared_role_filter(input, role)
      fitted <- if (role == "exposure") block$fit_x$sample_ids else block$fit_y$sample_ids
      values <- if (role == "exposure") x else y
      missing_pheno <- !requested %in% names(values)
      nonfinite_pheno <- !missing_pheno & !is.finite(values[requested])
      missing_covar <- if (is.null(covar)) rep(FALSE, length(requested)) else !requested %in% rownames(covar)
      incomplete_covar <- if (is.null(covar) || !ncol(covar)) rep(FALSE, length(requested)) else {
        !missing_covar & !qtl2mr_complete_rows(covar[match(requested, rownames(covar)), , drop = FALSE])
      }
      chr_ids <- dimnames(input$allele_probs[[as.character(block$chr)]])[[1L]]
      data.frame(marker = block$marker, role = role, n_requested = length(requested),
                 n_fitted = length(fitted), n_dropped = length(setdiff(requested, fitted)),
                 n_missing_phenotype_id = sum(missing_pheno), n_nonfinite_phenotype = sum(nonfinite_pheno),
                 n_missing_covar_id = sum(missing_covar), n_incomplete_covar = sum(incomplete_covar),
                 n_missing_probability_id = sum(!requested %in% chr_ids),
                 dropped_ids = paste(setdiff(requested, fitted), collapse = ";"),
                 stringsAsFactors = FALSE)
    }))
  }))
  covar_excluded <- sample_alignment$n_incomplete_covar + sample_alignment$n_missing_covar_id
  if (any(covar_excluded > 0L)) {
    counts <- tapply(covar_excluded, sample_alignment$role, max)
    warning("Missing covariate rows or non-finite covariate values excluded samples from fitting (maximum per locus: ",
      paste(paste0(names(counts), "=", counts), collapse = ", "),
      "); inspect sample_alignment (n_missing_covar_id, n_incomplete_covar) and the fitted sample sets.", call. = FALSE)
  }
  sd_x <- phenotype_scaling$sd[phenotype_scaling$role == "exposure"]
  sd_y <- phenotype_scaling$sd[phenotype_scaling$role == "outcome"]
  common_sd <- function(z) length(z) && all(is.finite(z)) && all(abs(z - z[[1L]]) <= 1e-10 * max(abs(z)))
  raw_conversion_available <- !isTRUE(standardize_phenotypes) || (common_sd(sd_x) && common_sd(sd_y))
  raw_multiplier <- if (!isTRUE(standardize_phenotypes)) 1 else if (raw_conversion_available) sd_y[[1L]] / sd_x[[1L]] else NA_real_
  summary <- data.frame(
    exposure = exposure,
    outcome = outcome,
    lmm_backend = backend_provenance$backend,
    n_loci = length(blocks),
    n_loci_selected = n_loci_selected,
    n_loci_dropped_nonestimable = n_loci_dropped,
    theta_so = fit$theta,
    se_so = fit$se,
    ci_lower_so = fit$profile_ci_lower,
    ci_upper_so = fit$profile_ci_upper,
    theta = fit$theta,
    se = fit$se,
    profile_ci_lower = fit$profile_ci_lower,
    profile_ci_upper = fit$profile_ci_upper,
    profile_p = fit$profile_p,
    primary_optimizer_status = fit$optimizer_status %||% NA_character_,
    primary_bracket_status = fit$bracket_status %||% NA_character_,
    primary_bracket_n_expand = fit$bracket_n_expand %||% NA_integer_,
    primary_profile_ci_finite = isTRUE(fit$profile_ci_finite),
    primary_reliable = isTRUE(fit$primary_reliable),
    numerical_ok = isTRUE(fit$numerical_ok),
    confidence_set_type = fit$confidence_set_type %||% "unresolved",
    inference_status = fit$inference_status %||% "not_calibrated",
    effect_units = if (isTRUE(standardize_phenotypes)) "SD_y per SD_x" else "outcome units per exposure unit",
    exposure_sd = if (common_sd(sd_x)) sd_x[[1L]] else NA_real_,
    outcome_sd = if (common_sd(sd_y)) sd_y[[1L]] else NA_real_,
    raw_conversion_available = raw_conversion_available,
    theta_raw_units = fit$theta * raw_multiplier,
    se_raw_units = fit$se * abs(raw_multiplier),
    profile_ci_lower_raw_units = fit$profile_ci_lower * raw_multiplier,
    profile_ci_upper_raw_units = fit$profile_ci_upper * raw_multiplier,
    first_order_theta = fit$first_order_theta %||% NA_real_,
    first_order_se = fit$first_order_se %||% NA_real_,
    theta_reported = reported$theta,
    se_reported = reported$se,
    reported_ci_lower = reported$ci_lower,
    reported_ci_upper = reported$ci_upper,
    reported_ci_type = reported$ci_type,
    reported_estimator = reported$estimator,
    reported_estimator_source = reported$source,
    reported_estimator_policy = reported$policy,
    reported_estimator_reason = reported$reason,
    theta_independence = fit$independence$theta,
    theta_overlap_delta = fit$theta_overlap_delta,
    theta_two_step_ivw = get_est("Two-step IVW", "theta"),
    se_two_step_ivw = get_est("Two-step IVW", "se"),
    theta_weighted_median = get_est("Weighted median", "theta"),
    se_weighted_median = get_est("Weighted median", "se"),
    ci_lower_weighted_median = get_est("Weighted median", "ci_lower"),
    ci_upper_weighted_median = get_est("Weighted median", "ci_upper"),
    ci_type_weighted_median = get_est_chr("Weighted median", "ci_type"),
    theta_weighted_mode = get_est("Weighted mode", "theta"),
    se_weighted_mode = get_est("Weighted mode", "se"),
    ci_lower_weighted_mode = get_est("Weighted mode", "ci_lower"),
    ci_upper_weighted_mode = get_est("Weighted mode", "ci_upper"),
    ci_type_weighted_mode = get_est_chr("Weighted mode", "ci_type"),
    diagnostic_bootstrap_B = as.integer(diagnostic_bootstrap_B %||% 0L),
    diagnostic_bootstrap_seed = diagnostic_bootstrap_seed %||% NA_integer_,
    Q_between = diagnostic_layer$heterogeneity$Q_between[[1L]],
    Q_between_df = diagnostic_layer$heterogeneity$df_between[[1L]],
    Q_between_p = diagnostic_layer$heterogeneity$Q_between_p[[1L]],
    Q = fit$Q,
    Q_chisq_df = fit$Q_chisq_df %||% NA_integer_,
    Q_p_chisq = fit$Q_p_chisq %||% NA_real_,
    Q_p_sim = fit$Q_p_sim %||% NA_real_,
    Q_sim_n = fit$Q_sim_n %||% NA_integer_,
    Q_sim_n_requested = fit$Q_sim_n_requested %||% NA_integer_,
    Q_sim_seed = fit$Q_sim_seed %||% NA_integer_,
    Q_sim_exceedances = fit$Q_sim_exceedances %||% NA_integer_,
    Q_sim_min_resolvable_p = fit$Q_sim_min_resolvable_p %||% NA_real_,
    Q_sim_null_mean = fit$Q_sim_null_mean %||% NA_real_,
    Q_sim_null_sd = fit$Q_sim_null_sd %||% NA_real_,
    Q_sim_method = fit$Q_sim_method %||% NA_character_,
    Q_p_reported = fit$Q_p_reported %||% NA_real_,
    Q_p_method = fit$Q_p_method %||% NA_character_,
    Q_p_note = fit$Q_p_note %||% NA_character_,
    cross_cov_model = cross_cov_model,
    n_loci_psd_adjusted = sum(fit$locus_diagnostics$psd_adjusted, na.rm = TRUE),
    n_loci_psd_invalid = sum(!fit$locus_diagnostics$joint_psd_valid, na.rm = TRUE),
    n_loci_raw_space_psd_adjusted = sum(fit$locus_diagnostics$raw_space_psd_adjusted, na.rm = TRUE),
    n_loci_raw_space_psd_invalid = sum(!fit$locus_diagnostics$raw_space_psd_valid, na.rm = TRUE),
    n_loci_star_space_psd_adjusted = sum(fit$locus_diagnostics$star_space_psd_adjusted, na.rm = TRUE),
    n_loci_star_space_psd_invalid = sum(!fit$locus_diagnostics$star_space_psd_valid, na.rm = TRUE),
    rho_xy_median = stats::median(fit$locus_diagnostics$rho_xy, na.rm = TRUE),
    rho_xy_raw_median = stats::median(fit$locus_diagnostics$rho_xy_raw, na.rm = TRUE),
    rho_xy_star_median = stats::median(fit$locus_diagnostics$rho_xy_star, na.rm = TRUE),
    max_cross_singular_value_median = stats::median(fit$locus_diagnostics$joint_max_cross_singular_value, na.rm = TRUE),
    raw_space_max_cross_singular_value_median = stats::median(fit$locus_diagnostics$raw_space_max_cross_singular_value, na.rm = TRUE),
    star_space_max_cross_singular_value_median = stats::median(fit$locus_diagnostics$star_space_max_cross_singular_value, na.rm = TRUE),
    standardize_phenotypes = isTRUE(standardize_phenotypes),
    stringsAsFactors = FALSE
  )
  summary <- qtl2mr_add_q_fields(summary, fit = fit)
  overlap_summary <- qtl2mr_locus_overlap_table(blocks)
  if (nrow(overlap_summary)) {
    summary$n_x_median <- stats::median(overlap_summary$n_x, na.rm = TRUE)
    summary$n_y_median <- stats::median(overlap_summary$n_y, na.rm = TRUE)
    summary$n_overlap_median <- stats::median(overlap_summary$n_overlap, na.rm = TRUE)
    summary$exposure_overlap_fraction_median <- stats::median(overlap_summary$exposure_overlap_fraction, na.rm = TRUE)
    summary$outcome_overlap_fraction_median <- stats::median(overlap_summary$outcome_overlap_fraction, na.rm = TRUE)
    summary$jaccard_overlap_fraction_median <- stats::median(overlap_summary$jaccard_overlap_fraction, na.rm = TRUE)
    summary$n_loci_independent_samples <- sum(overlap_summary$independent_samples, na.rm = TRUE)
    summary$n_loci_complete_overlap <- sum(overlap_summary$complete_overlap, na.rm = TRUE)
  }
  summary$sample_policy <- "role_specific_complete"
  summary$confidence_set_bounded <- qtl2mr_confidence_set_bounded(fit)
  hsq_table <- qtl2mr_hsq_table(blocks, exposure = exposure, outcome = outcome)
  founder_residuals <- do.call(rbind, lapply(blocks, founder_residual_scores, theta = fit$theta))
  founder_ratios <- do.call(rbind, lapply(blocks, founder_coordinate_ratios, theta_ref = fit$theta))
  out <- list(
	    summary = summary,
	    fit = fit,
	    locus_blocks = blocks,
	    selected_loci = selected_retained,
	    selected_loci_original = selected,
	    dropped_loci = dropped_loci,
	    locus_diagnostics = fit$locus_diagnostics,
    diagnostic_layer = diagnostic_layer,
    overlap_summary = overlap_summary,
    sample_alignment = sample_alignment,
    input_alignment = qtl2mr_input_alignment(input),
    phenotype_scaling = phenotype_scaling,
    founder_residual_scores_global = founder_residuals,
    founder_coordinate_ratios = founder_ratios,
    hsq = hsq_table,
    backend_provenance = backend_provenance,
    sample_ids = unique(c(
      unlist(lapply(blocks, function(b) b$fit_x$sample_ids), use.names = FALSE),
      unlist(lapply(blocks, function(b) b$fit_y$sample_ids), use.names = FALSE)
    )),
    overlap_sample_ids = Reduce(intersect, lapply(blocks, function(b) intersect(b$fit_x$sample_ids, b$fit_y$sample_ids))),
    prepared_metadata = input$metadata
  )
  class(out) <- c("qtl2mr_prepared_pilot", "list")
  if (isTRUE(write_outputs) && !is.null(output_dir)) {
	    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	    utils::write.csv(summary, file.path(output_dir, "prepared_multilocus_founder_mr_summary.csv"), row.names = FALSE)
	    utils::write.csv(selected_retained, file.path(output_dir, "prepared_selected_loci.csv"), row.names = FALSE)
	    utils::write.csv(selected, file.path(output_dir, "prepared_selected_loci_original.csv"), row.names = FALSE)
	    utils::write.csv(dropped_loci, file.path(output_dir, "prepared_dropped_loci.csv"), row.names = FALSE)
	    utils::write.csv(overlap_summary, file.path(output_dir, "prepared_overlap_summary.csv"), row.names = FALSE)
    saveRDS(fit, file.path(output_dir, "prepared_multilocus_founder_mr.rds"))
  }
  out
}

fit_mr_locus_block_overlap <- function(genoprobs,
                                       exposure,
                                       outcome,
                                       kinship_x = NULL,
                                       kinship_y = NULL,
                                       addcovar_x = NULL,
                                       addcovar_y = NULL,
                                       contrasts = NULL,
                                       hsq_x = NULL,
                                       hsq_y = NULL,
                                       marker = NULL,
                                       sample_ids_x = NULL,
                                       sample_ids_y = NULL,
                                       overlap_ids = NULL,
                                       O_xy = NULL,
                                       standardize_phenotypes = FALSE,
                                       reml = TRUE,
                                       fc_cols = paste0("FC", 1:7),
                                       ridge = 1e-10,
                                       cross_cov_model = c("star", "raw", "zero"),
                                       enforce_psd = TRUE,
                                       stop_on_invalid_psd = TRUE,
                                       tol = 1e-12) {
  cross_cov_model <- match.arg(cross_cov_model)
  selected <- qtl2mr_select_marker_genoprobs(genoprobs, marker = marker)
  fit_x <- fit_founder_lmm_block(
    genoprobs = genoprobs,
    pheno = exposure,
    kinship = kinship_x,
    addcovar = addcovar_x,
    contrasts = contrasts,
    hsq = hsq_x,
    marker = selected$marker,
    sample_ids = sample_ids_x,
    standardize_pheno = standardize_phenotypes,
    reml = reml,
    fc_cols = fc_cols,
    tol = tol
  )
  fit_y <- fit_founder_lmm_block(
    genoprobs = genoprobs,
    pheno = outcome,
    kinship = kinship_y,
    addcovar = addcovar_y,
    contrasts = contrasts,
    hsq = hsq_y,
    marker = selected$marker,
    sample_ids = sample_ids_y,
    standardize_pheno = standardize_phenotypes,
    reml = reml,
    fc_cols = fc_cols,
    tol = tol
  )
  if (is.null(O_xy)) {
    O_xy <- make_overlap_matrix(fit_x$sample_ids, fit_y$sample_ids, overlap_ids = overlap_ids)
  } else {
    O_xy <- as.matrix(O_xy)
    if (nrow(O_xy) != length(fit_x$sample_ids) || ncol(O_xy) != length(fit_y$sample_ids)) {
      stop("O_xy dimensions do not match exposure/outcome fitted sample IDs", call. = FALSE)
    }
    if (is.null(rownames(O_xy))) rownames(O_xy) <- fit_x$sample_ids
    if (is.null(colnames(O_xy))) colnames(O_xy) <- fit_y$sample_ids
  }
  cross_raw <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "raw", O_xy = O_xy)
  cross_star <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "star", O_xy = O_xy)
  cross_zero <- qtl2mr_lmm_cross_covariance(fit_x, fit_y, fc_cols, space = "zero", O_xy = O_xy)
  selected_cross <- switch(cross_cov_model, raw = cross_raw, star = cross_star, zero = cross_zero)
  Sigma_x <- fit_x$vcov_fc
  Sigma_y <- fit_y$vcov_fc
  psd_raw <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_raw$Sigma_xy, ridge = ridge)
  psd_star <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_star$Sigma_xy, ridge = ridge)
  psd_zero <- qtl2mr_enforce_joint_psd(Sigma_x, Sigma_y, cross_zero$Sigma_xy, ridge = ridge)
  psd <- switch(cross_cov_model, raw = psd_raw, star = psd_star, zero = psd_zero)
  Sigma_xy <- selected_cross$Sigma_xy
  if (isTRUE(enforce_psd)) {
    if (!isTRUE(psd$psd_valid) && isTRUE(stop_on_invalid_psd)) {
      stop("Selected ", cross_cov_model, " Sigma_xy fails joint PSD validation at marker ", selected$marker,
           " (status=", psd$status, ")", call. = FALSE)
    }
    Sigma_x <- psd$Sigma_x
    Sigma_y <- psd$Sigma_y
    Sigma_xy <- psd$Sigma_xy
  }
  list(
    marker = selected$marker,
    chr = selected$chr,
    beta_x = fit_x$beta_fc,
    beta_y = fit_y$beta_fc,
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    Sigma_xy = Sigma_xy,
    Sigma_x_raw = fit_x$vcov_fc,
    Sigma_y_raw = fit_y$vcov_fc,
    Sigma_xy_raw = selected_cross$Sigma_xy,
    Sigma_xy_raw_space = cross_raw$Sigma_xy,
    Sigma_xy_star_space = cross_star$Sigma_xy,
    Sigma_xy_zero = cross_zero$Sigma_xy,
    rho_xy = selected_cross$rho_xy,
    rho_star = selected_cross$rho_star,
    rho_xy_raw = cross_raw$rho_xy,
    rho_xy_star = cross_star$rho_xy,
    n_overlap = selected_cross$n_overlap,
    n_overlap_raw = cross_raw$n_overlap,
    n_overlap_star = cross_star$n_overlap,
    O_xy_summary = selected_cross$O_xy_summary,
    fit_x = fit_x,
    fit_y = fit_y,
    founders = fit_x$founders,
    contrasts = fit_x$contrasts,
    founder_H = fit_x$founder_H,
    hsq_x = fit_x$hsq,
    hsq_y = fit_y$hsq,
    cross_cov_model = cross_cov_model,
    psd_adjusted = isTRUE(enforce_psd) && isTRUE(psd$adjusted),
    joint_psd_candidate_adjusted = isTRUE(psd$adjusted),
    joint_psd_valid = isTRUE(psd$psd_valid),
    joint_psd_status = psd$status,
    joint_psd_large_violation = isTRUE(psd$large_violation),
    joint_min_eigenvalue = psd$min_eigenvalue,
    joint_max_cross_singular_value = psd$max_cross_singular_value %||% NA_real_,
    raw_space_psd_adjusted = isTRUE(psd_raw$adjusted),
    raw_space_psd_valid = isTRUE(psd_raw$psd_valid),
    raw_space_psd_status = psd_raw$status,
    raw_space_psd_large_violation = isTRUE(psd_raw$large_violation),
    raw_space_joint_min_eigenvalue = psd_raw$min_eigenvalue,
    raw_space_max_cross_singular_value = psd_raw$max_cross_singular_value %||% NA_real_,
    star_space_psd_adjusted = isTRUE(psd_star$adjusted),
    star_space_psd_valid = isTRUE(psd_star$psd_valid),
    star_space_psd_status = psd_star$status,
    star_space_psd_large_violation = isTRUE(psd_star$large_violation),
    star_space_joint_min_eigenvalue = psd_star$min_eigenvalue,
    star_space_max_cross_singular_value = psd_star$max_cross_singular_value %||% NA_real_,
    standardize_phenotypes = isTRUE(standardize_phenotypes)
  )
}
