# Opt-in founder-design QC (min_founder_copies). A locus is dropped for a trait
# when the rarest founder's expected allele count on the samples that trait is
# fitted on, E = min_k sum_i 2 p_ik, is below the threshold. The probabilities
# are those the package fits: after the probability policy's normalisation.
# With min_founder_copies = NULL (default) nothing here runs.

qtl2mr_check_min_founder_copies <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || is.logical(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop("min_founder_copies must be NULL or one finite number >= 0", call. = FALSE)
  }
  as.numeric(x)
}

# The samples prepare_lmm_context() fits a trait on (same filters, in the same
# order), without fitting anything: in the requested set, with a phenotype
# value, covariate row and kinship row, then finite phenotype and complete
# covariates. Signals nothing; an empty set is returned where the context would
# fail for lack of samples.
qtl2mr_fitted_sample_ids <- function(pheno, kinship = NULL, addcovar = NULL, sample_ids = NULL,
                                     trait_name = NULL, tol = 1e-12) {
  qtl2mr_quietly({
    pheno0 <- qtl2mr_as_named_numeric(pheno, trait_name %||% "phenotype")
    base_ids <- if (is.null(sample_ids)) names(pheno0) else as.character(sample_ids)
    ids <- base_ids[base_ids %in% names(pheno0)]
    if (!is.null(addcovar) && nrow(as.matrix(addcovar))) {
      ids <- ids[ids %in% rownames(as.matrix(addcovar))]
    }
    if (!is.null(kinship)) ids <- ids[ids %in% rownames(as.matrix(kinship))]
    if (!length(ids)) return(character())
    addcovar0 <- qtl2mr_as_addcovar(addcovar, ids = ids, tol = tol)
    keep <- is.finite(pheno0[ids])
    if (ncol(addcovar0)) keep <- keep & qtl2mr_complete_rows(addcovar0)
    ids[keep]
  })
}

# E_min = min over founders of sum_i 2 p_ik over the given samples at one marker.
qtl2mr_founder_copies <- function(genoprobs_chr, marker, sample_ids) {
  idx <- match(marker, dimnames(genoprobs_chr)[[3L]])
  if (is.na(idx)) stop("Marker not found in chromosome allele probabilities: ", marker, call. = FALSE)
  P <- matrix(genoprobs_chr[sample_ids, , idx], nrow = length(sample_ids))
  min(2 * colSums(P))
}

qtl2mr_chr_kinship <- function(kinship, chr) {
  if (is.list(kinship) && !is.matrix(kinship) && !is.data.frame(kinship)) kinship[[chr]] else kinship
}

# E_min of each locus for one trait of a prepared input whose probabilities are
# already normalised at these markers, on the samples the trait is fitted with
# in that role (its own complete samples, or the role samples of the input).
qtl2mr_trait_founder_copies <- function(input, trait, role, loci) {
  out <- rep(NA_real_, nrow(loci))
  if (!nrow(loci)) return(out)
  values <- qtl2mr_prepared_trait(input, trait)
  covar <- qtl2mr_prepared_covar(input)
  requested <- qtl2mr_prepared_role_filter(input, role)
  for (chr in unique(as.character(loci$chr))) {
    rows <- which(as.character(loci$chr) == chr)
    ids <- qtl2mr_fitted_sample_ids(values, qtl2mr_chr_kinship(input$kinship, chr), covar, requested, trait)
    for (j in rows) out[[j]] <- qtl2mr_founder_copies(input$allele_probs[[chr]], loci$marker[[j]], ids)
  }
  out
}

# Outcome-side QC of one pair: E_min of every instrument of the exposure on the
# samples the outcome is fitted with (own complete, or the pair's common
# complete samples), from a prepared input normalised at these markers.
# Returns the instrument rows with E_min_y and pass.
qtl2mr_outcome_founder_qc <- function(input, outcome, rows, threshold) {
  E <- qtl2mr_trait_founder_copies(input, outcome, "outcome", rows)
  data.frame(marker = as.character(rows$marker), chr = as.character(rows$chr),
    pos = if ("pos" %in% names(rows)) rows$pos else NA_real_, E_min_y = E, pass = !(E < threshold),
    stringsAsFactors = FALSE)
}

# The same QC for a pair built from resources, as qtl2mr_pair_input() and
# run_qtl2mr_batch() apply it: the pair input for all instruments, normalised as
# the fit would normalise it, without signalling anything.
qtl2mr_pair_founder_qc <- function(resources, rows, exposure, outcome, threshold) {
  input <- qtl2mr_quietly(qtl2mr_prepared_probability_copy(
    qtl2mr_pair_prepared(resources, rows, exposure, outcome), rows$marker))
  qtl2mr_outcome_founder_qc(input, outcome, rows, threshold)
}
