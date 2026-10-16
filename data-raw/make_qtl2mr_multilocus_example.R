# Extract the frozen, wholly synthetic Case A X1-to-Y teaching input.
# No phenotype generation, resampling, locus search or parameter tuning occurs.
# Usage: Rscript data-raw/make_qtl2mr_multilocus_example.R INPUT.rds OUTPUT.rda
#
# Source: CTC tutorial construction, 2026-10-15, v0.4 background-QTL Case A.
# Generation used master-seed provenance 2026101501: synthetic probabilities
# and family-AR kinship, X-only first-stage construction, then one algebraic
# observational-target solve before MR. The full workshop history is not a
# runtime dependency. The supplied frozen input is needed only for re-extraction.

extract_qtl2mr_multilocus_example <- function(input_file, output_file) {
  expected_md5 <- "3737665d82440de2cc9fed5b2c64b5a9"
  source_sha256 <- "23cf8f6a2c38637bfea6a0c9855f234247c6345ed3d46e8ff21271cd79ffb309"
  stopifnot(file.exists(input_file),
            identical(unname(tools::md5sum(input_file)), expected_md5))
  original <- readRDS(input_file)
  ids <- sprintf("synthetic_mouse_%03d", seq_len(400L))
  stopifnot(identical(rownames(original$pheno), ids),
            identical(rownames(original$covar), ids),
            identical(colnames(original$covar), c("batch", "age")),
            identical(dim(original$kinship), c(400L, 400L)),
            identical(dimnames(original$kinship), list(ids, ids)))

  selected <- original$marker_map$assigned_exposure == "x1" &
    original$marker_map$role == "iv"
  marker_map <- original$marker_map[selected, c("marker", "chr", "pos"), drop = FALSE]
  rownames(marker_map) <- NULL
  stopifnot(nrow(marker_map) == 40L, !anyDuplicated(marker_map$marker))
  chromosomes <- names(original$allele_probs)
  allele_probs <- lapply(chromosomes, function(chr) {
    markers <- marker_map$marker[marker_map$chr == chr]
    original$allele_probs[[chr]][, , markers, drop = FALSE]
  })
  names(allele_probs) <- chromosomes
  class(allele_probs) <- class(original$allele_probs)
  pheno <- original$pheno[, c("x1", "y_A"), drop = FALSE]
  names(pheno) <- c("x1", "y")
  stopifnot(identical(pheno$x1, original$pheno$x1),
            identical(pheno$y, original$pheno$y_A),
            all(vapply(allele_probs, function(a) {
              identical(dimnames(a)[[1L]], ids) &&
                identical(dimnames(a)[[2L]], LETTERS[1:8]) &&
                all(is.finite(a) & a >= 0 & a <= 1) &&
                max(abs(apply(a, c(1L, 3L), sum) - 1)) < 1e-12
            }, logical(1))))

  qtl2mr_multilocus_example <- list(
    allele_probs = allele_probs,
    pheno = pheno,
    covar = original$covar,
    kinship = original$kinship,
    marker_map = marker_map,
    selected_loci = marker_map$marker,
    exposure = "x1",
    outcome = "y",
    true_theta = -0.25,
    metadata = list(
      purpose = "Wholly synthetic multi-locus illustration; not calibration evidence.",
      source_id = "CTC tutorial 2026-10-15, v0.4 background-QTL Case A",
      source_sha256 = source_sha256,
      master_seed = 2026101501L,
      n_samples = 400L,
      n_selected_loci = 40L,
      sample_overlap = "complete",
      standardize_phenotypes = FALSE,
      effect_units = "outcome units per exposure unit",
      kinship_note = "Synthetic family-AR matrix M/2: families of four, AR coefficient 0.3; not LOCO.",
      selection_note = "Forty prespecified X1 loci; no QTL scan or post-MR selection.",
      construction_note = "First-stage amplitudes were set using X only. The sample covariate-adjusted OLS slope was then targeted at 0.05 before MR.",
      subset_note = "X1 and the original Case A outcome are unchanged; omitted traits and background loci still contribute to the stored outcome.",
      limitations = "Engineered teaching fixture, not a realistic DO genome, independent-IV validation or confidence-interval calibration."
    )
  )
  if (file.exists(output_file)) {
    previous <- new.env(parent = emptyenv())
    stopifnot(identical(load(output_file, previous), "qtl2mr_multilocus_example"),
              identical(previous$qtl2mr_multilocus_example, qtl2mr_multilocus_example))
    message("Existing dataset is identical; no file overwritten.")
  } else {
    stopifnot(dir.exists(dirname(output_file)))
    save(qtl2mr_multilocus_example, file = output_file, compress = "xz", version = 2)
  }
  invisible(qtl2mr_multilocus_example)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Supply the frozen synthetic input and destination .rda paths.")
  extract_qtl2mr_multilocus_example(args[[1L]], args[[2L]])
}
