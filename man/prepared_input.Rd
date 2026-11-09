\name{prepared_input}
\alias{prepared_input}
\alias{validate_qtl2mr_input}
\alias{fit_mr_locus_block_overlap}
\title{Validate the input for one analysis}
\description{
\if{latex}{\out{\sloppy}}
Validate a single-pair input made by \code{\link{qtl2mr_pair_input}()}, and fit
a single locus when exposure and outcome sample sets differ.
}
\usage{
validate_qtl2mr_input(input, founder_dim = 8L,
  require_selected_loci = FALSE, exposure = NULL, outcome = NULL)

fit_mr_locus_block_overlap(genoprobs, exposure, outcome,
  kinship_x = NULL, kinship_y = NULL, addcovar_x = NULL,
  addcovar_y = NULL, contrasts = NULL, hsq_x = NULL, hsq_y = NULL,
  marker = NULL, sample_ids_x = NULL, sample_ids_y = NULL,
  overlap_ids = NULL, O_xy = NULL, standardize_phenotypes = FALSE,
  reml = TRUE, fc_cols = paste0("FC", 1:7), ridge = 1e-10,
  cross_cov_model = c("star", "raw", "zero"), enforce_psd = TRUE,
  stop_on_invalid_psd = TRUE, tol = 1e-12)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{input}{A single-pair input from \code{\link{qtl2mr_pair_input}()}.}
\item{founder_dim}{Expected founder dimension; defaults to eight founders.}
\item{require_selected_loci}{Whether validation should require
\code{selected_loci}.}
\item{exposure, outcome}{For \code{validate_qtl2mr_input()}, optional exposure
and outcome trait names to check against the input; they must be two different
traits. For \code{fit_mr_locus_block_overlap()}, named exposure and outcome
phenotype vectors.}
\item{genoprobs}{qtl2-style allele probability object: a named list of sample by
founder by marker arrays, for example \code{qtl2::genoprob_to_alleleprob()}
output.}
\item{kinship_x, kinship_y}{Optional kinship matrices for the exposure and the
outcome fit.}
\item{addcovar_x, addcovar_y}{Optional exposure/outcome covariates.}
\item{contrasts}{Founder contrast matrix.}
\item{hsq_x, hsq_y}{Optional finite heritability values in [0, 1]; out-of-range
values are errors, never clamped. If omitted and kinship is supplied, the
internal profile ML/REML estimator is used.}
\item{marker}{Marker ID.}
\item{sample_ids_x, sample_ids_y}{Optional sample order for the exposure and
the outcome fit.}
\item{overlap_ids}{IDs present in both exposure and outcome samples.}
\item{O_xy}{Explicit exposure-by-outcome overlap selector matrix for low-level
experiments only. The high-level workflow requires harmonized sample IDs and
derives overlap internally.}
\item{standardize_phenotypes}{Whether to standardize exposure and outcome before
fitting.}
\item{reml}{Whether internal heritability estimation uses REML.}
\item{fc_cols}{Founder contrast columns.}
\item{ridge}{Numerical ridge used in PSD and inverse operations.}
\item{cross_cov_model}{Cross-covariance model. The high-level analyses use
\code{"star"}; \code{"raw"} is diagnostic only.}
\item{enforce_psd, stop_on_invalid_psd}{Joint PSD enforcement controls.}
\item{tol}{QR/rank tolerance.}
}
\details{
\if{latex}{\out{\sloppy}}
A single-pair input is built in two steps:
\code{\link{qtl2mr_prepare_resources}()} holds the shared allele
probabilities, kinship, phenotypes, covariates and marker map, and
\code{\link{qtl2mr_pair_input}()} adds the instruments of one exposure and the
pair's trait names, after the outcome founder QC. One input represents one MR
analysis. The high-level fit always uses the general partial-overlap engine;
full overlap and no overlap are special cases. Users must harmonize sample IDs
before calling qtl2mr.

The validator checks sample IDs, covariates, kinship names and coverage, marker
IDs, selected loci, and the eight-founder autosomal additive assumption.
Unknown IDs and incomplete matching are not silently treated as intentional
sample exclusions. Partial exposure/outcome overlap is supported and recorded.
Chromosome X is unsupported and cannot be fitted as an autosome.

The input retains the original probability arrays unchanged. Its
\code{metadata$probability_normalization} records the tolerances and the
per-marker adjustment counts and maximum original sum deviations of the
rounding normalisation described in \code{\link{qtl2mr_prepare_resources}}. An
adjustment emits a warning and is applied only to an analysis-local copy; the
result retains this record in \code{prepared_metadata}.

The analysis standardizes phenotypes: the effect is SD of outcome per SD of
exposure, and raw-unit columns are added when one conversion applies. Scaling
and sample-accounting tables are retained with the result. Per-locus
differences in fitting samples can produce different SDs, so there may not be
one valid pooled raw-unit conversion.

\code{fit_mr_locus_block_overlap()} is a low-level helper for partial or nested
sample-overlap experiments. It fits exposure and outcome founder LMMs on their
own sample sets and constructs star-space \code{Sigma_xy} through the explicit
\code{O_xy} bilinear route.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
d <- qtl2mr_example
resources <- qtl2mr_prepare_resources(d$allele_probs, d$kinship, d$pheno,
  d$covar, d$marker_map)
input <- qtl2mr_pair_input(resources,
  data.frame(exposure = "x", marker = d$selected_loci), "x", "y")
validate_qtl2mr_input(input, require_selected_loci = TRUE)
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
\code{validate_qtl2mr_input()} expects analysis-ready components;
\code{vignette("from-qtl2", package = "qtl2mr")} shows how to produce them
from a qtl2 cross. The current scope is eight-founder populations, autosomal
additive loci, continuous traits, additive covariates, selected markers and
one kinship source. \code{fit_mr_locus_block_overlap()} is part of the
advanced low-level interface.
}
\value{
\if{latex}{\out{\sloppy}}
\code{validate_qtl2mr_input()} returns \code{TRUE} invisibly.
\code{fit_mr_locus_block_overlap()} returns a locus block compatible with
\code{fit_multilocus_founder_mr()}.
}
\keyword{methods}
