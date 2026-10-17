\name{prepared_input}
\alias{prepared_input}
\alias{as_qtl2mr_prepared_input}
\alias{validate_qtl2mr_input}
\alias{fit_mr_locus_block_overlap}
\title{Prepare and validate the input for one analysis}
\description{
\if{latex}{\out{\sloppy}}
Create and validate an explicit qtl2MR prepared-input object, and fit a
single locus when exposure and outcome sample sets differ.
}
\usage{
as_qtl2mr_prepared_input(allele_probs, pheno, covar = NULL,
  kinship = NULL, marker_map, selected_loci = NULL, exposure = NULL,
  outcome = NULL, exposure_samples = NULL, outcome_samples = NULL,
  hsq = NULL, sample_ids = NULL, metadata = list(), use_all_loci = FALSE,
  probability_policy = c("normalize_rounding", "strict"))

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
\item{allele_probs, genoprobs}{qtl2-style allele probability object: a named
list of sample by founder by marker arrays, for example
\code{qtl2::genoprob_to_alleleprob()} output. Checked values must be finite and
within [0, 1] up to floating-point error of 1e-8; genotype-state probabilities
(36 states) and allele dosages must be converted first.}
\item{probability_policy}{Prepared-input sum policy. The default
\code{"normalize_rounding"} accepts a maximum absolute founder-sum deviation
of \code{n_founders * 5e-5} (0.0004 for eight founders), the worst-case bound
for rounding each probability to four decimal places. Sample-marker vectors
whose deviation exceeds the legacy strict tolerance of \code{1e-6} are divided
by their sum in an analysis-local copy. \code{"strict"} rejects deviations above
\code{1e-6}. Neither policy clips values or repairs negative/nonfinite entries.}
\item{pheno}{Phenotype data frame with sample IDs as row names.}
\item{covar}{Optional covariate data frame or matrix with sample IDs as row
names.}
\item{kinship, kinship_x, kinship_y}{Kinship matrix, or chromosome-keyed list of
LOCO kinship matrices for prepared inputs. The public prepared adapter uses one
union-sample kinship source and subsets it separately for exposure and outcome
fits. Separate exposure/outcome kinships are not supported by the prepared adapter.}
\item{marker_map}{Marker map as a data frame, qtl2-style list, or named vector.
A data frame must contain \code{marker} and \code{chr}, or marker row names and
a chromosome column.}
\item{selected_loci}{Marker vector or data frame containing selected loci.
May be omitted while preparing resources, but analysis requires an explicit
selection or deliberate \code{use_all_loci = TRUE}.}
\item{use_all_loci}{Explicit opt-in to use every mapped marker when no
selected-locus list is supplied. Defaults to false. This is not a recommended
instrument-selection rule.}
\item{exposure, outcome}{For prepared-input construction and validation,
optional exposure and outcome trait names stored with the prepared object or
supplied to \code{run_qtl2mr_analysis()}; they must be two different traits.
An empty \code{selected_loci} is an error. For
\code{fit_mr_locus_block_overlap()}, named exposure and outcome phenotype
vectors.}
\item{exposure_samples, outcome_samples}{Optional role-specific sample filters.
If omitted, finite exposure/outcome values, allele probabilities, covariates,
and kinship availability define each LMM sample set. Overlap is inferred from
shared sample IDs.}
\item{hsq, hsq_x, hsq_y}{Optional finite numeric heritability values in [0,1];
out-of-range values are errors, never clamped. Low-level \code{hsq_x} and
\code{hsq_y} are scalars. Prepared \code{hsq} accepts an unnamed scalar (both
traits, all chromosomes), a uniquely named numeric vector or list keyed by
chromosome, phenotype name, or \code{exposure}/\code{outcome} role, or a two-level
named structure in either chromosome-to-role/trait or role/trait-to-chromosome
order. Examples are \code{0.3}, \code{c(exposure=0.2, outcome=0.3)},
\code{list("1"=0.2)}, and \code{list(exposure=c("1"=0.2, "2"=0.3))}.
Applicable entries must agree; conflicting entries, unknown keys, unnamed
multi-value inputs, and unsupported structures are errors. Omitted or NULL
entries request internal profile ML/REML estimation when kinship is supplied.
A nonzero supplied value requires kinship.}
\item{sample_ids, sample_ids_x, sample_ids_y}{Optional legacy/global sample
filter or low-level sample order for fitting. New prepared workflows should
prefer \code{exposure_samples} and \code{outcome_samples}.}
\item{metadata}{Optional user metadata stored with the prepared input.}
\item{input}{Prepared input object.}
\item{founder_dim}{Expected founder dimension; defaults to eight founders.}
\item{require_selected_loci}{Whether validation should require
\code{selected_loci}.}
\item{addcovar_x, addcovar_y}{Optional exposure/outcome covariates.}
\item{contrasts}{Founder contrast matrix.}
\item{marker}{Marker ID.}
\item{overlap_ids}{IDs present in both exposure and outcome samples.}
\item{O_xy}{Explicit exposure-by-outcome overlap selector matrix for low-level
experiments only. Public prepared workflows require harmonized sample IDs and
derive overlap internally.}
\item{standardize_phenotypes}{Whether to standardize exposure and outcome before
fitting.}
\item{reml}{Whether internal heritability estimation uses REML.}
\item{fc_cols}{Founder contrast columns.}
\item{ridge}{Numerical ridge used in PSD and inverse operations.}
\item{cross_cov_model}{Cross-covariance model. Public analyses should use
\code{"star"}; \code{"raw"} is diagnostic only.}
\item{enforce_psd, stop_on_invalid_psd}{Joint PSD enforcement controls.}
\item{tol}{QR/rank tolerance.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{as_qtl2mr_prepared_input()} is the recommended explicit input scaffold
for qtl2mr. One prepared object represents one MR analysis: a shared
genetic scaffold, exposure/outcome phenotypes, optional shared covariates, one
union-sample kinship source, marker map, and selected loci. The high-level
prepared adapter always uses the generalized partial-overlap engine; full
overlap and no overlap are special cases. Users must harmonize sample IDs before
calling qtl2MR.

The validator checks sample IDs, covariates, kinship names and coverage, marker
IDs, selected loci, and the eight-founder autosomal additive assumption.
Unknown IDs and incomplete matching are not silently treated as intentional
sample exclusions. Partial exposure/outcome overlap is supported and recorded.
Chromosome X is unsupported and cannot be fitted as an autosome.

The prepared object retains the original probability arrays unchanged. Its
\code{metadata$probability_normalization} records the policy, tolerances, and
per-marker adjustment counts and maximum original sum deviations. An adjustment
emits a warning and is repeated only on an analysis-local copy; the prepared
result retains this disclosure in \code{prepared_metadata}. The bound permits
ordinary four-to-six-decimal exported probabilities, not arbitrary rescaling.
As before, only selected markers are value-checked when a selection is supplied;
the full array dimensions and marker/founder identifiers are always checked.

High-level analysis standardizes phenotypes by default: the effect is SD of
outcome per SD of exposure. Scaling and sample-accounting tables are retained
with the result. Per-locus differences in fitting samples can produce different
SDs, so there may not be one valid pooled raw-unit conversion.

\code{fit_mr_locus_block_overlap()} is a low-level helper for partial or nested
sample-overlap experiments. It fits exposure and outcome founder LMMs on their
own sample sets and constructs star-space \code{Sigma_xy} through the explicit
\code{O_xy} bilinear route.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
d <- qtl2mr_example
input <- as_qtl2mr_prepared_input(
  d$allele_probs, d$pheno, covar = d$covar, kinship = d$kinship,
  marker_map = d$marker_map, selected_loci = d$selected_loci,
  exposure = "x", outcome = "y")
validate_qtl2mr_input(input, require_selected_loci = TRUE)
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
\code{as_qtl2mr_prepared_input()} and \code{validate_qtl2mr_input()} expect
analysis-ready components; \code{vignette("from-qtl2", package = "qtl2mr")}
shows how to produce them from a qtl2 cross. The current scope is eight-founder
populations, autosomal additive loci, continuous traits, additive covariates,
selected markers and one kinship source.
\code{fit_mr_locus_block_overlap()} is part of the advanced low-level
interface.
}
\value{
\if{latex}{\out{\sloppy}}
\code{as_qtl2mr_prepared_input()} returns an object of class
\code{qtl2mr_prepared_input}. \code{validate_qtl2mr_input()} returns
\code{TRUE} invisibly. \code{fit_mr_locus_block_overlap()} returns a locus block
compatible with \code{fit_multilocus_founder_mr()}.
}
\keyword{methods}
