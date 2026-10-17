\name{run_qtl2mr_batch}
\alias{run_qtl2mr_batch}
\alias{qtl2mr_prepare_resources}
\alias{qtl2mr_select_loci_from_scan}
\alias{qtl2mr_select_loci_for_traits}
\title{Run qtl2mr for many exposure-outcome pairs}
\description{
\if{latex}{\out{\sloppy}}
Prepare shared inputs once, select instruments for each exposure from a qtl2
scan, and fit many exposure-outcome pairs with compact tabular output. This is
the standard multi-pair workflow; \code{run_qtl2mr_analysis()} is the
single-pair workflow with a full report.
}
\usage{
qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar,
  qtl_scan = NULL, marker_map = NULL, metadata = list(),
  probability_policy = c("normalize_rounding", "strict"))

qtl2mr_select_loci_from_scan(qtl_scan, marker_map, instrument_trait,
  lod_thresh = 4, exclude_chr = "X", max_loci = NULL)

qtl2mr_select_loci_for_traits(qtl_scan, marker_map, instrument_traits,
  lod_thresh = 4, exclude_chr = "X", max_loci = NULL)

run_qtl2mr_batch(resources, pairs, instrument_rule = list(),
  sample_policy = c("pair_common_complete", "batch_common_complete"),
  on_nonestimable_locus = c("fail", "drop"), min_sample_retention = 0.85,
  fit_options = list(), n_workers = 1L, fail_fast = FALSE, resume = FALSE,
  force = FALSE, output_dir = NULL, progress = TRUE, ...)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{allele_probs}{Founder allele probabilities: a chromosome-named list of
sample by founder by marker arrays, for example
\code{qtl2::genoprob_to_alleleprob()} output.}
\item{kinship}{One kinship matrix, a chromosome-named LOCO kinship list such as
\code{qtl2::calc_kinship(probs, type = "loco")}, or \code{NULL} for ordinary
least squares.}
\item{pheno}{Numeric phenotype data frame with sample IDs as row names.}
\item{covar}{Numeric covariate data frame or matrix with sample IDs as row
names, or \code{NULL}.}
\item{qtl_scan}{LOD scores used for instrument selection: a marker by trait
matrix such as \code{qtl2::scan1()} output, or a trait-named list of
single-column matrices.}
\item{marker_map}{Marker map with \code{marker}, \code{chr} and optional
\code{pos} columns, or a chromosome-named qtl2 map list. If omitted it is built
from allele-probability dimnames.}
\item{metadata}{Optional list stored with the resources.}
\item{probability_policy}{\code{"normalize_rounding"} (default) accepts small
founder-probability sum errors from rounding and normalizes an analysis-local
copy, with a warning and an adjustment record; \code{"strict"} rejects them.}
\item{instrument_trait, instrument_traits}{Trait name(s) whose LOD column is
used to select instruments; names must match exactly.}
\item{lod_thresh}{LOD threshold for contiguous-island peak picking within each
chromosome. It is a selection rule, not a significance calibration.}
\item{exclude_chr}{Chromosomes excluded from selection (default \code{"X"};
only autosomal loci are supported).}
\item{max_loci}{Optional cap on the number of selected loci, keeping the
highest LOD peaks.}
\item{resources}{Object returned by \code{qtl2mr_prepare_resources()}.}
\item{pairs}{Data frame or CSV path with \code{exposure} and \code{outcome}
columns. Optional per-pair columns: \code{pair_id}, \code{instrument_trait},
\code{lod_thresh}, \code{exclude_chr}, \code{max_loci} and
\code{selected_loci_file} (a CSV with \code{marker} and \code{chr} columns that
replaces scan-based selection for that pair).}
\item{instrument_rule}{List of default selection options: \code{lod_thresh},
\code{exclude_chr} and \code{max_loci}.}
\item{sample_policy}{Which samples each pair uses. \code{"pair_common_complete"}
(default) fits each pair on the samples complete for that exposure, that
outcome, the covariates and the selected chromosomes.
\code{"batch_common_complete"} fits every pair on one sample set complete for
all traits in the batch; this can change estimates when traits have different
missingness and is intended for screening.}
\item{on_nonestimable_locus}{What to do when a selected locus has a
non-estimable founder contrast in the fitted samples. \code{"fail"} (default)
records the pair as an error; \code{"drop"} removes the locus, refits, and
records the dropped locus and reason in \code{selected_loci} and the warnings.}
\item{min_sample_retention}{For \code{"batch_common_complete"} only: the minimum
fraction of eligible samples that the shared sample set must keep. The run
stops before writing anything if retention is lower; set \code{0} to accept any
loss.}
\item{fit_options}{List of fit settings. Defaults:
\code{standardize_phenotypes = TRUE}, \code{cross_cov_model = "star"},
\code{heterogeneity_method = "none"}, \code{q_sim_n = 0},
\code{diagnostic_bootstrap_B = 0},
\code{reported_estimator_policy = "diagnostic_only"}, \code{min_samples = 50},
\code{min_loci = 1}, \code{cache_contexts = TRUE} and
\code{cache_lmm_fits = TRUE}. Unknown names are an error.}
\item{n_workers}{Number of processes. Pairs are grouped by exposure; on Unix
groups run in parallel with \code{parallel::mclapply()}, elsewhere they run
sequentially. Results do not depend on \code{n_workers}.}
\item{fail_fast}{Stop at the first failing pair instead of recording an error
row and continuing.}
\item{resume}{Continue an earlier run in \code{output_dir}; see Details.}
\item{force}{Replace existing batch outputs in \code{output_dir}. Only files
written by \code{run_qtl2mr_batch()} are removed.}
\item{output_dir}{Optional output directory for CSV tables and resume records.}
\item{progress}{Print one progress line per pair (or per exposure group when
running in parallel).}
\item{...}{Deprecated arguments \code{sample_rule} and \code{output_mode}, which
are accepted with a warning until qtl2mr 0.6.0. Any other name is an error.}
}
\details{
\if{latex}{\out{\sloppy}}
Each pair is fitted exactly as a single prepared analysis restricted to the
pair's samples: per-locus founder LMMs with the kinship for that locus's
chromosome, overlap-aware star-space covariance, and the second-order
multi-locus estimator. Fits that share a trait, chromosome, marker, sample set
and covariates are computed once and reused across pairs.

With \code{"pair_common_complete"}, exposure and outcome use the same samples.
A single \code{run_qtl2mr_analysis()} instead lets the exposure use all samples
with an exposure value even when their outcome is missing. Restricting the
single analysis to the pair's samples reproduces the batch estimate exactly.

When \code{output_dir} is given, every completed pair is saved as one record in
\code{batch_pair_records/}. With \code{resume = TRUE} the run continues from
those records: a run stopped at any point, including during a parallel run,
continues where it left off; the directory may be moved; and
\code{n_workers} may change. A resume is refused, and nothing in
\code{output_dir} is changed, when the inputs, pairs, options or package code
differ from the original run, or when a record or a published table has been
modified. Pairs that ended in an error are fitted again on resume.

Batch defaults favour throughput: no Q simulation, no bootstrap intervals for
weighted median/mode, no figures or HTML. Use \code{run_qtl2mr_analysis()} for
a full report on an individual pair.
}
\value{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_prepare_resources()} returns a \code{qtl2mr_resources} object.
\code{qtl2mr_select_loci_from_scan()} and \code{qtl2mr_select_loci_for_traits()}
return data frames with \code{instrument_trait}, \code{marker}, \code{chr},
\code{pos}, \code{exposure_lod}, \code{locus_rank} and the selection settings.

\code{run_qtl2mr_batch()} returns a \code{qtl2mr_batch_result} list of data
frames, written as \code{batch_<name>.csv} when \code{output_dir} is given:
\code{main_results} (one row per pair: status, error message, sample counts,
the reported and raw second-order estimates with profile confidence limits,
\code{primary_reliable}, \code{confidence_set_type}, \code{effect_units}, Q
statistics, comparator estimates and \code{sample_policy}),
\code{method_comparison}, \code{selected_loci} (with per-locus retention
status), \code{overlap_qc}, \code{first_stage_summary} (per-locus block F
statistics and null-model hsq), \code{warnings}, \code{errors},
\code{sample_exclusions}, \code{sample_qc}, \code{probability_normalization}
and \code{run_manifest}. The run configuration is written to
\code{batch_config.json}.
}
\seealso{
\code{\link{run_qtl2mr_analysis}} for one pair with a full report;
\code{\link{run_qtl2mr_batch_fast}} for the deprecated entry points.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
scan <- matrix(ex$marker_map$exposure_lod, ncol = 1,
  dimnames = list(ex$marker_map$marker, ex$exposure))
resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno,
  ex$covar, qtl_scan = scan, marker_map = ex$marker_map)
qtl2mr_select_loci_from_scan(scan, ex$marker_map, "x", lod_thresh = 3)
batch <- run_qtl2mr_batch(resources, data.frame(exposure = "x", outcome = "y"),
  instrument_rule = list(lod_thresh = 3), fit_options = list(min_samples = 20L),
  progress = FALSE)
batch$main_results[, c("pair_id", "status", "theta", "profile_ci_lower",
  "profile_ci_upper", "sample_policy")]
}
\keyword{methods}
