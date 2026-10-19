\name{run_qtl2mr_batch}
\alias{run_qtl2mr_batch}
\alias{qtl2mr_prepare_resources}
\alias{qtl2mr_select_instruments}
\alias{qtl2mr_pair_input}
\title{Select instruments and run qtl2mr for many exposure-outcome pairs}
\description{
\if{latex}{\out{\sloppy}}
Prepare shared inputs once, select one instrument set per exposure, and fit
many exposure-outcome pairs with compact tabular output. Every batch row is
the same fit that \code{run_qtl2mr_analysis()} gives for that pair;
\code{qtl2mr_pair_input()} builds the single-pair input for a full report.
}
\usage{
qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar,
  marker_map = NULL, metadata = list(),
  probability_policy = c("normalize_rounding", "strict"))

qtl2mr_select_instruments(resources, qtl_scan, exposures = NULL,
  lod_thresh = 4, exclude_chr = "X", max_loci = NULL, F_X_min = NULL,
  n_workers = 1L, progress = TRUE)

qtl2mr_pair_input(resources, instruments, exposure, outcome)

run_qtl2mr_batch(resources, pairs, instruments,
  on_nonestimable_locus = c("fail", "drop"),
  reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
  heterogeneity_method = c("chisq", "simulation", "auto", "none"),
  heterogeneity_nsim = 1000L, heterogeneity_seed = 1L,
  diagnostic_bootstrap_B = 0L, diagnostic_bootstrap_seed = 1L,
  standardize_phenotypes = TRUE, cross_cov_model = c("star", "raw", "zero"),
  min_samples = 50L, min_loci = 1L, n_workers = 1L, fail_fast = FALSE,
  resume = FALSE, overwrite = FALSE, output_dir = NULL, progress = TRUE)
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
\item{marker_map}{Marker map with \code{marker}, \code{chr} and optional
\code{pos} columns, or a chromosome-named qtl2 map list. If omitted it is built
from allele-probability dimnames.}
\item{metadata}{Optional list stored with the resources.}
\item{probability_policy}{\code{"normalize_rounding"} (default) accepts small
founder-probability sum errors from rounding and normalizes an analysis-local
copy, with a warning and an adjustment record; \code{"strict"} rejects them.}
\item{resources}{Object returned by \code{qtl2mr_prepare_resources()}.}
\item{qtl_scan}{LOD scores: a marker by trait matrix such as
\code{qtl2::scan1()} output, a trait-named list of single-column matrices, or
a marker-named LOD vector for one trait.}
\item{exposures}{Traits to select instruments for; \code{NULL} uses every trait
in \code{qtl_scan}. Each must be a numeric column of \code{pheno}.}
\item{lod_thresh}{LOD threshold for peak picking: within each chromosome, the
maximum of each contiguous region above the threshold, in map order (marker
positions, or map rows when positions are missing), is a candidate. It is a
selection rule, not a significance calibration.}
\item{exclude_chr}{Chromosomes left out of the candidate table (default
\code{"X"}; only autosomal loci can be fitted).}
\item{max_loci}{Optional cap on the number of selected loci per exposure,
keeping the highest LOD peaks among those that pass the other rules.}
\item{F_X_min}{Optional minimum block F for the exposure. \code{NULL} (default)
does not filter on F.}
\item{n_workers}{Number of processes. On Unix, exposures (selection) or groups
of pairs sharing an exposure (batch) run in parallel with
\code{parallel::mclapply()}; elsewhere they run sequentially. Results do not
depend on \code{n_workers}.}
\item{progress}{Print one progress line per exposure (selection) or per pair
(batch; per exposure group when running in parallel).}
\item{instruments}{Instrument table, or the path of a CSV file holding one,
with one row per exposure and marker: usually
\code{qtl2mr_select_instruments()} output. Only \code{exposure} and
\code{marker} are required. When a logical \code{selected} column is present
only rows with \code{selected = TRUE} are used. \code{chr} must agree with the
marker map when given; other columns are ignored.}
\item{exposure, outcome}{Trait names of one pair.}
\item{pairs}{Data frame or CSV path with \code{exposure} and \code{outcome}
columns and an optional unique \code{pair_id}. Each pair uses the instruments
of its exposure.}
\item{on_nonestimable_locus, reported_estimator_policy, heterogeneity_method,
heterogeneity_nsim, heterogeneity_seed, diagnostic_bootstrap_B,
diagnostic_bootstrap_seed, standardize_phenotypes, cross_cov_model}{Fit
settings, with the same meaning as in \code{\link{run_qtl2mr_analysis}}. The
batch defaults differ only in \code{heterogeneity_method = "chisq"} and
\code{diagnostic_bootstrap_B = 0}.}
\item{min_samples}{Minimum number of complete samples required for both the
exposure and the outcome of a pair.}
\item{min_loci}{Minimum number of instruments required for a pair, also after
non-estimable loci are dropped.}
\item{fail_fast}{Stop at the first failing pair instead of recording an error
row and continuing. A pair that fails the checks made before fitting (missing
traits or instruments, too few loci or samples) stops the run before anything
is written; when a pair fails during fitting, completed pairs and the failing
pair are written to \code{output_dir} before the error is raised. With \code{n_workers > 1}, each
worker stops at its first failure, but exposure groups already running in other
workers finish first.}
\item{resume}{Continue an earlier run in \code{output_dir}; see Details.}
\item{overwrite}{Replace existing batch outputs in \code{output_dir}. Only files
written by \code{run_qtl2mr_batch()} are removed.}
\item{output_dir}{Optional output directory for CSV tables and resume records.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_select_instruments()} works on the exposure only. For every LOD
candidate it fits the exposure founder LMM with the same samples, covariates,
LOCO kinship and null-model hsq as the MR analysis and records the block F,
\code{F_X_block}: the joint Wald statistic of the seven founder contrasts
divided by 7 (see the README for its scale). It equals the block F of every pair
using that exposure; when the pairs are fitted with
\code{standardize_phenotypes = FALSE} the two agree up to rounding. A scan
without trait names (a plain LOD vector or an unnamed column) can only be used
for one exposure. A candidate is not selected, in
this order, when it is not autosomal (\code{drop_reason = "non_autosomal"}),
when the exposure founder effects cannot be estimated (\code{"nonestimable_x"}
or \code{"x_fit_error"}), when \code{F_X_block < F_X_min} (\code{"F"}), or when
\code{max_loci} better candidates by LOD were already kept
(\code{"max_loci"}).

\code{run_qtl2mr_batch()} fits each pair exactly as
\code{run_qtl2mr_analysis(qtl2mr_pair_input(resources, instruments, exposure,
outcome))} with the same settings: the exposure and the outcome each use their
own complete samples and their overlap enters the star-space covariance. The
columns of \code{main_results} other than \code{pair_id}, \code{status},
\code{error_message}, \code{locus_fit_status} and \code{runtime_sec} are
identical to the same columns of the single-pair \code{main_result}. Locus fits
that share a trait, chromosome, marker and sample set are computed once and
reused across pairs; with \code{n_workers > 1} this reuse happens within each
worker only. Pairs whose exposure has no selected instrument are recorded as
\code{no_instruments} errors.

When \code{output_dir} is given, every completed pair is saved as one record in
\code{batch_pair_records/}. With \code{resume = TRUE} the run continues from
those records: a run stopped at any point, including during a parallel run,
continues where it left off; the directory may be moved; and
\code{n_workers} may change. A resume is refused, and nothing in
\code{output_dir} is changed, when the inputs, pairs, instruments of the
exposures in \code{pairs} (including their order), options or package code
differ from the original run, or when a record or a published table has been modified. An output folder
whose package-owned files or records are symbolic links is refused. Pairs that
ended in an error are fitted again on resume.

Batch results are for screening many pairs: no figures, HTML, leave-one-out or
instrument-set sensitivity. Rerun a pair of interest with
\code{qtl2mr_pair_input()} and \code{run_qtl2mr_analysis()} for the full
report, including the Steiger direction check.
}
\value{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_prepare_resources()} returns a \code{qtl2mr_resources} object.

\code{qtl2mr_select_instruments()} returns a data frame with one row per LOD
candidate: \code{exposure}, \code{marker}, \code{chr}, \code{pos}, \code{lod},
\code{F_X_block}, \code{n_x}, \code{hsq_x}, \code{selected},
\code{drop_reason}, \code{locus_rank} (rank among selected loci by LOD),
\code{selected_by} and the settings \code{lod_thresh}, \code{F_X_min},
\code{max_loci} and \code{exclude_chr}.

\code{qtl2mr_pair_input()} returns a \code{qtl2mr_prepared_input} for
\code{\link{run_qtl2mr_analysis}}.

\code{run_qtl2mr_batch()} returns a \code{qtl2mr_batch_result} list of data
frames, written as \code{batch_<name>.csv} when \code{output_dir} is given:
\code{main_results} (one row per pair), \code{method_comparison},
\code{locus_results} (one row per pair and instrument: \code{kept},
\code{drop_reason}, exposure and outcome block F, partial R2,
\code{steiger_pass}, sample counts and null-model hsq), \code{overlap_qc},
\code{warnings}, \code{errors}, \code{sample_exclusions} (one row per trait
and exclusion reason), \code{probability_normalization},
\code{input_alignment} and \code{run_manifest}. Every table has its documented
columns even when it has no rows. The run configuration is written to
\code{batch_config.json}.
}
\seealso{
\code{\link{run_qtl2mr_analysis}} for one pair with a full report.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
scan <- matrix(ex$marker_map$exposure_lod, ncol = 1,
  dimnames = list(ex$marker_map$marker, ex$exposure))
resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno,
  ex$covar, marker_map = ex$marker_map)
instruments <- qtl2mr_select_instruments(resources, scan, lod_thresh = 3,
  progress = FALSE)
instruments[, c("exposure", "marker", "lod", "F_X_block", "selected")]
batch <- run_qtl2mr_batch(resources, data.frame(exposure = "x", outcome = "y"),
  instruments, min_samples = 20L, progress = FALSE)
batch$main_results[, c("pair_id", "status", "theta", "profile_ci_lower",
  "profile_ci_upper", "n_loci_steiger_fail")]
}
\keyword{methods}
