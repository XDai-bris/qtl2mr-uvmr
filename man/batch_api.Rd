\name{qtl2mr_batch}
\alias{qtl2mr_batch}
\alias{qtl2mr_prepare_resources}
\alias{qtl2mr_select_loci_from_scan}
\alias{qtl2mr_select_loci_for_traits}
\alias{qtl2mr_subset_allele_probs_to_loci}
\alias{qtl2mr_make_common_complete_samples}
\alias{plan_qtl2mr_batch_exact_groups}
\alias{run_qtl2mr_batch}
\alias{run_qtl2mr_batch_exact_grouped}
\alias{run_qtl2mr_batch_fast}
\alias{run_qtl2mr_batch_parallel_by_exposure}
\alias{write_qtl2mr_batch_outputs}
\alias{qtl2mr_read_mr_block_cache}
\alias{refit_qtl2mr_mr_block_cache}
\alias{refit_qtl2mr_mr_block_cache_grid}
\title{Compact batch qtl2MR workflow}
\description{
\if{latex}{\out{\sloppy}}
Prepare shared qtl2MR resources, select exposure instruments from qtl2 scan
results, run many exposure/outcome pairs with the batch LMM engine, and write
compact batch CSV outputs.
}
\usage{
qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar,
  qtl_scan = NULL, marker_map = NULL, metadata = list(),
  probability_policy = c("normalize_rounding", "strict"))

qtl2mr_select_loci_from_scan(qtl_scan, marker_map, instrument_trait,
  lod_thresh = 4, exclude_chr = "X", max_loci = NULL)

qtl2mr_select_loci_for_traits(qtl_scan, marker_map, instrument_traits,
  lod_thresh = 4, exclude_chr = "X", max_loci = NULL)

qtl2mr_subset_allele_probs_to_loci(allele_probs, selected_loci,
  marker_col = "marker", chr_col = "chr", drop_empty_chr = TRUE)

qtl2mr_make_common_complete_samples(resources, exposure, outcome,
  selected_chr = NULL, sample_rule = "common_complete")

plan_qtl2mr_batch_exact_groups(resources, pairs, instrument_rule = list(),
  min_samples = 50L, min_loci = 1L)

run_qtl2mr_batch(resources, pairs, instrument_rule = list(),
  sample_rule = list(mode = "common_complete"), fit_options = list(),
  output_mode = c("compact", "qc"), fail_fast = FALSE, resume = FALSE,
  force = FALSE, output_dir = NULL, progress = TRUE)

run_qtl2mr_batch_exact_grouped(resources, pairs, plan = NULL,
  instrument_rule = list(), fit_options = list(),
  output_mode = c("compact", "qc"), fail_fast = FALSE, n_workers = 1L,
  output_dir = NULL, progress = TRUE, force = FALSE)

run_qtl2mr_batch_fast(resources, pairs,
  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
  sample_retention_min = 0.85, allow_sample_loss = FALSE,
  fit_options = list(), output_dir = NULL, progress = TRUE,
  fail_fast = FALSE, resume = FALSE, force = FALSE)

run_qtl2mr_batch_parallel_by_exposure(resources, pairs, fast = TRUE,
  n_workers = 2L,
  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
  sample_retention_min = 0.85, allow_sample_loss = FALSE,
  fit_options = list(), output_dir = NULL, progress = TRUE,
  fail_fast = FALSE, resume = FALSE, force = FALSE)

write_qtl2mr_batch_outputs(batch, output_dir)

qtl2mr_read_mr_block_cache(cache)

refit_qtl2mr_mr_block_cache(cache, F_X_min = NA_real_,
  steiger_filter = FALSE, min_loci = 1L, allow_cross_cov = TRUE,
  ridge = 1e-10, heterogeneity_method = "none", q_sim_n = 0L,
  q_sim_seed = 1L, cache_file = NA_character_)

refit_qtl2mr_mr_block_cache_grid(cache,
  F_X_min_grid = c(NA_real_, 1, 2, 3, 5, 10),
  steiger_filter = FALSE, min_loci = 1L, cache_file = NULL, ...)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{allele_probs}{qtl2-style founder allele probabilities: a chromosome-keyed
list of sample by founder by marker arrays.}
\item{kinship}{One kinship matrix or a chromosome-keyed LOCO kinship list.}
\item{pheno}{Phenotype data frame with sample IDs as row names.}
\item{covar}{Numeric covariate data frame or matrix with sample IDs as row
names.}
\item{qtl_scan}{qtl2 scan object or list of trait-specific scan matrices.}
\item{marker_map}{Marker map containing \code{marker}, \code{chr}, and
optionally \code{pos}, or a chromosome-named qtl2 map list. If omitted it is constructed from allele-probability
dimnames.}
\item{metadata}{Optional metadata stored with the resource object.}
\item{probability_policy}{Founder-probability sum policy, propagated to each
prepared pair. \code{"normalize_rounding"} permits small rounding deviations
within the prepared-input tolerance and normalizes an analysis-local copy;
\code{"strict"} rejects sums outside the strict tolerance. Original resource
arrays are preserved. See \code{as_qtl2mr_prepared_input()}.}
\item{instrument_trait}{Trait used to select instruments from \code{qtl_scan}.
For a multi-column scan, this must exactly match one unique, non-blank column
name; qtl2mr never silently uses another trait's column.}
\item{instrument_traits}{Character vector of traits used to pre-select
instruments for a many-pair batch run.}
\item{lod_thresh}{One finite numeric value strictly greater than zero, used for
contiguous-island peak picking independently within each chromosome. It is
not a universal significance threshold or LD-pruning procedure.}
\item{exclude_chr}{Chromosomes excluded from selected loci, default
\code{"X"}.}
\item{max_loci}{Optional selected-locus cap for smoke tests.}
\item{selected_loci}{Selected-locus table containing marker and chromosome
columns. For \code{qtl2mr_subset_allele_probs_to_loci()}, this may contain
duplicates across traits; the function keeps the unique marker/chromosome
pairs.}
\item{marker_col, chr_col}{Column names identifying selected marker IDs and
chromosomes.}
\item{drop_empty_chr}{Whether to remove chromosomes with no selected markers
from the subset allele-probability object.}
\item{resources}{Object returned by \code{qtl2mr_prepare_resources()}.}
\item{exposure, outcome}{Trait names in \code{resources$pheno}.}
\item{selected_chr}{Chromosomes represented by selected loci, used to verify
LOCO kinship sample coverage.}
\item{min_samples}{Minimum complete sample count for a planned exact pair.}
\item{min_loci}{Minimum selected-locus count for each instrument request.}
\item{sample_rule}{Sample rule string or list. \code{"common_complete"} and
\code{"pair_common_complete"} use the pair-specific complete exposure/outcome
sample set. \code{"batch_common_complete"} uses one common sample set across
all traits in the batch, which can improve cross-pair context-cache reuse.
These batch policies are not the role-specific complete-case policy of an
unrestricted single prepared analysis. For an exact cross-API comparison,
use identical selected loci, options, and explicit exposure/outcome sample
filters in the prepared analysis.}
\item{sample_retention_min}{Minimum allowed ratio of batch-common complete
samples to genotype/covariate/kinship eligible samples for
\code{run_qtl2mr_batch_fast()}. Default is \code{0.85}.}
\item{allow_sample_loss}{Whether \code{run_qtl2mr_batch_fast()} may continue
when sample retention is below \code{sample_retention_min}. Default stops.}
\item{fast}{For \code{run_qtl2mr_batch_parallel_by_exposure()}, whether to use
the batch-common fast screening path. If \code{FALSE}, pair-specific complete
samples are used within each exposure group.}
\item{n_workers}{Number of exposure groups to run concurrently on Unix-like
systems for \code{run_qtl2mr_batch_parallel_by_exposure()}. For
\code{run_qtl2mr_batch_exact_grouped()}, only \code{1} is currently supported;
use external exposure-block jobs for parallel exact execution.}
\item{plan}{Optional object returned by
\code{plan_qtl2mr_batch_exact_groups()}. If omitted, the exact grouped runner
creates one before fitting.}
\item{pairs}{Data frame or CSV path with \code{exposure} and \code{outcome}
columns. Optional columns include \code{pair_id}, \code{instrument_trait},
\code{lod_thresh}, \code{exclude_chr}, \code{max_loci}, and
\code{selected_loci_file}.}
\item{instrument_rule}{List of default locus-selection options.}
\item{fit_options}{List of fit options. Defaults are compact:
\code{standardize_phenotypes = TRUE}, \code{cross_cov_model = "star"},
\code{heterogeneity_method = "none"}, \code{q_sim_n = 0},
\code{diagnostic_bootstrap_B = 0},
\code{reported_estimator_policy = "diagnostic_only"},
\code{engine = "batch"}, and
\code{cache_contexts = TRUE}, \code{cache_lmm_fits = TRUE}. Optional cache
settings include \code{save_mr_block_cache = TRUE},
\code{run_first_stage_filter_grid = TRUE}, \code{first_stage_filter_grid},
and \code{steiger_filter}. These write stripped per-locus MR block caches and
post-hoc F/Steiger-filter refits without saving sample-level residuals or
rerunning LMMs.}
\item{output_mode}{Compact output mode. \code{"qc"} currently writes the same
batch-level tables with diagnostic summaries.}
\item{fail_fast}{Whether to stop at the first failed pair or parallel exposure
group. Default records failures and continues with successful work; parallel
group failures are retained with their exposure-group key.}
\item{resume}{Opt into reuse of completed pairs only when the saved full
input/configuration/implementation fingerprint and output checksums agree.
Defaults to false. Legacy or mismatched checkpoints are rejected, not reused.
The parallel-by-exposure runner checks its parent request and all existing
worker checkpoints before any worker starts. Integrity failures are fatal
even when \code{fail_fast = FALSE}; existing outputs are not replaced by an
empty or partial merge. Older parallel outputs without a parent fingerprint
must use a new output directory or an explicit \code{force = TRUE} rerun.}
\item{force}{Explicitly authorize replacement of existing batch outputs.
Without \code{force = TRUE} or a validated resume, an occupied output directory
is refused. This does not authorize deleting unrelated files.}
\item{output_dir}{Optional directory. When supplied, outputs are updated after
each pair.}
\item{progress}{Whether to print one-line-per-pair progress.}
\item{batch}{Batch result object returned by \code{run_qtl2mr_batch()} or
\code{run_qtl2mr_batch_fast()}.}
\item{cache}{A \code{qtl2mr_mr_block_cache} object or path to a cache RDS.}
\item{F_X_min}{Minimum exposure block first-stage F statistic. \code{NA}
keeps all loci.}
\item{F_X_min_grid}{Numeric vector of \code{F_X_min} thresholds; \code{NA}
represents all loci.}
\item{steiger_filter}{Whether to keep only loci with
\code{partial_R2_X_block >= partial_R2_Y_block}. This uses cached block-level
R2 approximations and does not rerun LMMs.}
\item{allow_cross_cov}{Whether to use cached \code{Sigma_xy} in refits.}
\item{ridge}{Small numerical ridge used for covariance inversions.}
\item{heterogeneity_method, q_sim_n, q_sim_seed}{Heterogeneity p-value options
passed to \code{fit_multilocus_founder_mr()} during cached refits.}
\item{cache_file}{Optional path recorded in refit output for provenance.}
\item{...}{Additional arguments passed from
\code{refit_qtl2mr_mr_block_cache_grid()} to
\code{refit_qtl2mr_mr_block_cache()}.}
}
\details{
\if{latex}{\out{\sloppy}}
The batch engine does not fit a joint multi-marker LMM. It preserves the
qtl2MR per-locus marginal founder-block estimand, but caches the
trait/chromosome LMM transform so each selected marker only recomputes the
marker-specific founder block. Direct and batch paths share the internal
base-R kinship decomposition and profile ML/REML heritability backend; qtl2 is
optional and used only for explicitly requested parity validation.

For large many-pair runs, first select instruments for all possible exposure
traits with \code{qtl2mr_select_loci_for_traits()}, then subset the full qtl2
allele probabilities with \code{qtl2mr_subset_allele_probs_to_loci()}. The
statistical estimator is unchanged because qtl2MR fits selected-locus
marginal founder blocks; non-selected markers are not used after instrument
selection is frozen.

The compact batch runner also caches LMM contexts and per-marker LMM fits
across pairs. Cache reuse is highest when analyses share the same trait,
chromosome, marker, sample set, covariate set, heritability setting, and
kinship source. Pair-specific complete-case sample sets preserve maximum
sample size per pair; batch-wide common samples can be faster but may change
estimates if the common sample set differs.

\code{plan_qtl2mr_batch_exact_groups()} is an exact planning layer. It selects
loci once per unique instrument request, computes each pair's exact
pair-specific complete sample set, and groups pairs only when the sorted sample
IDs are identical. \code{run_qtl2mr_batch_exact_grouped()} executes this plan
without using batch-common samples. It writes \code{batch_exact_group_plan.csv},
\code{batch_sample_groups.csv}, and \code{batch_workload_summary.csv} in
addition to the compact batch tables. Exact and pair-specific workflows fail
the affected pair when a selected founder contrast is non-estimable; they do
not silently change the retained instrument set.

\code{run_qtl2mr_batch_fast()} is a screening convenience wrapper around the
batch-common sample mode. It forces compact output, enables both caches, writes
\code{batch_sample_qc.csv}, and stops by default if the retained batch-common
sample fraction is below 85 percent. As an explicit screening path it permits
audited dropping of non-estimable loci by default. The original and retained
locus sets, per-locus status and drop reason, pair-level drop count/status,
structured warning, and run-level drop policy/count are written to the compact
tables and manifest.

\code{run_qtl2mr_batch_parallel_by_exposure()} splits a compact batch by
exposure trait. This is the recommended parallelization unit for large
many-pair screening because all outcomes for the same exposure share the same
instrument selection and can reuse exposure-side LMM fits within a worker. In
fast mode, the global batch-common sample set is computed once and passed to
all exposure groups so group results remain comparable.

Batch configs, main results, and manifests identify \code{sample_policy} as
\code{"pair_common_complete"} or \code{"batch_common_complete"}. ``Exact''
describes preservation of pair-complete batching, not equivalence to a single
prepared analysis that uses all eligible exposure samples independently of
outcome missingness.

\code{batch_sample_exclusions.csv} records input-value eligibility by pair and
reason, including absent covariate IDs and NA/NaN/Inf covariate values. Its
denominator is the resource sample-ID set, before complete-case selection;
reasons may overlap and must not be summed as distinct people. The ledger
does not replace per-locus fitted/overlap counts or audit optional user sample
filters and chromosome-specific kinship/probability coverage. Covariate losses
also generate structured rows in \code{batch_warnings.csv}. The fast retention
ratio continues to use covariate-eligible samples as its denominator.
\code{batch_probability_normalization.csv} records per-pair selected-marker
rounding adjustments; the config records the probability policy. Neither
ledger changes sample selection or imputes missing values.

An ordinary parallel worker failure creates an error status for each affected
pair and retains the exposure-group root cause. Run-manifest error counts
therefore agree with the main-result statuses. A resume-integrity failure is
not an ordinary recoverable worker failure.

Default batch settings are designed for many-pair analysis and omit HTML,
figures, large RDS files, Q simulation calibration, and weighted median/mode
bootstrap intervals. Use the full reporting workflow for paper-style
single-pair reports.

When \code{save_mr_block_cache = TRUE}, the batch runner writes one stripped
\code{mr_block_cache/*.rds} file per successful pair plus
\code{batch_mr_block_cache_manifest.csv}. Each cache retains the selected-locus
founder effects and covariance matrices
(\code{beta_x}, \code{beta_y}, \code{Sigma_x}, \code{Sigma_y},
\code{Sigma_xy}), overlap metadata, block first-stage F statistics, and
block-level X/Y partial-R2 approximations. It intentionally omits large
sample-level residual and influence matrices. The cache can therefore be used
to rapidly test F-threshold and Steiger-style instrument filters with
\code{refit_qtl2mr_mr_block_cache()} or
\code{refit_qtl2mr_mr_block_cache_grid()} without rerunning the LMM layer, as
long as the phenotype, covariates, kinship, selected loci, and sample set are
unchanged.
}
\value{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_prepare_resources()} returns a \code{qtl2mr_resources} object.
\code{plan_qtl2mr_batch_exact_groups()} returns a \code{qtl2mr_exact_group_plan}
object. \code{run_qtl2mr_batch()}, \code{run_qtl2mr_batch_exact_grouped()},
\code{run_qtl2mr_batch_fast()}, and
\code{run_qtl2mr_batch_parallel_by_exposure()} return a
\code{qtl2mr_batch_result} list containing batch-level result tables.
\code{write_qtl2mr_batch_outputs()} writes compact CSV/config files and returns
the output directory invisibly. \code{refit_qtl2mr_mr_block_cache()} and
\code{refit_qtl2mr_mr_block_cache_grid()} return compact data frames with
filtered qtl2MR estimates and the number of loci retained by each filter.
}
\keyword{methods}
