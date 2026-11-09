\name{run_qtl2mr_batch}
\alias{run_qtl2mr_batch}
\alias{qtl2mr_prepare_resources}
\alias{qtl2mr_select_instruments}
\alias{qtl2mr_instrument_summary}
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
  marker_map = NULL, metadata = list())

qtl2mr_select_instruments(resources, qtl_scan, exposures = NULL,
  lod_thresh = 4, peak_rule = c("peakdrop", "spacing"),
  peakdrop = 1.5, min_spacing = 1, exclude_chr = "X",
  min_founder_copies = 5, n_workers = 1L, progress = TRUE)

qtl2mr_instrument_summary(instruments, F_breaks = c(1.5, 2, 5, 10))

qtl2mr_pair_input(resources, instruments, exposure, outcome,
  min_founder_copies = 5)

run_qtl2mr_batch(resources, pairs, instruments, min_samples = 50L,
  min_loci = 1L, min_founder_copies = 5, subsets = NULL, n_workers = 1L,
  fail_fast = FALSE, resume = FALSE, overwrite = FALSE, output_dir = NULL,
  progress = TRUE, write_locus_blocks = FALSE)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{allele_probs}{Founder allele probabilities: a chromosome-named list of
sample by founder by marker arrays, for example
\code{qtl2::genoprob_to_alleleprob()} output. Every chromosome must hold the
same samples.}
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
\item{resources}{Object returned by \code{qtl2mr_prepare_resources()}.}
\item{qtl_scan}{LOD scores: a marker by trait matrix such as
\code{qtl2::scan1()} output, a trait-named list of single-column matrices, or
a marker-named LOD vector for one trait.}
\item{exposures}{Traits to select instruments for; \code{NULL} uses every trait
in \code{qtl_scan}. Each must be a numeric column of \code{pheno}.}
\item{lod_thresh}{LOD threshold for peak picking: only markers with
LOD \eqn{\ge} \code{lod_thresh} can be peaks. It is a selection rule, not a
significance calibration.}
\item{peak_rule}{How LOD peaks become candidates; see Details.
\code{"peakdrop"} (default) is the qtl2 peak rule; \code{"spacing"} is a
sensitivity rule that keeps the maximum of each contiguous LOD island and then
applies \code{min_spacing}.}
\item{peakdrop}{For \code{peak_rule = "peakdrop"}: the LOD drop that separates
two peaks on one chromosome (default 1.5). \code{Inf} keeps one peak per
chromosome.}
\item{min_spacing}{For \code{peak_rule = "spacing"}: the minimum distance
between kept candidates on one chromosome, in the units of the marker map
positions (Mb for a qtl2 physical map; default 1).}
\item{exclude_chr}{Chromosomes left out of the candidate table (default
\code{"X"}; only autosomal loci can be fitted).}
\item{n_workers}{Number of processes. On Unix, exposures (selection) run in
parallel with \code{parallel::mclapply()}. In a batch, the trait LMM fits are
computed in parallel by trait; with several exposures, the groups of pairs
sharing an exposure then run in parallel workers, as before (warnings raised
inside a group worker are not returned to the caller); with one exposure the
pairs are fitted in order in the calling process. Elsewhere everything runs
sequentially. Results do not depend on \code{n_workers}.}
\item{progress}{Print one progress line per exposure (selection) or per pair
(batch; per exposure group when running in parallel).}
\item{instruments}{Instrument table, or the path of a CSV file holding one,
with one row per exposure and marker: usually
\code{qtl2mr_select_instruments()} output. Only \code{exposure} and
\code{marker} are required. When a \code{selected} column is present it must
be logical without missing values, and only rows with \code{selected = TRUE}
are used. \code{chr} must agree with the marker map when given; other columns
are ignored. A hand-written table such as
\code{data.frame(exposure = "x", marker = c("m1", "m2"))} is enough for
\code{qtl2mr_pair_input()}. \code{qtl2mr_instrument_summary()}
needs the data frame returned by \code{qtl2mr_select_instruments()} (at least
\code{exposure}, \code{selected}, \code{drop_reason} and \code{F_X_block}).}
\item{F_breaks}{Block F values at which \code{qtl2mr_instrument_summary()}
counts the fitted candidates reaching them (one column
\code{n_F_ge_<break>} each); distinct finite numbers.}
\item{exposure, outcome}{Trait names of one pair.}
\item{pairs}{Data frame or CSV path with \code{exposure} and \code{outcome}
columns and an optional unique \code{pair_id}. Each pair uses the instruments
of its exposure.}
\item{min_samples}{Minimum number of complete samples required for both the
exposure and the outcome of a pair.}
\item{min_loci}{Minimum number of instruments required for a pair, also after
the outcome founder QC.}
\item{fail_fast}{Stop at the first failing pair instead of recording an error
row and continuing. A pair that fails the checks made before fitting (missing
traits or instruments, too few loci or samples) stops the run before anything
is written; when a pair fails during fitting, completed pairs and the failing
pair are written to \code{output_dir} before the error is raised. With
\code{n_workers > 1}, each exposure group stops at its own first failing pair,
but every exposure group is still run before the error is raised.}
\item{resume}{Continue an earlier run in \code{output_dir}; see Details.}
\item{overwrite}{Replace existing batch outputs in \code{output_dir}. Only files
written by \code{run_qtl2mr_batch()} are removed.}
\item{output_dir}{Optional output directory for CSV tables and resume records.}
\item{write_locus_blocks}{Also write the per-locus founder blocks of every
completed pair (requires \code{output_dir}); see Details. The default writes
and returns nothing extra.}
\item{min_founder_copies}{Threshold of the founder QC (default 5); \code{NULL}
turns the QC off. A locus is dropped when the rarest founder's expected allele
count, \eqn{E = \min_k \sum_i 2 p_{ik}}, is below \code{min_founder_copies},
where the sum runs over the samples the trait's LMM is fitted on and \eqn{p}
are the founder probabilities as fitted (after rounding normalisation). In
\code{qtl2mr_select_instruments()} it applies to the exposure; in
\code{run_qtl2mr_batch()} and \code{qtl2mr_pair_input()} to the outcome of each
pair, so a locus can be dropped for one pair and kept for another. See
Details.}
\item{subsets}{Optional instrument subsets; \code{NULL} (default) is off. A
named list with one element per subset, each a list of rules:
\code{F_X_min} (one number; keeps loci with a finite block F of the exposure
\eqn{\ge} the value) and \code{steiger = TRUE} (keeps loci with
\code{steiger_pass = TRUE}, so \code{FALSE} and \code{NA} are dropped;
\code{FALSE} means no Steiger rule). The rules of a subset are combined by
intersection. Names must be unique, use only letters, digits, \code{.} and
\code{_}, and not be \code{"all"}; an unknown rule or an empty rule list is an
error. For example \code{list(F1.5 = list(F_X_min = 1.5), F1.5_SF =
list(F_X_min = 1.5, steiger = TRUE))}. See Details.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_select_instruments()} works on the exposure only. LOD candidates
are found separately on each chromosome, in map order (marker positions, then
marker map rows; map rows alone when positions are missing), and never across
a chromosome boundary. Markers with a missing or non-finite LOD are skipped.
With \code{peak_rule = "peakdrop"}, as in \code{qtl2::find_peaks()} without
intervals, a new peak starts only after the LOD has fallen by more than
\code{peakdrop} below the current peak and risen again by more than
\code{peakdrop}; otherwise a higher LOD moves the current peak. Each peak is
represented by its maximum-LOD marker (the first in map order when tied).
Unlike qtl2, a marker exactly at \code{lod_thresh} counts as above the
threshold. With \code{peak_rule = "spacing"} the candidates start from the
maximum of each contiguous run (island) of markers with LOD \eqn{\ge}
\code{lod_thresh}, the only rule of qtl2mr 0.6.1 and earlier; they are then
taken in decreasing LOD, and any within \code{min_spacing} of one already kept
on the same chromosome is dropped, so the kept peaks are more than
\code{min_spacing} apart. No rule is LD pruning: nearby peaks can still carry
correlated founder haplotypes.

For every LOD candidate that passes the founder QC, it fits the exposure
founder LMM with the same samples, covariates, LOCO kinship and null-model hsq
as the MR analysis and records the block F, \code{F_X_block}: the joint Wald
statistic of the seven founder contrasts divided by 7 (see the README for its
scale). It equals the block F of every pair using that exposure. A scan without
trait names (a plain LOD vector or an unnamed column) can only be used for one
exposure. A candidate is not selected, in this order, when it is not autosomal
(\code{drop_reason = "non_autosomal"}), when it fails the founder QC
(\code{"founder_qc"}; such a candidate is not fitted, so it has no block F and
cannot be \code{"nonestimable_x"}), or when the exposure founder effects cannot
be estimated (\code{"nonestimable_x"} or \code{"x_fit_error"}). Every other
candidate is selected: there is no F filter and no cap on the number of loci.
To use only the stronger loci, pass \code{subsets} (for example
\code{list(F2 = list(F_X_min = 2))}) to \code{run_qtl2mr_batch()} or
\code{\link{run_qtl2mr_analysis}()}, or edit the table.

\strong{Founder QC.} A locus where one founder is almost absent from the fitted
samples gives a near-singular founder design: its founder contrasts are poorly
estimated or not estimable. The QC removes such loci before they are fitted.
For each locus it computes \eqn{E}, the expected number of copies of the
rarest founder allele, as defined under \code{min_founder_copies}, and drops
the locus when \eqn{E} is below the threshold (default 5).
\itemize{
\item Exposure side: \code{qtl2mr_select_instruments()} checks each candidate on
the exposure's samples after \code{"non_autosomal"} and before fitting
(\code{drop_reason = "founder_qc"}; the table gains \code{E_min_x}).
\item Outcome side: \code{qtl2mr_pair_input()} and \code{run_qtl2mr_batch()}
check each instrument of a pair on the outcome's samples before fitting. A
failing locus is not fitted for that pair, is listed in \code{locus_results}
with \code{kept = FALSE} and \code{drop_reason = "founder_qc_outcome"}, and is
counted in \code{main_results$n_loci_dropped_founder_qc};
\code{n_loci_selected} counts the loci that pass. A pair left with fewer than
\code{min_loci} loci ends as a \code{too_few_loci} error
(\code{qtl2mr_pair_input()} stops when no locus is left). The QC result of a
single-pair input is recorded in its \code{metadata$founder_qc_outcome}, and a
batch row still equals the single analysis.
\item The order of the drops is \code{non_autosomal}, then the founder QC, then
non-estimable loci. In selection a non-estimable locus is dropped; in a pair, a
locus that passes the QC but is not estimable stops the pair with an error.
\item \code{qtl2mr_pair_input()} and \code{run_qtl2mr_batch()} do not repeat the
exposure-side check, so an instrument table written by hand, or selected with
\code{min_founder_copies = NULL}, is used as given for the exposure.
}
Typical effect: in a Diversity Outbred cohort of about 800 mice, with about
3,000 candidate loci across 80 traits, a threshold of 5 removed under 1\% of
the exposure loci, mostly loci whose design was near-singular because a founder
was almost absent, and removed an outcome locus in about 2\% of the pairs.
\code{min_founder_copies = NULL} turns the QC off, for example to reproduce
results made without it.

\code{run_qtl2mr_batch()} fits each pair exactly as
\code{run_qtl2mr_analysis(qtl2mr_pair_input(resources, instruments, exposure,
outcome))}: the exposure and the outcome each use their own complete samples,
their overlap enters the star-space covariance, and the fit settings are the
fixed settings described in \code{\link{run_qtl2mr_analysis}} (non-estimable
loci fail the pair; chi-square Q; weighted median and mode without bootstrap
standard errors; star-space covariance; standardised phenotypes). They are
still recorded in \code{batch_config.json} and the result tables. The batch
does not run the leave-one-locus-out or instrument-set sensitivity refits. The
columns of \code{main_results} other than \code{pair_id}, \code{status},
\code{error_message}, \code{locus_fit_status} and \code{runtime_sec} are
identical to the same columns of the single-pair \code{main_result}. To fit
both traits of every pair on one common sample set, build the resources with
\code{\link{qtl2mr_combine_resources}(common_only = TRUE)} or subset the inputs
first.
Locus fits that share a trait, chromosome, marker, sample set and phenotype
values are computed once and reused across pairs: the batch first fits, trait by
trait, every locus its pairs need, then fits the pairs in order. Fits are kept
for the R session (also those of \code{qtl2mr_select_instruments()} and of
single analyses) within a memory budget, by default 4 GiB as measured by
\code{object.size()} (option \code{qtl2mr.cache_budget}; see Value). A batch
works in waves that each hold up to half of that budget. The least recently used
fits are released first, a batch keeps the fits it is using, and fits that no
remaining pair needs become the first to be released. Reuse never changes a
result. Pairs whose exposure has no selected instrument are recorded as
\code{no_instruments} errors.

With \code{write_locus_blocks = TRUE}, \code{batch_locus_blocks/pair_<index>.rds}
holds, for each completed pair, a named list of per-locus blocks with
\code{marker}, \code{chr}, \code{pos}, \code{beta_x}, \code{beta_y} (founder
contrasts), \code{Sigma_x}, \code{Sigma_y}, \code{Sigma_xy} (as used by the
estimator), \code{F_X_block}, \code{F_Y_block}, \code{partial_R2_X_block},
\code{partial_R2_Y_block}, \code{steiger_pass}, \code{df_x}, \code{df_y},
\code{n_x}, \code{n_y}, \code{n_overlap}, \code{hsq_x}, \code{hsq_y},
\code{rho_xy}, \code{cross_cov_model}, \code{joint_psd_status} and short
\code{fit_x}/\code{fit_y} lists (\code{df}, \code{n}, \code{hsq}); the
attribute \code{qtl2mr_locus_blocks} records the pair and settings. The list can
be passed to \code{\link{fit_multilocus_founder_mr}()}:
\code{fit_multilocus_founder_mr(blocks, allow_cross_cov = TRUE, ridge = 1e-10,
heterogeneity_method = "chisq")} reproduces the pair's \code{theta}, \code{se},
\code{Q} and \code{Q_p_chisq}, and a subset of the list refits the pair on those
loci.
\code{batch_locus_blocks_manifest.csv} lists the files in pair order. The files
belong to the batch output: they are removed by \code{overwrite = TRUE}, and on
resume they are written again from the pair records.

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

With \code{subsets}, the instruments a pair is fitted on are its maximal set
and the main row is its fit, unchanged. For each subset, the loci of the pair
that were fitted (after the outcome founder QC) and pass the subset's rules are
judged on that pair's own \code{F_X_block} and \code{steiger_pass} (the values
in \code{locus_results}), and the multilocus combination is refitted on their
blocks with the same fixed settings as the main fit. The LMMs are not refitted, so a
subset row equals the main row of a separate run whose instruments are only
that subset's loci, provided the fitted samples of each locus do not depend on
which chromosomes carry instruments (they do not when every sample has
probabilities and kinship on every chromosome). A subset whose loci are those
of the maximal set or of an earlier subset copies that row and records it in
\code{same_loci_as}. A subset with fewer than \code{min_loci} loci gets status
\code{"too_few_loci"} and no estimate; a subset whose refit fails gets status
\code{"error"} and an \code{errors} row with stage \code{"subset"}, and does
not stop the run under \code{fail_fast} or change the pair's other rows. A
failure elsewhere in a pair's subset step (selecting the loci or assembling the
rows) gives every subset of that pair status \code{"error"} in the same way,
and the main row and locus rows stay as without subsets. Warnings raised by the
estimator itself are raised again by each subset refit, so they can appear
several times. Pairs
without a fit get status \code{"pair_not_fitted"}. The specification is part of
the configuration and the resume fingerprint; subset rows are saved in the pair
records, so they are restored on resume. The locus blocks written by
\code{write_locus_blocks} are those of the maximal set.

Batch results are for screening many pairs: no figures, HTML, leave-one-out or
instrument-set sensitivity. Rerun a pair of interest with
\code{qtl2mr_pair_input()} and \code{run_qtl2mr_analysis()} for the full
report, including the Steiger direction check, and then
\code{\link{qtl2mr_diagnostics}()} for its simulated Q p-value and bootstrap
standard errors.

\code{qtl2mr_prepare_resources()} checks the founder probabilities once. Sums
that differ from one only by rounding (at most \code{n_founders * 5e-5}, 0.0004
for eight founders, the worst case for probabilities rounded to four decimal
places) are accepted: sample-marker vectors off by more than \code{1e-6} are
divided by their sum in an analysis-local copy, with a warning and a record in
\code{metadata$probability_normalization}; the original arrays are kept. Larger deviations are an error; values are
never clipped, and negative or non-finite entries are never repaired.
}
\value{
\if{latex}{\out{\sloppy}}
\code{qtl2mr_prepare_resources()} returns a \code{qtl2mr_resources} object.

\code{qtl2mr_select_instruments()} returns a data frame with one row per LOD
candidate: \code{exposure}, \code{marker}, \code{chr}, \code{pos}, \code{lod},
\code{F_X_block}, \code{n_x}, \code{hsq_x}, \code{selected},
\code{drop_reason}, \code{locus_rank} (rank among selected loci by LOD),
\code{selected_by} (\code{"lod_peakdrop_max"} or
\code{"lod_island_max_min_spacing"}) and
the settings \code{lod_thresh}, \code{peak_rule}, \code{peakdrop} and
\code{min_spacing} (each \code{NA} when the rule does not use it),
\code{rule_version} (\code{"tophit_v1"}), \code{F_X_min},
\code{max_loci} and \code{exclude_chr}. \code{F_X_min} and \code{max_loci} are
always \code{NA}: selection no longer has these rules, and the columns are kept
so that the table format is unchanged. With the founder QC on (the default) it
also has \code{E_min_x} (after \code{hsq_x}; \code{NA} for loci that are not
autosomal) and the setting \code{min_founder_copies}.

\code{qtl2mr_instrument_summary()} returns one row per exposure, in table
order: \code{exposure}; \code{n_candidates}; the number of candidates dropped
for each reason, \code{n_dropped_non_autosomal}, \code{n_dropped_founder_qc}
(only when the table was selected with the founder QC on),
\code{n_dropped_nonestimable_x}, \code{n_dropped_x_fit_error},
\code{n_dropped_F} and \code{n_dropped_max_loci} (0 for tables selected by
this version) and one column for any other reason present; \code{n_selected}; \code{F_min}, \code{F_q25},
\code{F_median}, \code{F_mean}, \code{F_q75} and \code{F_max} of
\code{F_X_block} over the selected loci (type-7 quantiles, as in the single
analysis's \code{first_stage_summary.csv}; \code{NA} when no locus is
selected); \code{n_F_ge_<break>}, the number of fitted candidates (finite
\code{F_X_block}, selected or not) with \code{F_X_block} at least each break;
and the settings columns of the table (\code{selected_by}, \code{lod_thresh},
\code{peak_rule}, \code{peakdrop}, \code{min_spacing}, \code{rule_version},
\code{F_X_min}, \code{max_loci}, \code{exclude_chr} and, when present,
\code{min_founder_copies}). Because selection has no F filter, the fitted
candidates are the selected loci, so \code{n_F_ge_<b>} is the number of
selected loci that a subset \code{list(F_X_min = b)} keeps.

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
\code{input_alignment} and \code{run_manifest} (one row; \code{status} is
\code{"complete"} or \code{"stopped_fail_fast"}). In \code{run_manifest},
\code{context_cache_hits} and \code{lmm_fit_cache_hits} count requests for a
key this run had asked for before, whether or not its value was still kept
(the misses count first requests, which may reuse a value kept from an earlier
call); \code{context_computed} and \code{lmm_fit_computed} count the contexts
and per-locus fits actually computed in the run, including recomputations of
released values and those made in workers; \code{cache_released_entries} and
\code{cache_released_bytes} count values released to keep within
\code{cache_budget_bytes}. The budget is set with
\code{options(qtl2mr.cache_budget = bytes)} (default 4 GiB, by
\code{object.size()}) and kinship eigendecompositions are kept up to
\code{options(qtl2mr.decomposition_budget = bytes)} (default 1 GiB, per
process, outside the cache budget); neither changes a result, but both raise
peak memory (a wave holds up to half the cache budget). A later call in the
same session reuses only what is still kept. Each table has the same
columns whether or not it has rows. With \code{min_founder_copies},
\code{main_results} also has \code{n_loci_dropped_founder_qc} and
\code{locus_results} also has \code{E_min_y} (the outcome's \eqn{E} at every
instrument of the pair); both columns are present whenever the QC is on, with
\code{NA} for pairs that stopped before it. The threshold is recorded in
\code{config}. With \code{subsets}, the list also holds \code{subset_results}
(written as \code{batch_subset_results.csv}): one row per pair of
\code{main_results} and subset, in that order, with \code{pair_id},
\code{exposure}, \code{outcome}, \code{subset}, \code{rules} (for example
\code{"F_X_min>=1.5;steiger"}), \code{n_loci}, \code{status} (\code{"ok"},
\code{"too_few_loci"}, \code{"error"} or \code{"pair_not_fitted"}),
\code{same_loci_as} (\code{"all"} or an earlier subset, otherwise \code{NA}),
\code{error_message} and the estimate and test columns of \code{main_results}
(those from \code{lmm_backend} to \code{standardize_phenotypes}, without
\code{n_loci_selected} and \code{n_loci_dropped_nonestimable}, which describe
the pair's LMM stage); \code{locus_results} ends with one logical column
\code{in_subset_<name>} per subset, \code{NA} on rows of loci that were not
fitted. The specification is recorded in \code{config$subsets}. The list also holds \code{config},
\code{resources_metadata} and \code{pairs}; the run configuration is written to
\code{batch_config.json}. With \code{write_locus_blocks = TRUE} it also holds
\code{locus_blocks_manifest}.
}
\seealso{
\code{\link{run_qtl2mr_analysis}} for one pair with a full report;
\code{\link{qtl2mr_combine_resources}} to combine two resources objects;
\code{\link{qtl2mr_diagnostics}} for the simulated Q p-value and bootstrap
standard errors of one pair.
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
qtl2mr_instrument_summary(instruments)[, c("exposure", "n_candidates",
  "n_selected", "F_median", "n_F_ge_1.5", "n_F_ge_2")]
batch <- run_qtl2mr_batch(resources, data.frame(exposure = "x", outcome = "y"),
  instruments, min_samples = 20L, progress = FALSE)
batch$main_results[, c("pair_id", "status", "theta", "profile_ci_lower",
  "profile_ci_upper", "n_loci_steiger_fail", "n_loci_dropped_founder_qc")]

# The same pair in the single-pair workflow gives the same estimate.
single <- run_qtl2mr_analysis(qtl2mr_pair_input(resources, instruments, "x", "y"),
  output_dir = tempfile("qtl2mr-"), make_html = FALSE, make_figures = FALSE,
  progress = FALSE)
identical(batch$main_results$theta, qtl2mr_primary_result(single)$theta)
}
\keyword{methods}
