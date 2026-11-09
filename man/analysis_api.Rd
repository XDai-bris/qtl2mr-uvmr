\name{run_qtl2mr_analysis}
\alias{run_qtl2mr_analysis}
\alias{report_qtl2mr_locus}
\alias{validate_qtl2_scan1coef}
\alias{qtl2mr_primary_result}
\alias{qtl2mr_warnings}
\alias{qtl2mr_output_paths}
\alias{qtl2mr_method_comparison}
\alias{export_qtl2mr_lite_result}
\alias{summary.qtl2mr_analysis_result}
\alias{print.qtl2mr_analysis_result}
\title{Run and report one exposure-outcome analysis}
\description{
\if{latex}{\out{\sloppy}}
Fit the founder-aware MR analysis for one exposure-outcome pair, write a
result folder with tables, figures and an HTML report, and read the main
parts of the result.
}
\usage{
run_qtl2mr_analysis(input, exposure = NULL, outcome = NULL,
  output_dir = NULL, overwrite = FALSE, make_html = TRUE,
  make_figures = TRUE, progress = TRUE, n_workers = 1L, subsets = NULL)

\method{summary}{qtl2mr_analysis_result}(object, ...)

\method{print}{qtl2mr_analysis_result}(x, ...)

qtl2mr_primary_result(result)

qtl2mr_warnings(result)

qtl2mr_output_paths(result, detail = FALSE)

qtl2mr_method_comparison(result)

report_qtl2mr_locus(result, locus_id, output_dir = NULL,
  make_html = TRUE, make_figures = TRUE)

export_qtl2mr_lite_result(result_dir, dest_zip = NULL,
  max_file_mb = 10, include_rds = FALSE)

validate_qtl2_scan1coef(genoprobs, pheno, kinship = NULL,
  addcovar = NULL, contrasts = NULL, hsq = NULL, marker = NULL,
  sample_ids = NULL, standardize_pheno = FALSE, reml = TRUE,
  tolerance = 1e-6)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{object, x}{A \code{qtl2mr_analysis_result}.}
\item{...}{Ignored.}
\item{input}{An object created by \code{\link{qtl2mr_pair_input}()}.}
\item{exposure, outcome}{Trait names; default to the names stored in
\code{input}. They must differ.}
\item{output_dir}{Destination directory. If omitted, a new directory is created
under \code{tempdir()}; set option \code{qtl2mr.output_root} or environment
variable \code{QTL2MR_OUTPUT_ROOT} to use another root.}
\item{overwrite}{Allow replacing an existing qtl2mr result directory. Only a
directory previously written by qtl2mr can be replaced; other non-empty
directories are refused.}
\item{make_html, make_figures}{Write \code{report.html} and the figure files.}
\item{progress}{\code{TRUE} (default) prints a message as each stage starts and
ends, and the result summary at the end; \code{FALSE} prints neither.
\code{logs/progress.log} is written either way.}
\item{n_workers}{Number of processes for the leave-one-locus-out and
instrument-set sensitivity refits. On Unix they run in parallel with
\code{parallel::mclapply()}; elsewhere sequentially. Results and warnings do
not depend on \code{n_workers}. Under \code{options(warn = 2)} or higher the
refits run in the calling process. A worker that is killed (for example out of
memory) stops the analysis.}
\item{subsets}{Optional instrument subsets, specified as in
\code{\link{run_qtl2mr_batch}}; \code{NULL} (default) is off. Each subset
keeps the fitted loci that pass its rules (judged on \code{F_X_block} and
\code{steiger_pass} of \code{block_first_stage.csv}) and refits the multilocus
combination with the fixed settings of the analysis; the LMMs are not
refitted. A subset needs at least one locus. A subset that cannot be fitted, or
a failure anywhere in this step, gives rows with status \code{"error"} and the
message. The main result and the existing sensitivity tables
(\code{instrument_sensitivity.csv}, \code{steiger_sensitivity.csv}) are not
changed. Use \code{subsets}, or edit the instrument table, to fit a pair on
fewer loci.}
\item{result}{A \code{qtl2mr_analysis_result}; for \code{report_qtl2mr_locus()}
also a path to its \code{result_object.rds}.}
\item{detail}{For \code{qtl2mr_output_paths()}, return the full file manifest
instead of the main files.}
\item{locus_id}{Marker ID for a single-locus report.}
\item{result_dir}{A result directory or \code{qtl2mr_analysis_result} to export.}
\item{dest_zip}{Destination zip path for the lite export.}
\item{max_file_mb}{Largest file included in the lite export.}
\item{include_rds}{Include RDS objects in the lite export.}
\item{genoprobs, pheno, kinship, addcovar, contrasts, hsq, marker, sample_ids, standardize_pheno, reml}{Inputs
for a one-marker founder LMM fitted by qtl2mr and by
\code{qtl2::scan1coef()} for comparison.}
\item{tolerance}{Largest accepted absolute difference in coefficients and
standard errors.}
}
\details{
\if{latex}{\out{\sloppy}}
At each selected locus the exposure and the outcome are regressed on the
founder probabilities (as seven founder contrasts) with a linear mixed model.
The polygenic background is the kinship for that locus's chromosome, weighted
by the null-model heritability (hsq) and then treated as fixed. hsq is always
estimated by REML with the intercept and covariates but no locus; without
kinship the founder effects are fitted by ordinary least squares.
\code{data/hsq.csv} records the values and their status. The locus effects are
combined by minimising a second-order objective that accounts for the
uncertainty of the exposure effects and for the covariance between exposure
and outcome effects that arises when the two traits share samples
(star-space \code{Sigma_xy}). The confidence interval is the profile-likelihood
set of the objective.

The exposure and the outcome each use their own complete samples; overlap is
derived from the sample IDs. \code{\link{run_qtl2mr_batch}} fits pairs the same
way, so a batch row equals this analysis of the pair with the same instruments.
To fit both traits on one common sample set, build the resources with
\code{\link{qtl2mr_combine_resources}(common_only = TRUE)} or subset the
inputs before \code{\link{qtl2mr_prepare_resources}()}.

The statistical settings are fixed, here and in \code{\link{run_qtl2mr_batch}}:
\itemize{
\item A selected locus whose founder contrasts cannot be estimated in the
exposure or outcome samples stops the analysis with an error. Remove it from
the instrument table; the founder QC of \code{\link{qtl2mr_pair_input}()}
already drops outcome loci where a founder is almost absent.
\item The reported estimate is the second-order estimate (\code{reported_estimator
= "Star-space qtl2mr"}), reported only when its profile confidence set is one
bounded interval and \code{NA} otherwise. The Two-step IVW, weighted median,
weighted mode and \code{Sigma_xy = 0} estimates are diagnostics only; no other
estimate replaces a missing one.
\item The Q heterogeneity p-value uses the chi-square reference
(\code{Q_p_method = "chisq"}).
\item The weighted median and weighted mode are point estimates, without a
bootstrap standard error or interval.
\item The exposure-outcome covariance from shared samples uses the star-space
model (\code{cross_cov_model = "star"}).
\item Phenotypes are standardised, so the effect is in SD of outcome per SD of
exposure (\code{effect_units}). \code{theta_raw_units}, \code{se_raw_units} and
\code{profile_ci_lower_raw_units}/\code{profile_ci_upper_raw_units} convert it to
original units when \code{raw_conversion_available} is \code{TRUE} (both traits
have one fitted SD).
\item The leave-one-locus-out refits and the instrument-set sensitivity,
including the Steiger direction check, always run.
}
These settings are still recorded in the result tables (for example
\code{cross_cov_model}, \code{standardize_phenotypes} and \code{Q_p_method}), so
the table format is unchanged. \code{\link{qtl2mr_diagnostics}()} adds the
simulated Q p-value and bootstrap standard errors for the weighted median and
mode afterwards, from the stored locus blocks.

\code{primary_reliable = TRUE} means the numerical fit found a single bounded
confidence interval. It does not guarantee coverage: with weak instruments
(block F_X below 2, nominal p > 0.05 for the exposure) profile intervals can be too narrow, and
cross-locus dependence or selecting instruments in the same samples can also
invalidate nominal uncertainty. When the confidence set is unbounded or has
several pieces, the reported estimate is \code{NA} and the full set is in
\code{data/confidence_set.csv}; see the warnings for what to inspect.

The Steiger direction check, run with the instrument-set sensitivity, removes
loci whose founder effects explain more outcome than exposure variance (partial
R2 of the outcome above that of the exposure, \code{steiger_pass = FALSE} in
\code{block_first_stage.csv}) and repeats the multilocus fit without refitting
the LMMs. \code{steiger_sensitivity.csv} compares it with the primary estimate,
including the shift in units of the primary confidence-interval width; an
estimate outside the primary interval adds the
\code{steiger_estimate_outside_ci} notice. Its \code{status} is \code{"ok"},
\code{"no_loci_removed"} (every locus passes; the primary fit is shown),
\code{"no_loci_left"} (every locus fails; no estimate) or \code{"failed"}. It is
a sensitivity analysis, not a filter of the primary estimate.

Progress has five stages: input preparation, locus fitting, leave-one-locus-out
refits, instrument-set refits and output writing. The refits re-solve the
multilocus fit without refitting the locus LMMs. On a laptop, the bundled
400-mouse, 40-locus example (\code{\link{qtl2mr_multilocus_example}}) takes
about half a minute with figures and the HTML report.

\code{validate_qtl2_scan1coef()} is an optional check that the qtl2mr LMM
reproduces \code{qtl2::scan1coef()} coefficients and standard errors for one
marker; it requires the qtl2 package. \code{export_qtl2mr_lite_result()} writes a
small zip without RDS objects or large tables for sharing; it does not
anonymise results.
}
\value{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()} invisibly returns a \code{qtl2mr_analysis_result}
and writes \code{report.html}, \code{summary.md}, \code{output_manifest.html},
\code{result_object.rds},
\code{figures/}, \code{logs/} and \code{data/}. The main data files are
\code{main_result.csv} (one row; see below), \code{warnings.csv} (warning ID,
severity, message, consequence and recommended action),
\code{method_comparison.csv}, \code{per_locus_results.csv},
\code{confidence_set.csv}, \code{hsq.csv} (null-model heritability per role and
chromosome with its status and source), \code{sample_alignment.csv},
\code{phenotype_scaling.csv}, \code{overlap_qc.csv}, \code{loo_mr.csv},
\code{block_first_stage.csv} (per-locus block F, partial R2 and
\code{steiger_pass}), \code{instrument_sensitivity.csv},
\code{steiger_sensitivity.csv} and \code{output_manifest.csv}, which lists
every file the analysis wrote (a later \code{qtl2mr_diagnostics()} file is not
added). With \code{subsets}, the result also holds \code{subset_results},
written as \code{subset_results.csv}: one row per subset with \code{exposure},
\code{outcome}, \code{subset}, \code{rules}, \code{n_loci}, \code{status},
\code{same_loci_as}, \code{error_message} and the columns of
\code{main_result.csv} other than \code{n_loci_selected} and
\code{n_loci_dropped_nonestimable}; see \code{\link{run_qtl2mr_batch}}. For
the same pair its rows equal the batch rows. \code{\link{qtl2mr_diagnostics}()}
later adds \code{data/diagnostics.csv} to the folder.

\code{qtl2mr_primary_result()} returns the one-row summary. Key columns:
\code{theta_reported}, \code{se_reported}, \code{reported_ci_lower},
\code{reported_ci_upper}, \code{reported_estimator} and
\code{reported_estimator_reason} (the reported result);
\code{theta}, \code{se} and \code{profile_ci_lower}/\code{upper} (the raw
second-order fit); \code{primary_reliable}, \code{confidence_set_type},
\code{confidence_set_bounded} and \code{inference_status}; \code{effect_units} and \code{theta_raw_units};
\code{n_loci}, \code{n_x_median}, \code{n_y_median} and
\code{n_overlap_median}; \code{Q}, \code{Q_p_reported} and \code{Q_p_method}
(\code{"chisq"}; the simulation columns such as \code{Q_p_sim} are \code{NA});
comparator estimates (\code{theta_two_step_ivw},
\code{theta_weighted_median}, \code{theta_weighted_mode},
\code{theta_independence}; \code{se_weighted_median} and
\code{se_weighted_mode} are \code{NA}); and \code{n_loci_steiger_fail}, the
number of loci whose exposure partial R2 is below the outcome partial R2.

\code{summary()} returns a one-row overview with \code{exposure},
\code{outcome}, \code{n_samples}, \code{n_loci}, the reported \code{estimate},
\code{se}, \code{ci_lower} and \code{ci_upper},
\code{effect_units}, \code{reported_estimator}, \code{primary_reliable},
\code{confidence_set_type} and \code{confidence_set_bounded} (which describe
the second-order profile set), \code{Q}, \code{Q_p}, \code{Q_p_method},
and \code{warnings} (IDs of warnings shown in the report).
\code{qtl2mr_warnings()}, \code{qtl2mr_method_comparison()} and
\code{qtl2mr_output_paths()} return data frames. \code{report_qtl2mr_locus()}
invisibly returns the locus report directory, \code{export_qtl2mr_lite_result()}
the zip path, and \code{validate_qtl2_scan1coef()} a one-row pass/fail summary.
}
\seealso{
\code{\link{qtl2mr_prepare_resources}} and \code{\link{qtl2mr_pair_input}} to
build the input; \code{\link{qtl2mr_select_instruments}} to choose instruments;
\code{\link{qtl2mr_combine_resources}} to combine two data sets;
\code{\link{qtl2mr_diagnostics}} for the simulated Q p-value and bootstrap
standard errors; \code{\link{run_qtl2mr_batch}} for many pairs;
\code{\link{refit_qtl2mr_subset}} for follow-up refits.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno,
  ex$covar, ex$marker_map)
input <- qtl2mr_pair_input(resources,
  data.frame(exposure = "x", marker = c("m1", "m2")), "x", "y")
result <- run_qtl2mr_analysis(input, output_dir = tempfile("qtl2mr-"),
  make_figures = FALSE, progress = FALSE)
qtl2mr_primary_result(result)[, c("theta_reported", "reported_ci_lower",
  "reported_ci_upper", "effect_units", "primary_reliable")]
qtl2mr_warnings(result)[, c("warning_id", "severity")]
}
\keyword{methods}
