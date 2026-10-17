\name{run_qtl2mr_analysis}
\alias{run_qtl2mr_analysis}
\alias{report_qtl2mr_locus}
\alias{validate_qtl2_scan1coef}
\alias{qtl2mr_primary_result}
\alias{qtl2mr_warnings}
\alias{qtl2mr_output_paths}
\alias{qtl2mr_method_comparison}
\alias{export_qtl2mr_lite_result}
\alias{plot_method_comparison}
\title{High-level qtl2MR analysis and reporting API}
\description{
\if{latex}{\out{\sloppy}}
Run the packaged qtl2MR workflow, write a standardized result directory, and
create whole-analysis or single-locus reports.
}
\usage{
run_qtl2mr_analysis(input, exposure = NULL, outcome = NULL,
  adapter = c("auto", "prepared"), output_dir = NULL, loci = NULL,
  overwrite = FALSE, make_html = TRUE, make_figures = TRUE,
  verbose = TRUE, run_qtl2_validation = FALSE,
  run_heteroskedasticity = TRUE, heteroskedasticity_B = 300L,
  diagnostic_bootstrap_B = 500L, diagnostic_bootstrap_seed = 1L,
  reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
  heterogeneity_method = c("auto", "simulation", "chisq", "none"),
  heterogeneity_nsim = 1000L, heterogeneity_seed = 1L,
  run_loo_mr = TRUE, run_instrument_sensitivity = TRUE,
  progress = NULL, ...)

report_qtl2mr_locus(result, locus_id, output_dir = NULL,
  make_html = TRUE, make_figures = TRUE)

validate_qtl2_scan1coef(genoprobs, pheno, kinship = NULL,
  addcovar = NULL, contrasts = NULL, hsq = NULL, marker = NULL,
  sample_ids = NULL, standardize_pheno = FALSE, reml = TRUE,
  tolerance = 1e-6)

qtl2mr_primary_result(result)

qtl2mr_warnings(result)

qtl2mr_output_paths(result, detail = FALSE)

qtl2mr_method_comparison(result)

export_qtl2mr_lite_result(result_dir, dest_zip = NULL,
  max_file_mb = 10, include_rds = FALSE)

plot_method_comparison(methods, main = "qtl2MR method comparison",
  xlab = "MR effect", show_value_labels = TRUE, primary_method = NULL, ...)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{input}{A \code{qtl2mr_prepared_input} object.}
\item{adapter}{High-level input adapter. \code{"auto"} accepts a
\code{qtl2mr_prepared_input}; \code{"prepared"} requests it explicitly.}
\item{exposure, outcome}{Trait names. Prepared inputs may store these names in
the prepared object or receive them here.}
\item{output_dir}{Optional destination directory. If omitted,
\code{run_qtl2mr_analysis()} creates a unique result directory under
\code{tempdir()/results}. Set option \code{qtl2mr.output_root} or environment
variable \code{QTL2MR_OUTPUT_ROOT} to opt into another root.}
\item{loci}{Optional selected loci. A single numeric value is interpreted as a
top-\eqn{n} smoke-test cap; a character vector is interpreted as marker IDs.}
\item{overwrite}{Whether a non-empty \code{output_dir} may be reused.}
\item{make_html, make_figures}{Whether to create report HTML and figure files.}
\item{verbose}{Whether to print routine run messages and the standard result
summary. If \code{progress = NULL}, it also controls whether progress is shown.}
\item{progress}{Progress display mode. \code{NULL} follows \code{verbose};
\code{TRUE} prints progress and writes \code{logs/progress.log};
\code{FALSE} disables progress output and progress logging; \code{"log"}
writes only \code{logs/progress.log}. Console progress uses messages, so
\code{suppressMessages()} hides it without disabling the log.}
\item{run_qtl2_validation}{Whether to record that optional qtl2 parity
validation should be performed separately. Defaults to \code{FALSE}. Primary
prepared analysis never invokes qtl2. When \code{TRUE}, the result records an
informational not-available status with directions to
\code{validate_qtl2_scan1coef()} and explicit generic inputs.}
\item{run_heteroskedasticity, heteroskedasticity_B}{Retained for compatibility;
adapter-specific heteroskedasticity diagnostics are not run automatically by
the prepared-input workflow.}
\item{diagnostic_bootstrap_B}{Parametric-bootstrap draws for weighted
median/mode diagnostics. Intervals are the original estimate plus or minus
1.96 times bootstrap SD, using independent normal draws of per-locus effects
and fixed loci/weights. They are approximate 95\% \code{bootstrap-SE Wald}
intervals, not profile CIs. Set to \code{0L} for point estimates only or an
integer of at least two; the default is \code{500L}.}
\item{diagnostic_bootstrap_seed}{Random seed for weighted median/mode
parametric-bootstrap confidence intervals.}
\item{reported_estimator_policy}{How to choose the user-facing reported
estimate when the star-space second-order optimizer is not well bracketed or
has a non-finite profile interval. The default, \code{"diagnostic_only"},
reports the star-space SO estimate for a numerically acceptable, bounded fit
and otherwise returns no reported estimate (estimate, SE, and CI are \code{NA}) while retaining
the raw SO and diagnostic comparator fields. \code{"strict_so"} always
reports the star-space SO estimate with reliability warnings.
\code{"auto_fallback"} is an explicit compatibility opt-in that reports the
Two-step IVW diagnostic estimate when SO is unreliable and a finite diagnostic
estimate is available.}
\item{heterogeneity_method}{How to report the aggregate Q heterogeneity
p-value. \code{"auto"} uses a chi-square reference for single-locus analyses
and conditional residual simulation for multi-locus analyses.}
\item{heterogeneity_nsim}{Number of conditional residual simulations used to
calibrate the aggregate multi-locus Q heterogeneity p-value. Set to \code{0L}
to skip calibration in timing-only runs.}
\item{heterogeneity_seed}{Random seed for heterogeneity simulation
calibration.}
\item{run_loo_mr}{Whether to run leave-one-locus-out global MR sensitivity.}
\item{run_instrument_sensitivity}{Whether to run selected-locus instrument-set
sensitivity refits.}
\item{show_value_labels}{Whether \code{plot_method_comparison()} should draw
estimate/interval text labels beside plotted points.}
\item{result}{A \code{qtl2mr_analysis_result} object, a pilot result where
documented, or a path to \code{result_object.rds} for locus reports.}
\item{detail}{For \code{qtl2mr_output_paths()}, whether to return the full
file-by-file output manifest. The default \code{FALSE} returns a compact index
of the main report, summary, result table, manifest, and RDS object.}
\item{result_dir}{A qtl2MR result directory or \code{qtl2mr_analysis_result}
object for lite export.}
\item{dest_zip}{Optional destination zip path for lite export.}
\item{max_file_mb}{Maximum file size included in the lite export.}
\item{include_rds}{Whether complete RDS objects are included in the lite
export. Defaults to \code{FALSE}.}
\item{locus_id}{Marker ID for a single-locus report.}
\item{genoprobs, pheno, kinship, addcovar, contrasts, hsq, marker,
sample_ids, standardize_pheno, reml}{Explicit generic inputs passed to the
internal founder-LMM fit and optional \code{qtl2::scan1coef()} comparator.}
\item{tolerance}{Maximum accepted absolute beta and SE difference for the
optional qtl2 parity comparison.}
\item{methods}{Method comparison table with \code{method}, \code{estimate},
\code{ci_lower}, and \code{ci_upper} columns.}
\item{primary_method}{Optional named list used by
\code{plot_method_comparison()} to annotate the reported estimate. This is
normally filled from \code{theta_reported} and related summary fields by the
high-level report writer.}
\item{main}{Plot title.}
\item{xlab}{X-axis label.}
\item{...}{For analysis, only named \code{standardize_phenotypes},
\code{cross_cov_model}, \code{ridge}, \code{engine}, \code{context_cache}, and
\code{fit_cache} arguments are forwarded. Unknown names are errors.
For plotting, additional graphical arguments.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()} is the user-facing destination workflow for the
package. It fits the selected-locus founder-aware LMMs,
computes the overlap-aware second-order Founder GLS-IVW estimate using
star-space \code{Sigma_xy}, runs available diagnostics, and writes a fixed
folder structure:
\code{data/}, \code{figures/}, \code{logs/}, \code{report.html},
\code{summary.md}, and \code{result_object.rds}.

When \code{progress} is enabled, the run reports seven coarse stages:
input/output preparation, selected-locus LMM/MR fitting, qtl2 beta/SE
validation, local founder exposure score residual-variance QC,
leave-one-locus-out MR, instrument-set sensitivity, and output writing. The
same progress lines are written to \code{logs/progress.log}. During Q simulation,
completed draw counts and elapsed time are reported at roughly ten intervals.
Each draw refits the common-effect model; the 1,000-draw default can dominate
runtime. Progress does not affect estimates or random-number draws. Choosing
\code{heterogeneity_method = "chisq"} explicitly avoids this simulation but
changes the Q p-value's reference distribution, not the primary MR fit.

Use \code{adapter = "prepared"} with an object from
\code{\link{as_qtl2mr_prepared_input}()}. The prepared adapter uses
role-specific exposure/outcome complete samples and derives overlap from
harmonized sample IDs. Primary analysis never requires qtl2; the exported
validator is a standalone optional parity check.

Profile confidence intervals are the primary uncertainty summary. The printed
SE is a curvature/Wald summary. The summary contains both the raw star-space SO
estimate (\code{theta_so}, and the backward-compatible \code{theta}) and the
policy-selected reported estimate (\code{theta_reported}). By default,
numerically acceptable bounded SO fits remain the reported estimate; unresolved,
unbounded or disconnected SO fits have \code{NA} reported estimate, SE, and CI. Two-step IVW
remains diagnostic unless \code{reported_estimator_policy = "auto_fallback"}
is explicitly requested.
The legacy reliable flag is a numerical bounded-fit check, not a guarantee of
coverage or causal validity. The full confidence-set components are written to
\code{data/confidence_set.csv}. Default effects are SD of outcome per SD of
exposure; \code{standardize_phenotypes = FALSE} uses original units.
\code{data/phenotype_scaling.csv} and \code{data/sample_alignment.csv} record
scaling and sample accounting. Weak instruments, cross-locus dependence and
same-sample selection remain inferential limitations even for a bounded fit.
The result records \code{sample_policy = "role_specific_complete"}; exact
batch analyses instead use complete cases for each exposure-outcome pair.
For numerical comparisons, supply the same sample masks and fitting options.
Aggregate heterogeneity uses
\code{heterogeneity_method = "auto"} by default: single-locus analyses report
\code{Q_p_chisq}, multi-locus analyses report simulation-calibrated
\code{Q_p_sim}, and the value selected for reporting is stored as
\code{Q_p_reported} with label \code{Q_p_method}.

\code{qtl2mr_primary_result()}, \code{qtl2mr_warnings()},
\code{qtl2mr_output_paths()}, and \code{qtl2mr_method_comparison()} are stable
accessors for user-facing parts of the result object.
The method-comparison table records the interval method/level and requested,
completed and failed bootstrap counts. Its plots and HTML reports label the
median/mode bootstrap-SE Wald intervals separately from the primary profile CI.

\code{export_qtl2mr_lite_result()} creates a compact portable zip for external
discussion by excluding large RDS objects and large sample-level tables by
default. Inclusion and exclusion decisions are recorded in
\code{SMOKE_RESULT_SIZE_AUDIT.csv}.
This is a size/portability helper, not anonymization or permission to disclose
study results. Review every included file before sharing restricted studies.

}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()}, the result accessors, lite export, locus report,
and qtl2 validation are user-facing APIs. The \code{prepared} adapter is the
portable explicit input scaffold.
}
\value{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()} returns an object of class
\code{qtl2mr_analysis_result}. \code{write_qtl2mr_results()} returns the same
object after adding manifest paths. \code{report_qtl2mr_locus()} invisibly
returns the locus report directory. Plot helpers invisibly return their plotted
data.
}
\keyword{methods}
