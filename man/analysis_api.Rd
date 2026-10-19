\name{run_qtl2mr_analysis}
\alias{run_qtl2mr_analysis}
\alias{report_qtl2mr_locus}
\alias{validate_qtl2_scan1coef}
\alias{qtl2mr_primary_result}
\alias{qtl2mr_warnings}
\alias{qtl2mr_output_paths}
\alias{qtl2mr_method_comparison}
\alias{export_qtl2mr_lite_result}
\title{Run and report one exposure-outcome analysis}
\description{
\if{latex}{\out{\sloppy}}
Fit the founder-aware MR analysis for one exposure-outcome pair, write a
result folder with tables, figures and an HTML report, and read the main
parts of the result.
}
\usage{
run_qtl2mr_analysis(input, exposure = NULL, outcome = NULL,
  output_dir = NULL, loci = NULL, overwrite = FALSE, make_html = TRUE,
  make_figures = TRUE, verbose = TRUE, diagnostic_bootstrap_B = 500L,
  diagnostic_bootstrap_seed = 1L,
  reported_estimator_policy = c("diagnostic_only", "strict_so", "auto_fallback"),
  heterogeneity_method = c("auto", "simulation", "chisq", "none"),
  heterogeneity_nsim = 1000L, heterogeneity_seed = 1L,
  run_loo_mr = TRUE, run_instrument_sensitivity = TRUE,
  progress = NULL, ...)

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
\item{input}{An object created by \code{\link{as_qtl2mr_prepared_input}()}.}
\item{exposure, outcome}{Trait names; default to the names stored in
\code{input}. They must differ.}
\item{output_dir}{Destination directory. If omitted, a new directory is created
under \code{tempdir()}; set option \code{qtl2mr.output_root} or environment
variable \code{QTL2MR_OUTPUT_ROOT} to use another root.}
\item{loci}{Optional restriction of the selected loci: a positive integer keeps
the first \eqn{n} rows of \code{selected_loci}; a character vector keeps the
named markers.}
\item{overwrite}{Allow replacing an existing qtl2mr result directory. Only a
directory previously written by qtl2mr can be replaced; other non-empty
directories are refused.}
\item{make_html, make_figures}{Write \code{report.html} and the figure files.}
\item{verbose}{Print the result summary at the end; also controls progress
when \code{progress = NULL}.}
\item{diagnostic_bootstrap_B}{Parametric-bootstrap draws for the weighted
median and weighted mode comparators (\code{0} for point estimates only).
These comparators are diagnostics, not the reported estimate.}
\item{diagnostic_bootstrap_seed}{Seed for those draws; the global random-number
state is restored afterwards.}
\item{reported_estimator_policy}{\code{"diagnostic_only"} (default) reports the
second-order estimate only when it has a single bounded profile confidence
interval, and otherwise reports \code{NA}. \code{"strict_so"} always reports the
second-order estimate. \code{"auto_fallback"} reports the Two-step IVW
comparator when the second-order fit is not reportable.}
\item{heterogeneity_method}{Reference distribution for the Q heterogeneity
p-value. \code{"auto"} uses chi-square for one locus and conditional residual
simulation for several loci. \code{"chisq"} avoids the simulation.}
\item{heterogeneity_nsim, heterogeneity_seed}{Number of Q simulation draws and
their seed. Each draw refits the model; see Details for runtime.}
\item{run_loo_mr}{Run leave-one-locus-out refits.}
\item{run_instrument_sensitivity}{Run refits on standard locus subsets: the
strongest loci by block F, loci with block F of at least 2, one chromosome
left out at a time, and without the most influential leave-one-out loci.}
\item{progress}{\code{NULL} follows \code{verbose}; \code{TRUE} prints progress
messages and writes \code{logs/progress.log}; \code{"log"} only writes the log;
\code{FALSE} disables both.}
\item{...}{Advanced fit settings: \code{standardize_phenotypes} (default
\code{TRUE}), \code{cross_cov_model} (default \code{"star"}) and \code{ridge};
\code{engine}, \code{context_cache} and \code{fit_cache} are internal tuning
hooks. Other names are an error.}
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
by the null-model heritability (hsq) estimated by REML with the intercept and
covariates but no locus, and then treated as fixed. The locus effects are
combined by minimising a second-order objective that accounts for the
uncertainty of the exposure effects and for the covariance between exposure
and outcome effects that arises when the two traits share samples
(star-space \code{Sigma_xy}). The confidence interval is the profile-likelihood
set of the objective.

The exposure and the outcome each use their own complete samples
(\code{sample_policy = "role_specific_complete"}); overlap is derived from the
sample IDs. With default standardisation the effect is in SD of outcome per SD
of exposure (\code{effect_units}); \code{theta_raw_units} converts it to
original units when both traits have one fitted SD.

\code{primary_reliable = TRUE} means the numerical fit found a single bounded
confidence interval. It does not guarantee coverage: with weak instruments
(block F_X below 2, nominal p > 0.05 for the exposure) profile intervals can be too narrow, and
cross-locus dependence or selecting instruments in the same samples can also
invalidate nominal uncertainty. When the confidence set is unbounded or has
several pieces, the reported estimate is \code{NA} and the full set is in
\code{data/confidence_set.csv}; see the warnings for what to inspect.

Runtime is dominated by the Q heterogeneity simulation when there are several
loci: each of the 1,000 default draws refits the model. On a laptop, the
bundled 400-mouse, 40-locus example takes roughly 10 to 15 minutes with the
default and under 2 minutes with \code{heterogeneity_method = "chisq"}. Progress messages report completed
draws. Progress has five stages: input preparation, locus fitting (including
the Q simulation), leave-one-locus-out refits, instrument-set refits and output
writing.

\code{validate_qtl2_scan1coef()} is an optional check that the qtl2mr LMM
reproduces \code{qtl2::scan1coef()} coefficients and standard errors for one
marker; it requires the qtl2 package. \code{export_qtl2mr_lite_result()} writes a
small zip without RDS objects or large tables for sharing; it does not
anonymise results.
}
\value{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()} invisibly returns a \code{qtl2mr_analysis_result}
and writes \code{report.html}, \code{summary.md}, \code{result_object.rds},
\code{figures/}, \code{logs/} and \code{data/}. The main data files are
\code{main_result.csv} (one row; see below), \code{warnings.csv} (warning ID,
severity, message, consequence and recommended action),
\code{method_comparison.csv}, \code{per_locus_results.csv},
\code{confidence_set.csv}, \code{hsq.csv} (null-model heritability per role and
chromosome with its status and source), \code{sample_alignment.csv},
\code{phenotype_scaling.csv}, \code{overlap_qc.csv}, \code{loo_mr.csv},
\code{instrument_sensitivity.csv} and \code{output_manifest.csv}, which lists
every file.

\code{qtl2mr_primary_result()} returns the one-row summary. Key columns:
\code{theta_reported}, \code{se_reported}, \code{reported_ci_lower},
\code{reported_ci_upper}, \code{reported_estimator} and
\code{reported_estimator_reason} (the reported result);
\code{theta}, \code{se} and \code{profile_ci_lower}/\code{upper} (the raw
second-order fit); \code{primary_reliable}, \code{confidence_set_type},
\code{confidence_set_bounded} and \code{inference_status}; \code{effect_units} and \code{theta_raw_units};
\code{n_loci}, \code{n_x_median}, \code{n_y_median} and
\code{n_overlap_median}; \code{Q}, \code{Q_p_reported} and \code{Q_p_method};
comparator estimates (\code{theta_two_step_ivw},
\code{theta_weighted_median}, \code{theta_weighted_mode},
\code{theta_independence}); and \code{n_loci_steiger_fail}, the number of loci
whose exposure partial R2 is below the outcome partial R2.

\code{summary()} returns a one-row overview with \code{exposure},
\code{outcome}, \code{n_samples}, \code{n_loci}, the reported \code{estimate},
\code{se}, \code{ci_lower} and \code{ci_upper} (under
\code{"auto_fallback"} these can be the Two-step IVW comparator),
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
\code{\link{as_qtl2mr_prepared_input}} to build the input;
\code{\link{run_qtl2mr_batch}} for many pairs;
\code{\link{refit_qtl2mr_subset}} for follow-up refits.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
prepared <- as_qtl2mr_prepared_input(ex$allele_probs, ex$pheno, ex$covar,
  ex$kinship, ex$marker_map, ex$selected_loci, ex$exposure, ex$outcome)
result <- run_qtl2mr_analysis(prepared, output_dir = tempfile("qtl2mr-"),
  make_figures = FALSE, verbose = FALSE, progress = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L)
qtl2mr_primary_result(result)[, c("theta_reported", "reported_ci_lower",
  "reported_ci_upper", "effect_units", "primary_reliable")]
qtl2mr_warnings(result)[, c("warning_id", "severity")]
}
\keyword{methods}
