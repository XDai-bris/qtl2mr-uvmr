\name{diagnostics}
\alias{diagnostics}
\alias{fit_multilocus_diagnostic_layer}
\alias{founder_residual_scores}
\alias{founder_coordinate_ratios}
\alias{fit_xprs_heteroskedasticity_diagnostic}
\alias{plot_per_locus_so_forest}
\alias{plot_founder_residual_heatmap}
\alias{plot_founder_ratio_forest}
\alias{plot.qtl2mr_xprs_heteroskedasticity_diagnostic}
\title{qtl2MR diagnostic helpers}
\description{
\if{latex}{\out{\sloppy}}
Diagnostic helpers for per-locus consistency, founder-coordinate summaries,
local founder exposure score residual-variance QC, and plotting.
}
\usage{
fit_multilocus_diagnostic_layer(locus_blocks, global_fit = NULL,
  allow_cross_cov = TRUE, ridge = 1e-10, level = 0.95,
  bootstrap_B = 500L, bootstrap_seed = 1L)

founder_residual_scores(block, theta = NULL, H = NULL,
  allow_cross_cov = TRUE, ridge = 1e-10)

founder_coordinate_ratios(block, H = NULL, theta_ref = NULL,
  allow_cross_cov = TRUE, weak_f = 10, level = 0.95)

fit_xprs_heteroskedasticity_diagnostic(traits, xprs,
  covariates = NULL, sample_ids = NULL, trait_roles = NULL,
  per_locus_scores = NULL, marker_info = NULL, include_loo = TRUE,
  include_localization = !is.null(per_locus_scores), n_bins = 5L,
  B = 600L, seed = NULL, include_score2 = TRUE,
  n_loci = NA_integer_)

plot_per_locus_so_forest(per_locus, estimates = NULL,
  global_theta = NULL, max_loci = 40L, top_n = NULL,
  order_by = c("F", "position", "effect", "se", "weight"),
  main = "Per-locus SO estimates", xlab = "MR effect", ...)

plot_founder_residual_heatmap(residual_scores,
  value = c("t_one_founder_vs_rest", "e", "z_diag"),
  row = c("marker", "trait"), main = "Founder residual diagnostics",
  label_threshold = 2, zlim = NULL, palette = NULL, ...)

plot_founder_ratio_forest(ratio_table, group = NULL, xlim = NULL,
  main = "Founder diagnostic ratios",
  xlab = "Founder-coordinate Wald ratio", ...)

\method{plot}{qtl2mr_xprs_heteroskedasticity_diagnostic}(x,
  type = c("global", "dashboard", "loo", "bin", "localization"),
  top_n = 12L, ...)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{locus_blocks}{List of fitted locus blocks.}
\item{global_fit}{Optional result from \code{fit_multilocus_founder_mr()}.}
\item{allow_cross_cov}{Whether diagnostics use \code{Sigma_xy}.}
\item{ridge}{Numerical ridge for covariance inversions.}
\item{level}{Confidence level.}
\item{bootstrap_B}{Parametric-bootstrap draws for weighted median/mode
diagnostic intervals. Use \code{0L} for point estimates only or an integer of
at least two. The default is \code{500L}.}
\item{bootstrap_seed}{Random seed for weighted median/mode parametric bootstrap.}
\item{block}{A locus block returned by \code{fit_mr_locus_block()} or
\code{fit_mr_locus_block_overlap()}.}
\item{theta}{Candidate or fitted MR effect. If omitted where allowed, a local
second-order fit is computed.}
\item{H}{Founder contrast matrix mapping seven founder contrasts back to
founder coordinates. If omitted, \code{block$founder_H} is used when
available.}
\item{theta_ref}{Reference global or local effect for founder-coordinate ratio
tables.}
\item{weak_f}{Heuristic threshold for marking weak denominator coordinates.}
\item{traits}{Named list of phenotype vectors for local founder exposure score
residual-variance diagnostics.}
\item{xprs}{A named local founder exposure score vector or an object containing
a summed score and optional per-locus scores.}
\item{covariates}{Optional covariate matrix or data frame.}
\item{sample_ids}{Optional sample ID subset.}
\item{trait_roles}{Optional labels such as \code{"exposure"} and
\code{"outcome"}.}
\item{per_locus_scores}{Optional sample-by-locus matrix of local founder
exposure scores. When supplied, per-locus tests are returned as localization
diagnostics.}
\item{marker_info}{Optional marker metadata containing marker, chromosome, and
rank columns.}
\item{include_loo}{Whether to compute leave-one-locus-out local founder score
sensitivity tests.}
\item{include_localization}{Whether to compute per-locus localization tests.}
\item{n_bins}{Number of score bins for residual-variance summaries.}
\item{B}{Bootstrap replicates for score-bin variance intervals.}
\item{seed}{Optional random seed.}
\item{include_score2}{Include a quadratic score term in residual-variance
tests.}
\item{n_loci}{Number of loci represented by the summed local score.}
\item{per_locus}{Data frame returned by the diagnostic layer.}
\item{estimates}{Optional estimator summary table.}
\item{global_theta}{Optional full global estimate reference line.}
\item{max_loci, top_n}{Maximum number of loci shown in large-M plots.}
\item{order_by}{Selection and ordering rule for large-M per-locus displays.}
\item{residual_scores}{Founder residual score table.}
\item{ratio_table}{Founder-coordinate ratio table.}
\item{group}{Optional grouping variable for ratio plots.}
\item{xlim}{Optional x-axis limits.}
\item{value}{Column used for heatmap color.}
\item{row}{Heatmap row grouping.}
\item{main, xlab}{Plot labels.}
\item{label_threshold}{Absolute value threshold for heatmap labels.}
\item{zlim}{Optional color scale limits.}
\item{palette}{Optional color palette.}
\item{x}{A \code{qtl2mr_xprs_heteroskedasticity_diagnostic} object for the S3
plot method.}
\item{type}{Plot type for the S3 method.}
\item{...}{Additional graphical arguments.}
}
\details{
\if{latex}{\out{\sloppy}}
Weighted median/mode intervals are normal approximations: the original point
estimate plus or minus the normal quantile times the bootstrap standard
deviation (1.96 at 95\%). All included loci and inverse-variance weights are
held fixed. Each draw samples per-locus effects independently from normal
distributions centered on their fitted effects, with their fitted SEs, then
recomputes the median/mode. No mice or loci are resampled. These are
\code{bootstrap-SE Wald} diagnostic intervals, not profile-likelihood intervals;
nominal coverage is not established. Covariance-estimation uncertainty and
cross-locus dependence are not propagated. Failed draws are counted and
suppress the affected interval rather than being silently discarded.
The weighted mode uses the same Gaussian kernel and bandwidth rule, with
local grids and peak refinement to avoid losing a peak to a distant outlier.

The local founder exposure score is
\eqn{S_X = Z_m \hat b_X} for a locus, or the sum of selected-locus founder
scores for a multi-locus diagnostic. It is local to the selected
founder-block instruments and is not a genome-wide polygenic score.

The residual-variance diagnostic checks whether residual scale changes across
this local score after mean adjustment. It is a QC diagnostic, not an
alternative causal estimator.
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
Diagnostic functions. Outputs should be used for QC, localization, and
sensitivity review; they are not primary causal estimators.
}
\value{
\if{latex}{\out{\sloppy}}
Diagnostic functions return data frames or lists of diagnostic tables. Plot
functions draw base R plots and invisibly return the plotted data.
The diagnostic-layer \code{estimates} table includes \code{theta}, \code{se},
\code{ci_lower}, \code{ci_upper}, \code{ci_type}, \code{ci_level}, and
\code{bootstrap_B}, \code{bootstrap_seed}, \code{bootstrap_complete_B},
\code{bootstrap_failed_B}. An unavailable interval has missing bounds and an
explanation in \code{note}. \code{level} applies to the scalar diagnostics;
the primary profile interval remains at its fitted 95\% level.
}
\keyword{methods}
