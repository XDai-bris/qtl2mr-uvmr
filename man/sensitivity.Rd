\name{sensitivity}
\alias{sensitivity}
\alias{loo_multilocus_founder_mr}
\alias{plot_loo_multilocus_mr}
\alias{refit_qtl2mr_subset}
\alias{fit_instrument_set_sensitivity}
\alias{plot_instrument_set_sensitivity}
\alias{plot_block_first_stage_dashboard}
\title{qtl2MR sensitivity diagnostics}
\description{
\if{latex}{\out{\sloppy}}
Leave-one-locus-out global MR, selected-locus subset refits, and denominator
strength dashboard plots.
}
\usage{
loo_multilocus_founder_mr(fit, parallel = FALSE,
  allow_cross_cov = TRUE, ridge = 1e-10, ...)

plot_loo_multilocus_mr(loo, top_n = 25L,
  show = c("theta", "delta"),
  main = "Leave-one-locus-out MR sensitivity", xlab = NULL, ...)

refit_qtl2mr_subset(result, subset = NULL, reason = NULL,
  allow_cross_cov = TRUE, ridge = 1e-10, set_name = "custom")

fit_instrument_set_sensitivity(result, sets = NULL, loo = NULL,
  weak_f_threshold = 10, top_n = c(10L, 25L, 50L),
  allow_cross_cov = TRUE, ridge = 1e-10)

plot_instrument_set_sensitivity(sensitivity, global_theta = NULL,
  main = "Instrument-set sensitivity", xlab = "MR effect", ...)

plot_block_first_stage_dashboard(block_first_stage, threshold = 10,
  main = "Block first-stage denominator dashboard", ...)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{fit}{A \code{fit_multilocus_founder_mr()} result or
\code{qtl2mr_analysis_result}.}
\item{parallel}{Reserved for future parallel LOO refits.}
\item{allow_cross_cov}{Whether to use \code{Sigma_xy}.}
\item{ridge}{Numerical ridge for refits.}
\item{loo}{LOO table returned by \code{loo_multilocus_founder_mr()}.}
\item{top_n}{Maximum number of LOO rows to plot, or top-N block-F instrument
sets to construct.}
\item{show}{Whether the LOO plot displays \code{theta[-m]} or the delta from
the full estimate.}
\item{main, xlab}{Plot labels.}
\item{result}{A \code{qtl2mr_analysis_result} or compatible object containing
locus blocks.}
\item{subset}{Marker vector to refit. If omitted, all loci are used.}
\item{reason}{Optional text describing the subset.}
\item{set_name}{Name for a single subset refit.}
\item{sets}{Named list of marker vectors. If omitted, default sets are built
from selected loci, chromosome metadata, block-F values, and optional LOO
influence.}
\item{weak_f_threshold}{Heuristic block-F threshold for weak-locus exclusions.}
\item{sensitivity}{Instrument-set sensitivity table.}
\item{global_theta}{Optional full global estimate reference line.}
\item{block_first_stage}{Table containing \code{F_X_block}.}
\item{threshold}{Heuristic first-stage threshold shown on the dashboard.}
\item{...}{Additional graphical arguments.}
}
\details{
\if{latex}{\out{\sloppy}}
LOO MR removes one selected locus at a time and refits the same overlap-aware
second-order estimator. Instrument-set sensitivity refits selected-locus
subsets without rerunning LOD pruning. Both are diagnostics and should not be
reported as alternative causal estimands.
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
Diagnostic helpers. Outputs are intended for sensitivity review and package
reports, not as replacement causal estimators.
}
\value{
\if{latex}{\out{\sloppy}}
\code{loo_multilocus_founder_mr()} returns a data frame with one row per locus.
\code{fit_instrument_set_sensitivity()} returns a data frame with one row per
instrument set. Plot functions invisibly return plotted data.
}
\keyword{methods}
