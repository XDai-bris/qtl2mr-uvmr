\name{refit_qtl2mr_subset}
\alias{refit_qtl2mr_subset}
\alias{fit_instrument_set_sensitivity}
\title{Refit an analysis on chosen subsets of its loci}
\description{
\if{latex}{\out{\sloppy}}
Re-estimate the causal effect from an existing analysis using only some of its
selected loci, without refitting the locus-level LMMs. Use these to answer
questions such as "what if I drop the loci near this gene?".
}
\usage{
refit_qtl2mr_subset(result, subset = NULL, reason = NULL,
  allow_cross_cov = TRUE, ridge = 1e-10, set_name = "custom")

fit_instrument_set_sensitivity(result, sets = NULL, loo = NULL,
  weak_f_threshold = 2, top_n = c(10L, 25L, 50L),
  allow_cross_cov = TRUE, ridge = 1e-10)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{result}{A \code{qtl2mr_analysis_result} from
\code{\link{run_qtl2mr_analysis}()}.}
\item{subset}{Character vector of marker IDs to keep. \code{NULL} keeps all loci.}
\item{reason}{Optional note stored with the refit.}
\item{allow_cross_cov}{Use the exposure-outcome covariance from sample overlap
(\code{TRUE}, as in the main analysis).}
\item{ridge}{Small numerical ridge used in matrix inversions.}
\item{set_name}{Label for the refit.}
\item{sets}{Named list of marker vectors. If \code{NULL}, the standard sets are
built: all loci; the 10, 25 and 50 strongest loci by block F (when there are
that many); loci with block F of at least \code{weak_f_threshold}; one
chromosome left out at a time; and, when \code{loo} is given, all but the five
most influential loci. \code{run_qtl2mr_analysis()} uses these sets and adds a
\code{steiger_filtered} row, which this function does not compute.}
\item{loo}{Optional leave-one-locus-out table, \code{result$loo_mr}.}
\item{weak_f_threshold}{Block F threshold used for the standard sets. Block F
is the joint Wald statistic of the 7 founder contrasts divided by 7, so F = 2 is
nominal p of about 0.05. The default changed from 10 to 2 in 0.6.0; the set is
named \code{exclude_F_lt_<threshold>}.}
\item{top_n}{Sizes of the strongest-locus sets.}
}
\details{
\if{latex}{\out{\sloppy}}
Refits reuse the stored locus-level effects, covariances and null-model hsq,
so they are fast and consistent with the main analysis. Subsets selected after
looking at the results are sensitivity analyses, not new primary estimates.
}
\value{
\if{latex}{\out{\sloppy}}
\code{refit_qtl2mr_subset()} returns a one-row data frame with columns
\code{set_name}, \code{n_loci}, \code{theta}, \code{se},
\code{profile_ci_lower}, \code{profile_ci_upper}, \code{Qmin} (minimised Q
statistic), \code{overlap_delta}, \code{status}, \code{reason} and
\code{notes}. \code{fit_instrument_set_sensitivity()} returns one such row per
set as a \code{qtl2mr_instrument_sensitivity} data frame; a set whose refit fails
has \code{status = "failed"} and the error message in \code{notes}.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
prepared <- as_qtl2mr_prepared_input(ex$allele_probs, ex$pheno, ex$covar,
  ex$kinship, ex$marker_map, ex$selected_loci, ex$exposure, ex$outcome)
result <- run_qtl2mr_analysis(prepared, output_dir = tempfile("qtl2mr-"),
  make_figures = FALSE, make_html = FALSE, verbose = FALSE, progress = FALSE,
  heterogeneity_method = "chisq", diagnostic_bootstrap_B = 0L)
refit_qtl2mr_subset(result, subset = "m1", reason = "strongest locus only")
}
\keyword{methods}
