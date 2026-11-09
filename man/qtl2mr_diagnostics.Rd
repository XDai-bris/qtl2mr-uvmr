\name{qtl2mr_diagnostics}
\alias{qtl2mr_diagnostics}
\title{Simulated Q p-value and bootstrap standard errors for one analysis}
\description{
\if{latex}{\out{\sloppy}}
Add two simulation-based diagnostics to a finished
\code{\link{run_qtl2mr_analysis}()} result: the Q heterogeneity p-value
calibrated by simulation, and bootstrap standard errors and intervals for the
weighted median and weighted mode comparators.
}
\usage{
qtl2mr_diagnostics(result, nsim = 1000L, seed = 1L)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{result}{A \code{qtl2mr_analysis_result} from
\code{\link{run_qtl2mr_analysis}()}, its output directory, or the path of its
\code{result_object.rds}.}
\item{nsim}{Number of Q simulation draws, and also the number of bootstrap
draws (a whole number, at least 1).}
\item{seed}{Seed for both the simulation and the bootstrap (a whole number).
The global random-number state is restored afterwards.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{run_qtl2mr_analysis()} refers Q to the chi-square distribution and gives
the weighted median and mode as point estimates. This function works post hoc
from the locus blocks stored in the result; it does not refit the locus LMMs.
\itemize{
\item The Q p-value is calibrated by conditional residual simulation under the
fitted common-effect null. Each draw re-solves the multilocus fit, so this is
the slow part: roughly 0.3 seconds per draw on a laptop for the bundled
40-locus example, about 5 minutes with the default 1,000 draws.
\item The weighted median and mode standard errors come from a parametric
bootstrap: loci and weights stay fixed, and each per-locus effect is sampled
independently from its estimate and standard error. The intervals are the
estimate plus or minus 1.96 bootstrap standard errors
(\code{ci_type = "bootstrap-SE Wald"}). They need at least two loci.
}
These diagnostics are not profile intervals and carry no coverage guarantee.
The main result, its tables and \code{result_object.rds} are not changed. At
the same \code{nsim} and \code{seed} the numbers equal those of qtl2mr
0.6.1.9005 run with \code{heterogeneity_method = "simulation"},
\code{heterogeneity_nsim = nsim} and \code{diagnostic_bootstrap_B = nsim}, each
with that seed.

For a pair of a \code{\link{run_qtl2mr_batch}()} run, first rerun the pair with
\code{\link{qtl2mr_pair_input}()} and \code{run_qtl2mr_analysis()}, then call
this function on that result.
}
\value{
\if{latex}{\out{\sloppy}}
A one-row data frame with \code{exposure}, \code{outcome}, \code{n_loci},
\code{nsim} and \code{seed}; \code{Q}, \code{Q_chisq_df} and
\code{Q_p_chisq} (as in the main result); the simulation fields
\code{Q_p_sim}, \code{Q_sim_n}, \code{Q_sim_n_requested}, \code{Q_sim_seed},
\code{Q_sim_exceedances}, \code{Q_sim_min_resolvable_p} (the smallest p-value
the draws can show), \code{Q_sim_null_mean}, \code{Q_sim_null_sd} and
\code{Q_sim_method}; \code{theta_weighted_median}, \code{se_weighted_median},
\code{ci_lower_weighted_median}, \code{ci_upper_weighted_median} and
\code{ci_type_weighted_median}, the same five columns for the weighted mode,
and \code{diagnostic_bootstrap_B} and \code{diagnostic_bootstrap_seed}. The
table is also written to \code{data/diagnostics.csv} in the result folder: the
folder given as \code{result} (or the folder of the given
\code{result_object.rds}), or, for a result object, its recorded
\code{output_dir}, when that folder exists. The file is not added to
\code{output_manifest.csv}, and a second call replaces it.
}
\seealso{
\code{\link{run_qtl2mr_analysis}}, \code{\link{qtl2mr_method_comparison}}.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno,
  ex$covar, ex$marker_map)
input <- qtl2mr_pair_input(resources,
  data.frame(exposure = "x", marker = c("m1", "m2")), "x", "y")
result <- run_qtl2mr_analysis(input, output_dir = tempfile("qtl2mr-"),
  make_html = FALSE, make_figures = FALSE, progress = FALSE)
# A small nsim keeps the example fast; use the default 1,000 in an analysis.
diagnostics <- qtl2mr_diagnostics(result, nsim = 100L)
diagnostics[, c("Q", "Q_p_chisq", "Q_p_sim", "se_weighted_median",
  "se_weighted_mode")]
}
\keyword{methods}
