\name{qtl2mr_combine_resources}
\alias{qtl2mr_combine_resources}
\title{Combine two resources objects}
\description{
\if{latex}{\out{\sloppy}}
Combine two objects from \code{\link{qtl2mr_prepare_resources}()}, for example
two data sets measured on overlapping mice, into one resources object that
\code{\link{qtl2mr_select_instruments}()}, \code{\link{qtl2mr_pair_input}()} and
\code{\link{run_qtl2mr_batch}()} accept like any other.
}
\usage{
qtl2mr_combine_resources(x, y, kinship = NULL, common_only = FALSE)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{x, y}{Objects returned by \code{\link{qtl2mr_prepare_resources}()}.}
\item{kinship}{Kinship for the combined samples, in the form accepted by
\code{\link{qtl2mr_prepare_resources}()}, for example
\code{qtl2::calc_kinship(probs, type = "loco")} on the combined allele
probabilities. \code{NULL} (default) takes the kinship of \code{x} and
\code{y} when that is possible; see Details.}
\item{common_only}{Keep only the genotyped samples present in both \code{x}
and \code{y}.}
}
\details{
\if{latex}{\out{\sloppy}}
The inputs must agree; each failed check is an error, and nothing is repaired
silently:
\itemize{
\item The marker map (\code{marker}, \code{chr} and \code{pos}, in the same
order), the chromosomes and markers of the allele probabilities, and the
founders (in the same order) must be identical.
\item No trait name may appear in both inputs; rename traits before combining.
\item Both inputs must have covariates with the same columns, or both none.
\item A sample ID present in both inputs means the same mouse, so its allele
probabilities (on every chromosome) and its covariates must be identical.
}
The combined object holds the samples of \code{x} followed by those only in
\code{y}. A sample present in only one input has \code{NA} for the traits of
the other input, so each trait is fitted on its own complete samples. A pair
whose exposure comes from \code{x} and whose outcome comes from \code{y} is then
a two-sample MR with partial (or no) sample overlap, which the estimator models
from the shared sample IDs.

\strong{Kinship.} When \code{kinship = NULL} and both inputs have a kinship that
covers all combined samples with identical values (one cohort, as when both
data sets come from the same genotyped mice), it is taken from them. When
neither input has a kinship, the combined object has none. Otherwise supply
\code{kinship} computed on the combined probabilities; qtl2mr does not compute
kinship. This includes kinships that each cover only their own input's samples
(for example two cohorts genotyped separately): the kinship between a sample
only in \code{x} and one only in \code{y} is then unknown. With qtl2, the
combined probabilities can be formed with \code{rbind()} of the two
\code{calc_genoprob} objects (or taken from the combined resources'
\code{allele_probs}) and passed to \code{qtl2::calc_kinship(probs, "loco")}.

Numbers are compared by value, so the same covariate value or marker position
stored as integer on one side and double on the other is accepted. Each input
must hold the same samples, in the same order, on every chromosome. The
combined object holds the
two inputs and the combined probabilities in memory at once, which for a full
Diversity Outbred cohort is a few times the size of one probability array.

\strong{One common sample set.} With \code{common_only = TRUE} only the
genotyped samples present in both inputs are kept. Each trait is still fitted
on its own complete samples within that set, so a sample missing one trait is
still used for the other. In a batch, every pair is restricted to the same
common set; to fit each pair on the samples complete for both of its traits,
combine or subset the inputs for that pair before
\code{\link{qtl2mr_pair_input}()}.

\strong{Covariates.} No cohort or data-set covariate is added. Whether one is
needed is a scientific decision; add it to \code{covar} before
\code{\link{qtl2mr_prepare_resources}()}.

The combined object is built with \code{\link{qtl2mr_prepare_resources}()},
which checks it as any other input (including the founder-probability
rounding check).
}
\value{
\if{latex}{\out{\sloppy}}
A \code{qtl2mr_resources} object. Its \code{metadata$combined_resources}
records \code{n_samples_x_only}, \code{n_samples_y_only} and
\code{n_samples_overlap} (genotyped samples of the inputs),
\code{n_samples_combined} (genotyped samples kept), \code{common_only},
\code{traits_x}, \code{traits_y} and \code{kinship_source} (\code{"x_and_y"},
\code{"supplied"} or \code{"none"}); \code{metadata$x_metadata} and
\code{metadata$y_metadata} keep the metadata of the inputs.
}
\seealso{
\code{\link{qtl2mr_prepare_resources}}, \code{\link{qtl2mr_pair_input}},
\code{\link{run_qtl2mr_batch}}.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example
ids <- rownames(ex$pheno)
# Two synthetic data sets on overlapping mice: x measured on the first 64
# mice, y on the last 64, so 32 mice have both traits.
part <- function(keep, trait) {
  qtl2mr_prepare_resources(
    lapply(ex$allele_probs, function(a) a[keep, , , drop = FALSE]),
    kinship = NULL, pheno = ex$pheno[keep, trait, drop = FALSE],
    covar = ex$covar[keep, , drop = FALSE], marker_map = ex$marker_map)
}
both <- qtl2mr_combine_resources(part(ids[1:64], "x"), part(ids[33:96], "y"))
unlist(both$metadata$combined_resources[c("n_samples_x_only",
  "n_samples_y_only", "n_samples_overlap", "n_samples_combined")])
input <- qtl2mr_pair_input(both,
  data.frame(exposure = "x", marker = c("m1", "m2")), "x", "y")
result <- run_qtl2mr_analysis(input, output_dir = tempfile("qtl2mr-"),
  make_html = FALSE, make_figures = FALSE, progress = FALSE)
qtl2mr_primary_result(result)[, c("theta_reported", "reported_ci_lower",
  "reported_ci_upper", "n_x_median", "n_y_median", "n_overlap_median")]
}
\keyword{methods}
