\name{qtl2mr_multilocus_example}
\alias{qtl2mr_multilocus_example}
\docType{data}
\title{Synthetic forty-locus founder-MR example}
\description{
\if{latex}{\out{\sloppy}}
A compact, wholly synthetic X1-to-Y example with kinship and complete sample
overlap. It illustrates data preparation, founder LMMs, multi-locus MR and
diagnostic reporting. It is not method-calibration or biological evidence.
}
\usage{
data("qtl2mr_multilocus_example", package = "qtl2mr")
}
\format{
A list containing:
\describe{
\item{allele_probs}{Chromosome-named \code{calc_genoprob} list of sample by
founder by marker arrays: 400 samples, eight founders A--H, 40 markers across
19 synthetic chromosome labels. These are allele, not diplotype, probabilities.}
\item{pheno}{A 400-row data frame with \code{x1} and \code{y}. Row names are
synthetic sample IDs shared by every input.}
\item{covar}{The same samples' \code{batch} and \code{age} covariates.}
\item{kinship}{A fixed 400 by 400 synthetic family-AR matrix M/2; families
have four members and AR coefficient 0.3. This is not LOCO kinship.}
\item{marker_map}{Forty rows with \code{marker}, \code{chr} and \code{pos}.}
\item{selected_loci}{The forty prespecified X1 marker names.}
\item{exposure, outcome}{Trait names \code{"x1"} and \code{"y"}.}
\item{true_theta}{The original structural coefficient, -0.25 in raw units.}
\item{metadata}{Synthetic provenance, source checksum, construction and unit notes.}
}
}
\details{
\if{latex}{\out{\sloppy}}
Use \code{standardize_phenotypes = FALSE} in \code{run_qtl2mr_analysis()} to
report the original units. The metadata records this choice but does not set
analysis options automatically. Heritabilities are estimated by the package;
no fixed hsq or precomputed effect estimates are supplied.

The example is an unchanged subset of the CTC tutorial's 2026-10-15
background-QTL Case A construction (master-seed provenance 2026101501).
It retains all samples and X1 loci. Other traits and background loci were
omitted from the input, not removed from the already generated phenotype.
The original Case A outcome was renamed from \code{y_A} to \code{y} without
altering its values.

Founder probabilities were generated from normalized exponential draws,
not a DO meiosis, HMM or LD model. The forty instruments were specified
before analysis, not selected by a QTL scan. X-only amplitude construction
put all founder-block F statistics above two; this is a teaching criterion,
not a universal strong-instrument threshold. A subsequent algebraic solve
targeted the sample covariate-adjusted observational slope at +0.05 before
MR. This outcome adaptation and the large background effects preclude using
the example as independent calibration or proof of IV exogeneity. The
observational interval includes zero; its point-estimate contrast with MR
is not evidence of a significant opposite observational association.
}
\examples{
data("qtl2mr_multilocus_example", package = "qtl2mr")
d <- qtl2mr_multilocus_example
input <- as_qtl2mr_prepared_input(
  d$allele_probs, d$pheno, covar = d$covar, kinship = d$kinship,
  marker_map = d$marker_map, selected_loci = d$selected_loci,
  exposure = d$exposure, outcome = d$outcome)
validate_qtl2mr_input(input, require_selected_loci = TRUE)
}
\seealso{
\if{latex}{\out{\sloppy}}
\code{\link{qtl2mr_example}}, \code{\link{run_qtl2mr_analysis}},
\code{vignette("multilocus-example", package = "qtl2mr")}}
\keyword{datasets}
