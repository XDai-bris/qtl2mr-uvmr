\name{qtl2mr_example}
\alias{qtl2mr_example}
\docType{data}
\title{Small deterministic synthetic qtl2mr teaching fixture}
\description{
\if{latex}{\out{\sloppy}}
A compact, deterministic synthetic prepared-input fixture for the getting-started
vignette and integration examples. It contains no animal, application, or
biological data and is not method calibration or validation evidence.
}
\usage{
data("qtl2mr_example", package = "qtl2mr")
}
\format{
A list with \code{allele_probs}, \code{pheno}, \code{covar}, \code{kinship},
\code{marker_map}, \code{selected_loci}, \code{exposure}, \code{outcome},
\code{true_theta}, and \code{metadata}. \code{allele_probs} is a
chromosome-named \code{calc_genoprob} list holding a 96 by 8 by 4
sample-by-founder-by-marker array (markers \code{m1} to \code{m4} on chromosome
1). \code{pheno} has the exposure \code{x} and outcome \code{y}; \code{covar}
has \code{sex} and \code{age}; \code{marker_map} has \code{marker}, \code{chr},
\code{pos} and \code{exposure_lod}; \code{selected_loci} is \code{c("m1", "m2")}. \code{kinship} is deliberately \code{NULL}
only to keep this teaching fixture self-contained.
}
\details{
\if{latex}{\out{\sloppy}}
The fixture uses seed 301 and normalized exponential founder probabilities.
The exposure is generated from founder contrasts at the first three markers;
the outcome uses a structural coefficient of 0.75 and an additional fourth
marker founder contrast. The listed \code{exposure_lod} values are illustrative
marker labels, not qtl2-derived calibration.

For real Diversity Outbred analyses, users should harmonize sample IDs,
select instruments, for example with \code{\link{qtl2mr_select_instruments}()}
from a qtl2 scan, and provide harmonized chromosome-specific LOCO kinship. The fixture's \code{kinship = NULL} is not a recommended real-data
analysis setting.
}
\examples{
data("qtl2mr_example", package = "qtl2mr")
str(qtl2mr_example, max.level = 1)
}
\keyword{datasets}
