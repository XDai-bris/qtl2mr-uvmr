\name{qtl2mr}
\alias{qtl2mr}
\docType{package}
\title{qtl2MR: founder-aware Mendelian randomization for qtl2-style mouse crosses}
\description{
\if{latex}{\out{\sloppy}}
The qtl2MR package fits selected-locus founder-aware linear mixed models,
constructs full founder-contrast covariance blocks, uses star-space overlap
covariance for one-sample or partial-overlap analyses, validates joint positive
semidefiniteness, and estimates an overlap-aware second-order Founder GLS-IVW
effect.
}
\details{
\if{latex}{\out{\sloppy}}
The public API is organized around three layers:
\enumerate{
\item core selected-locus LMM and MR functions;
\item high-level result/report functions for standard local analyses;
\item diagnostic helpers for validation, sensitivity, and QC.
}

Use \code{\link{run_qtl2mr_analysis}()} for the standard prepared-input
analysis entrypoint, or \code{\link{as_qtl2mr_prepared_input}()} to construct
an explicit analysis-ready input object.
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
Development research package for selected-locus multiparent analyses with
additive covariates and harmonized sample identifiers.
}
\keyword{package}
