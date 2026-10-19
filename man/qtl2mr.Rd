\name{qtl2mr}
\alias{qtl2mr}
\alias{qtl2mr-package}
\docType{package}
\title{qtl2mr: founder-aware Mendelian randomization in multiparent mice}
\description{
\if{latex}{\out{\sloppy}}
qtl2mr estimates the causal effect of one trait on another in Diversity
Outbred and related eight-founder mouse populations, using selected QTL as
instruments represented by founder haplotype probabilities.
}
\details{
\if{latex}{\out{\sloppy}}
Three main functions cover ordinary use:
\describe{
\item{Instrument selection}{\code{\link{qtl2mr_select_instruments}()} picks
one instrument set per exposure from a qtl2 scan, with the exposure block F of
every candidate, after \code{\link{qtl2mr_prepare_resources}()}.}
\item{One exposure-outcome pair}{\code{\link{as_qtl2mr_prepared_input}()} (or
\code{\link{qtl2mr_pair_input}()} from selected instruments) then
\code{\link{run_qtl2mr_analysis}()}, which writes tables, figures and an HTML
report with sensitivity analyses. Read results with
\code{\link{qtl2mr_primary_result}()} and \code{\link{qtl2mr_warnings}()}.}
\item{Many pairs}{\code{\link{run_qtl2mr_batch}()} fits every pair exactly as
the single-pair analysis would and writes compact tables, with optional
parallel workers and resume.}
}
Follow-up refits on locus subsets use \code{\link{refit_qtl2mr_subset}()} and
\code{\link{fit_instrument_set_sensitivity}()}. The low-level fitting functions
(\code{\link{fit_founder_lmm_block}()}, \code{\link{fit_mr_locus_block}()},
\code{\link{fit_multilocus_founder_mr}()}, \code{\link{second_order_gls_ivw}()})
are an advanced interface for method work.

The vignettes \code{vignette("getting-started", package = "qtl2mr")},
\code{vignette("from-qtl2", package = "qtl2mr")} and
\code{vignette("multilocus-example", package = "qtl2mr")} walk through
installation, preparing your own qtl2 data, and reading a full result.

The method is described in a manuscript in preparation. Numerical success of a
fit is not evidence of calibrated coverage or causal validity; see the Details
of \code{\link{run_qtl2mr_analysis}()}.
}
\keyword{package}
