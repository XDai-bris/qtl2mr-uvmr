\name{run_qtl2mr_batch_fast}
\alias{run_qtl2mr_batch_fast}
\alias{run_qtl2mr_batch_exact_grouped}
\alias{run_qtl2mr_batch_parallel_by_exposure}
\alias{qtl2mr-deprecated}
\title{Deprecated batch entry points}
\description{
\if{latex}{\out{\sloppy}}
These functions are deprecated in qtl2mr 0.5.0 and will be removed in 0.6.0.
Each is now a thin wrapper that calls \code{\link{run_qtl2mr_batch}()} with
the equivalent settings and gives the same numerical results as qtl2mr 0.4.x.
Output folders use the \code{run_qtl2mr_batch()} layout.
}
\usage{
run_qtl2mr_batch_fast(resources, pairs,
  instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
  sample_retention_min = 0.85, allow_sample_loss = FALSE,
  fit_options = list(), output_dir = NULL, progress = TRUE,
  fail_fast = FALSE, resume = FALSE, force = FALSE)

run_qtl2mr_batch_exact_grouped(resources, pairs, plan = NULL,
  instrument_rule = list(), fit_options = list(), output_mode = "compact",
  fail_fast = FALSE, n_workers = 1L, output_dir = NULL, progress = TRUE,
  force = FALSE)

run_qtl2mr_batch_parallel_by_exposure(resources, pairs, fast = TRUE,
  n_workers = 2L, instrument_rule = list(lod_thresh = 4, exclude_chr = "X"),
  sample_retention_min = 0.85, allow_sample_loss = FALSE,
  fit_options = list(), output_dir = NULL, progress = TRUE,
  fail_fast = FALSE, resume = FALSE, force = FALSE)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{resources, pairs, instrument_rule, fit_options, output_dir, progress, fail_fast, resume, force, n_workers}{As
in \code{\link{run_qtl2mr_batch}()}.}
\item{sample_retention_min}{Mapped to \code{min_sample_retention}.}
\item{allow_sample_loss}{\code{TRUE} maps to \code{min_sample_retention = 0}.}
\item{fast}{\code{TRUE} maps to \code{sample_policy = "batch_common_complete"}
and \code{on_nonestimable_locus = "drop"}; \code{FALSE} to the defaults.}
\item{plan}{Ignored with a warning; shared fits are reused automatically.}
\item{output_mode}{Ignored.}
}
\details{
\if{latex}{\out{\sloppy}}
Replacement calls:
\tabular{ll}{
\code{run_qtl2mr_batch_fast(r, p)} \tab \code{run_qtl2mr_batch(r, p, sample_policy = "batch_common_complete", on_nonestimable_locus = "drop")} \cr
\code{run_qtl2mr_batch_exact_grouped(r, p)} \tab \code{run_qtl2mr_batch(r, p)} \cr
\code{run_qtl2mr_batch_parallel_by_exposure(r, p, fast = TRUE, n_workers = k)} \tab \code{run_qtl2mr_batch(r, p, sample_policy = "batch_common_complete", on_nonestimable_locus = "drop", n_workers = k)} \cr
\code{run_qtl2mr_batch_parallel_by_exposure(r, p, fast = FALSE, n_workers = k)} \tab \code{run_qtl2mr_batch(r, p, n_workers = k)} \cr
}
A \code{drop_nonestimable_loci} entry in \code{fit_options} keeps its earlier
meaning and is mapped to \code{on_nonestimable_locus}.
}
\value{A \code{qtl2mr_batch_result}, as returned by \code{run_qtl2mr_batch()}.}
\keyword{internal}
