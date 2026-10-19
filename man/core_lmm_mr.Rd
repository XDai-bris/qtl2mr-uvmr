\name{core_lmm_mr}
\alias{core_lmm_mr}
\alias{make_founder_contrasts}
\alias{fit_founder_lmm_block}
\alias{fit_mr_locus_block}
\alias{fit_multilocus_founder_mr}
\alias{second_order_gls_ivw}
\title{Advanced interface: founder LMM blocks and second-order MR}
\description{
\if{latex}{\out{\sloppy}}
Low-level founder-aware LMM blocks at single markers and the overlap-aware
second-order estimator. They are intended for method development and checks;
most analyses should use the high-level functions
(\code{\link{qtl2mr_select_instruments}()}, \code{\link{run_qtl2mr_analysis}()}
and \code{\link{run_qtl2mr_batch}()}), which add input validation, sample
accounting, diagnostics and reporting.
}
\usage{
make_founder_contrasts(founders = LETTERS[1:8])

fit_founder_lmm_block(genoprobs, pheno, kinship = NULL,
  addcovar = NULL, contrasts = NULL, hsq = NULL, marker = NULL,
  sample_ids = NULL, standardize_pheno = FALSE, reml = TRUE,
  fc_cols = paste0("FC", 1:7), tol = 1e-12)

fit_mr_locus_block(genoprobs, exposure, outcome, kinship = NULL,
  addcovar = NULL, contrasts = NULL, hsq_x = NULL, hsq_y = NULL,
  marker = NULL, sample_ids = NULL, standardize_phenotypes = FALSE,
  reml = TRUE, fc_cols = paste0("FC", 1:7), ridge = 1e-10,
  cross_cov_model = c("star", "raw", "zero"), enforce_psd = TRUE,
  stop_on_invalid_psd = TRUE, tol = 1e-12)

fit_multilocus_founder_mr(locus_blocks, allow_cross_cov = TRUE,
  ridge = 1e-10, theta_grid = NULL, compute_independence = TRUE,
  heterogeneity_method = c("auto", "simulation", "chisq", "none"),
  q_sim_n = 0L, q_sim_seed = 1L, keep_q_sim = FALSE,
  progress_callback = NULL)

second_order_gls_ivw(beta_x, beta_y, Sigma_x, Sigma_y,
  theta_grid = NULL, allow_cross_cov = FALSE, Sigma_xy = NULL,
  ridge = 1e-10)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{founders}{Founder labels. The high-level workflow assumes eight founders for the
standard autosomal additive model.}
\item{genoprobs}{qtl2-style founder allele probabilities: a named list of
sample-by-founder-by-marker arrays.}
\item{pheno}{Named phenotype vector for one trait.}
\item{exposure, outcome}{Named exposure and outcome phenotype vectors.}
\item{kinship}{Optional LOCO kinship matrix for the selected chromosome.}
\item{addcovar}{Optional additive covariate matrix or data frame with sample
row names.}
\item{contrasts}{Founder contrast matrix. If omitted,
\code{make_founder_contrasts()} is used.}
\item{hsq, hsq_x, hsq_y}{Optional trait-specific heritability values used to
construct the LMM residual covariance. If omitted and \code{kinship} is
supplied, the internal qtl2mr profile ML/REML estimator is used.}
\item{marker}{Selected marker ID. Required when \code{genoprobs} contains more
than one marker.}
\item{sample_ids}{Optional sample IDs and order for fitting.}
\item{standardize_pheno, standardize_phenotypes}{Whether to standardize traits
before fitting.}
\item{reml}{Whether internal heritability estimation should use REML.}
\item{fc_cols}{Founder contrast coefficient columns used as the MR block.}
\item{tol}{QR rank tolerance.}
\item{ridge}{Numerical ridge used in covariance inversions and PSD checks.}
\item{cross_cov_model}{Cross-covariance model. Public LMM analyses should use
\code{"star"}. \code{"raw"} is retained as technical QC only, and
\code{"zero"} gives the independence sensitivity.}
\item{enforce_psd}{Whether to validate the assembled joint covariance block.}
\item{stop_on_invalid_psd}{Whether to stop when the selected cross-covariance
model has a large joint PSD violation.}
\item{locus_blocks}{List of locus blocks returned by
\code{fit_mr_locus_block()} or \code{fit_mr_locus_block_overlap()}.}
\item{allow_cross_cov}{Whether to include \code{Sigma_xy} in the second-order
objective.}
\item{theta_grid}{Optional finite effect values added as search hints. These
do not restrict the search domain, which includes both infinite tails.}
\item{compute_independence}{Whether to compute a \code{Sigma_xy = 0}
sensitivity estimate.}
\item{heterogeneity_method}{How to choose the reported aggregate Q
heterogeneity p-value. \code{"auto"} reports a chi-square reference for
single-locus fits and simulation calibration for multi-locus fits.}
\item{q_sim_n}{Number of conditional residual simulations for the aggregate
Q heterogeneity p-value. The low-level default is \code{0L};
\code{run_qtl2mr_analysis()} runs 1,000 draws by default for several loci
(its \code{heterogeneity_nsim}), while \code{run_qtl2mr_batch()} defaults to the
chi-square reference. If simulation is selected with
\code{q_sim_n = 0L}, the low-level fit warns and returns \code{Q_p_reported =
NA}.}
\item{q_sim_seed}{Random seed for Q simulation calibration.}
\item{keep_q_sim}{Whether to retain the simulated Q values in the returned
object.}
\item{progress_callback}{Optional function accepting \code{stage},
\code{label} and \code{detail}. Receives coarse Q-simulation counts and
elapsed time. Callbacks must not generate random numbers or change
\code{RNGkind()}; some generators retain state outside \code{.Random.seed}.
The built-in analysis reporter is RNG-free and leaves simulated values unchanged.}
\item{beta_x, beta_y}{Founder-contrast exposure and outcome effect vectors.}
\item{Sigma_x, Sigma_y}{Within-trait covariance matrices.}
\item{Sigma_xy}{Cross-trait covariance matrix. Required when
\code{allow_cross_cov = TRUE}.}
}
\details{
\if{latex}{\out{\sloppy}}
\code{fit_founder_lmm_block()} fits one selected marker in founder-contrast
space and returns founder effects, covariance matrices, residuals, and
influence operators. Kinship is decomposed internally after forming
\eqn{2K}; \pkg{qtl2} is not required for primary fitting.

\code{fit_mr_locus_block()} fits exposure and outcome on joint complete
samples and constructs \code{Sigma_x}, \code{Sigma_y}, and one-sample
overlap-aware \code{Sigma_xy}. The high-level functions do not use it: they fit
each trait on its own complete samples and model the partial overlap, so their
estimates differ from this function when exposure and outcome have different
missing samples. The default \code{cross_cov_model = "star"}
matches the manuscript method for LMM analyses. Raw-space covariance is
reported only as technical QC and is not a public estimator.

\code{fit_multilocus_founder_mr()} minimizes
\deqn{
Q(\theta)=\sum_m d_m(\theta)^T\Omega_m(\theta)^{-1}d_m(\theta)
}
across selected locus blocks, where
\deqn{
\Omega_m(\theta)=\Sigma_{Y,m}+\theta^2\Sigma_{X,m}
-\theta(\Sigma_{XY,m}+\Sigma_{XY,m}^T).
}
The chi-square reference field is always available as \code{Q_p_chisq}. When
\code{q_sim_n > 0}, the simulation field \code{Q_p_sim} is calibrated by
conditional residual simulation under the fitted common-effect null. The
reported p-value is \code{Q_p_reported}, with method label \code{Q_p_method}.

Optimization uses trait-normalized coordinates and a compact angular search
with refined candidate extrema and explicit tail limits. This finite numerical
search is not a proof of the global optimum; \code{search_certified} is false.
The profile confidence set uses a nominal chi-square cutoff. It can be
unbounded, disconnected or the whole real line. These cases must not be
interpreted as an ordinary finite interval.

The legacy \code{primary_reliable} field only denotes a numerically acceptable
fit with a single bounded confidence-set component. It does not establish
nominal coverage, adequate instrument strength or causal validity.
\code{inference_status = "not_calibrated"} deliberately keeps that distinction.
Weak instruments and same-sample selection can invalidate nominal profile
coverage. The summed objective assumes no cross-locus covariance; overlap
correction within loci does not remove that assumption.
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
Core functions for the supported autosomal additive estimator.
}
\value{
\if{latex}{\out{\sloppy}}
\code{make_founder_contrasts()} returns a founder contrast matrix.
\code{fit_founder_lmm_block()} and \code{fit_mr_locus_block()} return fitted
blocks. \code{fit_multilocus_founder_mr()} and
\code{second_order_gls_ivw()} return MR estimate lists with effect, standard
error, profile interval, objective diagnostics, optimizer status, and, when
requested, simulation-calibrated Q heterogeneity fields.
\code{confidence_set} is a data frame of lower/upper endpoints, including
infinities. \code{confidence_set_type}, \code{numerical_ok},
\code{objective_tail} and \code{search_method} describe numerical behavior.
The legacy \code{profile_ci_lower} and \code{profile_ci_upper} do not replace
the component table for a disconnected set. Curvature \code{se} and Wald
\code{p} are different from profile uncertainty.
}
\examples{
bx <- c(0.8, -0.5, 0.6)
by <- 0.5 * bx + c(0.01, -0.02, 0.01)
fit <- second_order_gls_ivw(bx, by, diag(0.01, 3), diag(0.02, 3))
fit[c("theta", "confidence_set_type", "numerical_ok", "inference_status")]
fit$confidence_set
}
\keyword{methods}
