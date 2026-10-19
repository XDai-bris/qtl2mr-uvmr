\name{validate_joint_psd}
\alias{validate_joint_psd}
\title{Advanced interface: validate joint positive semidefiniteness of MR covariance blocks}
\description{
\if{latex}{\out{\sloppy}}
Check whether \code{Sigma_x}, \code{Sigma_y}, and \code{Sigma_xy} form a valid
joint covariance block for overlap-aware second-order Mendelian randomization.
}
\usage{
validate_joint_psd(Sigma_x, Sigma_y, Sigma_xy, ridge = 1e-10,
  tol = 1e-8, pass_tol = 1e-8)
}
\arguments{
\if{latex}{\out{\sloppy}}
\item{Sigma_x}{Within-exposure covariance matrix for founder-contrast effects.}
\item{Sigma_y}{Within-outcome covariance matrix for founder-contrast effects.}
\item{Sigma_xy}{Cross-trait covariance matrix between exposure and outcome
founder-contrast effects.}
\item{ridge}{Numerical ridge used when forming inverse square-root covariance
matrices for the canonical-correlation check.}
\item{tol}{Eigenvalue tolerance for the full joint covariance matrix.}
\item{pass_tol}{Tolerance for the maximum whitened cross-covariance singular
value.}
}
\details{
\if{latex}{\out{\sloppy}}
The validator applies two equivalent numerical checks: the minimum eigenvalue
of the assembled block covariance matrix and the largest singular value of the
whitened cross-covariance. The public function reports diagnostics only; it
does not repair or clip \code{Sigma_xy}.
}
\section{Lifecycle}{
\if{latex}{\out{\sloppy}}
Validation helper. Use it to audit covariance blocks before or after fitting;
it is not an estimator.
}
\value{
\if{latex}{\out{\sloppy}}
A list with logical \code{psd_valid}, eigenvalue diagnostics, the maximum
whitened cross-covariance singular value, and the singular values used in the
canonical-correlation check.
}
\examples{
# A valid block: whitened cross-covariance singular values below one.
validate_joint_psd(diag(2), diag(2), diag(0.5, 2))$psd_valid
# An impossible block: cross-covariance larger than both variances allow.
validate_joint_psd(diag(2), diag(2), diag(1.5, 2))$psd_valid
}
\keyword{methods}
