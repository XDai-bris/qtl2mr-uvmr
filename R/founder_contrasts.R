make_orthonormal_founder_contrasts <- function(founders = LETTERS[1:8]) {
  founders <- as.character(founders)
  if (length(founders) != 8L) stop("Expected exactly 8 founders", call. = FALSE)
  H <- stats::contr.helmert(length(founders))
  H <- apply(H, 2, function(v) v / sqrt(sum(v^2)))
  A <- cbind(mean = rep(1, length(founders)) / sqrt(length(founders)), H)
  rownames(A) <- founders
  colnames(A) <- c("mean", paste0("FC", seq_len(ncol(H))))
  validate_founder_contrasts(A, founders)
  A
}

make_founder_contrasts <- make_orthonormal_founder_contrasts

validate_founder_contrasts <- function(A, founders = rownames(A), tol = 1e-10) {
  if (!is.matrix(A)) A <- as.matrix(A)
  founders <- as.character(founders)
  if (nrow(A) != 8L || ncol(A) != 8L) stop("Founder contrast matrix must be 8 x 8", call. = FALSE)
  if (!identical(rownames(A), founders)) stop("Founder contrast row names must match founders", call. = FALSE)
  if (!identical(colnames(A), c("mean", paste0("FC", 1:7)))) {
    stop("Founder contrast columns must be mean, FC1, ..., FC7", call. = FALSE)
  }
  err <- max(abs(t(A) %*% A - diag(8)))
  if (!is.finite(err) || err > tol) stop("Founder contrast matrix is not orthonormal", call. = FALSE)
  fc_sum <- max(abs(colSums(A[, -1L, drop = FALSE])))
  if (!is.finite(fc_sum) || fc_sum > tol) stop("Founder contrasts are not mean-zero", call. = FALSE)
  invisible(TRUE)
}
