test_that("founder contrast matrix is orthonormal", {
  A <- make_founder_contrasts(LETTERS[1:8])
  expect_equal(dim(A), c(8L, 8L))
  expect_equal(colnames(A), c("mean", paste0("FC", 1:7)))
  expect_lt(max(abs(t(A) %*% A - diag(8))), 1e-12)
  expect_lt(max(abs(colSums(A[, -1, drop = FALSE]))), 1e-12)
})
