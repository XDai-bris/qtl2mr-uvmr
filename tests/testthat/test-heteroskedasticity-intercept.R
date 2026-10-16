test_that("score-adjusted residual mean model includes an intercept by default", {
  ids <- paste0("id", 1:40)
  score <- stats::setNames(seq(-1, 1, length.out = length(ids)), ids)
  set.seed(531)
  y <- stats::setNames(2.5 + 0.4 * score - 0.15 * score^2 + rnorm(length(ids), sd = 0.02), ids)

  resid <- score_adjusted_residuals(y, score, covariates = NULL)
  reference <- stats::lm(y ~ score + I(score^2))

  expect_lt(abs(mean(resid$residual)), 1e-10)
  expect_equal(as.numeric(resid$residual), as.numeric(stats::resid(reference)), tolerance = 1e-10)
  expect_match(unique(resid$mean_model_used), "intercept")
})
