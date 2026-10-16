test_that("scalar-normal bootstrap keeps loci and weights fixed and uses original points", {
  x <- c(-0.4, -0.1, 0.15, 0.3, 0.7)
  se <- c(.08, .3, .12, .2, .25)
  w <- 1 / se^2
  B <- 80L
  seed <- 19L
  med <- function(x) {
    i <- order(x)
    x[i][which(cumsum(w[i] / sum(w)) >= .5)[1L]]
  }
  mode <- qtl2mr:::qtl2mr_weighted_mode
  draws <- qtl2mr:::qtl2mr_with_seed(seed, {
    vapply(seq_len(B), function(i) {
      y <- x + se * rnorm(length(x))
      c(med(y), mode(y, w))
    }, c(0, 0))
  })
  point <- c(med(x), mode(x, w))
  set.seed(704)
  before <- .Random.seed
  actual <- qtl2mr:::qtl2mr_bootstrap_scalar_sensitivities(x, se, B, seed)
  expect_identical(.Random.seed, before)
  expect_equal(actual$se, apply(draws, 1, sd), tolerance = 1e-12)
  expect_equal(actual$ci_lower, point - qnorm(.975) * actual$se)
  expect_equal(actual$ci_upper, point + qnorm(.975) * actual$se)
  expect_equal((actual$ci_lower + actual$ci_upper) / 2, point)
  expect_true(all(actual$ci_type == "bootstrap-SE Wald"))
  expect_true(all(actual$bootstrap_complete_B == B & actual$bootstrap_failed_B == 0L))
  narrower <- qtl2mr:::qtl2mr_bootstrap_scalar_sensitivities(x, se, B, seed, level = .9)
  expect_equal(narrower$se, actual$se)
  expect_equal(narrower$ci_lower, point - qnorm(.95) * actual$se)
  expect_equal(narrower$ci_upper, point + qnorm(.95) * actual$se)
})

test_that("mode search resolves narrow peaks despite a distant effect", {
  x <- c(seq(-.5, .5, length.out = 39), 10000)
  w <- c(seq(1, 2, length.out = 39), 1e-5)
  h <- max(bw.nrd0(x), .03)
  density <- function(t) sum(w / sum(w) * dnorm((t - x) / h)) / h
  reference <- optimize(function(t) -density(t), c(-1, 1), tol = 1e-10)$minimum
  mode <- qtl2mr:::qtl2mr_weighted_mode
  actual <- mode(x, w)
  expect_equal(actual, reference, tolerance = 1e-7)
  expect_equal(actual, mode(x, w, grid_step = 1 / 16), tolerance = 1e-7)
  expect_equal(actual, mode(rev(x), rev(w)), tolerance = 1e-7)
  expect_equal(mode(rep(.7, 4), c(1, 2, 3, 4)), .7, tolerance = 1e-7)
  expect_equal(mode(c(-2, -2.1, 2, 2.1), c(1, 1, 3, 3), bandwidth = .1),
               2.05, tolerance = 1e-7)
  expect_error(mode(c(0, 1), c(1, 1), bandwidth = 0), "bandwidth")
  expect_error(mode(c(0, 1), 1), "equal lengths")
})

test_that("bootstrap validation and failed draws never yield a truncated-sample CI", {
  f <- qtl2mr:::qtl2mr_bootstrap_scalar_sensitivities
  for (B in list(-1, 1, 2.5, NA_real_, c(2, 3), Inf)) expect_error(f(c(0, 1), c(.1, .2), B), "B must")
  expect_error(f(c(0, 1), .1), "matching")
  expect_error(f(c(0, 1), c(.1, 0)), "positive se")
  expect_error(f(c(0, 1), c(.1, .2), level = c(.9, .95)), "level")
  absent <- f(c(0, 1), c(.1, .2), B = 0)
  expect_true(all(is.na(absent$ci_lower)))
  expect_true(all(absent$bootstrap_complete_B == 0L))
  expect_true(all(is.na(f(0, .1, B = 10)$ci_lower)))
  original <- qtl2mr:::qtl2mr_weighted_mode
  calls <- 0L
  local_mocked_bindings(qtl2mr_weighted_mode = function(...) {
    calls <<- calls + 1L
    if (calls == 3L) stop("controlled numerical failure")
    original(...)
  }, .package = "qtl2mr")
  failed <- f(c(0, .2, .4), c(.1, .2, .15), B = 10)
  expect_equal(failed$bootstrap_complete_B, c(10L, 9L))
  expect_equal(failed$bootstrap_failed_B, c(0L, 1L))
  expect_true(is.finite(failed$se[1L]))
  expect_true(all(is.na(failed[2L, c("se", "ci_lower", "ci_upper", "ci_type")])))
  expect_match(failed$note[2L], "1/10", fixed = TRUE)
})

test_that("captions describe only available normal-bootstrap intervals", {
  d <- data.frame(method = c("Weighted median", "Weighted mode"),
    estimate = c(.1, .2), ci_lower = c(-.2, -.3), ci_upper = c(.4, .7),
    ci_type = "bootstrap-SE Wald", ci_level = .95, source = "current_package_diagnostic")
  note <- qtl2mr:::qtl2mr_diagnostic_ci_caption(d)
  expect_match(note, "95%", fixed = TRUE)
  expect_match(note, "1.96 x parametric-bootstrap SD", fixed = TRUE)
  detail <- qtl2mr:::qtl2mr_diagnostic_ci_caption(d, detailed = TRUE)
  expect_match(detail, "Loci and inverse-variance weights are fixed", fixed = TRUE)
  expect_match(detail, "not profile-likelihood", fixed = TRUE)
  d$ci_level <- .9
  expect_match(qtl2mr:::qtl2mr_diagnostic_ci_caption(d), "1.645", fixed = TRUE)
  d$ci_type <- "locus bootstrap percentile"
  expect_identical(qtl2mr:::qtl2mr_diagnostic_ci_caption(d), "")
  d$ci_type <- "bootstrap-SE Wald"
  d$ci_lower <- NA_real_
  expect_identical(qtl2mr:::qtl2mr_diagnostic_ci_caption(d), "")
})

test_that("diagnostic levels do not relabel primary profiles in HTML", {
  r <- make_fake_qtl2mr_analysis_result()
  d <- fit_multilocus_diagnostic_layer(r$pilot$locus_blocks, global_fit = r$pilot$fit,
    bootstrap_B = 20L, bootstrap_seed = 1L, level = .9)
  expect_equal(d$estimates$ci_level, c(.95, .9, .9, .9))
  r$pilot$diagnostic_layer <- d
  r$method_comparison <- qtl2mr:::qtl2mr_method_comparison_table(r$pilot)
  html_file <- tempfile(fileext = ".html")
  qtl2mr:::write_qtl2mr_html_report(r, html_file)
  text <- paste(readLines(html_file), collapse = "\n")
  expect_match(text, "90% [", fixed = TRUE)
  expect_match(text, "95% [", fixed = TRUE)
  expect_false(grepl("<th>95% CI</th>", text, fixed = TRUE))
  expect_match(text, "1.645 x parametric-bootstrap SD", fixed = TRUE)
})

test_that("bootstrap SD uses safe scaling and unrepresentable bounds stay unavailable", {
  f <- qtl2mr:::qtl2mr_bootstrap_scalar_sensitivities
  large <- f(c(-1e154, 1e154), c(1e154, 1e154), B = 10, seed = 1)
  expect_true(all(is.finite(large$se)))
  expect_true(all(is.finite(large$ci_lower) & is.finite(large$ci_upper)))
  calls <- 0L
  local_mocked_bindings(qtl2mr_weighted_median = function(...) {
    calls <<- calls + 1L
    if (calls %% 2L) -1e308 else 1e308
  }, .package = "qtl2mr")
  huge <- f(c(-1, 1), c(.1, .1), B = 10, seed = 1)
  expect_true(is.na(huge$ci_type[1L]))
  expect_true(is.na(huge$se[1L]) && is.na(huge$ci_lower[1L]) && is.na(huge$ci_upper[1L]))
  expect_equal(huge$bootstrap_complete_B[1L], 10L)
  expect_equal(huge$bootstrap_failed_B[1L], 0L)
  expect_match(huge$note[1L], "non-finite numerical SD or bounds", fixed = TRUE)
})
