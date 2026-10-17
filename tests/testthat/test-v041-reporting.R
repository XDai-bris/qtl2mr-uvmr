v041_reporting_input <- function() {
  data("qtl2mr_example", package = "qtl2mr", envir = environment())
  d <- qtl2mr_example
  as_qtl2mr_prepared_input(d$allele_probs, d$pheno, covar = d$covar,
    kinship = d$kinship, marker_map = d$marker_map,
    selected_loci = d$selected_loci, exposure = d$exposure, outcome = d$outcome)
}

v041_reporting_run <- function(input = v041_reporting_input(), make_html = TRUE,
                              progress = FALSE, method = "chisq", nsim = 0L) {
  run_qtl2mr_analysis(input, output_dir = tempfile("qtl2mr-v041-report-"),
    make_html = make_html, make_figures = FALSE, verbose = FALSE,
    diagnostic_bootstrap_B = 0L, heterogeneity_method = method,
    heterogeneity_nsim = nsim, run_loo_mr = FALSE,
    run_instrument_sensitivity = FALSE, progress = progress)
}

test_that("reports distinguish the Q statistic, p-value and reference", {
  result <- v041_reporting_run()
  s <- qtl2mr_primary_result(result)
  txt <- paste(readLines(file.path(result$output_dir, "summary.md")), collapse = "\n")
  html <- paste(readLines(file.path(result$output_dir, "report.html")), collapse = "\n")
  for (label in c("Q statistic", "Q p-value", "Q reference", "Sample policy")) {
    expect_match(txt, label, fixed = TRUE)
    expect_match(html, label, fixed = TRUE)
  }
  expect_false(grepl("Reported Q", html, fixed = TRUE))
  expect_match(txt, paste0("Q statistic: ", qtl2mr:::qtl2mr_fmt_num(s$Q, 6)), fixed = TRUE)
  expect_match(txt, paste0("Q p-value: ", qtl2mr:::qtl2mr_fmt_num(s$Q_p_reported, 6)), fixed = TRUE)
  expect_match(html, "chi-square; df=", fixed = TRUE)
  expect_identical(s$sample_policy, "role_specific_complete")
  expect_false("diagnostic_unavailable" %in% qtl2mr_warnings(result)$warning_id)
})

test_that("disabled HTML is absent from the advertised output paths", {
  result <- v041_reporting_run(make_html = FALSE)
  expect_false(file.exists(file.path(result$output_dir, "report.html")))
  expect_false("report_html" %in% qtl2mr_output_paths(result)$item)
  expect_false("report.html" %in% qtl2mr_output_paths(result, detail = TRUE)$path_relative)
  expect_false(any(grepl("Outputs: report.html", summarise_qtl2mr_result(result), fixed = TRUE)))
})

test_that("Q progress preserves random draws and can be suppressed without losing logs", {
  input <- v041_reporting_input()
  base <- qtl2mr:::run_prepared_qtl2mr(input, exposure = "x", outcome = "y",
    write_outputs = FALSE, verbose = FALSE, diagnostic_bootstrap_B = 0L,
    heterogeneity_method = "chisq")
  events <- character()
  callback <- function(stage, label, detail) events <<- c(events, detail)
  a <- qtl2mr:::calibrate_multilocus_q(base$locus_blocks, base$fit,
    nsim = 3L, seed = 41L, keep_simulated = TRUE)
  set.seed(54)
  state <- .Random.seed
  b <- qtl2mr:::calibrate_multilocus_q(base$locus_blocks, base$fit,
    nsim = 3L, seed = 41L, keep_simulated = TRUE, progress_callback = callback)
  expect_identical(a, b)
  expect_identical(.Random.seed, state)
  observer <- function(stage, label, detail) invisible(runif(4L))
  observed <- qtl2mr:::calibrate_multilocus_q(base$locus_blocks, base$fit,
    nsim = 3L, seed = 41L, keep_simulated = TRUE, progress_callback = observer)
  expect_identical(a, observed)
  expect_identical(.Random.seed, state)
  expect_match(events[1L], "0/3 draws", fixed = TRUE)
  expect_match(tail(events, 1L), "3/3 draws", fixed = TRUE)

  log <- tempfile("qtl2mr-v041-progress-", fileext = ".log")
  reporter <- qtl2mr:::qtl2mr_make_progress_reporter(TRUE, log_file = log)
  expect_message(qtl2mr:::qtl2mr_capture_messages(reporter(2L, "Q", "1/3 draws")), "1/3 draws")
  expect_message(suppressMessages(reporter(2L, "Q", "2/3 draws")), NA)
  expect_true(any(grepl("2/3 draws", readLines(log), fixed = TRUE)))
  integrated <- v041_reporting_run(progress = "log", method = "simulation", nsim = 2L)
  expect_true(any(grepl("2/2 draws", readLines(file.path(integrated$output_dir,
    "logs", "progress.log")), fixed = TRUE)))
  # Defaults and the primary objective have not changed to obtain faster examples.
  expect_identical(formals(run_qtl2mr_analysis)$heterogeneity_nsim, 1000L)
  expect_identical(eval(formals(run_qtl2mr_analysis)$heterogeneity_method)[1L], "auto")
  expect_equal(integrated$pilot$fit$theta, base$fit$theta, tolerance = 1e-12)
  expect_equal(integrated$pilot$fit$profile_ci_lower, base$fit$profile_ci_lower, tolerance = 1e-12)
})

test_that("input exclusions and rounding adjustments reach the report", {
  input <- v041_reporting_input()
  input$covar[1L, 1L] <- NA_real_
  expect_warning(result <- v041_reporting_run(input), "covariate")
  expect_true("incomplete_covariates" %in% qtl2mr_warnings(result)$warning_id)
  expect_true(any(result$pilot$sample_alignment$n_incomplete_covar > 0L))
  html <- paste(readLines(file.path(result$output_dir, "report.html")), collapse = "\n")
  expect_match(html, "incomplete_covariates", fixed = TRUE)

  input <- v041_reporting_input()
  input$allele_probs <- lapply(input$allele_probs, round, digits = 4L)
  expect_warning(rounded <- v041_reporting_run(input), "probabilit")
  expect_true("probability_rounding_normalized" %in% qtl2mr_warnings(rounded)$warning_id)
  ledger <- read.csv(file.path(rounded$output_dir, "data", "probability_normalization.csv"))
  expect_gt(nrow(ledger), 0L)
  expect_true(all(ledger$n_normalized > 0L))
})

test_that("low-level OLS paths do not ignore supplied heritability", {
  y <- setNames(seq_len(20L), paste0("s", seq_len(20L)))
  X <- cbind(Intercept = 1, x = seq_len(20L))
  for (h in c(-.1, 1.1, NA_real_, .4)) {
    expect_error(qtl2mr:::qtl2mr_lmm_transform(y, X, hsq = h), "hsq")
    expect_error(qtl2mr:::prepare_lmm_context(y, hsq = h), "hsq")
  }
  expect_identical(qtl2mr:::qtl2mr_lmm_transform(y, X, hsq = 0)$lmm_backend,
    "ordinary_least_squares")
})
