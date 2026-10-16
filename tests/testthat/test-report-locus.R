test_that("locus report writes Round 2 required files and handles bad locus IDs", {
  result <- make_fake_qtl2mr_analysis_result()
  main_dir <- tempfile("qtl2mr-main-")
  result <- write_qtl2mr_results(result, output_dir = main_dir, make_html = FALSE, make_figures = FALSE)
  loc_dir <- tempfile("qtl2mr-locus-")

  out <- report_qtl2mr_locus(result, "m1", output_dir = loc_dir, make_html = TRUE, make_figures = TRUE)

  expect_true(file.exists(file.path(out, "report.html")))
  expect_true(file.exists(file.path(out, "data", "locus_betas.csv")))
  expect_true(file.exists(file.path(out, "data", "locus_covariance_qc.csv")))
  expect_true(file.exists(file.path(out, "data", "founder_residuals.csv")))
  expect_true(file.exists(file.path(out, "data", "founder_coordinate_ratios.csv")))
  expect_true(file.exists(file.path(out, "data", "locus_summary.csv")))
  expect_true(file.exists(file.path(out, "figures", "founder_beta_comparison.png")))
  expect_true(file.exists(file.path(out, "figures", "founder_residual_scores.png")))
  expect_true(file.exists(file.path(out, "figures", "founder_coordinate_ratio_forest.png")))
  expect_true(file.exists(file.path(out, "figures", "covariance_psd_summary.png")))

  expect_error(report_qtl2mr_locus(result, "missing_marker", output_dir = tempfile()), "locus_id not found")
})
