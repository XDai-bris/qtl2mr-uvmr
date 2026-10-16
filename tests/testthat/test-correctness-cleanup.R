test_that("default analysis output stays below tempdir and explicit roots remain opt-in", {
  portable_path <- function(path) chartr("\\", "/", path)
  old_option <- getOption("qtl2mr.output_root")
  old_env <- Sys.getenv("QTL2MR_OUTPUT_ROOT", unset = NA_character_)
  old_wd <- getwd()
  artificial_cwd <- tempfile("qtl2mr-artificial-cwd-")
  dir.create(artificial_cwd)
  on.exit({
    options(qtl2mr.output_root = old_option)
    if (is.na(old_env)) Sys.unsetenv("QTL2MR_OUTPUT_ROOT") else Sys.setenv(QTL2MR_OUTPUT_ROOT = old_env)
    setwd(old_wd)
  }, add = TRUE)

  setwd(artificial_cwd)
  options(qtl2mr.output_root = NULL)
  Sys.unsetenv("QTL2MR_OUTPUT_ROOT")
  default_dir <- qtl2mr:::qtl2mr_analysis_output_dir("traitA", "traitB")
  expect_identical(
    normalizePath(qtl2mr:::qtl2mr_default_output_root(), winslash = "/"),
    normalizePath(tempdir(), winslash = "/")
  )
  expect_true(startsWith(
    portable_path(default_dir),
    paste0(portable_path(qtl2mr:::qtl2mr_default_output_root()), "/results/result_exp-")
  ))
  expect_false(startsWith(portable_path(default_dir), paste0(portable_path(artificial_cwd), "/")))
  expect_false(dir.exists(file.path(artificial_cwd, "outputs")))

  option_root <- tempfile("qtl2mr-option-root-")
  options(qtl2mr.output_root = option_root)
  option_dir <- qtl2mr:::qtl2mr_analysis_output_dir("traitA", "traitB")
  expect_true(startsWith(portable_path(option_dir), paste0(portable_path(option_root), "/results/result_exp-")))

  options(qtl2mr.output_root = NULL)
  env_root <- tempfile("qtl2mr-env-root-")
  Sys.setenv(QTL2MR_OUTPUT_ROOT = env_root)
  env_dir <- qtl2mr:::qtl2mr_analysis_output_dir("traitA", "traitB")
  expect_true(startsWith(portable_path(env_dir), paste0(portable_path(env_root), "/results/result_exp-")))

  explicit <- tempfile("qtl2mr-explicit-")
  expect_identical(
    portable_path(qtl2mr:::qtl2mr_analysis_output_dir("traitA", "traitB", explicit)),
    portable_path(normalizePath(explicit, winslash = "/", mustWork = FALSE))
  )
})

test_that("per-locus forest uses finite estimates when intervals are unavailable", {
  per_locus <- data.frame(
    marker = c("m1", "m2"), theta = c(0.2, 0.4),
    profile_ci_lower = c(NA_real_, NA_real_), profile_ci_upper = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  png <- tempfile(fileext = ".png")
  grDevices::png(png)
  expect_silent(plot_per_locus_so_forest(per_locus, global_theta = 0.3))
  grDevices::dev.off()
  expect_true(file.exists(png) && file.info(png)$size > 0)
  per_locus$theta[] <- NA_real_
  expect_error(plot_per_locus_so_forest(per_locus), "no finite theta")
})

test_that("Q simulation counts distinguish requested from completed simulations", {
  unavailable <- qtl2mr:::calibrate_multilocus_q(
    locus_blocks = list(), fit = list(theta = NA_real_, Q = NA_real_), nsim = 7L
  )
  expect_equal(unavailable$Q_sim_n_requested, 7L)
  expect_equal(unavailable$Q_sim_n, 0L)
  expect_false(is.finite(unavailable$Q_p_sim))
  expect_false(is.finite(unavailable$Q_sim_min_resolvable_p))
  not_run <- qtl2mr:::calibrate_multilocus_q(locus_blocks = list(), nsim = 0L)
  expect_equal(not_run$Q_sim_n_requested, 0L)
  expect_equal(not_run$Q_sim_n, 0L)
})

test_that("finite Q simulations use a plus-one empirical p-value", {
  block <- list(
    beta_x = stats::setNames(c(0.35, -0.2, 0.15), paste0("FC", 1:3)),
    beta_y = stats::setNames(c(0.29, -0.16, 0.12), paste0("FC", 1:3)),
    Sigma_x = diag(c(0.015, 0.012, 0.018)),
    Sigma_y = diag(c(0.02, 0.025, 0.022)),
    Sigma_xy = matrix(0, 3, 3), marker = "m1", chr = "1"
  )
  calibration <- qtl2mr:::calibrate_multilocus_q(
    list(m1 = block), fit = list(theta = 0.8, Q = 1e12), nsim = 5L, seed = 13L
  )
  expect_equal(calibration$Q_sim_n, 5L)
  expect_equal(calibration$Q_sim_exceedances, 0L)
  expect_equal(calibration$Q_p_sim, 1 / 6)
  expect_gt(calibration$Q_p_sim, 0)
})

test_that("simulation Q labels display the stored plus-one p-value", {
  zero <- qtl2mr:::qtl2mr_q_sim_label(1 / 6, exceedances = 0L, min_p = 1 / 6, digits = 6)
  expect_identical(zero, qtl2mr:::qtl2mr_fmt_num(1 / 6, digits = 6))
  expect_false(startsWith(zero, "<"))
  expect_identical(qtl2mr:::qtl2mr_q_sim_label(NA_real_), "NA")

  zero_row <- data.frame(
    Q_p_method = "simulation", Q_p_reported = 1 / 6,
    Q_sim_exceedances = 0L, Q_sim_min_resolvable_p = 1 / 6
  )
  ordinary_row <- data.frame(
    Q_p_method = "simulation", Q_p_reported = 3 / 6,
    Q_sim_exceedances = 2L, Q_sim_min_resolvable_p = 1 / 6
  )
  expect_identical(
    qtl2mr:::qtl2mr_q_reported_label(zero_row, digits = 6),
    qtl2mr:::qtl2mr_fmt_num(1 / 6, digits = 6)
  )
  expect_identical(
    qtl2mr:::qtl2mr_q_reported_label(ordinary_row, digits = 6),
    qtl2mr:::qtl2mr_fmt_num(3 / 6, digits = 6)
  )
})
