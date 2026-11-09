# Helpers added in 0.6.1.9006: combining two resources objects, and the
# simulation and bootstrap diagnostics of a finished single-pair analysis.

# ---------------------------------------------------------------------------
# qtl2mr_combine_resources()
# ---------------------------------------------------------------------------

qtl2mr_combine_ids <- function(resources) {
  list(genotyped = resources$sample_ids, pheno = rownames(resources$pheno),
    covar = if (is.null(resources$covar)) character() else rownames(resources$covar))
}

# The same values in two vectors: numbers compared as doubles (so integer and
# double storage of equal values agree), anything else as text; NA positions
# must match.
qtl2mr_same_values <- function(a, b) {
  if (length(a) != length(b)) return(FALSE)
  if (is.numeric(a) && is.numeric(b)) return(identical(as.numeric(a), as.numeric(b)))
  if (is.numeric(a) != is.numeric(b)) return(FALSE)
  identical(is.na(a), is.na(b)) && identical(as.character(a)[!is.na(a)], as.character(b)[!is.na(b)])
}

# The marker map columns both inputs must agree on, compared by value.
qtl2mr_combine_same_map <- function(mx, my) {
  cols <- intersect(c("marker", "chr", "pos"), names(mx))
  identical(cols, intersect(c("marker", "chr", "pos"), names(my))) &&
    all(vapply(cols, function(k) qtl2mr_same_values(mx[[k]], my[[k]]), logical(1)))
}

# allele_probs of x and y on one chromosome, samples of x first, then those only in y.
qtl2mr_combine_probs <- function(ax, ay, ids) {
  dn <- dimnames(ax)
  out <- array(NA_real_, dim = c(length(ids), dim(ax)[2L], dim(ax)[3L]),
    dimnames = list(ids, dn[[2L]], dn[[3L]]))
  in_x <- ids %in% dn[[1L]]
  out[in_x, , ] <- ax[ids[in_x], , , drop = FALSE]
  only_y <- !in_x
  if (any(only_y)) out[only_y, , ] <- ay[ids[only_y], , , drop = FALSE]
  out
}

# Rows of x, then the rows only in y; columns as in x.
qtl2mr_combine_rows <- function(dx, dy, ids) {
  only_y <- setdiff(ids, rownames(dx))
  out <- rbind(dx, dy[only_y, names(dx), drop = FALSE])
  out[ids[ids %in% rownames(out)], , drop = FALSE]
}

# Kinship of the combined samples taken from x and y when both cover all of them
# and agree exactly; NULL when they cannot be combined. Entries between a sample
# only in x and one only in y are known only when a kinship covers both, so a
# kinship of each cohort's own samples cannot be combined.
qtl2mr_combine_kinship <- function(kx, ky, ids) {
  as_list <- function(k) if (is.matrix(k) || is.data.frame(k)) list(`1` = as.matrix(k)) else lapply(k, as.matrix)
  lx <- as_list(kx)
  ly <- as_list(ky)
  if (!identical(names(lx), names(ly)) || (is.matrix(kx) || is.data.frame(kx)) != (is.matrix(ky) || is.data.frame(ky))) {
    stop("the kinships of x and y have different structures (one matrix vs a chromosome list, or different ",
      "chromosomes); supply kinship for the combined samples", call. = FALSE)
  }
  covers <- function(K) all(ids %in% rownames(K)) && all(ids %in% colnames(K))
  out <- lapply(names(lx), function(chr) {
    Kx <- lx[[chr]]
    Ky <- ly[[chr]]
    if (!covers(Kx) || !covers(Ky)) return(NULL)
    a <- Kx[ids, ids, drop = FALSE]
    if (!identical(as.numeric(a), as.numeric(Ky[ids, ids, drop = FALSE]))) return(NULL)
    a
  })
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  if (is.matrix(kx) || is.data.frame(kx)) return(out[[1L]])
  stats::setNames(out, names(lx))
}

qtl2mr_combine_resources <- function(x, y, kinship = NULL, common_only = FALSE) {
  if (!inherits(x, "qtl2mr_resources") || !inherits(y, "qtl2mr_resources")) {
    stop("x and y must both come from qtl2mr_prepare_resources()", call. = FALSE)
  }
  common_only <- qtl2mr_check_flag(common_only, "common_only")
  # 1. Markers, chromosomes, positions and founders.
  if (!qtl2mr_combine_same_map(x$marker_map, y$marker_map)) {
    stop("x and y must have the same marker map (marker, chr and pos values, in the same order)", call. = FALSE)
  }
  if (!identical(names(x$allele_probs), names(y$allele_probs))) {
    stop("x and y must have allele_probs for the same chromosomes, in the same order", call. = FALSE)
  }
  for (chr in names(x$allele_probs)) {
    dx <- dimnames(x$allele_probs[[chr]])
    dy <- dimnames(y$allele_probs[[chr]])
    # Each input must hold the same samples on every chromosome (as qtl2 genoprobs do).
    if (!identical(dx[[1L]], x$sample_ids) || !identical(dy[[1L]], y$sample_ids)) {
      stop("each input must have the same samples, in the same order, on every chromosome of allele_probs ",
        "(chromosome ", chr, " differs)", call. = FALSE)
    }
    if (!identical(dx[[2L]], dy[[2L]])) stop("x and y must have the same founders in the same order", call. = FALSE)
    if (!identical(dx[[3L]], dy[[3L]])) {
      stop("x and y must have the same markers in allele_probs on chromosome ", chr, call. = FALSE)
    }
  }
  # 2. Traits.
  clash <- intersect(names(x$pheno), names(y$pheno))
  if (length(clash)) {
    stop("x and y share trait name(s): ", paste(utils::head(clash, 5L), collapse = ", "),
      "; rename them before combining", call. = FALSE)
  }
  # 3. Covariates.
  if (is.null(x$covar) != is.null(y$covar)) {
    stop("x and y must both have covariates or both have none", call. = FALSE)
  }
  if (!is.null(x$covar) && !setequal(names(x$covar), names(y$covar))) {
    stop("x and y must have the same covariate columns", call. = FALSE)
  }
  ix <- qtl2mr_combine_ids(x)
  iy <- qtl2mr_combine_ids(y)
  # 4. A sample ID in both inputs is the same mouse: its probabilities and
  # covariates must be identical.
  shared_geno <- intersect(ix$genotyped, iy$genotyped)
  for (chr in names(x$allele_probs)) {
    if (length(shared_geno) && !identical(as.numeric(x$allele_probs[[chr]][shared_geno, , , drop = FALSE]),
                                          as.numeric(y$allele_probs[[chr]][shared_geno, , , drop = FALSE]))) {
      stop("samples in both x and y must have identical allele_probs (chromosome ", chr, "); ",
        "a shared sample ID means the same mouse", call. = FALSE)
    }
  }
  if (!is.null(x$covar)) {
    shared_covar <- intersect(ix$covar, iy$covar)
    a <- x$covar[shared_covar, , drop = FALSE]
    b <- y$covar[shared_covar, names(x$covar), drop = FALSE]
    same <- vapply(names(x$covar), function(k) qtl2mr_same_values(a[[k]], b[[k]]), logical(1))
    if (length(shared_covar) && !all(same)) {
      stop("samples in both x and y must have the same covariate values (column ", names(same)[!same][[1L]],
        " differs); a shared sample ID means the same mouse", call. = FALSE)
    }
  }
  # Samples: x first, then those only in y.
  geno <- union(ix$genotyped, iy$genotyped)
  pheno_ids <- union(ix$pheno, iy$pheno)
  covar_ids <- union(ix$covar, iy$covar)
  common <- shared_geno
  if (common_only) {
    if (!length(common)) stop("common_only = TRUE but x and y share no genotyped samples", call. = FALSE)
    geno <- common
    pheno_ids <- pheno_ids[pheno_ids %in% common]
    covar_ids <- covar_ids[covar_ids %in% common]
  }
  # The list keeps the attributes of x (for example the qtl2 class "calc_genoprob").
  allele_probs <- x$allele_probs
  for (chr in names(allele_probs)) {
    allele_probs[[chr]] <- qtl2mr_combine_probs(x$allele_probs[[chr]], y$allele_probs[[chr]], geno)
  }
  # 5. Traits of one input are NA for samples only in the other.
  pheno <- data.frame(row.names = pheno_ids)
  for (trait in names(x$pheno)) pheno[[trait]] <- x$pheno[match(pheno_ids, rownames(x$pheno)), trait]
  for (trait in names(y$pheno)) pheno[[trait]] <- y$pheno[match(pheno_ids, rownames(y$pheno)), trait]
  covar <- if (is.null(x$covar)) NULL else qtl2mr_combine_rows(x$covar, y$covar, covar_ids)
  # Kinship: supplied, or taken from x and y when they agree on all combined samples.
  kinship_source <- "supplied"
  if (is.null(kinship)) {
    if (is.null(x$kinship) && is.null(y$kinship)) {
      kinship_source <- "none"
    } else if (is.null(x$kinship) || is.null(y$kinship)) {
      stop("only one of x and y has a kinship; supply kinship for the combined samples", call. = FALSE)
    } else {
      kinship <- qtl2mr_combine_kinship(x$kinship, y$kinship, geno)
      if (is.null(kinship)) {
        stop("the kinships of x and y do not both cover all combined samples with identical values ",
          "(for example, two cohorts, or kinships of each input's own samples); supply kinship computed on the ",
          "combined probabilities, e.g. qtl2::calc_kinship(probs, type = \"loco\")", call. = FALSE)
      }
      kinship_source <- "x_and_y"
    }
  }
  drop_record <- function(m) m[setdiff(names(m), "probability_normalization")]
  record <- list(
    n_samples_x_only = length(setdiff(ix$genotyped, iy$genotyped)),
    n_samples_y_only = length(setdiff(iy$genotyped, ix$genotyped)),
    n_samples_overlap = length(shared_geno),
    n_samples_combined = length(geno),
    common_only = common_only,
    traits_x = names(x$pheno), traits_y = names(y$pheno),
    kinship_source = kinship_source
  )
  metadata <- list(combined_resources = record,
    x_metadata = drop_record(x$metadata %||% list()), y_metadata = drop_record(y$metadata %||% list()))
  qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar, marker_map = x$marker_map, metadata = metadata)
}

# ---------------------------------------------------------------------------
# qtl2mr_diagnostics()
# ---------------------------------------------------------------------------

qtl2mr_diagnostics_columns <- c(
  "exposure", "outcome", "n_loci", "nsim", "seed",
  "Q", "Q_chisq_df", "Q_p_chisq", "Q_p_sim", "Q_sim_n", "Q_sim_n_requested", "Q_sim_seed",
  "Q_sim_exceedances", "Q_sim_min_resolvable_p", "Q_sim_null_mean", "Q_sim_null_sd", "Q_sim_method",
  "theta_weighted_median", "se_weighted_median", "ci_lower_weighted_median", "ci_upper_weighted_median",
  "ci_type_weighted_median",
  "theta_weighted_mode", "se_weighted_mode", "ci_lower_weighted_mode", "ci_upper_weighted_mode",
  "ci_type_weighted_mode",
  "diagnostic_bootstrap_B", "diagnostic_bootstrap_seed")

# The result and the folder diagnostics.csv goes to: the folder given (or that
# of the given result_object.rds); for a result object, its recorded output_dir,
# but only when the result_object.rds there is the same run (otherwise the
# folder was overwritten by another analysis, or moved, and nothing is written).
qtl2mr_read_analysis_result <- function(result) {
  dir <- NULL
  if (is.character(result) && length(result) == 1L && !is.na(result)) {
    path <- if (dir.exists(result)) file.path(result, "result_object.rds") else result
    if (!file.exists(path)) stop("no result_object.rds found at ", result, call. = FALSE)
    dir <- dirname(normalizePath(path))
    result <- readRDS(path)
  } else if (inherits(result, "qtl2mr_analysis_result")) {
    recorded <- result$output_dir
    saved <- if (is.character(recorded) && length(recorded) == 1L && !is.na(recorded)) file.path(recorded, "result_object.rds") else NA_character_
    same_run <- !is.na(saved) && file.exists(saved) && tryCatch({
      s <- readRDS(saved)
      identical(s$created_at, result$created_at) && identical(s$exposure, result$exposure) &&
        identical(s$outcome, result$outcome)
    }, error = function(e) FALSE)
    if (same_run) {
      dir <- recorded
    } else {
      warning("the result's output folder no longer holds this run (moved, removed or overwritten); ",
        "diagnostics.csv is not written. Pass the result folder to write it there.", call. = FALSE)
    }
  }
  if (!inherits(result, "qtl2mr_analysis_result")) {
    stop("result must be a run_qtl2mr_analysis() result or its output directory", call. = FALSE)
  }
  if (!length(result$pilot$locus_blocks)) stop("result has no locus blocks", call. = FALSE)
  attr(result, "qtl2mr_diagnostics_dir") <- dir
  result
}

# Writes data/diagnostics.csv into a result folder: never through a link, and
# atomically (a temporary file in data/ renamed into place). It replaces an
# earlier diagnostics.csv; output_manifest.csv is not changed.
qtl2mr_write_diagnostics_csv <- function(table, dir) {
  data_dir <- file.path(dir, "data")
  if (!dir.exists(data_dir)) return(invisible(NULL))
  if (!identical(qtl2mr_assert_known_path_node(data_dir, "result data folder"), "directory")) {
    stop("the result's data folder is not a plain directory; refusing to write diagnostics.csv", call. = FALSE)
  }
  target <- file.path(data_dir, "diagnostics.csv")
  if (!qtl2mr_assert_known_path_node(target, "diagnostics.csv") %in% c("absent", "file")) {
    stop("data/diagnostics.csv is a link or not a regular file; refusing to write it", call. = FALSE)
  }
  tmp <- tempfile("diagnostics-", tmpdir = data_dir, fileext = ".csv.tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(table, tmp, row.names = FALSE)
  if (!file.rename(tmp, target)) stop("could not write data/diagnostics.csv", call. = FALSE)
  invisible(target)
}

qtl2mr_diagnostics <- function(result, nsim = 1000L, seed = 1L) {
  result <- qtl2mr_read_analysis_result(result)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim != round(nsim) || nsim < 2) {
    stop("nsim must be one whole number >= 2 (it is both the number of Q simulations and of bootstrap draws)",
      call. = FALSE)
  }
  nsim <- as.integer(nsim)
  seed <- qtl2mr_check_whole(seed, "seed")
  pilot <- result$pilot
  blocks <- pilot$locus_blocks
  # The multilocus stage of the analysis, refitted from its blocks with the
  # simulated Q p-value and the bootstrap: the same computation as an analysis
  # run with heterogeneity_method = "simulation" and diagnostic_bootstrap_B = nsim
  # (each seeded with seed) in 0.6.1.9005. The LMMs are not refitted, and the
  # main result is not changed.
  stage <- qtl2mr_multilocus_stage(blocks, ridge = 1e-10, heterogeneity_method = "simulation",
    q_sim_n = nsim, q_sim_seed = seed, diagnostic_bootstrap_B = nsim, diagnostic_bootstrap_seed = seed,
    reported_estimator_policy = "diagnostic_only")
  summary <- qtl2mr_multilocus_summary(blocks, stage, pilot$phenotype_scaling, pilot$locus_first_stage,
    exposure = result$exposure, outcome = result$outcome,
    n_loci_selected = pilot$summary$n_loci_selected[[1L]],
    n_loci_dropped = pilot$summary$n_loci_dropped_nonestimable[[1L]],
    standardize_phenotypes = pilot$summary$standardize_phenotypes[[1L]],
    cross_cov_model = pilot$summary$cross_cov_model[[1L]],
    diagnostic_bootstrap_B = nsim, diagnostic_bootstrap_seed = seed)$summary
  summary$nsim <- nsim
  summary$seed <- seed
  out <- summary[qtl2mr_diagnostics_columns]
  rownames(out) <- NULL
  dir <- attr(result, "qtl2mr_diagnostics_dir")
  if (is.character(dir) && length(dir) == 1L && !is.na(dir)) qtl2mr_write_diagnostics_csv(out, dir)
  out
}
