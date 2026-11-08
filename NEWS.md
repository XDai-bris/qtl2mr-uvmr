# qtl2mr 0.6.1.9004 (development)

Not released. Instrument subsets in one run and an F report of an instrument
table, on top of 0.6.1.9003; with the default `subsets = NULL` every result,
table, file and warning is the same as in 0.6.1.9003.

- `subsets` in `run_qtl2mr_batch()` and `run_qtl2mr_analysis()`: a named list
  of rule lists, for example
  `list(F1.5 = list(F_X_min = 1.5), F1.5_SF = list(F_X_min = 1.5, steiger = TRUE))`.
  `F_X_min` keeps loci with a finite exposure block F at least the value;
  `steiger = TRUE` keeps loci with `steiger_pass = TRUE` (FALSE and NA are
  dropped). Rules within a subset intersect. Names must be unique, use only
  letters, digits, `.` and `_`, and not be `all`; unknown rules, empty rule
  lists and invalid values are refused before anything is fitted.
- Each pair is fitted as before on its instruments (the maximal set), and the
  main row is unchanged. For each subset, the pair's fitted loci that pass the
  rules, judged on that pair's own `F_X_block` and `steiger_pass`, are refitted
  in the multilocus combination with the run's settings (heterogeneity method,
  simulation draws and seed, diagnostic bootstrap, reported estimator policy,
  cross-covariance model); the LMMs are not refitted. A subset row therefore
  equals the main row of a separate run on that subset's loci. The main fit and
  the subset refits share one internal function for the multilocus stage and
  its summary row.
- A subset whose loci equal those of the maximal set or of an earlier subset
  copies that row (`same_loci_as`). Fewer loci than `min_loci` (one in the
  single analysis) gives status `too_few_loci`; a failing refit gives status
  `error` and, in the batch, an `errors` row with stage `"subset"`, without
  stopping the run or changing the pair's main row; pairs without a fit get
  `pair_not_fitted`. A failure elsewhere in a pair's subset step (selecting the
  loci or assembling the rows) gives every subset of that pair status `error`
  in the same way, and the main row, the locus rows and the other outputs stay
  as without subsets.
- Warnings raised by the estimator itself (for example the warning that
  simulation was selected with `heterogeneity_nsim = 0`) are raised again by
  each subset refit, so they can appear several times on the console; they do
  not enter any table.
- Batch: new table `subset_results` (`batch_subset_results.csv`), one row per
  pair and subset with `pair_id`, `exposure`, `outcome`, `subset`, `rules`,
  `n_loci`, `status`, `same_loci_as`, `error_message` and the estimate and test
  columns of `main_results`; `locus_results` gains `in_subset_<name>` (NA for
  loci that were not fitted); `config$subsets` and the resume fingerprint record
  the specification, and subset rows are restored on resume. The locus blocks
  of `write_locus_blocks` hold the maximal set.
- Single analysis: `result$subset_results` and `data/subset_results.csv`, with
  the columns of `main_result.csv`. `instrument_sensitivity.csv`,
  `steiger_sensitivity.csv` and the default sensitivity sets are unchanged.
- New exported `qtl2mr_instrument_summary(instruments, F_breaks = c(1.5, 2, 5, 10))`:
  the F report of a `qtl2mr_select_instruments()` table, one row per exposure
  with `n_candidates`, the count for each drop reason, `n_selected`, the
  minimum, quartiles, median, mean and maximum of `F_X_block` over the selected
  loci, `n_F_ge_<break>` over all fitted candidates, and the selection
  settings. No existing function or output changes.

# qtl2mr 0.6.1.9003 (development)

Not released. An opt-in founder-design QC on top of 0.6.1.9002; with the
default `min_founder_copies = NULL` every result, table, file and warning is
the same as in 0.6.1.9002.

- `min_founder_copies` in `qtl2mr_select_instruments()`, `run_qtl2mr_batch()`
  and `qtl2mr_pair_input()`: a locus is dropped when the rarest founder's
  expected allele count, E = min over founders of the sum over samples of
  2 p, is below the threshold. E is computed on the samples the trait's LMM is
  fitted on, from the probabilities as fitted (after the probability policy's
  normalisation).
- In selection it applies to the exposure on its own complete samples, after
  `non_autosomal` and before fitting: a failing candidate is not fitted and
  gets `drop_reason = "founder_qc"`; `max_loci` counts only loci that pass. The
  table gains `E_min_x` and `min_founder_copies`.
- In the batch it applies to each pair's outcome, on the samples the outcome
  is fitted with, before fitting: a failing locus is dropped for that pair only
  (`locus_results`: `kept = FALSE`, `drop_reason = "founder_qc_outcome"`,
  `E_min_y`; `main_results`: `n_loci_dropped_founder_qc`), so it cannot stop the
  pair as non-estimable. `qtl2mr_pair_input()` applies the same QC, so a batch
  row equals the single analysis. When the QC is on, these columns are always
  present (0 or NA where no locus was checked).
- The batch applies only the outcome-side rule. The exposure side is applied by
  `qtl2mr_select_instruments()`, so an instrument table selected without the QC
  is not re-checked for the exposure.
- Because `max_loci` counts only loci that pass, a QC drop can bring in a
  lower-LOD candidate that the QC-off selection left out (`max_loci`), and
  `locus_rank` is renumbered over the passing loci. With `max_loci` unset, the
  QC-on selection is the QC-off selection minus the QC drops (apart from
  candidates that fail both the QC and estimation, now labelled `founder_qc`).

# qtl2mr 0.6.1.9002 (development)

Not released. A speed build: at the same settings and worker count, every
result, table and warning is the same as in 0.6.1.9001, and
`run_qtl2mr_analysis()` with `n_workers > 1` gives the same as with one
process. Defaults are unchanged; three arguments are added. Fits now stay
cached in the R session, within a fixed memory bound (see below).

- `qtl2mr_prepare_resources()` checks the founder probabilities once instead of
  twice, and the sum-to-one check no longer calls `sum()` once per sample and
  marker; the checks, tolerances and the normalization record are unchanged.
  Validation of a kinship matrix that already passed is not repeated, and when
  rounding is normalized only the markers in use are copied.
- Each fit of the second-order estimator reuses the objective values of its
  search grid instead of evaluating them again; heterogeneity-simulation draws
  skip the parts of the fit that they do not report.
- `run_qtl2mr_analysis(n_workers = 1L)`: on Unix, the refits of the
  heterogeneity simulation, the diagnostic bootstrap, leave-one-locus-out and
  the instrument-set sensitivity can run in several processes. Random numbers
  are drawn in the calling process in the same order, so results do not depend
  on `n_workers`. Under `options(warn = 2)` or higher the refits run in the
  calling process; a worker that is killed (for example out of memory) stops
  the analysis.
- LMM contexts and locus fits are reused across the R session: fits made by
  `qtl2mr_select_instruments()` serve the MR stage, and single analyses and
  batches share them. They are keyed by trait, chromosome, sample set,
  phenotype and covariate values, kinship and marker probabilities, so only
  identical fits are shared. The kinship decomposition is computed once per
  chromosome and sample set, and the overlap products once per pair. The cache
  holds at most 2 GiB as measured by `object.size()`; the least recently used
  entries are released first, and a running batch keeps those it is using.
  Small session memos (kinship checks and fingerprints, up to eight kinship
  decompositions) are outside this bound.
- `run_qtl2mr_batch()` first fits, trait by trait (in parallel by trait with
  `n_workers > 1`), every locus its pairs need, then fits the pairs as before:
  in order in the calling process for a single exposure, and in parallel
  exposure-group workers for several exposures, as in 0.6.1.9001 (warnings
  raised inside those workers are still not returned). Fits that no later pair
  needs become the first to be released and otherwise stay in the session for
  later calls.
- `run_qtl2mr_batch(write_locus_blocks = FALSE)`: with `TRUE`, the per-locus
  founder blocks of each pair are written to `batch_locus_blocks/` with a
  manifest; `fit_multilocus_founder_mr()` refits a pair from them.
- `sample_policy = c("own_complete", "pair_common_complete")` in
  `run_qtl2mr_batch()` and `qtl2mr_pair_input()`: the default keeps each trait
  on its own complete samples; `"pair_common_complete"` fits both traits on the
  pair's common complete samples, as the 0.5.0 batch did.

# qtl2mr 0.6.1.9001 (development)

Not released. A selection-only build: instrument selection changes for
ordinary scans; MR estimates for a given instrument set are unchanged.

- `qtl2mr_select_instruments()` now finds LOD peaks with the qtl2 peak rule
  (`peak_rule = "peakdrop"`, the new default): on each chromosome, a second
  peak needs the LOD to fall by more than `peakdrop` (default 1.5) below the
  current peak and rise again by more than `peakdrop`, as in
  `qtl2::find_peaks()`. Each peak is represented by its maximum-LOD marker. A
  marker exactly at `lod_thresh` still counts as above it. Earlier versions
  took the maximum of every contiguous run with LOD ≥ `lod_thresh`, so a QTL
  whose LOD hovered around the threshold was split into several nearby
  candidates with nearly identical founder effects, and a separate peak inside
  a long above-threshold run was missed.
- Two other rules are available: `peak_rule = "spacing"` keeps the old
  candidates in decreasing LOD and drops any within `min_spacing` (default
  1 Mb, in marker-map units) of one already kept, as a sensitivity analysis;
  `peak_rule = "legacy"` gives the 0.6.1 selection, for reproducing earlier
  results only. No rule is LD pruning.
- The instrument table records the rule: `selected_by`, `peak_rule`,
  `peakdrop`, `min_spacing` and `rule_version`.
- Markers with a missing or non-finite LOD are skipped by every rule instead of
  splitting a peak region. A peak region that reaches the first or last marker
  of a chromosome is handled by every rule.

# qtl2mr 0.6.1

A patch release from a post-release audit. Estimates for ordinary inputs are
unchanged from 0.6.0; no arguments or outputs were added or removed.

- Results no longer depend on the units of the traits when
  `standardize_phenotypes = FALSE`. Eigenvalue floors in matrix inverses and
  square roots, the overlap residual-variance check and the REML
  residual-sum-of-squares check were absolute and are now relative to the
  data. Previously, traits with very small values (variances around 1e-16)
  could change block F and the Steiger flags, or stop with "Residual variance
  on the overlap is zero". The block F reported by
  `qtl2mr_select_instruments()`, which standardizes internally, now also
  matches raw-unit analyses.
- Allele probabilities must hold the same samples on every chromosome, as
  qtl2 output does; other inputs are refused with a direct message. Before,
  a sample missing on an unused chromosome made the batch sample check
  undercount, and a sample missing on a used chromosome stopped the fit with
  "subscript out of bounds".
- A Steiger refit that fails inside the direction check now raises the
  `steiger_sensitivity_error` warning and writes
  `logs/steiger_sensitivity_error.log`, as documented; previously only the
  `failed` status in the table recorded it.

# qtl2mr 0.6.0

This release separates instrument selection from MR fitting and makes the batch
workflow the batch form of the single-pair analysis: for the same pair,
instruments and settings, every batch row equals `run_qtl2mr_analysis()`. The
single-pair estimates are unchanged from 0.5.0; batch estimates change only for
pairs whose exposure and outcome have different missing samples (see below).

## Three main functions

- New `qtl2mr_select_instruments(resources, qtl_scan, exposures, lod_thresh,
  exclude_chr, max_loci, F_X_min, n_workers, progress)` builds one instrument
  set per exposure (`exposures = NULL` means every trait in the scan). It keeps every LOD candidate with its exposure block F
  (`F_X_block`, fitted with the same samples, covariates, LOCO kinship and
  null-model hsq as the MR stage, so it equals the block F of every pair using
  that exposure, up to rounding when pairs are fitted with
  `standardize_phenotypes = FALSE`), `selected` and `drop_reason`
  (`non_autosomal`, `nonestimable_x`, `x_fit_error`, `F`, `max_loci`). `F_X_min` is off by
  default. A scan without trait names can be used for only one exposure.
- New `qtl2mr_pair_input(resources, instruments, exposure, outcome)` gives the
  single-pair input for any pair of a batch.
- `run_qtl2mr_batch(resources, pairs, instruments, ...)` takes the instrument
  table (a data frame or CSV path; only `exposure` and `marker` are required,
  and only rows with `selected = TRUE` are used when that column exists). Each
  pair is fitted through the same input and settings as the single analysis.
- `as_qtl2mr_prepared_input()` uses `selected_loci` as given and refuses an
  instruments table with unselected rows or rows of several or another
  exposure; use `qtl2mr_pair_input()` for such tables.

## Breaking changes and migration

- Removed without deprecation: `qtl2mr_select_loci_from_scan()`,
  `qtl2mr_select_loci_for_traits()`, the `qtl_scan` argument of
  `qtl2mr_prepare_resources()`, and the batch arguments `instrument_rule`,
  `sample_policy`, `min_sample_retention` and `fit_options`. Pair columns
  `instrument_trait`, `lod_thresh`, `exclude_chr`, `max_loci` and
  `selected_loci_file` are now an error rather than silently ignored.
- Batch samples: every pair now fits the exposure and the outcome on their own
  complete samples, as the single analysis always did. 0.5.0 fitted both on
  the samples complete for both traits, so pairs with different missingness in
  exposure and outcome get new (and now identical-to-single) estimates; other
  pairs are unchanged. `batch_sample_qc.csv` is gone and
  `batch_sample_exclusions.csv` has one row per trait and reason, now also
  counting samples missing from the genotype probabilities or kinship.
- Both entry points share flat argument names: `on_nonestimable_locus`
  (now also in `run_qtl2mr_analysis()`, default `"fail"`),
  `reported_estimator_policy`, `heterogeneity_method`, `heterogeneity_nsim`
  and `heterogeneity_seed` (formerly `q_sim_n`/`q_sim_seed` in batch),
  `diagnostic_bootstrap_B`/`_seed`, `standardize_phenotypes` and
  `cross_cov_model`. `run_qtl2mr_analysis()` no longer takes `...`
  (`ridge`, `engine` and the caches are internal). In batch, `force` is now
  `overwrite`, and `min_samples`/`min_loci` are top-level arguments.
- The batch Q p-value defaults to chi-square (`heterogeneity_method = "chisq"`)
  instead of none.
- Main results: the duplicate `theta_so`, `se_so`, `ci_lower_so` and
  `ci_upper_so` columns and the `sample_policy` column are removed (use
  `theta`, `se`, `profile_ci_lower`, `profile_ci_upper`);
  `n_loci_steiger_fail` is added. Batch `main_results` drops
  `instrument_trait`, `lod_thresh` and `exclude_chr`. Apart from `pair_id`,
  `status`, `error_message`, `locus_fit_status` and `runtime_sec`, its columns
  are a subset of the single `main_result`, with equal values for the same
  settings.
- `fit_instrument_set_sensitivity()` now defaults to `weak_f_threshold = 2`,
  so the standard set `exclude_F_lt_10` is now `exclude_F_lt_2`; pass
  `weak_f_threshold = 10` for the old set.
- Arguments of `run_qtl2mr_analysis()` after `verbose` were reordered to match
  `run_qtl2mr_batch()`; pass them by name.
- Batch warnings now use the single-analysis catalog (for example
  `prepared_overlap_engine` instead of `batch_shared_samples`), and the batch
  run manifest has a `status` column (`complete` or `stopped_fail_fast`) with
  the same layout in both cases.
- Batch tables: `batch_selected_loci.csv` and `batch_first_stage_summary.csv`
  are merged into `batch_locus_results.csv` (one row per pair and instrument:
  kept, drop reason, exposure and outcome block F, partial R2,
  `steiger_pass`, sample counts and hsq). Every batch table has its documented
  columns even when it has no rows (empty `overlap_qc` and
  `method_comparison` tables previously had different columns).
- Batch folders from 0.5.0 cannot be resumed.

## Steiger direction check

- `run_qtl2mr_analysis()` (with `run_instrument_sensitivity = TRUE`, the
  default) repeats the multilocus fit without loci whose partial R2 for the
  outcome exceeds that for the exposure. The result is a `steiger_filtered` row
  in `instrument_sensitivity.csv`, a comparison with the primary estimate in
  `data/steiger_sensitivity.csv` and the report (including the shift in
  primary-CI widths), and the notice `steiger_estimate_outside_ci` when the
  filtered estimate leaves the primary 95% CI. It does not filter the primary
  estimate. `block_first_stage.csv` gains `F_Y_block`, partial R2 and
  `steiger_pass`. A failed check is logged and raises
  `steiger_sensitivity_error`.

## Fixes

- LOD islands are now formed in marker-map order when the map has no
  positions; previously the scan's row order was used, so a scan whose rows
  were not in map order could merge or split islands.
- Counts and seeds (`heterogeneity_nsim`, `diagnostic_bootstrap_B` and the
  seeds) are checked when a run starts and stored as integers, so `1000` and
  `1000L` resume the same batch.
- A missing `pairs` CSV path now reports that the file does not exist.
- A single-locus analysis no longer fails while drawing the block-F dashboard
  figure (present since 0.4.x).
- `run_qtl2mr_batch()` refuses an output folder whose package-owned tables,
  checkpoint or pair records are symbolic links, before writing anything.
  Previously a resume or overwrite could write through such a link to a file
  outside `output_dir`.
- The weak-instrument guide is now block F_X < 2 instead of < 10. Block F_X is
  the joint Wald statistic of the seven founder contrasts divided by 7, so
  F_X = 2 is nominal p of about 0.05; the F > 10 rule was derived for single
  instruments and is not calibrated for founder blocks. This changes when
  `weak_block_first_stage` is raised, the `threshold` in
  `first_stage_summary`, the dashboard guide line and the standard
  instrument set, now `exclude_F_lt_2`. The warning now reads "k of m
  retained loci have block F_X < 2". The README explains the scale and its
  relation to LOD. The per-founder-coordinate ratio flag (one degree of
  freedom) keeps its F < 10 rule.
- The HTML report describes how hsq was obtained for the sources actually
  used (REML null model, user supplied, or no kinship) instead of always
  stating REML.
- Continuous integration also runs on pull requests, and the optional qtl2
  job now runs the from-qtl2 onboarding test.

# qtl2mr 0.5.0

This release slims the package to two workflows and folds in the unreleased
0.4.1 fixes. For the same inputs and options, reported estimates, standard
errors, profile confidence intervals and Q statistics are unchanged from 0.4.x
(checked to be bitwise identical on the bundled and synthetic batch fixtures).

## Breaking changes and migration

- `run_qtl2mr_batch()` is the single batch entry point. New arguments:
  `sample_policy` (`"pair_common_complete"`, default, or
  `"batch_common_complete"`), `on_nonestimable_locus` (`"fail"` or `"drop"`),
  `min_sample_retention` (batch-common only, default 0.85), and `n_workers`
  (parallel by exposure). The `sample_rule` and `output_mode` arguments and
  `fit_options$drop_nonestimable_loci` were removed.
- Removed `run_qtl2mr_batch_fast()`, `run_qtl2mr_batch_exact_grouped()` and
  `run_qtl2mr_batch_parallel_by_exposure()`. Each is an option combination of
  `run_qtl2mr_batch()`, which gives the same numbers; see the migration table
  in the [0.5.0 README](https://github.com/XDai-bris/qtl2mr-uvmr/blob/v0.5.0/README.md#upgrading-from-04x).
- Removed the `run_qtl2mr_analysis()` arguments `adapter`,
  `run_qtl2_validation`, `run_heteroskedasticity` and `heteroskedasticity_B`,
  which had no effect. `output_dir` is now the fourth argument.
- `summary()` of an analysis result returns a one-row overview (`exposure`,
  `outcome`, `n_samples`, `n_loci`, the reported `estimate`, `se`, `ci_lower`,
  `ci_upper`, `effect_units`, `reported_estimator`, `primary_reliable`,
  `confidence_set_type`, `confidence_set_bounded`, `Q`, `Q_p`, `Q_p_method`,
  `sample_policy` and the IDs of warnings shown in the report) instead of the
  full 99-column table; `qtl2mr_primary_result()` still returns every column.
- New column `confidence_set_bounded` in single-analysis and batch results,
  appended after `sample_policy`: TRUE when every piece of the second-order
  confidence set is finite, FALSE when any piece reaches minus or plus infinity
  (including the whole real line), NA when the set is unresolved or the pair
  failed.
- The package and method name is written `qtl2mr` throughout; result labels
  change accordingly (for example `"Star-space qtl2mr"` in
  `reported_estimator` and the method-comparison table).
- Result folders no longer include `data/validation_tests.csv` (a package test
  inventory, not run evidence). `logs/pilot_messages.log` is now
  `logs/fit_messages.log`; the lite export writes `README_LITE.md` and
  `LITE_EXPORT_AUDIT.csv`.
- No longer exported, but still used internally (their tables and figures are
  still written by `run_qtl2mr_analysis()` or `run_qtl2mr_batch()`):
  `loo_multilocus_founder_mr()`, `founder_residual_scores()`,
  `founder_coordinate_ratios()`, `fit_multilocus_diagnostic_layer()`,
  `qtl2mr_make_common_complete_samples()`, `write_qtl2mr_batch_outputs()` and
  the plotting helpers `plot_method_comparison()`, `plot_per_locus_so_forest()`,
  `plot_founder_residual_heatmap()`, `plot_founder_ratio_forest()`,
  `plot_loo_multilocus_mr()`, `plot_instrument_set_sensitivity()` and
  `plot_block_first_stage_dashboard()`.
- Removed: `plan_qtl2mr_batch_exact_groups()`,
  `qtl2mr_subset_allele_probs_to_loci()`, `qtl2mr_read_mr_block_cache()`,
  `refit_qtl2mr_mr_block_cache()`, `refit_qtl2mr_mr_block_cache_grid()` and
  `fit_xprs_heteroskedasticity_diagnostic()` (with its print/plot methods),
  together with the xPRS/score heteroskedasticity diagnostics, the
  founder-residual null simulation and other code that no workflow used.
- Batch output layout: completed pairs are stored in `batch_pair_records/`; the
  exact-group planning tables and the per-exposure `by_exposure/` folders are
  no longer written. `batch_main_results.csv` drops the redundant
  `sample_rule` column and ends with `sample_policy` and
  `confidence_set_bounded`. 0.4.x batch folders cannot be resumed.
- Removed, without a deprecation period, the batch `fit_options` documented in
  0.4.x for saving per-pair block caches (`save_mr_block_cache`,
  `mr_block_cache_level`, `mr_block_cache_compress`), the post-hoc F-threshold
  refit grid (`run_first_stage_filter_grid`, `first_stage_filter_grid`,
  `first_stage_filter_min_loci`) and the Steiger-style filter
  (`steiger_filter`); they now give "Unknown option(s)". Also removed:
  `batch_mr_block_cache_manifest.csv`, `batch_first_stage_filter_results.csv`
  and `mr_block_cache/`. Per-locus block F statistics remain in
  `batch_first_stage_summary.csv`; F-threshold subsets of one analysis can be
  refitted with `fit_instrument_set_sensitivity()`.
- `batch_run_manifest.csv` has one layout for all batch runs: columns that
  belonged to the removed entry points or caches (for example `output_mode`,
  `sample_mode`, `sample_approximation`, `parallel_strategy`,
  `execution_strategy` and the block-cache and filter-grid columns) are gone,
  and `n_pairs_resumed` is added. `batch_config.json` records
  `sample_policy`, `on_nonestimable_locus` and `min_sample_retention` instead
  of `sample_rule`, and no cache settings.
- Single-analysis outputs no longer include the always-empty
  `qtl2_validation*.csv` and `heteroskedasticity_qc.csv`. Use
  `validate_qtl2_scan1coef()` for the optional qtl2 comparison.

## Batch resume

- Resume continues from per-pair records: an interrupted run (including a
  parallel one) continues where it stopped, the output folder may be moved,
  and `n_workers` may change. Changed inputs, options or package code, or
  modified records or tables, are refused before anything is written. Pairs
  that failed are refitted on resume.
- Invalid pairs become error rows under every sample policy; previously a
  missing trait aborted a whole `batch_common_complete` run. The batch-common
  sample set and `batch_sample_exclusions.csv` use only the valid pairs'
  traits.
- With `fail_fast = TRUE` the completed pairs and the failing pair are written
  before the error. If a parallel worker dies, the pairs it had finished are
  kept and only its unfinished pairs are reported as errors.
- The resume fingerprint covers founder probabilities only at the markers the
  run uses, so starting a batch no longer serialises genome-wide probability
  arrays; only package-written tables are checksummed.
- The batch-common retention gate is checked before any output is written.
- Exact-group planning was removed: shared fits were already reused through
  the fit cache, and planning saved about 5% of runtime in benchmarks, within
  run-to-run variation.

## Reporting and checks

- Null-model heritability is reported per role and chromosome
  (`data/hsq.csv`, summary line, HTML section, and `hsq_x`/`hsq_y` with status
  in the batch first-stage table), with notices for boundary estimates and
  kinship-free OLS fits. How hsq is estimated is unchanged.
- Identical exposure and outcome, an empty selected-locus set and malformed
  trait names fail at validation with direct messages.
- Founder probabilities may exceed [0, 1] by floating-point error up to 1e-8,
  as in 0.4.0; larger excursions report the observed range.
- The single-analysis covariate notice also counts samples without a
  covariate row.
- Batch warnings describe the shared batch sample set instead of the
  single-analysis role-specific note.
- Progress shows five real stages instead of seven.
- Continuous integration runs on every push (pull requests were added in
  0.6.0).

## Fixes carried over from the unreleased 0.4.1

- Corrected chromosome subsetting in the qtl2 onboarding recipe and added a
  runnable `cross2`-to-prepared-input example with regression coverage.
- Reports distinguish the Q statistic, Q p-value and reference distribution;
  Q simulation reports completed draws and elapsed time. The default method
  and draw count are unchanged.
- Founder-probability rounding has an explicit policy:
  `probability_policy = "normalize_rounding"` (default) normalizes small
  rounding errors in an analysis-local copy and records them; `"strict"`
  rejects them. Previously such inputs were rejected.
- User-supplied `hsq` values are validated and used. Named vectors (for
  example `qtl2::est_herit()` output) and chromosome-keyed lists were
  previously ignored and hsq was re-estimated, so analyses that pass hsq this
  way can now give different estimates. Conflicting entries are an error.
- Weighted-mode bandwidth no longer depends on measurement units; the
  diagnostic weighted-mode estimate can change slightly where the old absolute
  bandwidth floor applied.
- Degenerate per-locus diagnostics return unavailable values with reasons
  instead of failing.
- Covariate exclusions and sample policies are recorded in the outputs
  (`sample_alignment.csv`, `batch_sample_exclusions.csv`, `sample_policy`).
- Added the full GPL-3 licence text and declared graphics, grDevices and
  parallel imports. Runs without HTML no longer list an absent report.

# qtl2mr 0.4.0

- Added a compact synthetic 400-sample, 40-locus dataset and a worked
  multi-locus vignette with kinship, diagnostics and HTML reporting.
  The original two-locus quick start remains available.
- Weighted median/mode now use fixed-locus, fixed-weight normal parametric
  draws and bootstrap-SE Wald intervals (500 draws by default), replacing
  locus-resampling percentile intervals. Figures, HTML and tables identify
  these approximate diagnostic intervals separately from the primary profile CI.
- Weighted-mode peak finding uses local grids and numerical refinement;
  extreme effects no longer spread a fixed global grid too thinly. The
  kernel/bandwidth definition is unchanged; point estimates can change slightly
  from the former grid approximation. Failed bootstrap draws suppress the
  affected interval and are recorded, not silently removed.
- Fixed figure links in transactionally written HTML reports and saved
  result objects so they point to published files, not removed staging paths.
- Sensitivity plots reserve room for long instrument-set labels.
- Scale-normalized second-order fitting and whole-real-line candidate search
  replace unit-dependent local bracketing. Confidence sets retain unbounded
  and disconnected components; numerical status is separate from calibration.
- Instrument selection, sample alignment, unknown options and batch resume
  now fail explicitly instead of changing analyses silently.
- Result reports state phenotype units, numerical limitations and the full
  confidence-set representation. Synthetic-only onboarding includes a
  nonzero-kinship LMM and a qtl2-to-prepared-input recipe.
- Fixed synthetic regression inputs make the whole-line optimizer tests
  reproducible across BLAS/LAPACK builds without relaxing their assertions.
- First public source release from a new repository history. Release assets
  include the source package with all three prebuilt HTML vignettes and a
  separate SHA-256 checksum file.

# qtl2mr 0.3.0

- Kinship decomposition and null-model ML/REML heritability estimation now use
  one internal base-R LMM backend in direct and batch workflows. `qtl2` is an
  optional dependency used only for explicitly requested numerical parity
  checks; its absence does not prevent primary analysis.

# qtl2mr 0.2.3

- Instrument scans now require exact `instrument_trait` matching and positive,
  finite LOD thresholds; batch selection no longer falls back silently to a
  first scan column.
- Default results use a safe temporary root while explicit roots remain
  supported; parallel failures and marker-fit errors retain structured
  provenance.
- Improved vignette-test, forest-plot, and Q-simulation robustness, including
  plus-one empirical p-values and requested-versus-completed simulation counts.

# qtl2mr 0.2.2

- Finite-ridge second-order MR and diagnostics now use basis-invariant
  trace/mean-diagonal scaling across founder-contrast rotations, while keeping
  raw sampling covariance distinct from the ridged working covariance where
  appropriate.
- Result-directory and lite-ZIP publication is transactional and fail-closed,
  with link/junction, nested-parent, canonical-ancestor, rollback/cleanup,
  overwrite/preflight, and alias/inside-result protections.
- Added coverage for dense/permutation/sign rotations, adversarial covariance
  scales, output races/failures, Windows junction/path behavior, and
  ancestor-link return identity.

# qtl2mr 0.2.1

- Added a buildable, self-contained getting-started vignette and the small
  deterministic `qtl2mr_example` synthetic prepared-input dataset for
  installation and integration teaching. These are not method calibration or
  biological/application evidence.
- Expanded pre-install onboarding guidance, including private-repository PAT
  handling, `v0.2.1` release/tag guidance, optional `build_vignettes = TRUE`
  installation for the formal vignette, and tested conda/qtl2 setup notes.

# qtl2mr 0.2.0

- **Compatibility-visible default:** high-level reported estimates now use
  `diagnostic_only`. A reliable star-space second-order estimate remains
  reported; an unreliable SO fit reports `NA` estimate/SE/CI while retaining
  raw SO and Two-step diagnostic fields. The previous reported fallback is
  available only with explicit `reported_estimator_policy = "auto_fallback"`.
- Confirmatory, canonical, and exact-grouped workflows now fail closed for a
  non-estimable locus rather than silently changing the instrument set. The
  explicit fast/screening path may drop such loci only with structured locus,
  warning, and run-manifest audit fields.
- Added generic prepared-input, batch, exact-grouped, fast, and
  exposure-parallel APIs with compact result tables and structured manifests.
- Added regression coverage for reported-estimator policy, fail-closed exact
  behavior, audited screening drops, and deterministic public unreliable-SO
  behavior.

# qtl2mr 0.1.0

- Initial local research package.
- qtl2-compatible full-covariance founder LMM fitter.
- Star-space overlap covariance and PSD validation.
- Overlap-aware second-order Founder GLS-IVW.
- Multi-locus diagnostic layer, LOO diagnostics, and instrument-set sensitivity.
