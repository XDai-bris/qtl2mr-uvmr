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
  in the README.
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
- Continuous integration runs on every push and pull request.

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
