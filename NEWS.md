# qtl2mr 0.5.0 (unreleased)

This release slims the package to two workflows and folds in the unreleased
0.4.1 fixes. For the same inputs and options, reported estimates, standard
errors, profile confidence intervals and Q statistics are unchanged from 0.4.x
(checked to be bitwise identical on the bundled and synthetic batch fixtures).

## Breaking changes and migration

- `run_qtl2mr_batch()` is the single batch entry point. New arguments:
  `sample_policy` (`"pair_common_complete"`, default, or
  `"batch_common_complete"`), `on_nonestimable_locus` (`"fail"` or `"drop"`),
  `min_sample_retention` (batch-common only), and `n_workers` (parallel by
  exposure). `sample_rule`, `output_mode` and
  `fit_options$drop_nonestimable_loci` are deprecated and mapped with a
  warning.
- `run_qtl2mr_batch_fast()`, `run_qtl2mr_batch_exact_grouped()` and
  `run_qtl2mr_batch_parallel_by_exposure()` are deprecated thin wrappers that
  give the same numbers as before; they will be removed in 0.6.0. See the
  migration table in the README or `?run_qtl2mr_batch_fast`.
- `run_qtl2mr_analysis()` no longer has `adapter`, `run_qtl2_validation`,
  `run_heteroskedasticity` or `heteroskedasticity_B`; they had no effect and
  are now ignored with a warning when supplied.
- No longer exported (internal or removed): `plan_qtl2mr_batch_exact_groups()`,
  `qtl2mr_subset_allele_probs_to_loci()`,
  `qtl2mr_make_common_complete_samples()`, `write_qtl2mr_batch_outputs()`,
  the block-cache refit functions, `loo_multilocus_founder_mr()`,
  `founder_residual_scores()`, `founder_coordinate_ratios()`,
  `fit_multilocus_diagnostic_layer()`, `fit_xprs_heteroskedasticity_diagnostic()`
  and the `plot_*()` helpers. Their tables and figures are still produced by
  `run_qtl2mr_analysis()`.
- Removed the xPRS/score heteroskedasticity diagnostics, the founder-residual
  null simulation and other code that no workflow used (about 1,700 lines).
- Batch output layout: completed pairs are stored in `batch_pair_records/`; the
  exact-group planning tables and the per-exposure `by_exposure/` folders are
  no longer written. `batch_main_results.csv` drops the redundant
  `sample_rule` column and ends with `sample_policy`. 0.4.x batch folders
  cannot be resumed.
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
  missing trait aborted a whole `batch_common_complete` run.
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

# qtl2MR 0.3.0

- Kinship decomposition and null-model ML/REML heritability estimation now use
  one internal base-R LMM backend in direct and batch workflows. `qtl2` is an
  optional dependency used only for explicitly requested numerical parity
  checks; its absence does not prevent primary analysis.

# qtl2MR 0.2.3

- Instrument scans now require exact `instrument_trait` matching and positive,
  finite LOD thresholds; batch selection no longer falls back silently to a
  first scan column.
- Default results use a safe temporary root while explicit roots remain
  supported; parallel failures and marker-fit errors retain structured
  provenance.
- Improved vignette-test, forest-plot, and Q-simulation robustness, including
  plus-one empirical p-values and requested-versus-completed simulation counts.

# qtl2MR 0.2.2

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

# qtl2MR 0.2.1

- Added a buildable, self-contained getting-started vignette and the small
  deterministic `qtl2mr_example` synthetic prepared-input dataset for
  installation and integration teaching. These are not method calibration or
  biological/application evidence.
- Expanded pre-install onboarding guidance, including private-repository PAT
  handling, `v0.2.1` release/tag guidance, optional `build_vignettes = TRUE`
  installation for the formal vignette, and tested conda/qtl2 setup notes.

# qtl2MR 0.2.0

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

# qtl2MR 0.1.0

- Initial local research package.
- qtl2-compatible full-covariance founder LMM fitter.
- Star-space overlap covariance and PSD validation.
- Overlap-aware second-order Founder GLS-IVW.
- Multi-locus diagnostic layer, LOO diagnostics, and instrument-set sensitivity.
