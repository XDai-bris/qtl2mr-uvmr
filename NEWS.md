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
