# qtl2mr

Founder-aware Mendelian randomization for Diversity Outbred and related
eight-founder mouse populations.

`qtl2mr` estimates the causal effect of an exposure trait on an outcome trait
using selected QTL as instruments. Each instrument is represented by founder
haplotype probabilities from `qtl2`, not by a single SNP dosage, and the
estimator accounts for exposure-outcome sample overlap.

**Current release:** [v0.6.0](https://github.com/XDai-bris/qtl2mr-uvmr/releases/tag/v0.6.0).
See [NEWS.md](NEWS.md) for what changed and the
[migration table](#upgrading-from-05x) if you used 0.5.x.

## Installation

Download `qtl2mr_0.6.0.tar.gz` and `SHA256SUMS` from the release assets. Use
the package file, not GitHub's automatic "Source code" archive, because the
package file contains the built vignettes. Check the SHA-256 against
`SHA256SUMS` (for example `shasum -a 256 -c SHA256SUMS` on macOS or
`sha256sum -c SHA256SUMS` on Linux), then in R:

```r
install.packages("fs")  # with conda: conda install -c conda-forge r-fs
install.packages("qtl2mr_0.6.0.tar.gz", repos = NULL, type = "source")
library(qtl2mr)
vignette("getting-started", package = "qtl2mr")
```

`qtl2` is not needed to run analyses. You will normally use it upstream to
compute allele probabilities, kinship and QTL scans:

```r
install.packages("qtl2", repos = "https://cloud.r-project.org")
```

## How to use qtl2mr

| Step | Use | You get |
|---|---|---|
| Choose instruments | `qtl2mr_prepare_resources()` then `qtl2mr_select_instruments()` | One instrument set per exposure, with the exposure block F of every LOD candidate |
| Study one exposure → outcome pair | `as_qtl2mr_prepared_input()` or `qtl2mr_pair_input()`, then `run_qtl2mr_analysis()` | A result folder with tables, figures, sensitivity analyses and an HTML report |
| Screen many pairs | `run_qtl2mr_batch()` | Compact tables with one row per pair; optional parallel workers and resume |

A batch row is the same fit as `run_qtl2mr_analysis()` on that pair with the
same instruments and settings, so pairs found in a batch can be studied in
detail without any change in the estimate.

### One pair

The package includes `qtl2mr_example`, a small synthetic dataset (96 mice,
eight founders, two selected loci).

```r
library(qtl2mr)
data("qtl2mr_example", package = "qtl2mr")

prepared <- as_qtl2mr_prepared_input(
  allele_probs = qtl2mr_example$allele_probs,
  pheno = qtl2mr_example$pheno,
  covar = qtl2mr_example$covar,
  kinship = qtl2mr_example$kinship,
  marker_map = qtl2mr_example$marker_map,
  selected_loci = qtl2mr_example$selected_loci,
  exposure = qtl2mr_example$exposure,
  outcome = qtl2mr_example$outcome,
  metadata = qtl2mr_example$metadata
)

result_dir <- file.path(tempdir(), "qtl2mr-example")
result <- run_qtl2mr_analysis(
  prepared,
  output_dir = result_dir,
  overwrite = TRUE,
  make_figures = FALSE,
  make_html = TRUE,
  verbose = FALSE,
  diagnostic_bootstrap_B = 0L,
  heterogeneity_method = "chisq",
  run_loo_mr = FALSE,
  run_instrument_sensitivity = FALSE,
  progress = FALSE
)

summary(result)                # one-row overview; qtl2mr_primary_result() has every column
qtl2mr_warnings(result)
browseURL(file.path(result$output_dir, "report.html"))
```

These settings skip the optional diagnostics (the Q simulation, bootstrap
intervals, leave-one-out refits and the instrument-set sensitivity, including
the Steiger direction check) so the example finishes in about a second; they
are not the package defaults. The reported estimate is about 0.892 SD of
outcome per SD of exposure, 95% CI [0.685, 1.126]. In original units
(`theta_raw_units`) it is about 0.93 (95% CI 0.72 to 1.18), which covers the
simulated `true_theta = 0.75`. The example is teaching data, not evidence about
the method.

### Many pairs

```r
# Your own inputs (see "Using your own data"): allele_probs, LOCO kinship,
# pheno, covar and map from qtl2, and scan = qtl2::scan1() output (markers x traits).
resources <- qtl2mr_prepare_resources(allele_probs = allele_probs, kinship = kinship,
                                      pheno = pheno, covar = covar, marker_map = map)
instruments <- qtl2mr_select_instruments(resources, scan, exposures = c("x", "w"),
                                         lod_thresh = 4)
pairs <- data.frame(exposure = c("x", "x", "w"), outcome = c("y", "z", "z"))

batch <- run_qtl2mr_batch(resources, pairs, instruments,
                          output_dir = "batch-results", n_workers = 4)
batch$main_results[, c("pair_id", "status", "theta_reported",
                       "reported_ci_lower", "reported_ci_upper")]

# A full report for one pair, with the same estimate as its batch row. With the
# package defaults (Q simulation, bootstrap) this takes minutes, not seconds.
detail <- run_qtl2mr_analysis(qtl2mr_pair_input(resources, instruments, "x", "z"),
                              output_dir = "x-to-z")
```

`qtl2mr_select_instruments()` takes, for each exposure, one peak per
contiguous above-threshold region per chromosome and fits the exposure at each
candidate to report its block F (`F_X_block`). `F_X_min` optionally drops weak
candidates; it is off by default. The table keeps every candidate with
`selected` and `drop_reason`, so you can inspect it, edit it or write your own
(only `exposure` and `marker` are required).

Each pair is fitted with the exposure and the outcome on their own complete
samples, exactly as in `run_qtl2mr_analysis()`. If a run is interrupted, call
it again with `resume = TRUE`: completed pairs are kept, the output folder may
have been moved, and `n_workers` may differ. A resume is refused if the
inputs, instruments, options or qtl2mr code changed (for example after an
upgrade). `n_workers > 1` uses forked processes on macOS and Linux; on Windows
the pairs run one after another.

Batch defaults use the chi-square Q p-value and skip bootstrap intervals and
sensitivity analyses; `batch_locus_results.csv` flags loci that fail the
Steiger direction check (`steiger_pass`). Rerun pairs of interest with
`run_qtl2mr_analysis()` for the full report.

## Using your own data

`vignette("from-qtl2", package = "qtl2mr")` walks from a `qtl2` cross to a
prepared input, including LOCO kinship and instrument selection. The inputs
are:

- `allele_probs`: founder **allele** probabilities, a chromosome-named list of
  sample × 8 founder × marker arrays (`qtl2::genoprob_to_alleleprob()`), not the
  36-state genotype probabilities;
- `pheno`, `covar`: numeric data frames with sample IDs as row names;
- `kinship`: a chromosome-named LOCO kinship list
  (`qtl2::calc_kinship(probs, type = "loco")`); `NULL` fits ordinary least
  squares;
- `marker_map`: a data frame with `marker`, `chr`, `pos`, or a qtl2 map list;
- `selected_loci`: the instrument markers. For a pair from an instruments
  table, use `qtl2mr_pair_input(resources, instruments, "x", "y")` instead;
  `as_qtl2mr_prepared_input()` uses `selected_loci` as given and refuses a
  table with unselected rows or other exposures. Choose the threshold before
  seeing MR results; the package never uses every marker unless you set
  `use_all_loci = TRUE`.

Only autosomal loci, continuous traits and additive covariates are supported.

## Reading a result

- **Estimate and units.** `theta_reported` with `reported_ci_lower` and
  `reported_ci_upper`. By default phenotypes are standardised, so the effect is
  SD of outcome per SD of exposure (`effect_units`); `theta_raw_units` gives
  the original-unit value when available. Use
  `standardize_phenotypes = FALSE` to fit in original units.
- **When the estimate is `NA`.** The second-order confidence set was unbounded
  or split into several pieces, usually because the instruments are weak. The
  pieces are in `data/confidence_set.csv`; check the block F statistics,
  consider stronger or fewer instruments, and do not substitute another
  estimator silently. Comparators (Two-step IVW, weighted median and mode) are
  in `qtl2mr_method_comparison(result)` as diagnostics.
- **`primary_reliable`.** A numerical check that the fit gave one bounded
  interval. It does not guarantee 95% coverage: with weak instruments
  (see [Instrument strength](#instrument-strength)) profile intervals can be
  too narrow.
- **Heritability.** `data/hsq.csv` lists the null-model heritability used for
  each trait and chromosome. It sets how much the kinship term reweights the
  LMM; a value of 0 means kinship had no effect there.
- **Warnings.** `qtl2mr_warnings(result)` gives each warning's consequence and
  recommended action.
- **Overview.** `summary(result)` gives a one-row overview; `confidence_set_bounded`
  is `FALSE` when the confidence set reaches plus or minus infinity.

## Instrument strength

`F_X_block` is reported for every locus. It is the joint Wald statistic of the
locus's seven founder contrasts on the exposure divided by its degrees of
freedom, F = β̂ᵀΣ̂⁻¹β̂ / 7. When the locus has no effect on the exposure, 7·F
is approximately chi-square with 7 degrees of freedom.

- **F = 2** gives 7·F = 14, almost exactly the 5% critical value of χ²₇
  (14.07): the founder effects on the exposure just reach nominal p ≈ 0.05.
  This is a minimal signal, not the mark of a strong instrument. The warning
  `weak_block_first_stage` lists how many retained loci have F_X < 2.
- **Relation to LOD.** For a scan on the eight-founder allele probabilities
  (7 degrees of freedom), a LOD score converts to a likelihood-ratio
  chi-square as 2·ln(10)·LOD ≈ 4.6·LOD, and in large samples this is close to
  the Wald statistic, so F ≈ 4.6·LOD / 7 ≈ 0.66·LOD (for strong loci the Wald
  statistic, and so F, is somewhat larger). A LOD threshold of 4 corresponds
  to F ≈ 2.6, and F = 2 to LOD ≈ 3.0. This holds when the scan and the MR fit
  use the same samples, covariates and kinship; otherwise it is only a rough
  guide.
- **Why not F > 10.** That rule comes from linear instrumental-variable
  regression with one instrument and a chosen tolerance for bias. It was not
  derived for a seven-degree-of-freedom founder block and is not calibrated
  here. Weak blocks can make profile intervals too narrow; use the LOD
  threshold and F to limit that risk, not as a coverage guarantee.

When the scan and the MR fit use the same samples, covariates and kinship, a
LOD threshold of 4 already implies F of about 2.6, so `F_X_min = 2` in
`qtl2mr_select_instruments()` rarely removes a locus. It matters when the scan
was run differently (for example without kinship, or on other samples) or when
a lower LOD threshold is used with F as the real gate.

## Runtime

With several loci, `run_qtl2mr_analysis()` calibrates the Q heterogeneity
p-value by simulation (`heterogeneity_nsim = 1000`), re-solving the multilocus
fit once per draw (the locus LMMs are not refitted). On a laptop, the bundled
400-mouse, 40-locus example (`qtl2mr_multilocus_example`) takes roughly
10–15 minutes with the default and under 2 minutes with
`heterogeneity_method = "chisq"`. Progress messages report completed draws.

## Advanced use

- **Non-estimable loci.** In `run_qtl2mr_analysis()` and `run_qtl2mr_batch()`,
  `on_nonestimable_locus = "drop"`
  removes a locus whose founder contrast cannot be estimated in the exposure or
  outcome samples and records it, instead of failing the pair.
- **Steiger direction check.** `run_qtl2mr_analysis()` repeats the fit without
  loci that explain more outcome than exposure variance and compares it with
  the primary estimate (`data/steiger_sensitivity.csv`). It runs with the
  instrument-set sensitivity (`run_instrument_sensitivity = TRUE`, the
  default) and is a sensitivity analysis; the primary estimate is not
  filtered.
- **Parallel batches.** With `n_workers > 1`, pairs sharing an exposure run in
  one process; fits are reused within an exposure group but not across groups.
  See `?run_qtl2mr_batch`.
- **Follow-up refits.** `refit_qtl2mr_subset()` and
  `fit_instrument_set_sensitivity()` re-estimate an analysis on locus subsets
  without refitting the LMMs.
- **Low-level interface.** `fit_founder_lmm_block()`, `fit_mr_locus_block()`,
  `fit_multilocus_founder_mr()` and `second_order_gls_ivw()` expose the
  estimator for method work (`?fit_mr_locus_block`).
- **qtl2 parity check.** `validate_qtl2_scan1coef()` compares the qtl2mr LMM
  with `qtl2::scan1coef()` for one marker.

## Upgrading from 0.5.x

| 0.5.x | 0.6.0 |
|---|---|
| `qtl2mr_prepare_resources(probs, K, pheno, covar, qtl_scan = scan, marker_map = map)` | `qtl2mr_prepare_resources(probs, K, pheno, covar, marker_map = map)`; the scan goes to `qtl2mr_select_instruments()` |
| `qtl2mr_select_loci_from_scan(scan, map, "x", lod_thresh = 4)` | `qtl2mr_select_instruments(resources, scan, exposures = "x", lod_thresh = 4)` |
| `run_qtl2mr_batch(r, p, instrument_rule = list(lod_thresh = 4))` | `run_qtl2mr_batch(r, p, qtl2mr_select_instruments(r, scan, exposures = unique(p$exposure), lod_thresh = 4))`; without `exposures` every trait in the scan is used |
| pair columns `instrument_trait`, `lod_thresh`, `exclude_chr`, `max_loci`, `selected_loci_file` | one instrument set per exposure in the `instruments` table; these columns are now an error |
| `sample_policy`, `min_sample_retention` | removed: every pair uses role-specific samples, as in `run_qtl2mr_analysis()` |
| `run_qtl2mr_batch(..., fit_options = list(...), force = TRUE)` | top-level arguments and `overwrite = TRUE`; `fit_options$q_sim_n`/`q_sim_seed` are `heterogeneity_nsim`/`heterogeneity_seed`; `engine`, `cache_contexts`, `cache_lmm_fits` and `verbose` are internal (low-level `fit_multilocus_founder_mr()` keeps `q_sim_n`) |
| `run_qtl2mr_analysis(..., standardize_phenotypes = , cross_cov_model = )` through `...` | named arguments; `...` is gone, and `ridge`, `engine` and the caches are internal. Arguments after `verbose` were reordered, so pass them by name |
| `batch$selected_loci`, `batch$first_stage_summary`, `batch$sample_qc` (and their CSV files) | `batch$locus_results` (`batch_locus_results.csv`); per-trait `batch$sample_exclusions` |
| batch `main_results` columns `instrument_trait`, `lod_thresh`, `exclude_chr` | removed; the instruments table records the selection |
| `theta_so`, `se_so`, `ci_lower_so`, `ci_upper_so`, `sample_policy` columns | `theta`, `se`, `profile_ci_lower`, `profile_ci_upper`; no sample-policy column |
| `fit_instrument_set_sensitivity()` default `weak_f_threshold = 10`, set `exclude_F_lt_10` | default 2, set `exclude_F_lt_2`; pass `weak_f_threshold = 10` for the old set |

Batch estimates change for pairs whose exposure and outcome have different
missing samples: 0.5.x fitted both traits on their common samples, 0.6.0 fits
each on its own complete samples. Pairs without such differences give the same
numbers. The batch Q p-value now defaults to chi-square. Batch folders from
0.5.x cannot be resumed. For older versions, see [NEWS.md](NEWS.md).

## Method in brief

At each selected locus, exposure and outcome are regressed on the founder
probabilities with a linear mixed model whose polygenic term uses the
leave-one-chromosome-out kinship, weighted by a null-model heritability that
is then held fixed. The locus estimates are combined by a second-order
objective that includes uncertainty in the exposure effects and the
exposure-outcome covariance created by shared samples. The confidence interval
is the profile set of that objective at a nominal chi-square cutoff.

Known limitations: weak instruments can give poor coverage even when the fit
succeeds; the objective does not model covariance between loci; selecting
instruments in the same samples can bias the estimate; and heritability is
treated as known. A methods manuscript is in preparation.

## Documentation

- `vignette("getting-started", package = "qtl2mr")`
- `vignette("from-qtl2", package = "qtl2mr")`
- `vignette("multilocus-example", package = "qtl2mr")`
- `?run_qtl2mr_analysis`, `?run_qtl2mr_batch`, `help(package = "qtl2mr")`
- [NEWS.md](NEWS.md)

## Citation

Until the methods paper is available, cite this repository and the package
version you used (`citation("qtl2mr")`, `packageVersion("qtl2mr")`).
