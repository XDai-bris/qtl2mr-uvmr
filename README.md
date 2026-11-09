# qtl2mr

Founder-aware Mendelian randomization for Diversity Outbred and related
eight-founder mouse populations.

`qtl2mr` estimates the causal effect of an exposure trait on an outcome trait
using selected QTL as instruments. Each instrument is represented by founder
haplotype probabilities from `qtl2`, not by a single SNP dosage, and the
estimator accounts for exposure-outcome sample overlap.

**Current release:** [v0.6.1](https://github.com/XDai-bris/qtl2mr-uvmr/releases/tag/v0.6.1).
See [NEWS.md](NEWS.md) for what changed,
[Changes since 0.6.1](#changes-since-061) for the development version in this
repository, and the [migration table](#upgrading-from-05x) if you used 0.5.x.

## Installation

Download `qtl2mr_0.6.1.tar.gz` and `SHA256SUMS` from the release assets. Use
the package file, not GitHub's automatic "Source code" archive, because the
package file contains the built vignettes. Check the SHA-256 against
`SHA256SUMS` (for example `shasum -a 256 -c SHA256SUMS` on macOS or
`sha256sum -c SHA256SUMS` on Linux), then in R:

```r
install.packages("fs")  # with conda: conda install -c conda-forge r-fs
install.packages("qtl2mr_0.6.1.tar.gz", repos = NULL, type = "source")
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
| Study one exposure → outcome pair | `qtl2mr_pair_input()`, then `run_qtl2mr_analysis()` | A result folder with tables, figures, sensitivity analyses and an HTML report |
| Screen many pairs | `run_qtl2mr_batch()` | Compact tables with one row per pair; optional parallel workers and resume |

Two helpers complete these steps: `qtl2mr_combine_resources()` combines two
data sets, and `qtl2mr_diagnostics()` adds the slower simulation diagnostics to
a single-pair result. A batch row is the same fit as `run_qtl2mr_analysis()` on
that pair with the same instruments, so pairs found in a batch can be studied
in detail without any change in the estimate.

This README describes the development version in this repository
(0.6.1.9006). Its arguments differ from the 0.6.1 release; see
[Changes since 0.6.1](#changes-since-061).

### One pair

The package includes `qtl2mr_example`, a small synthetic dataset (96 mice,
eight founders, two selected loci).

```r
library(qtl2mr)
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example

resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno, ex$covar, ex$marker_map)
input <- qtl2mr_pair_input(resources, data.frame(exposure = "x", marker = c("m1", "m2")), "x", "y")
result <- run_qtl2mr_analysis(input, output_dir = tempfile())

summary(result)                # one-row overview; qtl2mr_primary_result() has every column
qtl2mr_warnings(result)
browseURL(file.path(result$output_dir, "report.html"))
```

`qtl2mr_prepare_resources()` holds the inputs shared by every pair, and
`qtl2mr_pair_input()` adds one pair and its instruments; a table with
`exposure` and `marker` columns is enough. The analysis prints its stages and
a summary (`progress = FALSE` turns this off) and finishes in about a second,
including the leave-one-out refits and the instrument-set sensitivity with the
Steiger direction check. The reported estimate is about 0.892 SD of outcome
per SD of exposure, 95% CI [0.685, 1.126]. In original units
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

# A full report for one pair, with the same estimate as its batch row, then
# its simulated Q p-value and bootstrap standard errors.
detail <- run_qtl2mr_analysis(qtl2mr_pair_input(resources, instruments, "x", "z"),
                              output_dir = "x-to-z")
qtl2mr_diagnostics(detail)
```

`qtl2mr_select_instruments()` finds the exposure's LOD peaks on each
chromosome with the qtl2 peak rule (LOD ≥ `lod_thresh`, default 4; a new peak
needs a drop of more than `peakdrop`, default 1.5, as in `qtl2::find_peaks()`)
and fits the exposure at each peak to report its block F (`F_X_block`).
`peak_rule = "spacing"` is a sensitivity rule: it keeps the maximum of each
contiguous island of markers with LOD ≥ `lod_thresh` (the selection of qtl2mr
0.6.1 and earlier) and then drops any candidate within `min_spacing` Mb of a
higher one. Every estimable peak that passes the founder QC is selected; there
is no F filter (use `subsets`, below). The table keeps every candidate with
`selected` and `drop_reason`, so you can inspect it, edit it or write your own
(only `exposure` and `marker` are required).

Each pair is fitted with the exposure and the outcome on their own complete
samples, exactly as in `run_qtl2mr_analysis()`; to fit both traits on one
common sample set, see [Combining two data sets](#combining-two-data-sets).
With `write_locus_blocks = TRUE` the batch also writes each pair's per-locus
founder blocks, which `fit_multilocus_founder_mr()` can refit. If a run is
interrupted, call it again with `resume = TRUE`: completed pairs are kept, the
output folder may have been moved, and `n_workers` may differ. A resume is
refused if the inputs, instruments, options or qtl2mr code changed (for example
after an upgrade). `n_workers > 1` uses forked processes on macOS and Linux;
on Windows the pairs run one after another.

The batch and the single analysis use the same fixed fit settings: the Q
p-value refers to chi-square and the weighted median and mode have no
bootstrap intervals. The batch skips the leave-one-out and sensitivity refits;
`batch_locus_results.csv` flags loci that fail the Steiger direction check
(`steiger_pass`). Rerun pairs of interest with `run_qtl2mr_analysis()` for the
full report.

## Founder QC

A locus where one founder is almost absent from the fitted mice has a
near-singular founder design: its founder effects are poorly estimated or not
estimable. For each locus qtl2mr computes E, the expected number of copies of
the rarest founder allele, E = min over founders of Σᵢ 2·pᵢₖ over the fitted
samples, and drops the locus when E is below `min_founder_copies` (default 5;
`NULL` turns the QC off).

- **Exposure side.** `qtl2mr_select_instruments()` checks each candidate
  before fitting it (`drop_reason = "founder_qc"`, column `E_min_x`).
- **Outcome side.** `qtl2mr_pair_input()` and `run_qtl2mr_batch()` check each
  instrument of a pair on the outcome's samples (`drop_reason =
  "founder_qc_outcome"`, counted in `n_loci_dropped_founder_qc`), so a locus
  can be dropped for one pair and kept for another. A hand-written instrument
  table is not checked on the exposure side.
- **Order.** Loci are dropped as non-autosomal first, then by the founder QC,
  then as non-estimable.

In one Diversity Outbred cohort of about 800 mice, with about 3,000 candidate
loci across 80 traits, a threshold of 5 removed under 1% of the exposure loci
(mostly loci whose design was near-singular because a founder was almost
absent) and removed an outcome locus in about 2% of the pairs. Every marker of
the two bundled examples passes.

## Combining two data sets

`qtl2mr_combine_resources(x, y, kinship = NULL, common_only = FALSE)` combines
two resources objects, for example two data sets measured on overlapping mice,
into one that every function above accepts.

```r
both <- qtl2mr_combine_resources(resources_a, resources_b)
both$metadata$combined_resources   # samples only in a, only in b, shared, kept
```

- The marker map (marker, chr, pos) and the founders must be identical, trait
  names must not clash, and the covariate columns must be the same.
- A sample ID in both inputs means the same mouse, so its probabilities and
  covariates must be identical. Every failed check is an error, never a silent
  repair.
- A mouse present in only one input is `NA` for the other input's traits, so
  each trait uses its own samples: an exposure from one data set and an outcome
  from the other give a two-sample MR with partial or no overlap.
- Kinship is taken from the inputs when both cover all combined mice with
  identical values (one cohort). For two separately genotyped cohorts, supply
  `kinship`, for example `qtl2::calc_kinship(probs, "loco")` on the combined
  probabilities; qtl2mr does not compute kinship.
- `common_only = TRUE` keeps only the mice present in both inputs. Each trait
  still uses its own complete samples within that set.
- No cohort covariate is added. Whether one is needed is your scientific
  decision; add it to `covar` before preparing the resources.

## Simulated Q p-value and bootstrap intervals

`run_qtl2mr_analysis()` refers the Q heterogeneity statistic to chi-square and
reports the weighted median and mode comparators as point estimates.
`qtl2mr_diagnostics(result, nsim = 1000L, seed = 1L)` adds, after the run, the
Q p-value calibrated by conditional residual simulation (`Q_p_sim`) and
bootstrap standard errors and intervals for the weighted median and mode. It
works from the locus blocks stored in the result (an object or a result
folder), without refitting the LMMs; `nsim` is both the number of simulations
and of bootstrap draws, and `seed` seeds both. It returns one row, writes
`data/diagnostics.csv` into the result folder and does not change the main
result. At the same `nsim` and seed its numbers equal those of 0.6.1.9005 with
`heterogeneity_method = "simulation"` and `diagnostic_bootstrap_B = nsim`. For
a batch pair, rerun the pair with `qtl2mr_pair_input()` and
`run_qtl2mr_analysis()` first.

## Using your own data

`vignette("from-qtl2", package = "qtl2mr")` walks from a `qtl2` cross to a
prepared input, including LOCO kinship and instrument selection. The inputs
are:

- `allele_probs`: founder **allele** probabilities, a chromosome-named list of
  sample × 8 founder × marker arrays (`qtl2::genoprob_to_alleleprob()`), not the
  36-state genotype probabilities; every chromosome must hold the same samples,
  as qtl2 output does. Sums that are off by rounding are normalised with a
  warning and a record;
- `pheno`, `covar`: numeric data frames with sample IDs as row names;
- `kinship`: a chromosome-named LOCO kinship list
  (`qtl2::calc_kinship(probs, type = "loco")`); `NULL` fits ordinary least
  squares;
- `marker_map`: a data frame with `marker`, `chr`, `pos`, or a qtl2 map list;
- instruments: a table with `exposure` and `marker` columns, from
  `qtl2mr_select_instruments()` or your own prespecified list. Choose the
  threshold before seeing MR results.

Only autosomal loci, continuous traits and additive covariates are supported.

## Reading a result

- **Estimate and units.** `theta_reported` with `reported_ci_lower` and
  `reported_ci_upper`. Phenotypes are standardised, so the effect is SD of
  outcome per SD of exposure (`effect_units`); `theta_raw_units`,
  `se_raw_units` and `profile_ci_lower_raw_units`/`profile_ci_upper_raw_units`
  give original units when `raw_conversion_available` is `TRUE`.
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
  each trait and chromosome, always estimated by REML. It sets how much the
  kinship term reweights the LMM; a value of 0 means kinship had no effect
  there.
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
LOD threshold of 4 already implies F of about 2.6, so a subset
`list(F_X_min = 2)` rarely removes a locus. It matters when the scan was run
differently (for example without kinship, or on other samples) or when a lower
LOD threshold is used with F as the real gate. `qtl2mr_instrument_summary()`
counts, per exposure, the selected loci at or above several F values.

## Runtime

`run_qtl2mr_analysis()` always runs the leave-one-locus-out and instrument-set
refits, which re-solve the multilocus fit without refitting the locus LMMs. On
a laptop, the bundled 400-mouse, 40-locus example
(`qtl2mr_multilocus_example`) takes about half a minute with figures and the
HTML report. `qtl2mr_diagnostics()` is slower: each draw re-solves the
multilocus fit, about 0.3 seconds per draw for that example, so the default
1,000 draws take about 5 minutes.

## Advanced use

- **Instrument subsets.** `subsets = list(F2 = list(F_X_min = 2), F2_SF =
  list(F_X_min = 2, steiger = TRUE))` in `run_qtl2mr_analysis()` or
  `run_qtl2mr_batch()` refits each pair on the fitted loci that pass each
  rule set, without refitting the LMMs, and adds `subset_results`; the main
  result is unchanged. See `?run_qtl2mr_batch`.
- **Non-estimable loci.** A selected locus whose founder contrast cannot be
  estimated in the exposure or outcome samples stops `run_qtl2mr_analysis()`
  and gives an error row in `run_qtl2mr_batch()`. Remove it from the instrument
  table. Selection already drops such exposure loci, and the founder QC drops
  outcome loci where a founder is almost absent.
- **Steiger direction check.** `run_qtl2mr_analysis()` repeats the fit without
  loci that explain more outcome than exposure variance and compares it with
  the primary estimate (`data/steiger_sensitivity.csv`). It is a sensitivity
  analysis; the primary estimate is not filtered.
- **Parallel runs and reuse.** In a batch, each trait's locus fits are computed
  once and shared by all its pairs, and within an R session fits made by
  `qtl2mr_select_instruments()` and earlier analyses are reused; results never
  depend on this. Cached fits stay in the session within a memory budget
  (option `qtl2mr.cache_budget`, default 4 GiB). With `n_workers > 1` the batch
  fits traits in parallel processes and, with several exposures, runs exposure
  groups in parallel workers; `run_qtl2mr_analysis(n_workers = )` runs its
  leave-one-out and sensitivity refits in parallel. Results are the same. See
  `?run_qtl2mr_batch`.
- **Follow-up refits.** `refit_qtl2mr_subset()` and
  `fit_instrument_set_sensitivity()` re-estimate an analysis on locus subsets
  without refitting the LMMs.
- **Low-level interface.** `fit_founder_lmm_block()`, `fit_mr_locus_block()`,
  `fit_multilocus_founder_mr()` and `second_order_gls_ivw()` expose the
  estimator for method work (`?fit_mr_locus_block`).
- **qtl2 parity check.** `validate_qtl2_scan1coef()` compares the qtl2mr LMM
  with `qtl2::scan1coef()` for one marker.

## Changes since 0.6.1

The development version 0.6.1.9006 simplifies the main functions: each
statistical choice has one setting, and the arguments below are gone. Calls
that pass a removed argument now fail with an "unused argument" error
(`peak_rule = "legacy"` fails in the check of `peak_rule`).

| 0.6.1 (and 0.6.1.900x) | 0.6.1.9006 |
|---|---|
| `as_qtl2mr_prepared_input(allele_probs, pheno, covar, kinship, marker_map, selected_loci, exposure, outcome)` | `qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar, marker_map)`, then `qtl2mr_pair_input(resources, data.frame(exposure = exposure, marker = selected_loci), exposure, outcome)`. The function is no longer exported |
| `hsq =` in `as_qtl2mr_prepared_input()` | removed; hsq is always estimated by the package's REML |
| `sample_policy = "pair_common_complete"` | `qtl2mr_combine_resources(..., common_only = TRUE)`, or subset the inputs; otherwise each trait uses its own complete samples (the former default) |
| `F_X_min`, `max_loci` in `qtl2mr_select_instruments()` | `subsets = list(F2 = list(F_X_min = 2))` in `run_qtl2mr_analysis()` or `run_qtl2mr_batch()`, or edit the instrument table. The table keeps the two columns, always `NA` |
| `peak_rule = "legacy"` | `peak_rule = "spacing"`, which starts from the same island maxima; qtl2mr 0.6.1 reproduces the old selection exactly |
| `loci =` in `run_qtl2mr_analysis()` | edit the instrument table, or use `subsets` |
| `heterogeneity_method = "simulation"`, `heterogeneity_nsim`, `heterogeneity_seed`, `diagnostic_bootstrap_B`, `diagnostic_bootstrap_seed` | `qtl2mr_diagnostics(result, nsim, seed)` after the run |
| `verbose` (with `progress`) | `progress` alone |
| `on_nonestimable_locus`, `reported_estimator_policy`, `cross_cov_model`, `standardize_phenotypes`, `run_loo_mr`, `run_instrument_sensitivity` | fixed: a non-estimable locus fails the pair; the second-order estimate is reported only when bounded; star-space covariance; standardised phenotypes, with the `_raw_units` columns for original units; leave-one-out and sensitivity always run in a single analysis and never in a batch |
| `probability_policy` in `qtl2mr_prepare_resources()` | fixed: rounding errors are normalised with a warning and a record |
| `min_founder_copies = NULL` by default (0.6.1.9003 to 0.6.1.9005; QC off) | `min_founder_copies = 5` by default; `NULL` turns the QC off |

The fixed settings are still recorded in the result tables and in
`batch_config.json`, so table formats are unchanged. They are the former batch
defaults, so a batch run with the same instruments and the
same `min_founder_copies` is designed to give results bit for bit identical to
0.6.1.9005. For a single pair, three defaults change: the founder QC is on;
with several loci the Q p-value refers to chi-square instead of simulation;
and the weighted median and mode have no bootstrap standard errors. The last
two do not change the estimate or its interval; the QC changes them only for
pairs where it drops a locus. To reproduce results that need a removed
option, use the version that had it, such as the 0.6.1 release.

## Upgrading from 0.5.x

| 0.5.x | 0.6.0 |
|---|---|
| `qtl2mr_prepare_resources(probs, K, pheno, covar, qtl_scan = scan, marker_map = map)` | `qtl2mr_prepare_resources(probs, K, pheno, covar, marker_map = map)`; the scan goes to `qtl2mr_select_instruments()` |
| `qtl2mr_select_loci_from_scan(scan, map, "x", lod_thresh = 4)` | `qtl2mr_select_instruments(resources, scan, exposures = "x", lod_thresh = 4)` |
| `run_qtl2mr_batch(r, p, instrument_rule = list(lod_thresh = 4))` | `run_qtl2mr_batch(r, p, qtl2mr_select_instruments(r, scan, exposures = unique(p$exposure), lod_thresh = 4))`; without `exposures` every trait in the scan is used |
| pair columns `instrument_trait`, `lod_thresh`, `exclude_chr`, `max_loci`, `selected_loci_file` | one instrument set per exposure in the `instruments` table; these columns are now an error |
| `sample_policy`, `min_sample_retention` | removed: every pair uses role-specific samples, as in `run_qtl2mr_analysis()` |
| `run_qtl2mr_batch(..., fit_options = list(...), force = TRUE)` | top-level arguments and `overwrite = TRUE`; `fit_options$q_sim_n`/`q_sim_seed` became `heterogeneity_nsim`/`heterogeneity_seed` (removed in 0.6.1.9006; see [Changes since 0.6.1](#changes-since-061)); `engine`, `cache_contexts`, `cache_lmm_fits` and `verbose` are internal (low-level `fit_multilocus_founder_mr()` keeps `q_sim_n`) |
| `run_qtl2mr_analysis(..., standardize_phenotypes = , cross_cov_model = )` through `...` | named arguments in 0.6.0 (both removed in 0.6.1.9006, with `verbose`; see [Changes since 0.6.1](#changes-since-061)); `...` is gone, and `ridge`, `engine` and the caches are internal. Pass arguments by name |
| `batch$selected_loci`, `batch$first_stage_summary`, `batch$sample_qc` (and their CSV files) | `batch$locus_results` (`batch_locus_results.csv`); per-trait `batch$sample_exclusions` |
| batch `main_results` columns `instrument_trait`, `lod_thresh`, `exclude_chr` | removed; the instruments table records the selection |
| `theta_so`, `se_so`, `ci_lower_so`, `ci_upper_so`, `sample_policy` columns | `theta`, `se`, `profile_ci_lower`, `profile_ci_upper`; no sample-policy column |
| `fit_instrument_set_sensitivity()` default `weak_f_threshold = 10`, set `exclude_F_lt_10` | default 2, set `exclude_F_lt_2`; pass `weak_f_threshold = 10` for the old set |

Batch estimates change for pairs whose exposure and outcome have different
missing samples: 0.5.x fitted both traits on their common samples, 0.6.0 fits
each on its own complete samples. Pairs without such differences give the same
numbers. The batch Q p-value now defaults to chi-square. Batch folders from
0.5.x cannot be resumed. Several 0.6.0 arguments named in this table were
removed in 0.6.1.9006; see [Changes since 0.6.1](#changes-since-061). For older
versions, see [NEWS.md](NEWS.md).

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
- `?run_qtl2mr_analysis`, `?run_qtl2mr_batch`, `?qtl2mr_combine_resources`,
  `?qtl2mr_diagnostics`, `help(package = "qtl2mr")`
- [NEWS.md](NEWS.md)

## Citation

Until the methods paper is available, cite this repository and the package
version you used (`citation("qtl2mr")`, `packageVersion("qtl2mr")`).
