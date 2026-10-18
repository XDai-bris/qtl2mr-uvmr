# qtl2mr

Founder-aware Mendelian randomization for Diversity Outbred and related
eight-founder mouse populations.

`qtl2mr` estimates the causal effect of an exposure trait on an outcome trait
using selected QTL as instruments. Each instrument is represented by founder
haplotype probabilities from `qtl2`, not by a single SNP dosage, and the
estimator accounts for exposure-outcome sample overlap.

**Current release:** [v0.5.0](https://github.com/XDai-bris/qtl2mr-uvmr/releases/tag/v0.5.0).
See [NEWS.md](NEWS.md) for what changed and the
[migration table](#upgrading-from-04x) if you used 0.4.x.

## Installation

Download `qtl2mr_0.5.0.tar.gz` and `SHA256SUMS` from the release assets. Use
the package file, not GitHub's automatic "Source code" archive, because the
package file contains the built vignettes. Check the SHA-256 against
`SHA256SUMS`, then in R:

```r
install.packages("fs")  # with conda: conda install -c conda-forge r-fs
install.packages("qtl2mr_0.5.0.tar.gz", repos = NULL, type = "source")
library(qtl2mr)
vignette("getting-started", package = "qtl2mr")
```

`qtl2` is not needed to run analyses. You will normally use it upstream to
compute allele probabilities, kinship and QTL scans:

```r
install.packages("qtl2", repos = "https://cloud.r-project.org")
```

## Two ways to use qtl2mr

| You want to analyse | Use | You get |
|---|---|---|
| One exposure → outcome pair | `as_qtl2mr_prepared_input()` then `run_qtl2mr_analysis()` | A result folder with tables, figures and an HTML report |
| Many pairs | `qtl2mr_prepare_resources()` then `run_qtl2mr_batch()` | Compact tables with one row per pair; optional parallel workers and resume |

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
  loci = 2,
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

qtl2mr_primary_result(result)
qtl2mr_warnings(result)
browseURL(file.path(result$output_dir, "report.html"))
```

These settings skip the optional diagnostics so the example finishes in about
a second; they are not the package defaults. The reported estimate is about
0.892 SD of outcome per SD of exposure, 95% CI [0.685, 1.126]. The example's
`true_theta = 0.75` is in original units, so the two numbers are not directly
comparable. The example is teaching data, not evidence about the method.

### Many pairs

```r
resources <- qtl2mr_prepare_resources(
  allele_probs, kinship, pheno, covar,
  qtl_scan = scan,          # e.g. qtl2::scan1() output, markers x traits
  marker_map = map
)
pairs <- data.frame(exposure = c("x", "x", "w"), outcome = c("y", "z", "z"))

batch <- run_qtl2mr_batch(
  resources, pairs,
  instrument_rule = list(lod_thresh = 4),
  output_dir = "batch-results",
  n_workers = 4
)
batch$main_results[, c("pair_id", "status", "theta_reported",
                       "reported_ci_lower", "reported_ci_upper")]
```

Instruments are selected for each exposure from the scan (one peak per
contiguous above-threshold region per chromosome). Each pair is fitted on the
samples complete for both of its traits. If a run is interrupted, call it
again with `resume = TRUE`: completed pairs are kept, the output folder may
have been moved, and `n_workers` may differ. A resume is refused if the
inputs, options or qtl2mr code changed (for example after an upgrade).

Batch defaults skip the Q simulation and bootstrap intervals; rerun pairs of
interest with `run_qtl2mr_analysis()` for a full report.

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
- `selected_loci`: the instrument markers, for example from
  `qtl2mr_select_loci_from_scan()`. Choose the threshold before seeing MR
  results; the package never uses every marker unless you set
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
  (block F below about 10, flagged as `weak_block_first_stage`) profile
  intervals can be too narrow.
- **Heritability.** `data/hsq.csv` lists the null-model heritability used for
  each trait and chromosome. It sets how much the kinship term reweights the
  LMM; a value of 0 means kinship had no effect there.
- **Warnings.** `qtl2mr_warnings(result)` gives each warning's consequence and
  recommended action.
- **Overview.** `summary(result)` gives a one-row overview; `confidence_set_bounded`
  is `FALSE` when the confidence set reaches plus or minus infinity.

## Runtime

With several loci, `run_qtl2mr_analysis()` calibrates the Q heterogeneity
p-value by simulation (`heterogeneity_nsim = 1000`), refitting the model once
per draw. On a laptop, the bundled 400-mouse, 40-locus example takes roughly
10–15 minutes with the default and under 2 minutes with
`heterogeneity_method = "chisq"`. Progress messages report completed draws.

## Advanced use

- **Batch options.** `sample_policy = "batch_common_complete"` fits every pair
  on one sample set complete for all batch traits (faster to share, but it can
  change estimates; `min_sample_retention` stops the run if too many samples
  are lost). `on_nonestimable_locus = "drop"` removes a locus whose founder
  contrast cannot be estimated and records it, instead of failing the pair.
  See `?run_qtl2mr_batch`.
- **Follow-up refits.** `refit_qtl2mr_subset()` and
  `fit_instrument_set_sensitivity()` re-estimate an analysis on locus subsets
  without refitting the LMMs.
- **Low-level interface.** `fit_founder_lmm_block()`, `fit_mr_locus_block()`,
  `fit_multilocus_founder_mr()` and `second_order_gls_ivw()` expose the
  estimator for method work (`?fit_mr_locus_block`).
- **qtl2 parity check.** `validate_qtl2_scan1coef()` compares the qtl2mr LMM
  with `qtl2::scan1coef()` for one marker.

## Upgrading from 0.4.x

| 0.4.x | 0.5.0 |
|---|---|
| `run_qtl2mr_batch(r, p)` | unchanged |
| `run_qtl2mr_batch(r, p, sample_rule = list(mode = "batch_common_complete"))` | `run_qtl2mr_batch(r, p, sample_policy = "batch_common_complete", min_sample_retention = 0)` (0.4.x had no retention gate here; the new default stops below 85% retention) |
| `run_qtl2mr_batch_fast(r, p)` | `run_qtl2mr_batch(r, p, sample_policy = "batch_common_complete", on_nonestimable_locus = "drop")` |
| `run_qtl2mr_batch_exact_grouped(r, p)` | `run_qtl2mr_batch(r, p)` |
| `run_qtl2mr_batch_parallel_by_exposure(r, p, n_workers = k)` | `run_qtl2mr_batch(r, p, sample_policy = "batch_common_complete", on_nonestimable_locus = "drop", n_workers = k)` |
| `run_qtl2mr_batch_parallel_by_exposure(r, p, fast = FALSE, n_workers = k)` | `run_qtl2mr_batch(r, p, n_workers = k)` |
| `fit_options = list(save_mr_block_cache = TRUE)`, `run_first_stage_filter_grid`, `steiger_filter`, `refit_qtl2mr_mr_block_cache()` | removed; for one pair, refit locus subsets with `fit_instrument_set_sensitivity()` |
| `run_qtl2mr_analysis(..., adapter = "prepared")` | drop `adapter` (and `run_qtl2_validation`, `run_heteroskedasticity`, `heteroskedasticity_B`) |
| `run_qtl2mr_analysis(input, "x", "y", "prepared", out)` (positional) | `run_qtl2mr_analysis(input, "x", "y", out)`: `output_dir` is now the fourth argument |
| `summary(result)` (99 columns) | `qtl2mr_primary_result(result)`; `summary()` is now a one-row overview |

The three old batch functions were removed in 0.5.0; the replacements give the
same numbers. Batch output folders from 0.4.x cannot be resumed by 0.5.0.
Several internal helpers are no longer exported, and result labels now read
`qtl2mr` (for example `"Star-space qtl2mr"`); see NEWS.

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
