# qtl2mr

Founder-aware Mendelian randomization for Diversity Outbred and related
eight-founder mouse populations. Selected QTL are the instruments, each
represented by its founder haplotype probabilities from `qtl2`, and the
estimator accounts for samples shared by the exposure and the outcome.

**Current release:** [v0.7.0](https://github.com/XDai-bris/qtl2mr-uvmr/releases/tag/v0.7.0).
What changed: [NEWS.md](NEWS.md).

## Installation

Download `qtl2mr_0.7.0.tar.gz` and `SHA256SUMS` from the release assets (not
GitHub's "Source code" archive, which lacks the built vignettes), check them
with `shasum -a 256 -c SHA256SUMS`, then in R:

```r
install.packages("fs")
install.packages("qtl2mr_0.7.0.tar.gz", repos = NULL, type = "source")
```

`qtl2` is not needed to run qtl2mr; you will normally use it beforehand for
allele probabilities, kinship and QTL scans.

## One pair

```r
library(qtl2mr)
data("qtl2mr_example", package = "qtl2mr")
ex <- qtl2mr_example

resources <- qtl2mr_prepare_resources(ex$allele_probs, ex$kinship, ex$pheno, ex$covar, ex$marker_map)
input <- qtl2mr_pair_input(resources, data.frame(exposure = "x", marker = c("m1", "m2")), "x", "y")
result <- run_qtl2mr_analysis(input, output_dir = tempfile())

summary(result)
browseURL(file.path(result$output_dir, "report.html"))
```

The result folder holds the tables, figures and an HTML report. Here the
estimate is 0.892 SD of outcome per SD of exposure (95% CI 0.685 to 1.126).

## Many pairs

```r
resources <- qtl2mr_prepare_resources(allele_probs, kinship, pheno, covar, map)
instruments <- qtl2mr_select_instruments(resources, scan, exposures = c("x", "w"))
pairs <- data.frame(exposure = c("x", "x", "w"), outcome = c("y", "z", "z"))

batch <- run_qtl2mr_batch(resources, pairs, instruments,
                          output_dir = "batch-results", n_workers = 4)
batch$main_results
```

A batch row equals `run_qtl2mr_analysis()` on that pair, so pairs of interest
can be rerun one at a time for the full report. An interrupted batch continues
with `resume = TRUE`.

## What qtl2mr does by default

- **Instruments:** one per LOD peak (LOD ≥ 4, qtl2's peak rule), fitted to
  report its first-stage strength `F_X_block`; nothing is filtered on F.
- **Founder QC:** a locus is dropped when its rarest founder is almost absent
  from the fitted mice (fewer than 5 expected copies; `min_founder_copies = NULL`
  turns this off).
- **Samples:** each trait uses its own complete samples.
- **Estimate:** the overlap-aware second-order estimate with its profile
  confidence interval, on standardised traits (original units are given as
  well). Q heterogeneity is referred to chi-square.

## Your own data

- `allele_probs`: founder allele probabilities from
  `qtl2::genoprob_to_alleleprob()`, a list by chromosome.
- `pheno`, `covar`: numeric data frames with sample IDs as row names.
- `kinship`: LOCO kinship from `qtl2::calc_kinship(probs, "loco")`.
- `marker_map`: `marker`, `chr`, `pos`, or a qtl2 map.
- `scan`: `qtl2::scan1()` output, or your own instrument table with `exposure`
  and `marker` columns.

Only autosomal loci, continuous traits and additive covariates are supported.
`vignette("from-qtl2", package = "qtl2mr")` walks through a qtl2 cross.

## Reading a result

- The estimate is `theta_reported` with `reported_ci_lower` and
  `reported_ci_upper`.
- `NA` means the confidence set was unbounded or split, usually because the
  instruments are weak; check `F_X_block` rather than switching estimator.
- Weak instruments can make intervals too narrow even when the fit succeeds.
- `qtl2mr_warnings(result)` explains each warning and what to do.

## More

- Refit on instrument subsets, for example F ≥ 2: `subsets` in
  `run_qtl2mr_analysis()` and `run_qtl2mr_batch()`.
- F summary of an instrument table: `qtl2mr_instrument_summary()`.
- Combine two data sets: `?qtl2mr_combine_resources`.
- Simulated Q p-value and bootstrap intervals: `?qtl2mr_diagnostics`.
- Vignettes: `getting-started`, `from-qtl2`, `multilocus-example`.

Upgrading from 0.6.x or earlier: many options were removed in 0.7.0; see
[NEWS.md](NEWS.md) for what replaces each one.

## Method in brief

At each locus, exposure and outcome are regressed on the founder
probabilities with a linear mixed model (LOCO kinship). The locus estimates
are combined by a second-order objective that includes the uncertainty in the
exposure effects and the covariance from shared samples; the interval is the
profile set of that objective. Limitations: weak instruments can give poor
coverage, covariance between loci is not modelled, selecting instruments in
the same samples can bias the estimate, and heritability is treated as known.

## Citation

A methods manuscript is in preparation. Until then, cite this repository and
the version you used (`citation("qtl2mr")`).
