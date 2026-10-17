# qtl2mr

Founder-aware Mendelian randomization for `qtl2`-style multiparent
mouse data.

`qtl2mr` is an R package for selected-locus MR analyses in Diversity
Outbred and related multiparent populations. Genetic instruments are
represented by founder haplotype probabilities rather than scalar SNP
dosages.

Current release: **v0.4.0**.
This development branch contains unreleased v0.4.1 fixes; v0.4.1 has not been
published. The installation links below intentionally point to v0.4.0.

## Installation

Download `qtl2mr_0.4.0.tar.gz` and `SHA256SUMS` from the
[v0.4.0 release](https://github.com/XDai-bris/qtl2mr-uvmr/releases/tag/v0.4.0).
Use the package file listed under Assets, not GitHub's automatically generated
source-code archive. Keep the downloaded files in your working directory.

Create a clean environment in a terminal:

```bash
conda create -n qtl2mr-demo -c conda-forge r-base=4.4 r-fs -y
conda activate qtl2mr-demo
R
```

In R:

```r
pkg_file <- "qtl2mr_0.4.0.tar.gz"
stopifnot(file.exists(pkg_file))
install.packages(pkg_file, repos = NULL, type = "source")
library(qtl2mr)
stopifnot(as.character(packageVersion("qtl2mr")) == "0.4.0")
vignette("getting-started", package = "qtl2mr")
```

The package includes three prebuilt HTML vignettes; you do not need to knit
them during installation. Compare the package file's SHA-256 with `SHA256SUMS`
before installing.

`qtl2` is not required for fitting prepared inputs. It is useful for upstream
genotype processing, QTL scans, or the optional comparison with `qtl2::scan1coef()`:

```r
install.packages("qtl2", repos = "https://cloud.r-project.org")
```

## Quick example

The package includes `qtl2mr_example`, a small synthetic dataset with
96 samples, eight founders, four markers, and two selected loci.

```r
library(qtl2mr)
data("qtl2mr_example", package = "qtl2mr")
```

Prepare the input:

```r
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

validate_qtl2mr_input(
  prepared,
  require_selected_loci = TRUE
)
```

Run the two-locus analysis:

```r
result_dir <- file.path(tempdir(), "qtl2mr-example")

result <- run_qtl2mr_analysis(
  prepared,
  adapter = "prepared",
  loci = 2,
  output_dir = result_dir,
  overwrite = TRUE,
  make_figures = FALSE,
  make_html = TRUE,
  verbose = FALSE,
  run_qtl2_validation = FALSE,
  run_heteroskedasticity = FALSE,
  diagnostic_bootstrap_B = 0L,
  heterogeneity_method = "chisq",
  run_loo_mr = FALSE,
  run_instrument_sensitivity = FALSE,
  progress = FALSE
)
```

These explicit fast-demo settings match the getting-started vignette; they
do not change package defaults. The default `heterogeneity_method = "auto"`
uses a chi-square Q reference for one locus and simulation for multiple loci
(`heterogeneity_nsim = 1000L`). Each simulated Q draw refits MR, so this can
dominate runtime even for two loci. Diagnostic bootstraps, leave-one-locus-out
fits and sensitivity analyses also add time. Use `progress = TRUE` to see
stage updates; the demo's approximate chi-square reference is not equivalent
to the default simulated reference or a calibration guarantee.

Read the reported estimate and any warnings:

```r
primary <- qtl2mr_primary_result(result)
warnings <- qtl2mr_warnings(result)

primary
warnings
```

Open the HTML report:

```r
browseURL(file.path(result$output_dir, "report.html"))
```

Read the estimate together with its confidence-set type and warnings. The
legacy `primary_reliable` field now means a numerically acceptable fit with
one bounded confidence-set component. It does not establish calibrated
coverage or a valid causal interpretation. A weak-founder-block warning
is a reason to investigate identification, not a calibrated pass/fail rule.

`qtl2mr_example` is synthetic teaching data. It is useful for checking an
installation and learning the API, but it is not calibration or biological
evidence.
Its `true_theta = 0.75` is in original outcome units per exposure unit;
the default output is in SD units. Do not compare those numbers directly.
With a single common scaling, standardized truth is raw truth times
SD(X)/SD(Y); the conversion is not a claim that this one fit recovers truth.

## Multi-locus example

For a larger worked example, load the bundled 400-sample, 40-locus synthetic
dataset and open its walkthrough:

```r
data("qtl2mr_multilocus_example", package = "qtl2mr")
vignette("multilocus-example", package = "qtl2mr")
```

It includes covariates and synthetic kinship, and walks through the MR
estimate, leave-one-locus-out checks, founder residuals and HTML report.
It is a constructed teaching example, not method-calibration evidence.
The shipped data reproduce the analysis, but the public extraction script
does not regenerate the upstream frozen simulation from its recorded seed.
It requires a separate frozen synthetic input that is not distributed here.
The omitted traits and background loci cannot be reconstructed from this
40-locus subset. See the walkthrough's provenance boundary.

## Using your own data

The recommended interface is a prepared input object:

```r
prepared <- as_qtl2mr_prepared_input(
  allele_probs = allele_probs,
  pheno = pheno,
  covar = covar,
  kinship = kinship,
  marker_map = marker_map,
  selected_loci = selected_loci,
  exposure = "trait_x",
  outcome = "trait_y"
)
```

The main inputs are:

- `allele_probs`: chromosome-named arrays of founder **allele probabilities**,
  with dimensions sample x eight founders x marker and probabilities summing
  to one across founders; these are not the full diplotype probabilities;
- `pheno`: exposure and outcome phenotypes, with sample IDs as row names;
- `covar`: additive covariates using the same sample IDs;
- `kinship`: preferably a chromosome-named LOCO kinship list;
- `marker_map`: marker chromosome and position information;
- `selected_loci`: the loci to use as instruments.

Sample IDs must be harmonized across these objects. The single-pair prepared
analysis fits exposure and outcome on their separate role-specific usable
samples; overlap is determined from the fitted sample IDs.

Locus selection takes place before the main MR analysis. The package does not
treat the bundled example thresholds as general instrument-selection rules.
Missing `selected_loci` is not permission to use every marker. Select loci
explicitly, or deliberately opt into `use_all_loci = TRUE`. For X-to-Y MR,
the exposure-selected loci are fitted in both X and Y. A union of X and Y
peaks is an explicit instrument-design choice, not an automatic requirement.

The default `standardize_phenotypes = TRUE` reports SD of outcome per SD of
exposure. The per-locus scaling and sample counts are recorded in
`data/phenotype_scaling.csv`; use `standardize_phenotypes = FALSE` for original
units. Do not back-transform a pooled result using a single SD ratio if the
underlying locus fits used different sample sets and SDs.

## Choosing a workflow

| Analysis | Function | Intended use |
|---|---|---|
| One exposure–outcome pair | `run_qtl2mr_analysis()` | Main analysis and full report |
| Several pairs, exact samples | `run_qtl2mr_batch()` | Pair-specific sample analysis |
| Large exact batches | `run_qtl2mr_batch_exact_grouped()` | Reuse work for identical sample sets |
| Fast multi-pair screen | `run_qtl2mr_batch_fast()` | Discovery-stage screening |
| Parallel exposure blocks | `run_qtl2mr_batch_parallel_by_exposure()` | Larger screening runs |

The exact batch workflows use **pair-common-complete** samples: both traits
and covariates must be usable for the pair. This is not the single-pair
analysis's role-specific sample policy. To compare the engines, restrict the
single analysis with `sample_ids` (and any explicit role-specific sample
filters) and the batch with `sample_rule$sample_ids` to the **same** usable
IDs, and match loci, kinship, hsq, scaling and fitting options. A different
sample mask is a different analysis, not by itself an engine discrepancy.

The fast workflow uses a common complete sample set across all requested
traits and may drop non-estimable loci, with those changes
recorded in its result tables and run manifest. Important findings from a fast
screen should be rerun with an exact workflow. Exact sample handling does not
by itself establish calibrated inference.
Batch resume is opt-in. A resumed run must match its saved input/configuration
fingerprint; changing the data or options requires a new run directory.

## Results

A standard analysis directory includes:

- `report.html`: the main review report;
- `summary.md`: a short text summary;
- `data/main_result.csv`: the reported estimate and uncertainty;
- `data/warnings.csv`: structured warnings;
- `data/output_manifest.csv`: files produced by the run;
- `result_object.rds`: the complete R result object.

The main programmatic accessors are:

```r
qtl2mr_primary_result(result)
qtl2mr_warnings(result)
qtl2mr_method_comparison(result)
qtl2mr_output_paths(result)
```

In the method-comparison plot, weighted median/mode use approximate 95%
intervals: the original estimate plus or minus 1.96 times parametric-bootstrap
SD (500 draws by default). Loci and weights stay fixed; per-locus effects are
drawn independently using their estimates and SEs. These are quick diagnostic
intervals, not profile CIs or a coverage guarantee. The plot and HTML report
label the distinction. Set `diagnostic_bootstrap_B = 0L` to omit these intervals.

A smaller package of files can be prepared for sharing with:

```r
export_qtl2mr_lite_result(result)
```

Lite export is not anonymization. Study results and identifiers may remain;
review the included files and obtain permission before sharing restricted data.

## Method summary

At each selected locus, `qtl2mr` estimates exposure and outcome effects in
founder-contrast space using kinship-aware GLS/LMM. The locus estimates are
combined with a second-order MR objective that includes founder-effect
uncertainty and cross-trait covariance arising from sample overlap.

The default reporting policy is `diagnostic_only`. A numerically acceptable,
bounded second-order estimate is reported. Unbounded, disconnected or unresolved
confidence sets do not become an ordinary finite CI. Inspect
`data/confidence_set.csv` and the raw fit when no estimate is reported.
Two-step IVW results remain diagnostic rather than automatic replacements.

The profile set uses a nominal chi-square cutoff. Weak instruments can give
poor coverage even when the numerical fit succeeds. The summed objective
does not include cross-locus covariance, and same-sample instrument selection
can introduce selection bias. Neither overlap correction nor a successful
optimizer removes these limitations. Whole-line numerical search is not a
mathematical certificate of the global optimum.

The package also provides leave-one-locus-out analyses, instrument-set
sensitivity, overlap summaries, founder residual diagnostics, first-stage
summaries, PSD checks, and optional comparison with `qtl2::scan1coef()`.
Kinship decomposition and ML/REML heritability estimation use the package's
internal base-R LMM backend; `qtl2` is used only when an explicit parity check
is requested and available.

## Current scope

The high-level workflow currently supports:

- eight-founder DO and related founder allele probabilities;
- autosomal additive selected loci;
- continuous exposure and outcome traits;
- additive covariates;
- one kinship source, preferably chromosome-specific LOCO kinship;
- full, partial, or zero sample overlap.

The following are not currently supported:

- chromosome X models;
- binary outcomes;
- automatic conversion of arbitrary `qtl2` cross objects;
- separate exposure and outcome kinship sources;
- a universally calibrated weak-instrument threshold for founder blocks.

## Documentation

- [Getting started](vignettes/getting-started.Rmd)
- [Multi-locus example](vignettes/multilocus-example.Rmd)
- [From qtl2 inputs](vignettes/from-qtl2.Rmd)
- Package help: `help(package = "qtl2mr")`
- Main workflow: `?run_qtl2mr_analysis`
- Batch workflows: `?run_qtl2mr_batch`
- Release history: [NEWS.md](NEWS.md)

## Citation

A methods manuscript is in preparation. Until it is available, cite this
repository and the installed package version used for the analysis.
