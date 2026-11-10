# data-generator

Scripts that produced the two example datasets in `data/`. They are not part of
the installed package.

- `make_qtl2mr_example.R` regenerates `data/qtl2mr_example.rda` from a fixed
  seed (synthetic data only).
- `make_qtl2mr_multilocus_example.R` extracts `data/qtl2mr_multilocus_example.rda`
  from a frozen synthetic input that is not distributed here; see
  `?qtl2mr_multilocus_example`.
