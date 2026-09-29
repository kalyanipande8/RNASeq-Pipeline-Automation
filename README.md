# RNA-Seq Pipeline Automation & Harmonization

A real, runnable RNA-seq processing pipeline written in **R**, run end-to-end
against a synthetic patient-derived-model compendium with a fixed random
seed. This project implements, honestly, the claims behind the resume
bullet:

> RNA-Seq Pipeline Automation & Harmonization | NGS Data Processing,
> Reproducible Research
> - Built automated RNA-seq processing pipeline in R; analyzed 320+ samples
>   from patient-derived models.
> - Harmonized expression data from multiple sequencing batches;
>   standardized metadata across 14 cohorts.
> - Developed QC scripts identifying quality issues in 8 datasets; improved
>   reproducibility and data accessibility.

## Disclaimer: synthetic data

**All data in this project is synthetically generated** by
`simulate_rnaseq_data.R` with a fixed random seed (`set.seed(42)`). There is
no real patient, tumor, PDX/organoid, or sequencing data anywhere in this
repository. Gene expression counts, cohort assignments, batch assignments,
and metadata (patient/model IDs, sequencing dates) are all simulated. The
"8 of 16 sequencing batches with quality issues" and "3 corrupted samples"
are **deliberately injected** ground truth, used to sanity-check that the QC
logic actually detects what it's supposed to detect (see `tests/`).

## Honest accounting: why R (not a Python substitution)

Before building anything, the environment was checked directly:

```
which Rscript      -> /usr/local/bin/Rscript (R 4.5.1, already installed)
df -h /             -> disk was low (single-digit GB free) at the time of
                        this build, on a machine that has been tight on
                        disk space across other recent projects in this
                        Career/projects tree
Rscript -e 'installed.packages()' -> testthat, ggplot2, jsonlite already
                        installed; no new package installation was needed
```

Because `Rscript` and every R package this pipeline needs (`testthat` for
tests, `ggplot2` for figures — both base-R-adjacent and already present)
were already available, **no new installation was required**, and the
resume bullet explicitly says "Built ... pipeline in R." So, unlike the
sibling `VariantAnalysisFramework` project (which substituted Python for
`bedtools`/`pysam` because those were genuinely missing and installing them
was judged unsafe given disk constraints), there was no reason to substitute
Python here — **every script in this pipeline is real, runnable R**, executed
with `Rscript`, and every number in the Results section below was read
directly off an actual run's output files. Nothing was hand-typed or
fabricated.

## What it does

1. **`simulate_rnaseq_data.R`** simulates a 600-gene x 322-sample raw count
   matrix (negative-binomial counts, kept to 600 genes so the whole pipeline
   runs in seconds) for patient-derived models spread across 14 cohorts and
   16 sequencing batches. Each cohort perturbs ~15% of genes (biological
   cohort effect); each batch applies a multiplicative depth/mean shift plus
   a dispersion inflation (technical batch effect). **8 of the 16 batches are
   deliberately degraded** (higher noise, lower/more variable depth) and
   **3 samples are deliberately corrupted** (near-empty libraries) as ground
   truth for the QC step. Metadata is written in **4 different raw formats**
   (`metadata_style_A/B/C/D.csv`) — different column names, different ID/
   cohort/batch spellings, and different date formats (ISO, US, `DD-Mon-YYYY`,
   and raw epoch-day) — mirroring how a real multi-cohort, multi-core
   sequencing compendium actually arrives.
2. **`harmonize.R`** performs the real harmonization:
   - **Metadata**: parses all 4 raw formats into one schema
     (`sample_id, cohort, batch, seq_date, raw_format`), standardizing cohort
     IDs (`Cohort_7` / `coh-07` / `COH-07` -> `COH-07`), batch IDs
     (`Batch-3` / `b3` -> `B03`), and dates (4 different input formats -> one
     `Date` type).
   - **Expression**: library-size (CPM) normalization followed by **quantile
     normalization** (implemented directly, no external package) on
     log2(CPM+1) values, correcting both per-sample sequencing depth and
     batch-to-batch distributional shift, producing one unified
     genes x samples matrix.
3. **`qc_checks.R`** runs three real QC checks across all 322 samples using
   the raw counts:
   - **Library-size outliers** (MAD rule + absolute floor on raw read count).
   - **Sample-to-global correlation**: each sample's log-expression profile
     is Spearman-correlated against the across-all-samples mean profile.
     A **batch** is flagged as a quality-issue batch when its *median*
     sample-to-global correlation falls below 0.87 — this is what actually
     surfaces systematically noisy sequencing runs, as opposed to one-off
     sample outliers.
   - **Mitochondrial-content proxy**: the top-expressed-gene fraction of a
     sample's total counts (a standard proxy QC metric when no explicit
     mitochondrial gene annotation exists).
   A sample is flagged `qc_fail` if it fails >=2 of the 3 checks, or if its
   library size is catastrophically low.
4. **`make_figures.R`** generates 2 real figures from the actual QC output.
5. **`run_all.R`** orchestrates all four steps.
6. **`tests/test_pipeline.R`** (testthat) makes real assertions against the
   actual output on disk: 320+ samples present, 14 cohorts present, one
   consistent harmonized schema, and — as a sanity check on the QC logic
   itself — that the deliberately-corrupted samples and deliberately-degraded
   batches are actually the ones the QC step flags.

## How to run

```bash
cd RNASeq-Pipeline-Automation
Rscript run_all.R          # simulate -> harmonize -> QC -> figures
Rscript tests/run_tests.R  # run the pipeline's testthat suite
```

Each script can also be run individually and in order
(`simulate_rnaseq_data.R`, `harmonize.R`, `qc_checks.R`, `make_figures.R`).

## Results (from an actual run on 2026-09-29)

```
=== 1/4 Simulating RNA-seq compendium ===
Simulated 322 samples x 600 genes across 14 cohorts and 16 batches
  Degraded batches (ground truth, 8 of 16): B01, B05, B16, B09, B10, B04, B02, B14
  Deliberately corrupted samples (ground truth): PDX-0116, PDX-0124, PDX-0229

=== 2/4 Harmonizing metadata + expression ===
Harmonized metadata for 322 samples across 14 cohorts (from 4 raw formats) into one schema
Harmonized expression: 600 genes x 322 samples -> CPM + quantile-normalized (log2 scale)

=== 3/4 Running QC checks ===
QC checks run on 322 samples across 16 batches
  14 samples flagged qc_fail (>=2 checks failed, or near-empty library)
  Quality issues concentrated in 8 of 16 batches: B14, B10, B16, B04, B09, B01, B02, B05
```

**QC correctly recovered all 8 injected degraded batches** (median
sample-to-global Spearman correlation, threshold 0.87):

| Batch | Median corr. | Flagged samples | Batch samples | % flagged |
|---|---|---|---|---|
| B14 | 0.708 | 2 | 19 | 10.5% |
| B10 | 0.721 | 4 | 21 | 19.0% |
| B16 | 0.726 | 1 | 21 | 4.8% |
| B04 | 0.765 | 3 | 20 | 15.0% |
| B09 | 0.767 | 2 | 30 | 6.7% |
| B01 | 0.776 | 0 | 17 | 0.0% |
| B02 | 0.815 | 1 | 19 | 5.3% |
| B05 | 0.818 | 0 | 23 | 0.0% |
| ...8 clean batches (B03, B06, B07, B08, B11, B12, B13, B15) | 0.912–0.924 | 0–1 | 16–26 | 0–6.2% |

All 3 deliberately corrupted samples (`PDX-0116`, `PDX-0124`, `PDX-0229`)
were flagged `qc_fail`, along with 11 additional samples whose profiles were
degraded by their batch's injected noise.

**Test suite**: `Rscript tests/run_tests.R` -> **`[ FAIL 0 | WARN 0 | SKIP 0 | PASS 18 ]`**
(18/18 assertions pass), including the sanity check that QC recovers exactly
the 8 injected degraded batches and all 3 injected corrupted samples.

## Figures

`output/figures/fig1_batch_qc_correlation.png` — per-batch median
sample-to-global correlation, with the 0.87 QC threshold and the 8
flagged batches highlighted.

`output/figures/fig2_library_size_by_batch.png` — raw library size per
sample by batch (log scale), colored by final QC pass/fail status.

## Project structure

```
RNASeq-Pipeline-Automation/
├── simulate_rnaseq_data.R      # step 1: simulate counts + 4-format raw metadata
├── harmonize.R                 # step 2: unify metadata schema + normalize expression
├── qc_checks.R                 # step 3: library/correlation/mito-proxy QC checks
├── make_figures.R              # step 4: 2 figures from real QC output
├── run_all.R                   # orchestrator
├── data/
│   ├── raw/                    # simulated raw counts + 4 raw-format metadata files
│   └── processed/              # unified expression matrix + unified metadata
├── output/
│   ├── qc_report_per_sample.csv
│   ├── qc_report_per_batch.csv
│   └── figures/
├── tests/
│   ├── test_pipeline.R         # testthat assertions
│   └── run_tests.R             # test runner
└── README.md
```

## Tech stack

- **R 4.5.1** (`Rscript`, base R + `stats`)
- **testthat** for unit tests
- **ggplot2** for figures
- Fixed-seed synthetic data generation (`set.seed(42)`) for full
  reproducibility
