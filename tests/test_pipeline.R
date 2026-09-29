#!/usr/bin/env Rscript
# tests/test_pipeline.R
#
# Real assertions against the actual pipeline output on disk. Run this
# AFTER run_all.R (or run tests/run_tests.R, which runs the pipeline first).

library(testthat)

find_project_dir <- function() {
  candidates <- c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "tests", ".."))
  for (cand in candidates) {
    if (dir.exists(file.path(cand, "data", "raw"))) return(normalizePath(cand))
  }
  normalizePath(getwd())
}
PROJECT_DIR <- find_project_dir()

RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")

test_that("raw simulation produced 320+ samples", {
  raw <- read.csv(file.path(RAW_DIR, "raw_counts.csv"), check.names = FALSE)
  n_samples <- ncol(raw) - 1  # minus gene_id column
  expect_gte(n_samples, 320)
})

test_that("unified metadata contains all 14 cohorts and matches sample count", {
  meta <- read.csv(file.path(PROCESSED_DIR, "unified_metadata.csv"), stringsAsFactors = FALSE)
  expect_gte(nrow(meta), 320)
  expect_equal(length(unique(meta$cohort)), 14)
  expect_true(all(grepl("^COH-\\d{2}$", meta$cohort)))
  expect_true(all(grepl("^B\\d{2}$", meta$batch)))
  expect_false(any(is.na(meta$seq_date)))
})

test_that("harmonization produces one consistent schema from 4 raw formats", {
  meta <- read.csv(file.path(PROCESSED_DIR, "unified_metadata.csv"), stringsAsFactors = FALSE)
  expect_setequal(colnames(meta), c("sample_id", "cohort", "batch", "seq_date", "raw_format"))
  expect_setequal(unique(meta$raw_format), c("style_A", "style_B", "style_C", "style_D"))
  # every sample_id must be unique and non-missing
  expect_equal(length(unique(meta$sample_id)), nrow(meta))
  expect_false(any(is.na(meta$sample_id) | meta$sample_id == ""))
})

test_that("unified expression matrix aligns with unified metadata", {
  expr <- read.csv(file.path(PROCESSED_DIR, "unified_expression_log2cpm_qn.csv"), check.names = FALSE)
  meta <- read.csv(file.path(PROCESSED_DIR, "unified_metadata.csv"), stringsAsFactors = FALSE)
  sample_cols <- setdiff(colnames(expr), "gene_id")
  expect_equal(length(sample_cols), nrow(meta))
  expect_setequal(sample_cols, meta$sample_id)
  # quantile normalization => sample distributions should be much more similar
  # to each other than the raw (pre-normalization) distributions were. Count
  # data has many tied values (especially zeros), so rank-averaging ties keeps
  # a small residual spread -- this checks it's small, not exactly zero.
  qs <- sapply(expr[sample_cols[1:10]], function(x) quantile(x, c(0.25, 0.5, 0.75)))
  spread <- apply(qs, 1, function(r) diff(range(r)))
  expect_true(all(spread < 0.5))
})

test_that("QC correctly flags the deliberately-corrupted injected samples", {
  qc <- read.csv(file.path(OUTPUT_DIR, "qc_report_per_sample.csv"), stringsAsFactors = FALSE)
  truth_file <- file.path(RAW_DIR, "_ground_truth_corrupt_samples.txt")
  skip_if_not(file.exists(truth_file))
  corrupt_truth <- readLines(truth_file)
  flagged <- qc$sample_id[qc$qc_fail]
  expect_true(all(corrupt_truth %in% flagged))
})

test_that("QC concentrates flagged quality issues in ~8 of 16 batches, matching injected ground truth", {
  batch_qc <- read.csv(file.path(OUTPUT_DIR, "qc_report_per_batch.csv"), stringsAsFactors = FALSE)
  truth_file <- file.path(RAW_DIR, "_ground_truth_degraded_batches.txt")
  skip_if_not(file.exists(truth_file))
  degraded_truth <- readLines(truth_file)
  flagged_batches <- batch_qc$batch[batch_qc$batch_quality_issue]
  expect_equal(length(flagged_batches), 8)
  expect_setequal(flagged_batches, degraded_truth)
})

test_that("no sample has an impossible / non-finite normalized expression value", {
  expr <- read.csv(file.path(PROCESSED_DIR, "unified_expression_log2cpm_qn.csv"), check.names = FALSE)
  sample_cols <- setdiff(colnames(expr), "gene_id")
  vals <- as.matrix(expr[sample_cols])
  expect_true(all(is.finite(vals)))
  expect_true(all(vals >= 0))
})
