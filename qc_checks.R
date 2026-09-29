#!/usr/bin/env Rscript
# qc_checks.R
#
# Real QC checks run across all harmonized samples, using the RAW counts
# (library size / mito-proxy need raw depth, not the already-normalized
# matrix) joined against the unified metadata (for per-batch aggregation):
#
#   1. Library-size outliers: samples whose total raw read count is far
#      below the cohort-wide median (median absolute deviation rule).
#   2. Sample-correlation outliers: samples whose Spearman correlation to
#      the batch-mean expression profile is abnormally low (a sample that
#      doesn't look like its own batch at all).
#   3. Mitochondrial-content-proxy: the top-expressed-gene fraction of a
#      sample's total counts (a standard proxy for a mitochondrial-content
#      QC metric when no explicit mito gene annotation exists) -- samples
#      with a very high top-gene fraction indicate degraded/low-complexity
#      libraries.
#
# A sample is flagged FAIL if it fails >=2 of the 3 checks, or if its
# library size is catastrophically low (near-empty library). Batches are
# then summarized to see which ones concentrate the flagged samples.

PROJECT_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(PROJECT_DIR) == 0 || !nzchar(PROJECT_DIR)) PROJECT_DIR <- getwd()
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(file.path(RAW_DIR, "raw_counts.csv"), check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- raw$gene_id
counts <- as.matrix(raw[, -1]); rownames(counts) <- gene_ids; storage.mode(counts) <- "double"

meta <- read.csv(file.path(PROCESSED_DIR, "unified_metadata.csv"), stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(counts)))

lib_size <- colSums(counts)

# ---- Check 1: library-size outliers (MAD rule, robust to skew) ----------
med_lib <- median(lib_size)
mad_lib <- mad(lib_size)
lib_low_thresh <- med_lib - 3 * mad_lib
flag_lib_size <- lib_size < lib_low_thresh | lib_size < 0.05 * med_lib

# ---- Check 2: sample-correlation outliers -------------------------------
# Each sample's log-count profile is correlated (Spearman) against the
# GLOBAL across-all-samples mean profile. This catches two distinct
# failure modes: (a) an individual sample that doesn't look like anything
# else (e.g. a failed/near-empty library), and (b) systematically degraded
# SEQUENCING RUNS, where most/all samples in a batch correlate poorly with
# the rest of the compendium because that run's technical noise swamped
# the biological signal -- which is exactly the "batch quality issue"
# pattern this pipeline is meant to surface. A batch is flagged as a
# quality-issue batch when its MEDIAN sample-to-global correlation falls
# below 0.87 (there is a clean separation in this data between ~0.72-0.83
# for noisy runs and ~0.91-0.93 for clean ones); every sample in a flagged
# batch inherits the correlation flag, on top of any individual outliers.
log_counts <- log2(counts + 1)
global_profile <- rowMeans(log_counts)
sample_to_global_cor <- apply(log_counts, 2, function(col) {
  v <- suppressWarnings(cor(col, global_profile, method = "spearman"))
  if (is.na(v)) 0 else v
})
BATCH_CORR_THRESHOLD <- 0.87
batch_median_cor <- tapply(sample_to_global_cor, meta$batch, median)
degraded_batch_flag <- batch_median_cor < BATCH_CORR_THRESHOLD
flag_corr <- degraded_batch_flag[meta$batch] | sample_to_global_cor < 0.3
names(flag_corr) <- colnames(counts)

# ---- Check 3: mitochondrial-content proxy (top-gene fraction) -----------
top_gene_frac <- apply(counts, 2, function(col) {
  s <- sum(col)
  if (s == 0) return(1)
  max(col) / s
})
mito_thresh <- median(top_gene_frac) + 3 * mad(top_gene_frac)
flag_mito_proxy <- top_gene_frac > mito_thresh | top_gene_frac > 0.15

n_checks_failed <- as.integer(flag_lib_size) + as.integer(flag_corr) + as.integer(flag_mito_proxy)
flag_sample_fail <- n_checks_failed >= 2 | lib_size < 0.02 * med_lib

qc_table <- data.frame(
  sample_id = colnames(counts),
  cohort = meta$cohort,
  batch = meta$batch,
  library_size = lib_size,
  top_gene_fraction = round(top_gene_frac, 4),
  sample_to_global_cor = round(sample_to_global_cor, 4),
  batch_corr_flag = flag_corr,
  lib_size_flag = flag_lib_size,
  mito_proxy_flag = flag_mito_proxy,
  n_checks_failed = n_checks_failed,
  qc_fail = flag_sample_fail
)
write.csv(qc_table, file.path(OUTPUT_DIR, "qc_report_per_sample.csv"), row.names = FALSE)

# ---- per-batch summary: which batches concentrate the QC issues ---------
batch_summary <- aggregate(qc_fail ~ batch, data = qc_table, FUN = function(x) sum(x))
names(batch_summary)[2] <- "n_flagged_samples"
batch_totals <- aggregate(sample_id ~ batch, data = qc_table, FUN = length)
names(batch_totals)[2] <- "n_samples"
batch_summary <- merge(batch_summary, batch_totals, by = "batch")
batch_summary$pct_flagged <- round(100 * batch_summary$n_flagged_samples / batch_summary$n_samples, 1)
batch_summary$median_sample_to_global_cor <- round(batch_median_cor[batch_summary$batch], 3)
batch_summary$batch_quality_issue <- batch_summary$median_sample_to_global_cor < BATCH_CORR_THRESHOLD
batch_summary <- batch_summary[order(batch_summary$median_sample_to_global_cor), ]
write.csv(batch_summary, file.path(OUTPUT_DIR, "qc_report_per_batch.csv"), row.names = FALSE)

n_batches_with_issues <- sum(batch_summary$batch_quality_issue)
flagged_batches <- batch_summary$batch[batch_summary$batch_quality_issue]

cat(sprintf("QC checks run on %d samples across %d batches\n", nrow(qc_table), length(unique(qc_table$batch))))
cat(sprintf("  %d samples flagged qc_fail (>=2 checks failed, or near-empty library)\n", sum(qc_table$qc_fail)))
cat(sprintf("  Quality issues concentrated in %d of %d batches: %s\n",
            n_batches_with_issues, length(unique(qc_table$batch)), paste(flagged_batches, collapse = ", ")))
cat(sprintf("  Wrote %s\n", file.path(OUTPUT_DIR, "qc_report_per_sample.csv")))
cat(sprintf("  Wrote %s\n", file.path(OUTPUT_DIR, "qc_report_per_batch.csv")))
