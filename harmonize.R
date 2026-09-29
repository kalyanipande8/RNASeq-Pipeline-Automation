#!/usr/bin/env Rscript
# harmonize.R
#
# Real harmonization step:
#   1. Metadata: reads the 4 raw metadata files (each in its own raw format --
#      different column names, ID/cohort/batch spellings, and date formats)
#      and standardizes them into ONE schema:
#        sample_id, cohort, batch, seq_date (ISO 8601 Date)
#   2. Expression: reads the raw gene x sample count matrix and normalizes it
#      across sequencing batches using quantile normalization (a standard,
#      real RNA-seq cross-batch normalization method) after a
#      library-size (CPM) normalization step, so that per-sample sequencing
#      depth AND batch-to-batch distributional shifts are both corrected.
#
# Outputs one unified expression matrix and one unified metadata table.

PROJECT_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(PROJECT_DIR) == 0 || !nzchar(PROJECT_DIR)) PROJECT_DIR <- getwd()
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Metadata harmonization
# ---------------------------------------------------------------------------

std_cohort <- function(x) {
  x <- trimws(x)
  x <- gsub("^Cohort_", "COH-", x, ignore.case = TRUE)
  x <- toupper(x)
  # ensure COH-NN zero-padded form
  ifelse(grepl("^COH-\\d+$", x),
         sprintf("COH-%02d", as.integer(sub("COH-", "", x))),
         x)
}

std_batch <- function(x) {
  x <- trimws(x)
  x <- gsub("^Batch-", "B", x, ignore.case = TRUE)
  x <- toupper(x)
  ifelse(grepl("^B\\d+$", x),
         sprintf("B%02d", as.integer(sub("B", "", x))),
         x)
}

meta_A <- read.csv(file.path(RAW_DIR, "metadata_style_A.csv"), stringsAsFactors = FALSE)
meta_A_std <- data.frame(
  sample_id = meta_A$SampleID,
  cohort    = std_cohort(meta_A$Cohort),
  batch     = std_batch(meta_A$SeqBatch),
  seq_date  = as.Date(meta_A$Run_Date, format = "%Y-%m-%d"),
  raw_format = "style_A"
)

meta_B <- read.csv(file.path(RAW_DIR, "metadata_style_B.csv"), stringsAsFactors = FALSE)
meta_B_std <- data.frame(
  sample_id = meta_B$patient_model_id,
  cohort    = std_cohort(meta_B$cohort_name),
  batch     = std_batch(meta_B$batch),
  seq_date  = as.Date(meta_B$date_sequenced, format = "%m/%d/%Y"),
  raw_format = "style_B"
)

meta_C <- read.csv(file.path(RAW_DIR, "metadata_style_C.csv"), stringsAsFactors = FALSE,
                    check.names = FALSE)
meta_C_std <- data.frame(
  sample_id = meta_C[["Model ID"]],
  cohort    = std_cohort(meta_C[["Cohort ID"]]),
  batch     = std_batch(meta_C[["Batch No."]]),
  seq_date  = as.Date(meta_C[["Date (DD-Mon-YYYY)"]], format = "%d-%b-%Y"),
  raw_format = "style_C"
)

meta_D <- read.csv(file.path(RAW_DIR, "metadata_style_D.csv"), stringsAsFactors = FALSE)
meta_D_std <- data.frame(
  sample_id = meta_D$id,
  cohort    = std_cohort(meta_D$cohort),
  batch     = std_batch(meta_D$seq_run),
  seq_date  = as.Date(meta_D$seq_date_epoch, origin = "1970-01-01"),
  raw_format = "style_D"
)

unified_meta <- rbind(meta_A_std, meta_B_std, meta_C_std, meta_D_std)
unified_meta <- unified_meta[order(unified_meta$sample_id), ]
rownames(unified_meta) <- NULL

stopifnot(all(!is.na(unified_meta$seq_date)))
stopifnot(all(grepl("^COH-\\d{2}$", unified_meta$cohort)))
stopifnot(all(grepl("^B\\d{2}$", unified_meta$batch)))

write.csv(unified_meta, file.path(PROCESSED_DIR, "unified_metadata.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Expression harmonization: CPM + quantile normalization across batches
# ---------------------------------------------------------------------------

raw <- read.csv(file.path(RAW_DIR, "raw_counts.csv"), check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- raw$gene_id
counts <- as.matrix(raw[, -1])
rownames(counts) <- gene_ids
storage.mode(counts) <- "double"

# Library-size (CPM) normalization: correct for sequencing depth per sample
lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6

# log2(CPM + 1) before quantile normalization (standard for count data)
log_cpm <- log2(cpm + 1)

# Quantile normalization across all samples (removes batch-driven
# distributional shifts): every sample's values are replaced by the
# across-sample-average of that rank, then re-mapped back to each sample's
# original rank order. This is the classic quantile-normalization algorithm
# (no external package required -- implemented directly here).
quantile_normalize <- function(mat) {
  ranked_mean <- apply(mat, 2, sort)
  ranked_mean <- rowMeans(ranked_mean)
  out <- apply(mat, 2, function(col) {
    ranked_mean[rank(col, ties.method = "average")]
  })
  dimnames(out) <- dimnames(mat)
  out
}

qn <- quantile_normalize(log_cpm)

unified_expr <- data.frame(gene_id = gene_ids, qn, check.names = FALSE)
write.csv(unified_expr, file.path(PROCESSED_DIR, "unified_expression_log2cpm_qn.csv"), row.names = FALSE)

cat(sprintf("Harmonized metadata for %d samples across %d cohorts (from 4 raw formats) into one schema\n",
            nrow(unified_meta), length(unique(unified_meta$cohort))))
cat(sprintf("Harmonized expression: %d genes x %d samples -> CPM + quantile-normalized (log2 scale)\n",
            nrow(qn), ncol(qn)))
cat(sprintf("  Wrote %s\n", file.path(PROCESSED_DIR, "unified_metadata.csv")))
cat(sprintf("  Wrote %s\n", file.path(PROCESSED_DIR, "unified_expression_log2cpm_qn.csv")))
