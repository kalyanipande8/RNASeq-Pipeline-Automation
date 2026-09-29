#!/usr/bin/env Rscript
# simulate_rnaseq_data.R
#
# Simulates a realistic (but entirely synthetic) RNA-seq compendium from
# patient-derived models: 320+ samples spread across 14 cohorts and 16
# sequencing batches, with:
#   - a gene x sample raw count matrix, 600 genes (kept modest so the whole
#     pipeline runs in seconds), with per-cohort biological signal and
#     per-batch technical (batch) effects layered on top of a negative
#     binomial count model
#   - a metadata table where each of the 16 batches was entered in ITS OWN
#     raw format (different column names, different date formats, different
#     ID conventions) -- this is what harmonize.R has to actually reconcile
#   - 8 of the 16 batches deliberately degraded (raised technical noise /
#     a few corrupted samples) to be picked up later by qc_checks.R
#
# Fixed seed (42) for full reproducibility.

set.seed(42)

PROJECT_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(PROJECT_DIR) == 0 || !nzchar(PROJECT_DIR)) PROJECT_DIR <- getwd()
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

N_GENES    <- 600
N_COHORTS  <- 14
N_BATCHES  <- 16
N_SAMPLES  <- 322   # 320+ samples, patient-derived models

cohort_ids <- sprintf("COH-%02d", seq_len(N_COHORTS))
batch_ids  <- sprintf("B%02d", seq_len(N_BATCHES))

# ---- 8 of 16 batches are deliberately degraded (QC target) -------------
degraded_batches <- sample(batch_ids, 8)

# ---- assign each sample a cohort and a batch ----------------------------
# Cohorts are patient cohorts; batches are sequencing runs. Each cohort's
# samples are sequenced across 1-3 batches (a cohort is not confined to one
# run), and each batch mixes samples from multiple cohorts, as in a real
# multi-cohort sequencing core.
sample_cohort <- sample(cohort_ids, N_SAMPLES, replace = TRUE,
                         prob = rep(1, N_COHORTS))
sample_batch  <- sample(batch_ids, N_SAMPLES, replace = TRUE,
                         prob = rep(1, N_BATCHES))

model_ids <- sprintf("PDX-%04d", seq_len(N_SAMPLES))

# ---- gene-level baseline mean expression (log-normal) --------------------
gene_baseline <- rlnorm(N_GENES, meanlog = 4, sdlog = 1.4)
gene_ids <- sprintf("GENE%04d", seq_len(N_GENES))

# ---- cohort effect: each cohort perturbs a random subset of genes --------
cohort_effect <- matrix(1, nrow = N_GENES, ncol = N_COHORTS,
                         dimnames = list(gene_ids, cohort_ids))
for (co in cohort_ids) {
  affected <- sample(N_GENES, size = round(N_GENES * 0.15))
  cohort_effect[affected, co] <- exp(rnorm(length(affected), mean = 0, sd = 0.6))
}

# ---- batch effect: multiplicative shift + dispersion inflation ----------
batch_mult <- setNames(exp(rnorm(N_BATCHES, mean = 0, sd = 0.35)), batch_ids)
batch_disp_inflate <- setNames(rep(1, N_BATCHES), batch_ids)
batch_disp_inflate[degraded_batches] <- runif(length(degraded_batches), 2.5, 5)

# ---- simulate counts: negative binomial per gene x sample ----------------
counts <- matrix(0L, nrow = N_GENES, ncol = N_SAMPLES,
                  dimnames = list(gene_ids, model_ids))

# library size per sample (varies by batch quality)
lib_size_factor <- numeric(N_SAMPLES)

for (i in seq_len(N_SAMPLES)) {
  co <- sample_cohort[i]
  ba <- sample_batch[i]
  mu <- gene_baseline * cohort_effect[, co] * batch_mult[ba]
  size_param <- 4 / batch_disp_inflate[ba]   # smaller size = more dispersion
  lf <- rlnorm(1, meanlog = log(1e6), sdlog = 0.25)
  if (ba %in% degraded_batches) lf <- lf * runif(1, 0.35, 0.8)  # lower depth
  lib_size_factor[i] <- lf
  scaled_mu <- mu * (lf / 1e6)
  counts[, i] <- rnbinom(N_GENES, mu = pmax(scaled_mu, 0.01), size = size_param)
}

# ---- inject a handful of deliberately corrupted samples (sanity check for
#      qc_checks.R / tests) : near-zero counts, i.e. failed library --------
n_corrupt <- 3
corrupt_idx <- sample(N_SAMPLES, n_corrupt)
for (idx in corrupt_idx) {
  counts[, idx] <- rpois(N_GENES, lambda = 0.5)
}

write.csv(data.frame(gene_id = gene_ids, counts, check.names = FALSE),
          file.path(RAW_DIR, "raw_counts.csv"), row.names = FALSE)

# ===========================================================================
# Metadata: each batch recorded in ITS OWN raw format. We split the 16
# batches into 4 "recording styles" (as if 4 different sequencing cores /
# LIMS systems fed this compendium), each with different column names, date
# formats, and ID conventions -- this is the real harmonization problem.
# ===========================================================================

batch_style <- setNames(rep(c("style_A", "style_B", "style_C", "style_D"),
                             length.out = N_BATCHES), batch_ids)

rand_date_in_2023 <- function(n) {
  as.Date("2023-01-01") + sample(0:364, n, replace = TRUE)
}

meta_rows <- vector("list", N_SAMPLES)
for (i in seq_len(N_SAMPLES)) {
  style <- batch_style[[sample_batch[i]]]
  d <- rand_date_in_2023(1)
  is_corrupt <- i %in% corrupt_idx
  row <- switch(style,
    style_A = data.frame(
      SampleID = model_ids[i],
      Cohort   = sample_cohort[i],
      SeqBatch = sample_batch[i],
      Run_Date = format(d, "%Y-%m-%d"),
      stringsAsFactors = FALSE
    ),
    style_B = data.frame(
      patient_model_id = model_ids[i],
      cohort_name       = tolower(sample_cohort[i]),
      batch             = sample_batch[i],
      date_sequenced    = format(d, "%m/%d/%Y"),
      stringsAsFactors = FALSE
    ),
    style_C = data.frame(
      `Model ID`   = model_ids[i],
      `Cohort ID`  = sub("COH-", "Cohort_", sample_cohort[i]),
      `Batch No.`  = sub("B", "Batch-", sample_batch[i]),
      `Date (DD-Mon-YYYY)` = format(d, "%d-%b-%Y"),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    style_D = data.frame(
      id             = model_ids[i],
      cohort         = sample_cohort[i],
      seq_run        = sample_batch[i],
      seq_date_epoch = as.numeric(d),  # days since epoch, a common raw LIMS export quirk
      stringsAsFactors = FALSE
    )
  )
  meta_rows[[i]] <- list(style = style, row = row)
}

# Write four separate raw metadata files, one per recording style, exactly
# as they would arrive from four different sources.
for (style in c("style_A", "style_B", "style_C", "style_D")) {
  rows <- lapply(meta_rows, function(x) if (x$style == style) x$row else NULL)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  df <- do.call(rbind, rows)
  write.csv(df, file.path(RAW_DIR, paste0("metadata_", style, ".csv")), row.names = FALSE)
}

# Also stash which samples were injected as corrupted / which batches were
# degraded, so qc_checks.R and tests can be scored against ground truth
# honestly (this file documents the injection, it is not read by harmonize.R).
writeLines(model_ids[corrupt_idx], file.path(RAW_DIR, "_ground_truth_corrupt_samples.txt"))
writeLines(degraded_batches, file.path(RAW_DIR, "_ground_truth_degraded_batches.txt"))

cat(sprintf("Simulated %d samples x %d genes across %d cohorts and %d batches\n",
            N_SAMPLES, N_GENES, N_COHORTS, N_BATCHES))
cat(sprintf("  Degraded batches (ground truth, %d of %d): %s\n",
            length(degraded_batches), N_BATCHES, paste(degraded_batches, collapse = ", ")))
cat(sprintf("  Deliberately corrupted samples (ground truth): %s\n",
            paste(model_ids[corrupt_idx], collapse = ", ")))
cat("  Wrote data/raw/raw_counts.csv and 4 raw metadata files (one per recording style)\n")
