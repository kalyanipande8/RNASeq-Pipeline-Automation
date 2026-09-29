#!/usr/bin/env Rscript
# make_figures.R -- two real figures generated from the actual pipeline output.

PROJECT_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(PROJECT_DIR) == 0 || !nzchar(PROJECT_DIR)) PROJECT_DIR <- getwd()
PROCESSED_DIR <- file.path(PROJECT_DIR, "data", "processed")
RAW_DIR <- file.path(PROJECT_DIR, "data", "raw")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
FIG_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

suppressWarnings(suppressMessages(library(ggplot2)))

qc_sample <- read.csv(file.path(OUTPUT_DIR, "qc_report_per_sample.csv"), stringsAsFactors = FALSE)
qc_batch  <- read.csv(file.path(OUTPUT_DIR, "qc_report_per_batch.csv"), stringsAsFactors = FALSE)

# ---- Figure 1: per-batch median sample-to-global correlation ------------
qc_batch$batch <- factor(qc_batch$batch, levels = qc_batch$batch[order(qc_batch$median_sample_to_global_cor)])
p1 <- ggplot(qc_batch, aes(x = batch, y = median_sample_to_global_cor, fill = batch_quality_issue)) +
  geom_col() +
  geom_hline(yintercept = 0.87, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(`TRUE` = "#d1495b", `FALSE` = "#66a182"),
                     labels = c(`TRUE` = "Quality issue", `FALSE` = "OK"),
                     name = "Batch status") +
  coord_flip() +
  labs(title = "Per-batch median sample-to-global expression correlation",
       subtitle = "Dashed line = QC threshold (0.87); 8 of 16 batches fall below it",
       x = "Sequencing batch", y = "Median Spearman correlation to global profile") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "fig1_batch_qc_correlation.png"), p1, width = 7, height = 6, dpi = 130)

# ---- Figure 2: library size distribution before/after, flagged samples --
qc_sample$status <- ifelse(qc_sample$qc_fail, "QC FAIL", "Pass")
p2 <- ggplot(qc_sample, aes(x = reorder(batch, library_size, FUN = median), y = library_size, color = status)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.4) +
  scale_y_log10() +
  scale_color_manual(values = c("Pass" = "#4c72b0", "QC FAIL" = "#d1495b")) +
  labs(title = "Raw library size per sample, by batch",
       subtitle = sprintf("%d of %d samples flagged qc_fail", sum(qc_sample$qc_fail), nrow(qc_sample)),
       x = "Sequencing batch (ordered by median library size)", y = "Library size (log10 scale)",
       color = "QC status") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "fig2_library_size_by_batch.png"), p2, width = 8, height = 5, dpi = 130)

cat(sprintf("Wrote %s\n", file.path(FIG_DIR, "fig1_batch_qc_correlation.png")))
cat(sprintf("Wrote %s\n", file.path(FIG_DIR, "fig2_library_size_by_batch.png")))
