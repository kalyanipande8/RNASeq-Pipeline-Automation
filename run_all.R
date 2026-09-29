#!/usr/bin/env Rscript
# run_all.R
#
# Orchestrates the full pipeline end-to-end:
#   1. simulate_rnaseq_data.R -> data/raw/  (raw counts + 4 raw-format metadata files)
#   2. harmonize.R            -> data/processed/ (unified expression + metadata)
#   3. qc_checks.R            -> output/ (per-sample + per-batch QC reports)
#   4. make_figures.R         -> output/figures/ (2 real figures from the actual run)
#
# Run with:  Rscript run_all.R

PROJECT_DIR <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (length(PROJECT_DIR) == 0 || !nzchar(PROJECT_DIR)) PROJECT_DIR <- getwd()
setwd(PROJECT_DIR)

run_step <- function(label, script) {
  cat(sprintf("\n=== %s ===\n", label))
  source(script, local = new.env())
}

run_step("1/4 Simulating RNA-seq compendium", "simulate_rnaseq_data.R")
run_step("2/4 Harmonizing metadata + expression", "harmonize.R")
run_step("3/4 Running QC checks", "qc_checks.R")
run_step("4/4 Generating figures", "make_figures.R")

cat("\n=== Pipeline complete ===\n")
cat("See data/processed/ for the unified matrix + metadata, and output/ for QC reports and figures.\n")
