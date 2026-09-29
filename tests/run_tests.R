#!/usr/bin/env Rscript
# tests/run_tests.R -- convenience wrapper: run the full pipeline, then the tests.
library(testthat)
PROJECT_DIR <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
test_file(file.path(PROJECT_DIR, "tests", "test_pipeline.R"))
