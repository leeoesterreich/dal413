#!/usr/bin/env Rscript

# Load required libraries
library(GSA)
library(glmnet)
library(caret)

# Source helper functions
source("Rscript/helper.R")

# Set input data paths
CN_data_path <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/DNA-based-predictors-of-non-genetic-cancer-phenotypes/Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt"
RNA_data_dir <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/DNA-based-predictors-of-non-genetic-cancer-phenotypes/HiSeqV2"

# Read copy number data
CN_data <- read.table(CN_data_path, header=TRUE, sep="\t", check.names=FALSE)

# Read RNA expression data
RNA_files <- list.files(RNA_data_dir, pattern=".txt$", full.names=TRUE)
RNA_data <- do.call(rbind, lapply(RNA_files, function(f) {
  read.table(f, header=TRUE, sep="\t", check.names=FALSE)
}))

# Calculate signature scores
source("Rscript/signature_score_and_segment_score_calculation.R")

# Calculate segment scores
segment_score <- calc_segments(CN_data, "CNA_segments.gmt", method="mean")
save(segment_score, file="segment_score.rda")

# Run association tests
source("Rscript/association_test.R")

# Build elastic net models
source("Rscript/Elastic_Net_modeling.R")

print("Analysis complete. Check the output files for results.") 