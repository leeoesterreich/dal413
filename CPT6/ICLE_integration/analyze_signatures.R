#!/usr/bin/env Rscript

# Script to analyze mutational signatures using maftools in R

# Load required libraries
suppressPackageStartupMessages({
  library(maftools)
  library(ggplot2)
})

# Set working directory
setwd("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/")

# Create output directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# File paths
merged_maf_file <- "results/merged_maf.maf"
output_dir <- "results"

# Check if merged MAF file exists
if (!file.exists(merged_maf_file)) {
  cat("Error: Merged MAF file not found at", merged_maf_file, "\n")
  cat("Please run the Python script first to generate the merged MAF file.\n")
  quit(status = 1)
}

# Read the merged MAF file
cat("Reading merged MAF file...\n")
tryCatch({
  maf_data <- read.maf(maf = merged_maf_file)
}, error = function(e) {
  cat("Error reading MAF file:", e$message, "\n")
  quit(status = 1)
})

# Summarize the MAF file
cat("Generating MAF summary...\n")
maf_summary <- getSampleSummary(maf_data)
write.csv(maf_summary, file.path(output_dir, "maf_summary.csv"), row.names = FALSE)

# Generate oncoplot for the top 20 genes
cat("Generating oncoplot...\n")
pdf(file.path(output_dir, "oncoplot_top20.pdf"), width = 12, height = 8)
oncoplot(maf = maf_data, top = 20)
dev.off()

# Extract and analyze mutational signatures
cat("Extracting mutational signatures...\n")
# Extract signatures using NMF
tryCatch({
  mat <- trinucleotideMatrix(maf_data, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
  signatures <- extractSignatures(mat = mat, nTry = 6, plotBestFitRes = TRUE)
  
  # Compare with COSMIC signatures
  cat("Comparing with COSMIC signatures...\n")
  cosmic_signatures <- getSignatures(reference = "cosmic_v3.2")
  comparison <- compareSignatures(signatures, cosmic_signatures)
  
  # Plot the signatures
  cat("Plotting signatures...\n")
  pdf(file.path(output_dir, "mutational_signatures.pdf"), width = 12, height = 8)
  plotSignatures(signatures, title_size = 0.8, cosmic_sigs = comparison)
  dev.off()
}, error = function(e) {
  cat("Error in signature analysis:", e$message, "\n")
  cat("Continuing with other analyses...\n")
})

# Generate rainfall plot for each sample
cat("Generating rainfall plots...\n")
samples <- unique(maf_data@data$Tumor_Sample_Barcode)
for (sample in samples[1:min(5, length(samples))]) {  # Limit to first 5 samples to avoid too many plots
  tryCatch({
    pdf(file.path(output_dir, paste0("rainfall_", sample, ".pdf")), width = 10, height = 6)
    rainfallPlot(maf = maf_data, detectChangePoints = TRUE, pointSize = 0.4, 
                 tsb = sample)
    dev.off()
  }, error = function(e) {
    cat("Error generating rainfall plot for sample", sample, ":", e$message, "\n")
  })
}

# Generate lollipop plots for key oncogenes from the list
cat("Generating lollipop plots for key oncogenes...\n")
if (file.exists("Oncogene_CPT6.csv")) {
  oncogenes <- read.csv("Oncogene_CPT6.csv", header = FALSE)$V1
  
  for (gene in oncogenes[1:min(10, length(oncogenes))]) {  # Limit to first 10 genes to avoid too many plots
    tryCatch({
      if (gene %in% maf_data@data$Hugo_Symbol) {
        pdf(file.path(output_dir, paste0("lollipop_", gene, ".pdf")), width = 10, height = 6)
        lollipopPlot(maf = maf_data, gene = gene, showMutationRate = TRUE)
        dev.off()
      }
    }, error = function(e) {
      cat("Error generating lollipop plot for gene", gene, ":", e$message, "\n")
    })
  }
} else {
  cat("Warning: Oncogene list file not found. Skipping lollipop plots.\n")
}

# Compare mutation load across samples
cat("Comparing mutation load across samples...\n")
pdf(file.path(output_dir, "mutation_load_comparison.pdf"), width = 10, height = 6)
plotmafSummary(maf = maf_data, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

cat("Analysis complete. Results saved to the 'results' directory.\n") 