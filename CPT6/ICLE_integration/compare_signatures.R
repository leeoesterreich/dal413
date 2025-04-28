#!/usr/bin/env Rscript

# Script to compare mutational signatures (SBS and COSMIC) between CPT6 and ICLE
# Author: Daisong

# Load required libraries
library(maftools)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome

# Set file paths
icle_maf_file <- "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
cpt6_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
output_dir <- "results/signatures"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to read MAF files
read_maf_data <- function(maf_file, source_name) {
  cat(sprintf("Reading %s MAF file...\n", source_name))
  maf <- read.maf(maf_file)
  return(maf)
}

# Function to extract trinucleotide context
extract_tricontext <- function(maf, source_name) {
  cat(sprintf("Extracting trinucleotide context for %s...\n", source_name))
  # Extract trinucleotide matrix
  trmat <- trinucleotideMatrix(maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
  return(trmat)
}

# Function to extract signatures using deconstructSigs
extract_signatures <- function(maf, source_name) {
  cat(sprintf("Extracting mutational signatures for %s...\n", source_name))
  
  # Extract signatures using maftools' built-in function
  # This uses NMF (Non-negative Matrix Factorization)
  sig_result <- extractSignatures(
    mat = maf, 
    nTry = 6,  # Try up to 6 signatures
    plotBestFitRes = TRUE,
    outDir = file.path(output_dir, paste0(source_name, "_signatures"))
  )
  
  return(sig_result)
}

# Function to compare with COSMIC signatures
compare_with_cosmic <- function(sig_result, source_name) {
  cat(sprintf("Comparing %s signatures with COSMIC signatures...\n", source_name))
  
  # Compare extracted signatures with COSMIC signatures
  cosmic_comparison <- compareSignatures(
    sig_result, 
    sig_db = "legacy", # Use legacy COSMIC signatures (v2)
    plotdata = TRUE
  )
  
  # Save the comparison plot
  pdf(file.path(output_dir, paste0(source_name, "_cosmic_comparison.pdf")), width = 12, height = 8)
  plotSignatures(cosmic_comparison)
  dev.off()
  
  return(cosmic_comparison)
}

# Function to perform signature enrichment analysis
signature_enrichment <- function(maf, source_name) {
  cat(sprintf("Performing signature enrichment analysis for %s...\n", source_name))
  
  # Perform signature enrichment analysis
  sig_enrich <- signatureEnrichment(
    maf = maf,
    minMut = 5, # Minimum mutations required
    useCNV = FALSE,
    plotEnrichment = TRUE
  )
  
  # Save the enrichment plot
  pdf(file.path(output_dir, paste0(source_name, "_signature_enrichment.pdf")), width = 12, height = 8)
  plotEnrichmentResults(sig_enrich)
  dev.off()
  
  return(sig_enrich)
}

# Function to perform SBS analysis using SigProfiler
sbs_analysis <- function(maf, source_name) {
  cat(sprintf("Performing SBS analysis for %s...\n", source_name))
  
  # Extract SBS (Single Base Substitution) signatures
  sbs_result <- plotMafSummary(
    maf = maf,
    rmOutlier = TRUE,
    addStat = 'median',
    dashboard = TRUE,
    titvRaw = TRUE
  )
  
  # Save the SBS plot
  pdf(file.path(output_dir, paste0(source_name, "_sbs_summary.pdf")), width = 12, height = 8)
  print(sbs_result)
  dev.off()
  
  return(sbs_result)
}

# Function to compare signature contributions between samples
compare_signature_contributions <- function(icle_cosmic, cpt6_cosmic) {
  cat("Comparing signature contributions between CPT6 and ICLE...\n")
  
  # Extract signature contributions
  icle_contribs <- icle_cosmic$cosine$contribution
  cpt6_contribs <- cpt6_cosmic$cosine$contribution
  
  # Combine contributions
  combined_contribs <- data.frame(
    Signature = rownames(icle_contribs),
    ICLE = icle_contribs[,1],
    CPT6 = cpt6_contribs[,1]
  )
  
  # Create a long format for plotting
  contribs_long <- tidyr::pivot_longer(
    combined_contribs,
    cols = c("ICLE", "CPT6"),
    names_to = "Source",
    values_to = "Contribution"
  )
  
  # Create a bar plot comparing contributions
  p <- ggplot(contribs_long, aes(x = Signature, y = Contribution, fill = Source)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Comparison of Mutational Signature Contributions",
      x = "COSMIC Signature",
      y = "Relative Contribution"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("ICLE" = "skyblue", "CPT6" = "coral"))
  
  # Save the comparison plot
  ggsave(file.path(output_dir, "signature_contribution_comparison.pdf"), p, width = 12, height = 8)
  ggsave(file.path(output_dir, "signature_contribution_comparison.png"), p, width = 12, height = 8, dpi = 300)
  
  return(p)
}

# Main function to run the signature comparison
compare_signatures_main <- function() {
  tryCatch({
    # Read MAF files
    icle_maf <- read_maf_data(icle_maf_file, "ICLE")
    cpt6_maf <- read_maf_data(cpt6_maf_file, "CPT6")
    
    # Extract trinucleotide context
    icle_tricontext <- extract_tricontext(icle_maf, "ICLE")
    cpt6_tricontext <- extract_tricontext(cpt6_maf, "CPT6")
    
    # Extract signatures
    icle_signatures <- extract_signatures(icle_maf, "ICLE")
    cpt6_signatures <- extract_signatures(cpt6_maf, "CPT6")
    
    # Compare with COSMIC signatures
    icle_cosmic <- compare_with_cosmic(icle_signatures, "ICLE")
    cpt6_cosmic <- compare_with_cosmic(cpt6_signatures, "CPT6")
    
    # Perform signature enrichment analysis
    icle_enrich <- signature_enrichment(icle_maf, "ICLE")
    cpt6_enrich <- signature_enrichment(cpt6_maf, "CPT6")
    
    # Perform SBS analysis
    icle_sbs <- sbs_analysis(icle_maf, "ICLE")
    cpt6_sbs <- sbs_analysis(cpt6_maf, "CPT6")
    
    # Compare signature contributions
    contrib_comparison <- compare_signature_contributions(icle_cosmic, cpt6_cosmic)
    
    # Generate a combined report
    cat("\nSignature Analysis Complete!\n")
    cat("Results saved to:", output_dir, "\n")
    
    # Print summary of top signatures for each dataset
    cat("\nTop COSMIC Signatures in ICLE:\n")
    print(head(sort(icle_cosmic$cosine$contribution[,1], decreasing = TRUE), 5))
    
    cat("\nTop COSMIC Signatures in CPT6:\n")
    print(head(sort(cpt6_cosmic$cosine$contribution[,1], decreasing = TRUE), 5))
    
  }, error = function(e) {
    cat("Error in signature comparison:", conditionMessage(e), "\n")
  })
}

# Run the main function
compare_signatures_main() 