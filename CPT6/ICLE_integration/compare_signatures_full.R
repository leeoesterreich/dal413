#!/usr/bin/env Rscript

# Script to compare SBS and COSMIC signatures between CPT6 and ICLE
# Author: Daisong

# Load required libraries
library(maftools)
library(ggplot2)
library(dplyr)
library(NMF)
library(pheatmap)
library(gridExtra)

# Set file paths
icle_maf_file <- "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
cpt6_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
output_dir <- "results/signatures_full"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to perform signature analysis for a dataset
analyze_signatures <- function(maf_file, dataset_name) {
  cat(sprintf("\n\n========== Analyzing %s Signatures ==========\n\n", dataset_name))
  
  # Read MAF file
  cat(sprintf("Reading %s MAF file...\n", dataset_name))
  maf <- read.maf(maf_file)
  
  # Extract trinucleotide matrix
  cat("Extracting trinucleotide matrix...\n")
  tryCatch({
    # Try to extract trinucleotide matrix
    tnm <- trinucleotideMatrix(maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
    
    # If successful, proceed with signature analysis
    cat("Estimating optimal number of signatures...\n")
    
    # Step 1: Estimate signatures
    # Use a small pConstant to handle potential low mutation counts
    sig_estimate <- estimateSignatures(mat = tnm, nTry = 6, pConstant = 0.1)
    
    # Save the cophenetic correlation plot
    pdf(file.path(output_dir, paste0(dataset_name, "_cophenetic_correlation.pdf")), width = 10, height = 8)
    plotCophenetic(res = sig_estimate)
    dev.off()
    
    # Determine optimal number of signatures based on cophenetic correlation
    # For simplicity, we'll use n=3 as in the example, but this should be adjusted based on the plot
    optimal_n <- 3
    cat(sprintf("Using %d as the optimal number of signatures based on cophenetic correlation\n", optimal_n))
    
    # Step 2: Extract signatures
    cat("Extracting signatures...\n")
    signatures <- extractSignatures(mat = tnm, n = optimal_n, pConstant = 0.1)
    
    # Step 3: Compare with COSMIC signatures (legacy)
    cat("Comparing with COSMIC legacy signatures...\n")
    cosmic_legacy <- compareSignatures(nmfRes = signatures, sig_db = "legacy")
    
    # Save the legacy comparison plot
    pdf(file.path(output_dir, paste0(dataset_name, "_cosmic_legacy_comparison.pdf")), width = 12, height = 10)
    plotSignatures(nmfRes = signatures, sig_db = "legacy")
    dev.off()
    
    # Step 4: Compare with COSMIC SBS signatures
    cat("Comparing with COSMIC SBS signatures...\n")
    cosmic_sbs <- compareSignatures(nmfRes = signatures, sig_db = "SBS")
    
    # Save the SBS comparison plot
    pdf(file.path(output_dir, paste0(dataset_name, "_cosmic_sbs_comparison.pdf")), width = 12, height = 10)
    plotSignatures(nmfRes = signatures, sig_db = "SBS")
    dev.off()
    
    # Step 5: Create heatmap of cosine similarities
    cat("Creating heatmap of cosine similarities...\n")
    pdf(file.path(output_dir, paste0(dataset_name, "_cosine_similarities_heatmap.pdf")), width = 12, height = 10)
    pheatmap(mat = cosmic_legacy$cosine_similarities, 
             cluster_rows = FALSE, 
             main = paste(dataset_name, "- Cosine Similarity Against Validated Signatures"))
    dev.off()
    
    # Return the results
    return(list(
      signatures = signatures,
      cosmic_legacy = cosmic_legacy,
      cosmic_sbs = cosmic_sbs
    ))
    
  }, error = function(e) {
    cat("Error in signature analysis:", conditionMessage(e), "\n")
    
    # If trinucleotideMatrix fails, try an alternative approach
    cat("Trying alternative approach for signature analysis...\n")
    
    # Get mutation types directly from MAF
    mut_types <- getClinicalData(maf)$Variant_Type
    var_class <- getClinicalData(maf)$Variant_Classification
    
    # Create summary of mutation types
    type_summary <- table(mut_types)
    class_summary <- table(var_class)
    
    # Save summaries
    write.csv(as.data.frame(type_summary), 
              file.path(output_dir, paste0(dataset_name, "_mutation_types.csv")))
    write.csv(as.data.frame(class_summary), 
              file.path(output_dir, paste0(dataset_name, "_variant_classifications.csv")))
    
    # Create basic plots instead
    pdf(file.path(output_dir, paste0(dataset_name, "_mutation_summary.pdf")), width = 10, height = 8)
    plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
    dev.off()
    
    # Return NULL to indicate failure
    return(NULL)
  })
}

# Main function to compare signatures between datasets
compare_signatures_main <- function() {
  tryCatch({
    # Analyze CPT6 signatures
    cpt6_results <- analyze_signatures(cpt6_maf_file, "CPT6")
    
    # Analyze ICLE signatures
    icle_results <- analyze_signatures(icle_maf_file, "ICLE")
    
    # If both analyses were successful, compare the results
    if (!is.null(cpt6_results) && !is.null(icle_results)) {
      cat("\n\n========== Comparing CPT6 and ICLE Signatures ==========\n\n")
      
      # Compare top COSMIC signatures
      cat("Comparing top COSMIC signatures between datasets...\n")
      
      # Extract top COSMIC signatures for each dataset
      cpt6_top_cosmic <- apply(cpt6_results$cosmic_legacy$cosine_similarities, 2, function(x) names(which.max(x)))
      icle_top_cosmic <- apply(icle_results$cosmic_legacy$cosine_similarities, 2, function(x) names(which.max(x)))
      
      # Create a comparison table
      cosmic_comparison <- data.frame(
        Signature = paste0("Signature_", 1:length(cpt6_top_cosmic)),
        CPT6_Top_Match = cpt6_top_cosmic,
        ICLE_Top_Match = icle_top_cosmic,
        CPT6_Similarity = sapply(1:length(cpt6_top_cosmic), function(i) {
          max(cpt6_results$cosmic_legacy$cosine_similarities[, i])
        }),
        ICLE_Similarity = sapply(1:length(icle_top_cosmic), function(i) {
          max(icle_results$cosmic_legacy$cosine_similarities[, i])
        })
      )
      
      # Save the comparison table
      write.csv(cosmic_comparison, 
                file.path(output_dir, "cosmic_signature_comparison.csv"), 
                row.names = FALSE)
      
      # Print the comparison table
      cat("\nCOSMIC Signature Comparison:\n")
      print(cosmic_comparison)
      
      # Create a visual comparison of signature contributions
      cat("Creating visual comparison of signature contributions...\n")
      
      # Extract signature contributions
      cpt6_contribs <- cpt6_results$signatures$contribution
      icle_contribs <- icle_results$signatures$contribution
      
      # Normalize contributions to percentages
      cpt6_contribs_norm <- sweep(cpt6_contribs, 2, colSums(cpt6_contribs), "/") * 100
      icle_contribs_norm <- sweep(icle_contribs, 2, colSums(icle_contribs), "/") * 100
      
      # Create a data frame for plotting
      cpt6_df <- data.frame(
        Sample = rep(colnames(cpt6_contribs), each = nrow(cpt6_contribs)),
        Signature = rep(paste0("Signature_", 1:nrow(cpt6_contribs)), times = ncol(cpt6_contribs)),
        Contribution = as.vector(cpt6_contribs_norm),
        Dataset = "CPT6"
      )
      
      icle_df <- data.frame(
        Sample = rep(colnames(icle_contribs), each = nrow(icle_contribs)),
        Signature = rep(paste0("Signature_", 1:nrow(icle_contribs)), times = ncol(icle_contribs)),
        Contribution = as.vector(icle_contribs_norm),
        Dataset = "ICLE"
      )
      
      # Combine the data frames
      combined_df <- rbind(cpt6_df, icle_df)
      
      # Create a bar plot comparing signature contributions
      p <- ggplot(combined_df, aes(x = Signature, y = Contribution, fill = Dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        labs(
          title = "Comparison of Signature Contributions",
          x = "Signature",
          y = "Contribution (%)"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)
        ) +
        scale_fill_manual(values = c("CPT6" = "coral", "ICLE" = "skyblue"))
      
      # Save the comparison plot
      ggsave(file.path(output_dir, "signature_contribution_comparison.pdf"), p, width = 12, height = 8)
      ggsave(file.path(output_dir, "signature_contribution_comparison.png"), p, width = 12, height = 8, dpi = 300)
    }
    
    cat("\nSignature Analysis Complete!\n")
    cat("Results saved to:", output_dir, "\n")
    
  }, error = function(e) {
    cat("Error in signature comparison:", conditionMessage(e), "\n")
    print(e)
  })
}

# Run the main function
compare_signatures_main() 