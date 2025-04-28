#!/usr/bin/env Rscript

# Script to calculate tumor mutation burden using maftools
# Comparing CPT6, ICLE, and TCGA_ILC samples

# Load required libraries
library(maftools)
library(ggplot2)
library(dplyr)

# Set file paths
icle_maf_file <- "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
cpt6_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
tcga_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/TCGA_ILC_BRCA_MAF.maf"
output_dir <- "results"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to calculate TMB using maftools' built-in function
calculate_tmb <- function(maf_file, source_name) {
  cat(sprintf("Calculating TMB for %s...\n", source_name))
  
  # Read MAF file
  maf <- read.maf(maf_file)
  
  # Calculate TMB using maftools' built-in function
  tmb_result <- tmb(maf = maf)
  
  # Add source column
  tmb_result$Source <- source_name
  
  return(tmb_result)
}

# Main function to generate TMB comparison plot
generate_tmb_comparison <- function() {
  # Calculate TMB for each source
  tryCatch({
    icle_tmb <- calculate_tmb(icle_maf_file, "ICLE")
    cpt6_tmb <- calculate_tmb(cpt6_maf_file, "CPT6")
    tcga_tmb <- calculate_tmb(tcga_maf_file, "TCGA_ILC")
    
    # Combine TMB data
    combined_tmb <- rbind(icle_tmb, cpt6_tmb, tcga_tmb)
    
    # Save TMB data
    write.csv(combined_tmb, file.path(output_dir, "tumor_mutation_burden_comparison.csv"), row.names = FALSE)
    
    # Calculate mean TMB values
    icle_mean_tmb <- mean(icle_tmb$total_perMB)
    tcga_mean_tmb <- mean(tcga_tmb$total_perMB)
    cpt6_value <- cpt6_tmb$total_perMB[1]
    
    # Perform t-test between ICLE and TCGA_ILC
    t_test_result <- t.test(icle_tmb$total_perMB, tcga_tmb$total_perMB)
    p_value <- t_test_result$p.value
    p_value_text <- ifelse(p_value < 0.001, "p < 0.001", 
                          ifelse(p_value < 0.01, "p < 0.01", 
                                ifelse(p_value < 0.05, "p < 0.05", 
                                      sprintf("p = %.3f", p_value))))
    
    # Set a fixed random seed for reproducibility
    set.seed(123)
    
    # Create a boxplot of TMB by source
    p <- ggplot(combined_tmb, aes(x = Source, y = total_perMB, fill = Source)) +
      # Add boxplots
      geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
      # Add jittered points
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
      # Add mean values as text
      annotate("text", x = 1, y = cpt6_value + 0.5,
               label = sprintf("Mean: %.2f", cpt6_value),
               color = "black", size = 4) +
      annotate("text", x = 2, y = icle_mean_tmb + 0.5,
               label = sprintf("Mean: %.2f", icle_mean_tmb),
               color = "black", size = 4) +
      annotate("text", x = 3, y = tcga_mean_tmb + 0.5,
               label = sprintf("Mean: %.2f", tcga_mean_tmb),
               color = "black", size = 4) +
      # Add significance bracket between ICLE and TCGA_ILC
      annotate("segment", x = 2, xend = 3, 
               y = max(combined_tmb$total_perMB) * 0.95, 
               yend = max(combined_tmb$total_perMB) * 0.95,
               color = "black", linewidth = 0.8) +
      # Add significance text
      annotate("text", x = 2.5, y = max(combined_tmb$total_perMB) * 0.98,
               label = p_value_text, size = 5) +
      # Customize appearance
      theme_classic() +
      scale_fill_manual(values = c("CPT6" = "black", "ICLE" = "skyblue", "TCGA_ILC" = "lightgreen"),
                        labels = c("CPT6" = "CPT6 (n=1)", 
                                  "ICLE" = "ICLE (n=17)", 
                                  "TCGA_ILC" = "TCGA_ILC (n=179)")) +
      labs(
        title = "Tumor Mutation Burden Comparison",
        y = "Mutations per Mb",
        x = ""
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12)
      )
    
    # Save the plot
    ggsave(file.path(output_dir, "tmb_boxplot_comparison.png"), p, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "tmb_boxplot_comparison.pdf"), p, width = 10, height = 8)
    
    cat("Tumor mutation burden analysis saved to", output_dir, "\n")
    
    # Print summary statistics
    cat("\nTMB Summary Statistics (mut/Mb):\n")
    for (source in unique(combined_tmb$Source)) {
      source_tmb <- combined_tmb$total_perMB[combined_tmb$Source == source]
      cat(sprintf("\n%s:\n", source))
      cat(sprintf("  Samples: %d\n", length(source_tmb)))
      cat(sprintf("  Mean: %.4f\n", mean(source_tmb)))
      cat(sprintf("  Median: %.4f\n", median(source_tmb)))
      cat(sprintf("  Min: %.4f\n", min(source_tmb)))
      cat(sprintf("  Max: %.4f\n", max(source_tmb)))
    }
    
    # Print t-test results
    cat("\nT-test Results (ICLE vs TCGA_ILC):\n")
    cat(sprintf("  t-statistic: %.4f\n", t_test_result$statistic))
    cat(sprintf("  p-value: %.6f\n", t_test_result$p.value))
    cat(sprintf("  95%% Confidence Interval: [%.4f, %.4f]\n", 
               t_test_result$conf.int[1], t_test_result$conf.int[2]))
    cat(sprintf("  Mean difference: %.4f\n", diff(t_test_result$estimate)))
    
  }, error = function(e) {
    cat("Error generating TMB comparison:", conditionMessage(e), "\n")
  })
}

# Run the main function
generate_tmb_comparison() 