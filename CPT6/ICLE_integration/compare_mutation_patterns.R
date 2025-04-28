#!/usr/bin/env Rscript

# Script to compare mutation patterns between CPT6 and ICLE
# Author: Daisong

# Load required libraries
library(maftools)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)

# Set file paths
icle_maf_file <- "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
cpt6_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
output_dir <- "results/mutation_patterns"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to analyze mutation patterns for a dataset
analyze_mutation_patterns <- function(maf_file, dataset_name) {
  cat(sprintf("\n\n========== Analyzing %s Mutation Patterns ==========\n\n", dataset_name))
  
  # Read MAF file
  cat(sprintf("Reading %s MAF file...\n", dataset_name))
  maf <- read.maf(maf_file)
  
  # 1. Basic mutation summary
  cat("Generating mutation summary...\n")
  pdf(file.path(output_dir, paste0(dataset_name, "_mutation_summary.pdf")), width = 10, height = 8)
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
  dev.off()
  
  # 2. Transition and Transversion analysis
  cat("Analyzing transition and transversion patterns...\n")
  pdf(file.path(output_dir, paste0(dataset_name, "_titv.pdf")), width = 10, height = 6)
  titv_result <- titv(maf = maf, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = titv_result)
  dev.off()
  
  # Save Ti/Tv summary
  write.csv(titv_result$fraction, 
            file.path(output_dir, paste0(dataset_name, "_titv_fraction.csv")), 
            row.names = TRUE)
  
  # 3. Oncoplot for top mutated genes
  cat("Generating oncoplot...\n")
  pdf(file.path(output_dir, paste0(dataset_name, "_oncoplot.pdf")), width = 12, height = 10)
  oncoplot(maf = maf, top = 20)
  dev.off()
  
  # 4. Lollipop plots for top genes
  cat("Generating lollipop plots for top genes...\n")
  top_genes <- getGeneSummary(maf)$Hugo_Symbol[1:5]
  
  for (gene in top_genes) {
    tryCatch({
      pdf(file.path(output_dir, paste0(dataset_name, "_", gene, "_lollipop.pdf")), width = 10, height = 6)
      lollipopPlot(maf = maf, gene = gene, AACol = "Protein_Change", labelPos = "all")
      dev.off()
    }, error = function(e) {
      cat("Error generating lollipop plot for gene", gene, ":", conditionMessage(e), "\n")
    })
  }
  
  # 5. Mutation spectrum
  cat("Analyzing mutation spectrum...\n")
  
  # Get variant classifications
  var_class <- getClinicalData(maf)$Variant_Classification
  var_type <- getClinicalData(maf)$Variant_Type
  
  # Count variant classifications and types
  class_counts <- table(var_class)
  type_counts <- table(var_type)
  
  # Save summaries
  write.csv(as.data.frame(class_counts), 
            file.path(output_dir, paste0(dataset_name, "_variant_classifications.csv")), 
            row.names = FALSE)
  write.csv(as.data.frame(type_counts), 
            file.path(output_dir, paste0(dataset_name, "_variant_types.csv")), 
            row.names = FALSE)
  
  # 6. Somatic interactions (co-occurrence and mutual exclusivity)
  cat("Analyzing somatic interactions...\n")
  tryCatch({
    # Check if there are enough genes for somatic interactions
    gene_count <- nrow(getGeneSummary(maf))
    if (gene_count >= 2) {
      pdf(file.path(output_dir, paste0(dataset_name, "_somatic_interactions.pdf")), width = 12, height = 10)
      somaticInteractions(maf = maf, top = min(25, gene_count), pvalue = 0.1)
      dev.off()
    } else {
      cat("Not enough genes for somatic interactions analysis (minimum 2 required).\n")
    }
  }, error = function(e) {
    cat("Error in somatic interactions analysis:", conditionMessage(e), "\n")
  })
  
  # 7. Rainfall plot (if applicable)
  cat("Generating rainfall plot...\n")
  tryCatch({
    if (dataset_name == "CPT6" || length(getSampleSummary(maf)$Tumor_Sample_Barcode) <= 5) {
      pdf(file.path(output_dir, paste0(dataset_name, "_rainfall.pdf")), width = 12, height = 6)
      rainfallPlot(maf = maf, detectChangePoints = TRUE, pointSize = 0.4)
      dev.off()
    } else {
      # For ICLE with many samples, just take a few
      cat("Generating rainfall plot for a subset of samples...\n")
      selected_samples <- getSampleSummary(maf)$Tumor_Sample_Barcode[1:5]
      maf_subset <- subsetMaf(maf, tsb = selected_samples)
      
      pdf(file.path(output_dir, paste0(dataset_name, "_rainfall_subset.pdf")), width = 12, height = 6)
      rainfallPlot(maf = maf_subset, detectChangePoints = TRUE, pointSize = 0.4)
      dev.off()
    }
  }, error = function(e) {
    cat("Error in rainfall plot generation:", conditionMessage(e), "\n")
  })
  
  # Return summary data
  return(list(
    maf = maf,
    titv = titv_result,
    class_counts = class_counts,
    type_counts = type_counts
  ))
}

# Function to compare mutation patterns between datasets
compare_mutation_patterns <- function(cpt6_results, icle_results) {
  cat("\n\n========== Comparing CPT6 and ICLE Mutation Patterns ==========\n\n")
  
  # 1. Compare variant classifications
  cat("Comparing variant classifications...\n")
  
  # Get all unique classifications
  all_classes <- unique(c(names(cpt6_results$class_counts), names(icle_results$class_counts)))
  
  # Create a comparison data frame
  class_comparison <- data.frame(
    Class = all_classes,
    CPT6_Count = sapply(all_classes, function(c) {
      if (c %in% names(cpt6_results$class_counts)) cpt6_results$class_counts[c] else 0
    }),
    ICLE_Count = sapply(all_classes, function(c) {
      if (c %in% names(icle_results$class_counts)) icle_results$class_counts[c] else 0
    })
  )
  
  # Calculate percentages
  class_comparison$CPT6_Percent <- class_comparison$CPT6_Count / sum(class_comparison$CPT6_Count) * 100
  class_comparison$ICLE_Percent <- class_comparison$ICLE_Count / sum(class_comparison$ICLE_Count) * 100
  
  # Save comparison
  write.csv(class_comparison, 
            file.path(output_dir, "variant_classification_comparison.csv"), 
            row.names = FALSE)
  
  # Create a plot
  class_long <- melt(
    class_comparison[, c("Class", "CPT6_Percent", "ICLE_Percent")],
    id.vars = "Class",
    variable.name = "Dataset",
    value.name = "Percentage"
  )
  
  # Clean up dataset names
  class_long$Dataset <- gsub("_Percent", "", class_long$Dataset)
  
  # Create a bar plot
  p1 <- ggplot(class_long, aes(x = Class, y = Percentage, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Comparison of Variant Classifications",
      x = "Variant Classification",
      y = "Percentage (%)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("CPT6" = "coral", "ICLE" = "skyblue"))
  
  # Save the plot
  ggsave(file.path(output_dir, "variant_classification_comparison.pdf"), p1, width = 12, height = 8)
  ggsave(file.path(output_dir, "variant_classification_comparison.png"), p1, width = 12, height = 8, dpi = 300)
  
  # 2. Compare variant types
  cat("Comparing variant types...\n")
  
  # Get all unique types
  all_types <- unique(c(names(cpt6_results$type_counts), names(icle_results$type_counts)))
  
  # Create a comparison data frame
  type_comparison <- data.frame(
    Type = all_types,
    CPT6_Count = sapply(all_types, function(t) {
      if (t %in% names(cpt6_results$type_counts)) cpt6_results$type_counts[t] else 0
    }),
    ICLE_Count = sapply(all_types, function(t) {
      if (t %in% names(icle_results$type_counts)) icle_results$type_counts[t] else 0
    })
  )
  
  # Calculate percentages
  type_comparison$CPT6_Percent <- type_comparison$CPT6_Count / sum(type_comparison$CPT6_Count) * 100
  type_comparison$ICLE_Percent <- type_comparison$ICLE_Count / sum(type_comparison$ICLE_Count) * 100
  
  # Save comparison
  write.csv(type_comparison, 
            file.path(output_dir, "variant_type_comparison.csv"), 
            row.names = FALSE)
  
  # Create a plot
  type_long <- melt(
    type_comparison[, c("Type", "CPT6_Percent", "ICLE_Percent")],
    id.vars = "Type",
    variable.name = "Dataset",
    value.name = "Percentage"
  )
  
  # Clean up dataset names
  type_long$Dataset <- gsub("_Percent", "", type_long$Dataset)
  
  # Create a bar plot
  p2 <- ggplot(type_long, aes(x = Type, y = Percentage, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Comparison of Variant Types",
      x = "Variant Type",
      y = "Percentage (%)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("CPT6" = "coral", "ICLE" = "skyblue"))
  
  # Save the plot
  ggsave(file.path(output_dir, "variant_type_comparison.pdf"), p2, width = 12, height = 8)
  ggsave(file.path(output_dir, "variant_type_comparison.png"), p2, width = 12, height = 8, dpi = 300)
  
  # 3. Compare Ti/Tv ratios
  cat("Comparing Ti/Tv ratios...\n")
  
  # Extract Ti/Tv data
  cpt6_titv <- cpt6_results$titv$fraction
  icle_titv <- icle_results$titv$fraction
  
  # Create a comparison data frame
  titv_comparison <- data.frame(
    Category = rownames(cpt6_titv),
    CPT6_Fraction = cpt6_titv$fraction,
    ICLE_Fraction = icle_titv$fraction
  )
  
  # Save comparison
  write.csv(titv_comparison, 
            file.path(output_dir, "titv_comparison.csv"), 
            row.names = FALSE)
  
  # Create a plot
  titv_long <- melt(
    titv_comparison,
    id.vars = "Category",
    variable.name = "Dataset",
    value.name = "Fraction"
  )
  
  # Clean up dataset names
  titv_long$Dataset <- gsub("_Fraction", "", titv_long$Dataset)
  
  # Create a bar plot
  p3 <- ggplot(titv_long, aes(x = Category, y = Fraction, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Comparison of Transition/Transversion Fractions",
      x = "Category",
      y = "Fraction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("CPT6" = "coral", "ICLE" = "skyblue"))
  
  # Save the plot
  ggsave(file.path(output_dir, "titv_comparison.pdf"), p3, width = 12, height = 8)
  ggsave(file.path(output_dir, "titv_comparison.png"), p3, width = 12, height = 8, dpi = 300)
  
  # 4. Compare commonly mutated genes
  cat("Comparing commonly mutated genes...\n")
  
  # Get top mutated genes
  cpt6_genes <- getGeneSummary(cpt6_results$maf)$Hugo_Symbol
  icle_genes <- getGeneSummary(icle_results$maf)$Hugo_Symbol
  
  # Limit to top 20 or less if fewer genes are available
  cpt6_top <- head(cpt6_genes, min(20, length(cpt6_genes)))
  icle_top <- head(icle_genes, min(20, length(icle_genes)))
  
  # Find common genes
  common_genes <- intersect(cpt6_top, icle_top)
  
  # Create a summary
  gene_summary <- data.frame(
    CPT6_Top_Genes = c(cpt6_top, rep(NA, max(0, length(icle_top) - length(cpt6_top)))),
    ICLE_Top_Genes = c(icle_top, rep(NA, max(0, length(cpt6_top) - length(icle_top))))
  )
  
  # Save summary
  write.csv(gene_summary, 
            file.path(output_dir, "top_genes_comparison.csv"), 
            row.names = FALSE)
  
  # Create a Venn diagram-like summary
  cat("\nTop Mutated Genes Comparison:\n")
  cat("CPT6 only:", setdiff(cpt6_top, icle_top), "\n")
  cat("ICLE only:", setdiff(icle_top, cpt6_top), "\n")
  cat("Common:", common_genes, "\n")
  
  # Save this summary as text
  sink(file.path(output_dir, "top_genes_venn_summary.txt"))
  cat("Top Mutated Genes Comparison:\n")
  cat("CPT6 only:", setdiff(cpt6_top, icle_top), "\n")
  cat("ICLE only:", setdiff(icle_top, cpt6_top), "\n")
  cat("Common:", common_genes, "\n")
  sink()
  
  # Return comparison results
  return(list(
    class_comparison = class_comparison,
    type_comparison = type_comparison,
    titv_comparison = titv_comparison,
    gene_summary = gene_summary
  ))
}

# Main function
main <- function() {
  tryCatch({
    # Analyze CPT6 mutation patterns
    cpt6_results <- analyze_mutation_patterns(cpt6_maf_file, "CPT6")
    
    # Analyze ICLE mutation patterns
    icle_results <- analyze_mutation_patterns(icle_maf_file, "ICLE")
    
    # Compare mutation patterns
    comparison_results <- compare_mutation_patterns(cpt6_results, icle_results)
    
    cat("\nMutation Pattern Analysis Complete!\n")
    cat("Results saved to:", output_dir, "\n")
    
  }, error = function(e) {
    cat("Error in mutation pattern analysis:", conditionMessage(e), "\n")
    print(e)
  })
}

# Run the main function
main() 