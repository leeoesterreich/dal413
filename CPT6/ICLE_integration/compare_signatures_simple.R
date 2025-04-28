#!/usr/bin/env Rscript

# Simple script to compare mutational signatures between CPT6 and ICLE
# Author: Daisong

# Load required libraries
library(maftools)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)  # For melt function

# Set file paths
icle_maf_file <- "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
cpt6_maf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
output_dir <- "results/signatures"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Main function to run the signature comparison
compare_signatures_main <- function() {
  tryCatch({
    # Read MAF files
    cat("Reading MAF files...\n")
    icle_maf <- read.maf(icle_maf_file)
    cpt6_maf <- read.maf(cpt6_maf_file)
    
    # 1. Basic mutation summary for both datasets
    cat("Generating mutation summaries...\n")
    
    # CPT6 summary
    pdf(file.path(output_dir, "CPT6_mutation_summary.pdf"), width = 10, height = 8)
    plotmafSummary(maf = cpt6_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
    dev.off()
    
    # ICLE summary
    pdf(file.path(output_dir, "ICLE_mutation_summary.pdf"), width = 10, height = 8)
    plotmafSummary(maf = icle_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
    dev.off()
    
    # 2. Transition and Transversion analysis
    cat("Analyzing transition and transversion patterns...\n")
    
    # CPT6 Ti/Tv
    pdf(file.path(output_dir, "CPT6_titv.pdf"), width = 10, height = 6)
    titv = titv(maf = cpt6_maf, plot = FALSE, useSyn = TRUE)
    plotTiTv(res = titv)
    dev.off()
    
    # ICLE Ti/Tv
    pdf(file.path(output_dir, "ICLE_titv.pdf"), width = 10, height = 6)
    titv = titv(maf = icle_maf, plot = FALSE, useSyn = TRUE)
    plotTiTv(res = titv)
    dev.off()
    
    # 3. Oncoplot for top mutated genes
    cat("Generating oncoplots...\n")
    
    # CPT6 oncoplot
    pdf(file.path(output_dir, "CPT6_oncoplot.pdf"), width = 10, height = 8)
    oncoplot(maf = cpt6_maf, top = 20)
    dev.off()
    
    # ICLE oncoplot
    pdf(file.path(output_dir, "ICLE_oncoplot.pdf"), width = 10, height = 8)
    oncoplot(maf = icle_maf, top = 20)
    dev.off()
    
    # 4. Lollipop plots for key genes
    cat("Generating lollipop plots for key genes...\n")
    
    # Get top mutated genes in both datasets
    cpt6_genes <- getGeneSummary(cpt6_maf)$Hugo_Symbol[1:5]
    icle_genes <- getGeneSummary(icle_maf)$Hugo_Symbol[1:5]
    
    # Combine and get unique genes
    key_genes <- unique(c(cpt6_genes, icle_genes))
    
    # Generate lollipop plots for each key gene
    for (gene in key_genes) {
      tryCatch({
        # CPT6 lollipop
        if (gene %in% cpt6_genes) {
          pdf(file.path(output_dir, paste0("CPT6_", gene, "_lollipop.pdf")), width = 10, height = 6)
          lollipopPlot(maf = cpt6_maf, gene = gene)
          dev.off()
        }
        
        # ICLE lollipop
        if (gene %in% icle_genes) {
          pdf(file.path(output_dir, paste0("ICLE_", gene, "_lollipop.pdf")), width = 10, height = 6)
          lollipopPlot(maf = icle_maf, gene = gene)
          dev.off()
        }
      }, error = function(e) {
        cat("Error generating lollipop plot for gene", gene, ":", conditionMessage(e), "\n")
      })
    }
    
    # 5. Compare mutation types between datasets
    cat("Comparing mutation types between datasets...\n")
    
    # Get variant classifications
    cpt6_var_class <- getClinicalData(cpt6_maf)$Variant_Classification
    icle_var_class <- getClinicalData(icle_maf)$Variant_Classification
    
    # Count variant classifications
    cpt6_class_counts <- table(cpt6_var_class)
    icle_class_counts <- table(icle_var_class)
    
    # Combine into a data frame
    all_classes <- unique(c(names(cpt6_class_counts), names(icle_class_counts)))
    
    class_comparison <- data.frame(
      Class = all_classes,
      CPT6 = sapply(all_classes, function(t) {
        if(t %in% names(cpt6_class_counts)) cpt6_class_counts[t] else 0
      }),
      ICLE = sapply(all_classes, function(t) {
        if(t %in% names(icle_class_counts)) icle_class_counts[t] else 0
      })
    )
    
    # Calculate percentages
    class_comparison$CPT6_Percent <- class_comparison$CPT6 / sum(class_comparison$CPT6) * 100
    class_comparison$ICLE_Percent <- class_comparison$ICLE / sum(class_comparison$ICLE) * 100
    
    # Create a long format for plotting
    class_long <- melt(
      class_comparison[, c("Class", "CPT6_Percent", "ICLE_Percent")],
      id.vars = "Class",
      variable.name = "Dataset",
      value.name = "Percentage"
    )
    
    # Clean up dataset names
    class_long$Dataset <- gsub("_Percent", "", class_long$Dataset)
    
    # Create a bar plot comparing variant classifications
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
    
    # Save the comparison plot
    ggsave(file.path(output_dir, "variant_classification_comparison.pdf"), p1, width = 12, height = 8)
    ggsave(file.path(output_dir, "variant_classification_comparison.png"), p1, width = 12, height = 8, dpi = 300)
    
    # 6. Compare mutation types (SNP, INS, DEL)
    cat("Comparing mutation types (SNP, INS, DEL)...\n")
    
    # Get variant types
    cpt6_var_type <- getClinicalData(cpt6_maf)$Variant_Type
    icle_var_type <- getClinicalData(icle_maf)$Variant_Type
    
    # Count variant types
    cpt6_type_counts <- table(cpt6_var_type)
    icle_type_counts <- table(icle_var_type)
    
    # Combine into a data frame
    all_types <- unique(c(names(cpt6_type_counts), names(icle_type_counts)))
    
    type_comparison <- data.frame(
      Type = all_types,
      CPT6 = sapply(all_types, function(t) {
        if(t %in% names(cpt6_type_counts)) cpt6_type_counts[t] else 0
      }),
      ICLE = sapply(all_types, function(t) {
        if(t %in% names(icle_type_counts)) icle_type_counts[t] else 0
      })
    )
    
    # Calculate percentages
    type_comparison$CPT6_Percent <- type_comparison$CPT6 / sum(type_comparison$CPT6) * 100
    type_comparison$ICLE_Percent <- type_comparison$ICLE / sum(type_comparison$ICLE) * 100
    
    # Create a long format for plotting
    type_long <- melt(
      type_comparison[, c("Type", "CPT6_Percent", "ICLE_Percent")],
      id.vars = "Type",
      variable.name = "Dataset",
      value.name = "Percentage"
    )
    
    # Clean up dataset names
    type_long$Dataset <- gsub("_Percent", "", type_long$Dataset)
    
    # Create a bar plot comparing variant types
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
    
    # Save the comparison plot
    ggsave(file.path(output_dir, "variant_type_comparison.pdf"), p2, width = 12, height = 8)
    ggsave(file.path(output_dir, "variant_type_comparison.png"), p2, width = 12, height = 8, dpi = 300)
    
    # 7. Compare SNV classes
    cat("Comparing SNV classes...\n")
    
    # Get SNV classes
    cpt6_snv <- getSampleSummary(cpt6_maf)$C>A + getSampleSummary(cpt6_maf)$C>G + 
                getSampleSummary(cpt6_maf)$C>T + getSampleSummary(cpt6_maf)$T>A + 
                getSampleSummary(cpt6_maf)$T>C + getSampleSummary(cpt6_maf)$T>G
    
    icle_snv <- getSampleSummary(icle_maf)$C>A + getSampleSummary(icle_maf)$C>G + 
                getSampleSummary(icle_maf)$C>T + getSampleSummary(icle_maf)$T>A + 
                getSampleSummary(icle_maf)$T>C + getSampleSummary(icle_maf)$T>G
    
    # Create a data frame for SNV classes
    snv_classes <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    
    snv_comparison <- data.frame(
      Class = snv_classes,
      CPT6 = c(
        sum(getSampleSummary(cpt6_maf)$`C>A`),
        sum(getSampleSummary(cpt6_maf)$`C>G`),
        sum(getSampleSummary(cpt6_maf)$`C>T`),
        sum(getSampleSummary(cpt6_maf)$`T>A`),
        sum(getSampleSummary(cpt6_maf)$`T>C`),
        sum(getSampleSummary(cpt6_maf)$`T>G`)
      ),
      ICLE = c(
        sum(getSampleSummary(icle_maf)$`C>A`),
        sum(getSampleSummary(icle_maf)$`C>G`),
        sum(getSampleSummary(icle_maf)$`C>T`),
        sum(getSampleSummary(icle_maf)$`T>A`),
        sum(getSampleSummary(icle_maf)$`T>C`),
        sum(getSampleSummary(icle_maf)$`T>G`)
      )
    )
    
    # Calculate percentages
    snv_comparison$CPT6_Percent <- snv_comparison$CPT6 / sum(snv_comparison$CPT6) * 100
    snv_comparison$ICLE_Percent <- snv_comparison$ICLE / sum(snv_comparison$ICLE) * 100
    
    # Create a long format for plotting
    snv_long <- melt(
      snv_comparison[, c("Class", "CPT6_Percent", "ICLE_Percent")],
      id.vars = "Class",
      variable.name = "Dataset",
      value.name = "Percentage"
    )
    
    # Clean up dataset names
    snv_long$Dataset <- gsub("_Percent", "", snv_long$Dataset)
    
    # Create a bar plot comparing SNV classes
    p3 <- ggplot(snv_long, aes(x = Class, y = Percentage, fill = Dataset)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(
        title = "Comparison of SNV Classes",
        x = "SNV Class",
        y = "Percentage (%)"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_fill_manual(values = c("CPT6" = "coral", "ICLE" = "skyblue"))
    
    # Save the comparison plot
    ggsave(file.path(output_dir, "snv_class_comparison.pdf"), p3, width = 12, height = 8)
    ggsave(file.path(output_dir, "snv_class_comparison.png"), p3, width = 12, height = 8, dpi = 300)
    
    # 8. Generate a summary report
    cat("\nMutation Analysis Complete!\n")
    cat("Results saved to:", output_dir, "\n")
    
    # Print summary of variant classifications
    cat("\nVariant Classification Comparison:\n")
    print(class_comparison[, c("Class", "CPT6_Percent", "ICLE_Percent")])
    
    # Print summary of variant types
    cat("\nVariant Type Comparison:\n")
    print(type_comparison[, c("Type", "CPT6_Percent", "ICLE_Percent")])
    
    # Print summary of SNV classes
    cat("\nSNV Class Comparison:\n")
    print(snv_comparison[, c("Class", "CPT6_Percent", "ICLE_Percent")])
    
  }, error = function(e) {
    cat("Error in mutation comparison:", conditionMessage(e), "\n")
    print(e)
  })
}

# Run the main function
compare_signatures_main() 