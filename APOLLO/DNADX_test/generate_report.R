#!/usr/bin/env Rscript

# Script to generate a comprehensive report of CNA signature analysis
# This script creates visualizations and summaries of the analysis results

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(knitr)
library(rmarkdown)

# Set paths
output_dir <- "./output"
dnadx_dir <- "./DNA-based-predictors-of-non-genetic-cancer-phenotypes"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check if required files exist
segment_scores_file <- file.path(output_dir, "segment_scores.rda")
gene_level_data_file <- file.path(output_dir, "gene_level_cn_data.rda")

if (!file.exists(segment_scores_file)) {
  cat("Error: Segment scores file not found:", segment_scores_file, "\n")
  cat("Please run the analyze_signatures.R script first.\n")
  quit(status = 1)
}

# Load the segment scores and priority list
cat("Loading analysis results...\n")
load(segment_scores_file)

# Function to create a bar plot of top signatures
create_signature_barplot <- function(signature_priority, top_n = 40) {
  # Check if there are any signatures
  if (nrow(signature_priority) == 0) {
    cat("Warning: No signatures available for barplot\n")
    return(NULL)
  }
  
  # Select top N signatures (or all if fewer than top_n)
  n_signatures <- min(top_n, nrow(signature_priority))
  top_signatures <- head(signature_priority, n_signatures)
  
  # Create a bar plot
  p <- ggplot(top_signatures, aes(x = reorder(segment, variance), y = variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste("Top", n_signatures, "CNA Signatures by Variance"),
      x = "Signature",
      y = "Variance"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    ) +
    # Add text labels with mutation counts, positioned to the right of the bars
    geom_text(aes(label = sprintf("%.3f", variance)), hjust = -0.3, size = 3)
  
  # Save the plot
  tryCatch({
    ggsave(file.path(output_dir, "top_signatures_barplot.pdf"), p, width = 10, height = 12)
    cat("Barplot saved to", file.path(output_dir, "top_signatures_barplot.pdf"), "\n")
  }, error = function(e) {
    cat("Error saving barplot:", e$message, "\n")
  })
  
  return(p)
}

# Function to create a correlation heatmap of top signatures
create_correlation_heatmap <- function(segment_scores, top_n = 40) {
  # Check if there are enough signatures
  if (nrow(segment_scores) < 2) {
    cat("Warning: Not enough signatures for correlation heatmap (need at least 2)\n")
    return(NULL)
  }
  
  # Remove NA values by replacing them with the column means
  for (i in 1:nrow(segment_scores)) {
    na_indices <- is.na(segment_scores[i, ])
    if (any(na_indices)) {
      non_na_values <- segment_scores[i, !na_indices]
      if (length(non_na_values) > 0) {
        segment_scores[i, na_indices] <- mean(non_na_values)
      } else {
        segment_scores[i, na_indices] <- 0
      }
    }
  }
  
  # Select top N signatures by variance (or all if fewer than top_n)
  segment_variance <- apply(segment_scores, 1, var)
  n_segments <- min(top_n, nrow(segment_scores))
  
  # Make sure we have valid row names
  if (is.null(rownames(segment_scores))) {
    rownames(segment_scores) <- paste0("Segment_", 1:nrow(segment_scores))
    top_segments <- paste0("Segment_", 1:n_segments)
  } else {
    top_segments <- names(sort(segment_variance, decreasing = TRUE)[1:n_segments])
  }
  
  # Ensure top_segments exist in the row names
  valid_segments <- intersect(top_segments, rownames(segment_scores))
  if (length(valid_segments) < 2) {
    cat("Warning: Not enough valid segments for correlation heatmap (need at least 2)\n")
    return(NULL)
  }
  
  cat("Creating correlation heatmap with", length(valid_segments), "segments\n")
  
  # Calculate correlation matrix
  tryCatch({
    cor_matrix <- cor(t(segment_scores[valid_segments, , drop = FALSE]), method = "pearson")
    
    # Create heatmap
    pheatmap(
      cor_matrix,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
      fontsize_row = 8,
      fontsize_col = 8,
      filename = file.path(output_dir, "signature_correlation_heatmap.pdf"),
      width = 12,
      height = 12
    )
    cat("Correlation heatmap saved to", file.path(output_dir, "signature_correlation_heatmap.pdf"), "\n")
  }, error = function(e) {
    cat("Error creating correlation heatmap:", e$message, "\n")
  })
}

# Function to create a PCA plot of samples based on segment scores
create_pca_plot <- function(segment_scores) {
  # Check if there are enough samples and signatures
  if (ncol(segment_scores) < 3 || nrow(segment_scores) < 3) {
    cat("Warning: Not enough samples or signatures for PCA (need at least 3 of each)\n")
    return(NULL)
  }
  
  # Perform PCA
  tryCatch({
    pca_result <- prcomp(t(segment_scores), scale. = TRUE)
    
    # Extract PC scores
    pca_data <- as.data.frame(pca_result$x)
    
    # Calculate variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    
    # Create PCA plot
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_text(aes(label = rownames(pca_data)), size = 2, vjust = -1) +
      labs(
        title = "PCA of Samples Based on CNA Signatures",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_minimal()
    
    # Save the plot
    ggsave(file.path(output_dir, "pca_plot.pdf"), p, width = 10, height = 8)
    cat("PCA plot saved to", file.path(output_dir, "pca_plot.pdf"), "\n")
    
    return(p)
  }, error = function(e) {
    cat("Error creating PCA plot:", e$message, "\n")
    return(NULL)
  })
}

# Function to create a sample clustering dendrogram
create_sample_dendrogram <- function(segment_scores) {
  # Check if there are enough samples
  if (ncol(segment_scores) < 2) {
    cat("Warning: Not enough samples for clustering (need at least 2)\n")
    return(NULL)
  }
  
  # Calculate distance matrix and perform hierarchical clustering
  tryCatch({
    # Calculate distance matrix
    dist_matrix <- dist(t(segment_scores), method = "euclidean")
    
    # Perform hierarchical clustering
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    # Create dendrogram plot
    pdf(file.path(output_dir, "sample_dendrogram.pdf"), width = 10, height = 8)
    plot(hc, main = "Sample Clustering Based on CNA Signatures", xlab = "", sub = "")
    dev.off()
    
    cat("Sample dendrogram saved to", file.path(output_dir, "sample_dendrogram.pdf"), "\n")
  }, error = function(e) {
    cat("Error creating sample dendrogram:", e$message, "\n")
  })
}

# Function to create an RMarkdown report
create_rmarkdown_report <- function(segment_scores, signature_priority) {
  # Create RMarkdown template
  rmd_content <- '---
title: "CNA Signature Analysis Report"
date: "`r format(Sys.time(), "%Y-%m-%d")`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(knitr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
```

## Introduction

This report presents the results of Copy Number Alteration (CNA) signature analysis. CNA signatures represent recurrent patterns of copy number changes across the genome that may be associated with specific cancer types, subtypes, or phenotypes.

## Sample Overview

The analysis included **`r ncol(segment_scores)`** samples. The following visualizations show how these samples relate to each other based on their CNA signature profiles.

### Sample Clustering

The dendrogram below shows hierarchical clustering of samples based on their CNA signature profiles. Samples that cluster together have similar CNA patterns.

```{r sample_clustering, fig.width=10, fig.height=8}
# Calculate distance matrix
dist_matrix <- dist(t(segment_scores), method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram
plot(hc, main = "Sample Clustering Based on CNA Signatures", xlab = "", sub = "")
```

### Principal Component Analysis

The PCA plot below shows the first two principal components of the CNA signature data. Samples that are close to each other in this plot have similar CNA signature profiles.

```{r pca_plot, fig.width=10, fig.height=8}
# Perform PCA
pca_result <- prcomp(t(segment_scores), scale. = TRUE)

# Extract PC scores
pca_data <- as.data.frame(pca_result$x)

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = rownames(pca_data)), size = 2, vjust = -1) +
  labs(
    title = "PCA of Samples Based on CNA Signatures",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal()
```

## CNA Signature Analysis

### Top CNA Signatures by Variance

The table below shows the top CNA signatures ranked by variance across samples. Signatures with higher variance are more informative for distinguishing between samples.

```{r top_signatures_table}
top_n <- min(20, nrow(signature_priority))
# Select only the simplified name, variance, and description columns
display_priority <- signature_priority[, c("segment", "variance", "description")]
colnames(display_priority) <- c("Signature", "Variance", "Description")
kable(head(display_priority, top_n),
      caption = paste("Top", top_n, "CNA Signatures by Variance"))
```

### Signature Barplot

The barplot below shows the top CNA signatures ranked by variance across samples.

```{r signature_barplot, fig.width=10, fig.height=10}
# Select top N signatures
top_n <- min(20, nrow(signature_priority))
top_signatures <- head(signature_priority, top_n)

# Create a bar plot
ggplot(top_signatures, aes(x = reorder(segment, variance), y = variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = paste("Top", top_n, "CNA Signatures by Variance"),
    x = "Signature",
    y = "Variance"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  ) +
  # Add text labels with mutation counts, positioned to the right of the bars
  geom_text(aes(label = sprintf("%.3f", variance)), hjust = -0.3, size = 3)
```

### Signature Correlation Heatmap

The heatmap below shows the correlation between the top CNA signatures. Signatures that cluster together may represent related biological processes.

```{r correlation_heatmap, fig.width=12, fig.height=12}
# Calculate variance of each segment score across samples
segment_variance <- apply(segment_scores, 1, var)

# Select top N signatures by variance
top_n <- min(30, nrow(segment_scores))
top_segments <- names(sort(segment_variance, decreasing = TRUE)[1:top_n])

# Get simplified names for the heatmap
simplified_names <- signature_priority$segment[match(top_segments, signature_priority$original_segment)]
if (length(simplified_names) == length(top_segments)) {
  # Create a copy of segment_scores with simplified row names
  segment_scores_simplified <- segment_scores[top_segments, , drop = FALSE]
  rownames(segment_scores_simplified) <- simplified_names
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(segment_scores_simplified), method = "pearson")
  
  # Create heatmap
  pheatmap(
    cor_matrix,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
    main = "Correlation Between Top CNA Signatures"
  )
} else {
  # Fallback to original names if there is a mismatch
  cor_matrix <- cor(t(segment_scores[top_segments, , drop = FALSE]), method = "pearson")
  
  # Create heatmap
  pheatmap(
    cor_matrix,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
    main = "Correlation Between Top CNA Signatures"
  )
}
```

### Segment Scores Heatmap

The heatmap below shows the segment scores for the top CNA signatures across all samples. Red indicates higher scores (amplification), while blue indicates lower scores (deletion).

```{r segment_scores_heatmap, fig.width=12, fig.height=12}
# Select top N signatures by variance
top_n <- min(20, nrow(segment_scores))
top_segments <- names(sort(segment_variance, decreasing = TRUE)[1:top_n])

# Get simplified names for the heatmap
simplified_names <- signature_priority$segment[match(top_segments, signature_priority$original_segment)]
if (length(simplified_names) == length(top_segments)) {
  # Create a copy of segment_scores with simplified row names
  segment_scores_simplified <- segment_scores[top_segments, , drop = FALSE]
  rownames(segment_scores_simplified) <- simplified_names
  
  # Create heatmap
  pheatmap(
    segment_scores_simplified,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = "Segment Scores Heatmap"
  )
} else {
  # Fallback to original names if there is a mismatch
  pheatmap(
    segment_scores[top_segments, ],
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = "Segment Scores Heatmap"
  )
}
```

## Interpretation of CNA Signatures

CNA signatures represent recurrent patterns of copy number alterations across the genome. Here\'s what the top signatures in this analysis mean:

```{r signature_interpretation}
# Select top signatures for interpretation
top_n <- min(5, nrow(signature_priority))
# Select only the simplified name, variance, and description columns
display_priority <- signature_priority[, c("segment", "variance", "description")]
colnames(display_priority) <- c("Signature", "Variance", "Description")
kable(head(display_priority, top_n),
      caption = paste("Top", top_n, "Recommended CNA Signatures"))
```

### What do these signatures mean?

CNA signatures with high variance across samples are the most informative for distinguishing between different cancer subtypes or phenotypes. Here\'s what some common signatures indicate:

- **Amplifications (Amp)**: Indicate extra copies of a gene or region, often leading to increased expression of oncogenes.
- **Deletions (Del)**: Indicate loss of a gene or region, often leading to decreased expression of tumor suppressor genes.
- **Chromosome arm alterations (e.g., 8.q)**: Indicate large-scale changes affecting entire chromosome arms.

### How to interpret enrichment:

- **High positive scores**: Indicate amplification or gain of the region.
- **High negative scores**: Indicate deletion or loss of the region.
- **Scores near zero**: Indicate normal copy number or balanced alterations.

Samples with similar patterns of CNA signatures may share biological characteristics or clinical outcomes.

## Conclusion

This analysis has identified the most informative CNA signatures in the dataset. These signatures can be used for further investigation of cancer subtypes, prediction of clinical outcomes, or identification of potential therapeutic targets.

'
  
  # Write the RMarkdown file
  writeLines(rmd_content, file.path(output_dir, "cna_signature_report.Rmd"))
  
  # Render the RMarkdown file to HTML
  tryCatch({
    cat("Rendering RMarkdown report...\n")
    render(file.path(output_dir, "cna_signature_report.Rmd"), output_format = "html_document", output_dir = output_dir)
    cat("Report generated successfully:", file.path(output_dir, "cna_signature_report.html"), "\n")
  }, error = function(e) {
    cat("Error rendering RMarkdown report:", e$message, "\n")
    cat("You may need to install pandoc to generate HTML reports.\n")
  })
}

# Main execution
cat("Creating visualizations...\n")

# Create barplot of top signatures
create_signature_barplot(signature_priority)

# Create correlation heatmap
create_correlation_heatmap(segment_scores)

# Create PCA plot
create_pca_plot(segment_scores)

# Create sample dendrogram
create_sample_dendrogram(segment_scores)

# Create RMarkdown report
create_rmarkdown_report(segment_scores, signature_priority)

cat("Report generation complete!\n") 