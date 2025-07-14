#!/usr/bin/env Rscript

# Correlation Analysis for Significant Signatures
# This script analyzes the correlations between the four significant signatures:
# - chr1_168438419_168880930__Amp_ (XCL1, XCL2, DPT)
# - chr14_56307190_56415610__Amp_ (No genes)
# - chr4_15313739_17632477__Del_ (C1QTNF7, FBXL5, CD38, FGFBP1, TAPT1, LDB2)
# - chr9_21489625_22474701__Del____CDKN2A (CDKN2A, CDKN2B)

# Load required libraries
library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)
library(gridExtra)

# Set paths
output_dir <- "../output"
survival_dir <- "."

# Load segment scores and signature priority
cat("Loading segment scores and signature priority...\n")
load(file.path(output_dir, "segment_scores.rda"))
signature_priority <- read.csv(file.path(output_dir, "signature_priority.csv"), stringsAsFactors = FALSE)

# Define the four significant signatures
significant_signatures <- c(
  "chr1_168438419_168880930__Amp_",
  "chr14_56307190_56415610__Amp_",
  "chr4_15313739_17632477__Del_",
  "chr9_21489625_22474701__Del____CDKN2A"
)

# Find indices of these signatures in segment_scores
sig_indices <- list()
for (sig_name in significant_signatures) {
  # Find the original segment name in signature_priority
  original_segment <- NULL
  for (i in 1:nrow(signature_priority)) {
    if (gsub("[^a-zA-Z0-9]", "_", signature_priority$segment[i]) == sig_name) {
      original_segment <- signature_priority$original_segment[i]
      break
    }
  }
  
  if (!is.null(original_segment)) {
    idx <- which(rownames(segment_scores) == original_segment)
    if (length(idx) > 0) {
      sig_indices[[sig_name]] <- idx
    }
  }
}

# Check if we found all signatures
if (length(sig_indices) < length(significant_signatures)) {
  missing_sigs <- setdiff(significant_signatures, names(sig_indices))
  cat("Warning: Could not find the following signatures:", paste(missing_sigs, collapse=", "), "\n")
}

# Extract scores for significant signatures
sig_scores <- data.frame(matrix(ncol = length(sig_indices), nrow = ncol(segment_scores)))
colnames(sig_scores) <- names(sig_indices)

for (i in 1:length(sig_indices)) {
  sig_name <- names(sig_indices)[i]
  sig_scores[, i] <- as.numeric(segment_scores[sig_indices[[sig_name]], ])
}

# Add sample names
sig_scores$Sample <- colnames(segment_scores)

# Create more readable labels for the signatures
sig_labels <- c(
  "chr1_168438419_168880930__Amp_" = "Chr1 Amp\n(XCL1, XCL2, DPT)",
  "chr14_56307190_56415610__Amp_" = "Chr14 Amp\n(No genes)",
  "chr4_15313739_17632477__Del_" = "Chr4 Del\n(C1QTNF7, FBXL5, CD38, FGFBP1, TAPT1, LDB2)",
  "chr9_21489625_22474701__Del____CDKN2A" = "Chr9 Del\n(CDKN2A, CDKN2B)"
)

# Calculate correlation matrix
cor_matrix <- cor(sig_scores[, names(sig_indices)], method = "pearson")
cat("Correlation matrix:\n")
print(cor_matrix)

# Find the actual range of correlations to set a more focused scale
min_cor <- min(cor_matrix[cor_matrix < 1])  # Exclude diagonal (1's)
max_cor <- max(cor_matrix)

# Use fixed scale from 0.6 to 1 as requested
min_cor_adj <- 0.6
max_cor_adj <- 1.0

cat("Actual correlation range:", min_cor, "to", max_cor, "\n")
cat("Using fixed correlation range for visualization:", min_cor_adj, "to", max_cor_adj, "\n")

# Create correlation plot
pdf(file.path(survival_dir, "four_signature_correlations.pdf"), width = 10, height = 8)

# Plot correlation matrix with corrplot - using the fixed range 0.6 to 1
corrplot(cor_matrix, method = "circle", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         col = colorRampPalette(c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"))(100),  # Blue-White-Red diverging palette
         title = paste0("Pearson Correlation Between Significant Signatures\nRange: 0.6 to 1.0"),
         mar = c(0, 0, 2, 0),
         is.corr = FALSE,  # Allow non-correlation matrix values
         zlim = c(0.6, 1.0),  # Fixed range from 0.6 to 1.0
         addCoef.col = "black",  # Add correlation coefficients in black
         number.cex = 0.9,  # Size of correlation text
         diag = TRUE)  # Show diagonal (self-correlations)

# Create a more detailed correlation plot with ggplot2
cor_melted <- melt(cor_matrix)
colnames(cor_melted) <- c("Signature1", "Signature2", "Correlation")

# Replace signature names with more readable labels
cor_melted$Signature1 <- factor(cor_melted$Signature1, 
                               levels = names(sig_indices),
                               labels = sig_labels[names(sig_indices)])
cor_melted$Signature2 <- factor(cor_melted$Signature2, 
                               levels = names(sig_indices),
                               labels = sig_labels[names(sig_indices)])

# Create a heatmap with ggplot2
p <- ggplot(cor_melted, aes(x = Signature1, y = Signature2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"),  # Blue-White-Red diverging palette
                     limits = c(0.6, 1.0),  # Fixed range from 0.6 to 1.0
                     name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),  # Rotate x-axis labels 90 degrees
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  coord_fixed() +
  geom_text(aes(label = sprintf("%.3f", Correlation), 
                color = ifelse(Correlation == 1, "black", "black")),  # Ensure diagonal values are visible
            size = 4) +
  scale_color_identity() +  # Use the colors directly
  labs(title = "Correlation Heatmap of Significant Signatures (Range: 0.6-1.0)",
       x = "", y = "")

print(p)

# Create scatter plots for each pair of signatures
plots <- list()
plot_idx <- 1

for (i in 1:(length(sig_indices)-1)) {
  for (j in (i+1):length(sig_indices)) {
    sig1 <- names(sig_indices)[i]
    sig2 <- names(sig_indices)[j]
    
    # Calculate correlation
    cor_val <- cor(sig_scores[, sig1], sig_scores[, sig2], method = "pearson")
    
    # Determine color for correlation text based on value (matching heatmap gradient)
    # Map correlation value from 0.6-1.0 range to 0-1 range for color interpolation
    color_scale <- colorRampPalette(c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"))(100)
    color_idx <- round((cor_val - 0.6) / 0.4 * 99) + 1  # Scale to 1-100 index
    color_idx <- max(1, min(100, color_idx))  # Ensure within bounds
    text_color <- color_scale[color_idx]
    
    # Create scatter plot with matching correlation annotation style
    plots[[plot_idx]] <- ggplot(sig_scores, aes_string(x = sig1, y = sig2)) +
      geom_point(color = "#0072B2", alpha = 0.7) +
      geom_smooth(method = "lm", color = "#D55E00", se = TRUE) +
      labs(
        title = paste0("Pearson Correlation: ", sprintf("%.3f", cor_val)),
        x = sig_labels[sig1],
        y = sig_labels[sig2]
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 12, face = "bold", color = text_color),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better visibility
      ) +
      # Add correlation range annotation
      annotate("text", x = min(sig_scores[, sig1], na.rm = TRUE) + 0.1 * diff(range(sig_scores[, sig1], na.rm = TRUE)), 
               y = max(sig_scores[, sig2], na.rm = TRUE) - 0.1 * diff(range(sig_scores[, sig2], na.rm = TRUE)),
               label = "Range: 0.6-1.0", size = 3, hjust = 0)
    
    plot_idx <- plot_idx + 1
  }
}

# Arrange all scatter plots in a grid
do.call(grid.arrange, c(plots, ncol = 2))

dev.off()

# Create a separate file with just the scatter plots
pdf(file.path(survival_dir, "signature_scatter_plots.pdf"), width = 10, height = 8)
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

# Create an additional visualization with a more focused range
# This will create a version with an even narrower range to highlight subtle differences
# Use a narrower range within the 0.6-1.0 scale
narrow_min <- 0.8  # Narrower minimum
narrow_max <- 1.0  # Keep maximum at 1.0

pdf(file.path(survival_dir, "four_signature_correlations_narrow_range.pdf"), width = 10, height = 8)

# Plot correlation matrix with corrplot - using a very narrow range
corrplot(cor_matrix, method = "circle", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         col = colorRampPalette(c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"))(100),
         title = paste0("Pearson Correlation Between Significant Signatures\nNarrow Range: 0.8 to 1.0",
                        " (for highlighting subtle differences)"),
         mar = c(0, 0, 2, 0),
         is.corr = FALSE,
         zlim = c(0.8, 1.0),  # Narrower fixed range from 0.8 to 1.0
         addCoef.col = "black",
         number.cex = 0.9,
         diag = TRUE)  # Show diagonal (self-correlations)

# Create a heatmap with ggplot2 using the narrow range
p_narrow <- ggplot(cor_melted, aes(x = Signature1, y = Signature2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"),
                     limits = c(0.8, 1.0),  # Narrower fixed range from 0.8 to 1.0
                     name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),  # Rotate x-axis labels 90 degrees
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  coord_fixed() +
  geom_text(aes(label = sprintf("%.3f", Correlation),
                color = ifelse(Correlation == 1, "black", "black")),  # Ensure diagonal values are visible
            size = 4) +
  scale_color_identity() +  # Use the colors directly
  labs(title = "Correlation Heatmap with Narrow Range (0.8-1.0)",
       subtitle = "Highlights subtle differences in correlation values",
       x = "", y = "")

print(p_narrow)

# Create a special version that highlights the diagonal values
pdf(file.path(survival_dir, "four_signature_correlations_with_diagonal.pdf"), width = 10, height = 8)

# Create a heatmap that specifically highlights the diagonal values
p_diagonal <- ggplot(cor_melted, aes(x = Signature1, y = Signature2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#2166AC", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#B2182B"),
                     limits = c(0.6, 1.0),  # Fixed range from 0.6 to 1.0
                     name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),  # Rotate x-axis labels 90 degrees
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  coord_fixed() +
  # Add text with special formatting for diagonal values
  geom_text(aes(label = sprintf("%.3f", Correlation),
                fontface = ifelse(Correlation == 1, "bold", "plain"),
                color = ifelse(Correlation == 1, "darkred", "black")),
            size = ifelse(cor_melted$Correlation == 1, 5, 4)) +  # Larger text for diagonal
  scale_color_identity() +  # Use the colors directly
  labs(title = "Correlation Heatmap with Highlighted Diagonal Values (Range: 0.6-1.0)",
       subtitle = "Diagonal values (1.000) represent self-correlations",
       x = "", y = "")

print(p_diagonal)
dev.off()

cat("Correlation analysis complete! Results are available in the survival directory.\n")
cat("Three correlation visualizations were created:\n")
cat("1. four_signature_correlations.pdf - with fixed range 0.6 to 1.0\n")
cat("2. signature_scatter_plots.pdf - scatter plots for each pair of signatures\n")
cat("3. four_signature_correlations_narrow_range.pdf - with a narrower range to highlight subtle differences\n")
cat("4. four_signature_correlations_with_diagonal.pdf - with highlighted diagonal values (self-correlations)\n") 