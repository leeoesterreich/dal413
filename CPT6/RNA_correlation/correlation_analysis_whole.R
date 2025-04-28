# Load required libraries
library(methods)
library(ggplot2)

# Load the filtered data
load("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/CPT6_filtered.RData")
genes_list <- read.csv("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/list_igv_corrected.csv")

# Calculate mean expression for all 12 columns
cat("Calculating mean expression for all 12 columns...\n")
mean_expression <- rowMeans(filtered_data)

# Convert to log2(CPM)
# Note: The data might already be log-transformed, but we'll ensure it's in log2(CPM)
# If the data is already log2-transformed, we don't need to apply log2 again
# Check if data appears to be log-transformed already
is_log_transformed <- all(mean_expression < 50) # Heuristic: most log-transformed expression values are < 50

if (!is_log_transformed) {
  cat("Converting to log2(CPM)...\n")
  # Convert to CPM if not already in that format
  # CPM = (count / total_count) * 1,000,000
  # We'll assume the data is already normalized if it's log-transformed
  log2_mean_expression <- log2(mean_expression + 1) # Add 1 to avoid log(0)
} else {
  cat("Data appears to be already log-transformed, using as is...\n")
  log2_mean_expression <- mean_expression
}

# Create a data frame with gene names and log2(CPM) mean expression
expression_df <- data.frame(
  gene = rownames(filtered_data),
  log2_mean_expression = log2_mean_expression
)

# Merge with the log2ratio from the CSV file
cat("Merging with log2ratio data...\n")
merged_data <- merge(expression_df, genes_list[, c("gene", "log2ratio")], by = "gene")

# Calculate Spearman correlation
correlation <- cor.test(merged_data$log2ratio, merged_data$log2_mean_expression, method = "spearman")
cat("Spearman correlation coefficient:", correlation$estimate, "\n")
cat("p-value:", correlation$p.value, "\n")

# Create a scatter plot
cat("Creating correlation plot...\n")
p <- ggplot(merged_data, aes(x = log2ratio, y = log2_mean_expression)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = paste("Spearman Correlation:", round(correlation$estimate, 3), 
                  "(p-value:", format.pval(correlation$p.value, digits = 3), ")"),
    x = "Gene log2ratio (CNV)",
    y = "log2(CPM) Mean Expression (All 12 samples)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

# Specifically annotate CDKN2A, CDKN2B, CDH1, and PTEN
genes_to_annotate <- c("CDKN2A", "CDKN2B", "CDH1", "PTEN")
annotate_data <- merged_data[merged_data$gene %in% genes_to_annotate, ]

# Make sure all requested genes are found
if (nrow(annotate_data) < length(genes_to_annotate)) {
  missing_genes <- setdiff(genes_to_annotate, annotate_data$gene)
  cat("Warning: The following genes were not found in the data:", paste(missing_genes, collapse=", "), "\n")
}

# Find other significant genes to annotate (top 5 by difference)
merged_data$diff <- abs(scale(merged_data$log2ratio) - scale(merged_data$log2_mean_expression))
# Exclude already annotated genes
other_genes <- merged_data[!merged_data$gene %in% genes_to_annotate, ]
top_genes <- other_genes[order(other_genes$diff, decreasing = TRUE), ][1:5, ]

# Create a custom annotation data frame with adjusted positions for CDKN2A and CDKN2B
custom_annotate_data <- annotate_data
# Find the indices of CDKN2A and CDKN2B
cdkn2a_idx <- which(custom_annotate_data$gene == "CDKN2A")
cdkn2b_idx <- which(custom_annotate_data$gene == "CDKN2B")

# Create columns for position adjustments
custom_annotate_data$y_adj <- custom_annotate_data$log2_mean_expression
custom_annotate_data$x_adj <- custom_annotate_data$log2ratio
custom_annotate_data$hjust_adj <- -0.2  # Default hjust for all labels

if (length(cdkn2a_idx) > 0) {
  # Move CDKN2A down and position horizontally
  custom_annotate_data$y_adj[cdkn2a_idx] <- custom_annotate_data$log2_mean_expression[cdkn2a_idx] - 0.3
  custom_annotate_data$x_adj[cdkn2a_idx] <- custom_annotate_data$log2ratio[cdkn2a_idx] + 0.05  # Now slightly to the right
  custom_annotate_data$hjust_adj[cdkn2a_idx] <- 0.5  # Center the text horizontally
}
if (length(cdkn2b_idx) > 0) {
  # Move CDKN2B up and position horizontally
  custom_annotate_data$y_adj[cdkn2b_idx] <- custom_annotate_data$log2_mean_expression[cdkn2b_idx] + 0.3
  custom_annotate_data$x_adj[cdkn2b_idx] <- custom_annotate_data$log2ratio[cdkn2b_idx] + 0.05  # Now slightly to the right
  custom_annotate_data$hjust_adj[cdkn2b_idx] <- 0.5  # Center the text horizontally
}

# Add gene labels for the specified genes (in black) with adjusted positions
p <- p + geom_text(data = custom_annotate_data, 
                  aes(x = x_adj, y = y_adj, label = gene), 
                  hjust = custom_annotate_data$hjust_adj, vjust = 0.5, 
                  size = 3.5, color = "black", fontface = "bold")

# Add gene labels for other significant genes
p <- p + geom_text(data = top_genes, aes(label = gene), 
                  hjust = -0.2, vjust = 0.5, size = 3, color = "black")

# Highlight the specifically annotated points
p <- p + geom_point(data = annotate_data, aes(x = log2ratio, y = log2_mean_expression), 
                   color = "red", size = 3, alpha = 0.7)

# Save the plot
pdf_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/expression_log2ratio_correlation_whole_annotated.pdf"
pdf(pdf_file, width = 8, height = 6)
print(p)
dev.off()

# Also save as PNG for easier viewing
png_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/expression_log2ratio_correlation_whole_annotated.png"
png(png_file, width = 800, height = 600)
print(p)
dev.off()

cat("Plot saved to:", pdf_file, "and", png_file, "\n")

# Save the correlation results
result_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/correlation_results_whole.csv"
write.csv(merged_data, result_file, row.names = FALSE)
cat("Correlation data saved to:", result_file, "\n")

# Print the data for the annotated genes
cat("\nData for specifically annotated genes:\n")
print(annotate_data[, c("gene", "log2ratio", "log2_mean_expression")])

# Print the data for other significant genes
cat("\nData for other significant genes:\n")
print(top_genes[, c("gene", "log2ratio", "log2_mean_expression")]) 