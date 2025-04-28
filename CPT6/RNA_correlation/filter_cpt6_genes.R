# Load required libraries
library(methods)

# Load the data
load("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/CPT6.RData")
genes_list <- read.csv("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/list_igv_corrected.csv")

# Get unique genes from the CSV file
unique_genes <- unique(genes_list$gene)
cat("Number of unique genes in CSV file:", length(unique_genes), "\n")

# Check which genes from the CSV are present in the combined_harmonized matrix
common_genes <- intersect(unique_genes, rownames(combined_harmonized))
cat("Number of genes found in both datasets:", length(common_genes), "\n")
cat("Genes found in both datasets:", paste(common_genes, collapse=", "), "\n")

# Filter the combined_harmonized matrix to keep only the common genes
filtered_data <- combined_harmonized[common_genes, ]
cat("Dimensions of filtered data:", dim(filtered_data)[1], "genes x", dim(filtered_data)[2], "samples\n")

# Save the filtered data
save(filtered_data, file="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/CPT6_filtered.RData")
cat("Filtered data saved to /bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/CPT6_filtered.RData\n")
