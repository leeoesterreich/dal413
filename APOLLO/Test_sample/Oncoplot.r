library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(grid)
 
file.choose()
> [1] "/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/CNVkit/tumor_sample1.annotated.cnr"
 
# Set the directory path based on your located file
cnr_files_path <- "/bgfs/alee/LO_LAB/Personal/Daisong/Test_sample/CNVkit"
 
# List all .cnr files in the directory
cnr_files <- list.files(path = cnr_files_path, pattern = "\\.cnr$", full.names = TRUE, ignore.case = TRUE)
cnr_files

# Load all .cnr files into a combined data frame with sample identifiers
cnr_data_list <- lapply(cnr_files, function(file) {
    sample_name <- tools::file_path_sans_ext(basename(file))
    fread(file, sep = "\t") %>%
        mutate(sample = sample_name)
})
combined_cnr_data <- bind_rows(cnr_data_list)

# Filter for High Amplification (HLAMP)
amplified_genes <- combined_cnr_data %>%
     filter(if_any(ends_with("Corrected_Call"), ~ . == "HLAMP")) %>% 
     separate_rows(gene, sep = ",")  # Expand rows if multiple genes are comma-separated

# Define genes of interest
genes_of_interest <- c("RAD21", "MYC", "PREX2", "PRDM14", "AGO2", "ELOC", "KMT2C", 
                       "NSD3", "FGFR1", "RECQL4", "SOX17", "NBN", "NTRK3", "DDR2", 
                       "HGF", "AURKA", "EGFR", "LYN", "RAD51D", "GNAS", "PTPRT", 
                       "MSI2", "PRKAR1A", "NUF2", "TERT", "RICTOR", "MCL1", "RIT1", 
                       "NTRK1", "PAK1", "VEGFA", "MDM4", "EIF4A2", "IGF1R", "BRCA1", 
                       "SMYD3", "CCND2", "GRIN2A", "CDK12", "ERBB2", "CDK6", "INHBA", 
                       "H3F3B", "CDK4", "MDM2", "RNF43", "BRIP1", "AXIN2", "SDHC", "ESR1")

# Filter for the genes of interest
filtered_amplified_genes <- amplified_genes %>%
    filter(gene %in% genes_of_interest) %>%
    distinct(gene, sample)  # Remove duplicates to get unique amplifications

# Create a binary matrix where rows are genes and columns are samples
oncoplot_data <- filtered_amplified_genes %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = sample, values_from = value, values_fill = 0)

# Convert to matrix format for the heatmap
oncoplot_matrix <- as.matrix(oncoplot_data[,-1])
rownames(oncoplot_matrix) <- oncoplot_data$gene


# Assuming `oncoplot_matrix` is a binary matrix (genes x samples) with 1 for amplified and 0 otherwise

# Calculate gene amplification frequencies in percentages
gene_amplifications <- rowSums(oncoplot_matrix)
gene_percentage <- (gene_amplifications / ncol(oncoplot_matrix)) * 100

# Sort genes based on amplification percentage (high to low)
sorted_genes <- names(sort(gene_percentage, decreasing = TRUE))
oncoplot_matrix <- oncoplot_matrix[sorted_genes, ]

# Define annotations
# Left annotation for amplification percentage
left_annotation <- rowAnnotation(
    `Frequency (%)` = anno_text(paste0(round(gene_percentage[sorted_genes], 0), "%"), 
                                just = "right", 
                                gp = gpar(fontsize = 8, col = "black")),
    annotation_name_side = "top"
)

# Right annotation for the number of amplified samples per gene (bar plot)
right_annotation <- rowAnnotation(
    `Amplified Samples` = anno_barplot(gene_amplifications[sorted_genes], 
                                       gp = gpar(fill = "red"), 
                                       border = FALSE),
    annotation_name_side = "bottom"
)

# Top annotation (optional) for sample-wise amplification counts
sample_amplifications <- colSums(oncoplot_matrix)
top_annotation <- HeatmapAnnotation(
    `Sample Amplifications` = anno_barplot(sample_amplifications, 
                                           gp = gpar(fill = "red"),
                                           border = FALSE)
)

# Define color scale for amplification
amplification_colors <- c("0" = "lightgrey", "1" = "red")

# For a PDF
pdf("oncoplot_output.pdf", width = 10, height = 12)

# For a PNG
png("oncoplot_output.png", width = 1000, height = 1200, res = 150)

# Assuming `oncoplot_matrix` and annotations are already set up as per previous code
Heatmap(
    oncoplot_matrix,
    name = "Amplifications",
    col = amplification_colors,
    show_row_names = TRUE,              # Show gene names on the right side
    row_names_side = "right",           # Place gene names on the right side
    show_column_names = TRUE,           # Show sample names at the bottom
    column_names_side = "bottom",
    cluster_rows = FALSE,               # Disable clustering for genes
    cluster_columns = FALSE,            # Disable clustering for samples
    left_annotation = left_annotation,  # Left side shows amplification frequency percentage
    right_annotation = right_annotation, # Right side shows number of amplified samples per gene
    top_annotation = top_annotation,    # Optional: Top shows sample frequency bar plot
    heatmap_legend_param = list(
        title = "Amplification", 
        at = c(0, 1),                   # Positions for 'No Amplification' and 'Amplification'
        labels = c("No Amplification", "Amplification")  # Matching labels
    ),
    row_title = "Genes",
    column_title = "Samples"
)

dev.off()
