##############################################################################
##########----------- Combine SET 1 and 3 for comparison -----------##########
##############################################################################
rm(list = ls())
set.seed(42)
library(Seurat)
library(tidyverse)
library(sva)
library(edgeR)
library(preprocessCore)
library(patchwork)
library(GSVA)
library(msigdbr)
library(Rtsne)
load("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/SET1.RData")
set1_data <- set1_data
set1_meta <- data.frame(cell = set1_meta,
                        passage =  c(rep("p5", 5), rep("p2", 6)))
load("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/SET3.RData")
set3_data <- set3_data
set3_meta <- data.frame(cell = factor(c(rep("E0771", 5), rep("B6.ILC", 6))),
                        passage = c(rep("p6", 3), rep("p7", 2), rep("p3", 3), rep("p5", 3)))

mouse_human_genes = read.csv("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/RNA_correlation/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0))){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return(output[1])
}
inhouse_gene <- lapply(rownames(set3_data), convert_mouse_to_human)
inhouse_gene[sapply(inhouse_gene, is.null)] <- NA
inhouse_gene <- unlist(inhouse_gene)
set3_data <- set3_data %>%
  mutate(gene = inhouse_gene)
set3_data <- set3_data %>% na.omit() 
id_tbl <- data.frame(ID = rownames(set3_data), 
                     RNA = set3_data$gene,
                     iqr = apply(set3_data[,-12], 1, IQR)) %>%
  group_by(RNA) %>%
  arrange(desc(iqr)) %>%
  slice_head(n = 1) %>%
  ungroup()
set3_data <- set3_data[id_tbl$ID, ]
rownames(set3_data) <- id_tbl$RNA
set3_data <- set3_data[, -12]

common_genes <- Reduce(intersect, list(rownames(set1_data), rownames(set3_data)))
set1_data <- set1_data[common_genes, ]
set3_data <- set3_data[common_genes, ]
icle_data <- icle_data[common_genes, ]

batch <- factor(c(rep(1, ncol(set1_data)), rep(2, ncol(set3_data))))
combined_cells <- as.matrix(cbind(set1_data, set3_data))
qn_normalized <- normalize.quantiles(combined_cells)
dimnames(qn_normalized) <- dimnames(combined_cells)
qn_normalized <- log2(qn_normalized + 1)
adjusted_cells <- sva::ComBat(
  dat = qn_normalized,
  batch = batch
)

icle_indices <- match(colnames(icle_data), colnames(adjusted_cells))
set1_indices <- match(colnames(set1_data), colnames(adjusted_cells))
set3_indices <- match(colnames(set3_data), colnames(adjusted_cells))

icle_harmonized <- adjusted_cells[, icle_indices]
set1_harmonized <- adjusted_cells[, set1_indices]
set3_harmonized <- adjusted_cells[, set3_indices]

# Prepare for t-SNE visualization
## Combine the harmonized datasets
combined_harmonized <- cbind(set1_harmonized, set3_harmonized)

## Create combined metadata
icle_meta <- data.frame(
  cell = icle_group$subtype_intrinsic,  # Assuming icle_group exists with subtype information
  passage = NA,
  dataset = "ICLE"
)

set1_meta$dataset <- "SET1"
set3_meta$dataset <- "SET3"

set1_meta <- set1_meta[, -2]  # Remove the second column
colnames(set1_meta)[1] <- "cell"  # Rename the first column to "cell"

combined_meta <- rbind(
  set1_meta,
  set3_meta
)

## Create labels with passage information for B6.ILC and E0771
plot_labels <- ifelse(
  combined_meta$cell %in% c("B6.ILC", "E0771"),
  paste0(combined_meta$cell, "_", combined_meta$passage),
  as.character(combined_meta$cell)
)

## Perform t-SNE
set.seed(42)  # for reproducibility
tsne_result <- Rtsne::Rtsne(t(combined_harmonized), perplexity = 2, check_duplicates = FALSE)

## Create data frame for plotting
plot_df <- data.frame(
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Group = plot_labels,
  Dataset = combined_meta$dataset
)

## Filter out NA and Normal cells
plot_df_filtered <- plot_df[!is.na(plot_df$Group) & plot_df$Group != "Normal", ]

## Get unique B6.ILC and E0771 passage combinations
b6ilc_passages <- unique(plot_df_filtered$Group[grepl("B6.ILC", plot_df_filtered$Group)])
e0771_passages <- unique(plot_df_filtered$Group[grepl("E0771", plot_df_filtered$Group)])

## Define color palette (extend as needed for all passage types)
color_values <- c(
  "Basal" = "#1f77b4",
  "CL" = "#1f77b4",
  "HER2" = "#ff7f0e",
  "LuminalA" = "#2ca02c",
  "LuminalB" = "#d62728",
  "E0771_p5" = "#BF40BF",
  "E0771_p6" = "#e377c2",
  "E0771_p7" = "#9467bd",
  "B6.ILC_p2" = "#8c564b",
  "B6.ILC_p3" = "#c49c94",
  "B6.ILC_p5" = "#7f7f7f", 
  "B6.ILC_p6" = "#bcbd22",
  "B6.ILC_p7" = "#17becf"
)

## Create the plot
ggplot(plot_df_filtered, aes(x = tSNE1, y = tSNE2, color = Group, shape = Dataset)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(data = subset(plot_df_filtered, grepl("B6.ILC|E0771", Group)), 
               aes(color = Group), 
               level = 0.95,
               linewidth = 1) +
  scale_color_manual(values = color_values) +
  scale_shape_manual(values = c(ICLE = 16, SET1 = 17, SET3 = 15)) +
  theme_minimal() +
  labs(
    title = "Cell samples by passage (t-SNE)",
    subtitle = "Combined ICLE, SET1, and SET3 datasets",
    x = "t-SNE1",
    y = "t-SNE2"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

## Compare SET2 with mice tumor database

# Prepare gene expression boxplots for CDH1 and ESR1
## Combine set1 and set3 harmonized data and metadata
combined_expr <- cbind(set1_harmonized, set3_harmonized)
combined_cell_meta <- rbind(
  transform(set1_meta, cell_passage = paste0(cell, "_", passage)),
  transform(set3_meta, cell_passage = paste0(cell, "_", passage))
)

## Check if the genes exist in our dataset
genes_of_interest <- c("CDH1", "ESR1")
missing_genes <- genes_of_interest[!genes_of_interest %in% rownames(combined_expr)]
if(length(missing_genes) > 0) {
  warning(paste("The following genes are not in the dataset:", paste(missing_genes, collapse=", ")))
}

## Extract expression data for these genes
gene_expr <- data.frame(t(combined_expr[genes_of_interest[genes_of_interest %in% rownames(combined_expr)], ]))
gene_expr$cell_type <- combined_cell_meta$cell
gene_expr$passage <- combined_cell_meta$passage
gene_expr$cell_passage <- combined_cell_meta$cell_passage

## Create a long format dataframe for plotting
gene_expr_long <- gene_expr %>%
  pivot_longer(cols = all_of(genes_of_interest[genes_of_interest %in% rownames(combined_expr)]),
               names_to = "gene",
               values_to = "expression")

## Plot CDH1 expression by passage for each cell type
cdh1_plot <- ggplot(subset(gene_expr_long, gene == "CDH1" & cell_type %in% c("B6.ILC", "E0771")), 
                    aes(x = passage, y = expression, color = cell_type)) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~cell_type, scales = "free_x") +
  labs(title = "CDH1 Expression by Passage", 
       y = "Expression (log2)", 
       x = "Passage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## Plot ESR1 expression by passage for each cell type
esr1_plot <- ggplot(subset(gene_expr_long, gene == "ESR1" & cell_type %in% c("B6.ILC", "E0771")), 
                   aes(x = passage, y = expression, color = cell_type)) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  facet_wrap(~cell_type, scales = "free_x") +
  labs(title = "ESR1 Expression by Passage", 
       y = "Expression (log2)", 
       x = "Passage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

## Display plots side by side
cdh1_plot / esr1_plot + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "",
                 theme = theme(plot.title = element_text(size = 16)))

# Create heatmaps showing all individual samples grouped by passage
## Get the HALLMARK_ESTROGEN_RESPONSE_EARLY gene set
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
estrogen_early_genes <- hallmark_sets %>%
  filter(gs_name %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE")) %>%
  pull(gene_symbol)

## Check which genes from the pathway are in our dataset
common_genes <- intersect(estrogen_early_genes, rownames(combined_expr))
if (length(common_genes) == 0) {
  stop("None of the genes from the pathway are in the dataset")
} else {
  message(paste("Found", length(common_genes), "out of", length(estrogen_early_genes), "estrogen response early genes in the dataset"))
}

## Limit to top variable genes if too many
if (length(common_genes) > 50) {
  gene_var <- apply(combined_expr[common_genes,], 1, var)
  common_genes <- names(sort(gene_var, decreasing = TRUE))[1:50]
  message("Limited to top 50 most variable genes")
}

## Extract expression data for the pathway genes
pathway_expr <- combined_expr[common_genes,]

## Create separate heatmaps for each cell line
# Define cell types to create heatmaps for
cell_types <- c("B6.ILC", "E0771")
heatmap_list <- list()

# Loop through each cell type and create a heatmap
for (cell_type in cell_types) {
  # Filter samples for this cell type
  sample_idx <- which(combined_cell_meta$cell == cell_type)
  
  # Extract expression data for these samples
  cell_expr <- pathway_expr[, sample_idx]
  cell_meta <- combined_cell_meta[sample_idx, ]
  
  # Order samples by passage
  passage_order <- order(as.numeric(gsub("p", "", cell_meta$passage)))
  cell_expr <- cell_expr[, passage_order]
  cell_meta <- cell_meta[passage_order, ]
  
  # Get annotation data for columns (samples)
  sample_anno <- data.frame(
    passage = cell_meta$passage,
    row.names = colnames(cell_expr)
  )
  
  # Create color palette for passages
  passage_colors <- setNames(
    colorRampPalette(c("#fde0dd", "#c51b8a"))(length(unique(sample_anno$passage))),
    sort(unique(sample_anno$passage))
  )
  
  anno_colors <- list(passage = passage_colors)
  
  # Create the heatmap
  heatmap_title <- paste0("Estrogen Response Early Genes in ", cell_type, " by Passage")

  pheatmap(
    cell_expr,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    scale = "row",  # Scale by row to see expression patterns clearly
    cluster_cols = FALSE,  # Keep samples grouped by passage
    cluster_rows = TRUE,   # Cluster genes by expression pattern
    show_rownames = TRUE,
    show_colnames = FALSE, # Too many samples to show column names
    annotation_col = sample_anno,
    annotation_colors = anno_colors,
    fontsize_row = 8,
    fontsize_col = 8,
    main = heatmap_title,
    width = 10,
    height = 12
  )
  
  # Store the heatmap in the list (without saving to file)
  heatmap_list[[cell_type]] <- pheatmap(
    cell_expr,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = sample_anno,
    annotation_colors = anno_colors,
    fontsize_row = 8,
    fontsize_col = 8,
    main = heatmap_title
  )
}

## Display the heatmaps side by side on screen (in addition to saving individual PDFs)
library(gridExtra)
library(grid)

grid.arrange(
  grid.grabExpr(print(heatmap_list[["B6.ILC"]])),
  grid.grabExpr(print(heatmap_list[["E0771"]])),
  ncol = 2,
  top = textGrob("Estrogen Response Early Genes Expression by Passage", 
                gp = gpar(fontsize = 16, font = 2))
)
