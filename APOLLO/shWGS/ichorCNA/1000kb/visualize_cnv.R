#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Read the data
gene_cnv <- read.table("results/gene_cnv_summary.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
segments <- read.table("results/all_segments.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Function to create chromosome-level plot
plotChromosomeCNV <- function(segments_data) {
    # Remove 'chr' prefix for ordering
    segments_data$chr_num <- gsub("chr", "", segments_data$chr)
    segments_data$chr_num <- factor(segments_data$chr_num, 
                                  levels = c(1:22, "X", "Y"))
    
    # Calculate frequency of CNV events per chromosome
    cnv_freq <- segments_data %>%
        group_by(chr, chr_num) %>%
        summarize(
            amp_freq = sum(median > 2.5) / n(),
            del_freq = sum(median < 1.5) / n()
        )
    
    # Create the plot
    p <- ggplot() +
        geom_bar(data=cnv_freq, aes(x=chr_num, y=amp_freq), stat="identity", fill="red", alpha=0.5) +
        geom_bar(data=cnv_freq, aes(x=chr_num, y=-del_freq), stat="identity", fill="blue", alpha=0.5) +
        theme_minimal() +
        coord_flip() +
        labs(x="Chromosome", y="Frequency", 
             title="Chromosome-level CNV Distribution",
             subtitle="Red: Amplification (CN>2.5), Blue: Deletion (CN<1.5)") +
        theme(axis.text.y = element_text(size=10))
    
    pdf("cnv_chromosome_plot.pdf", width=10, height=8)
    print(p)
    dev.off()
}

# Function to create bubble plot
plotCNVBubble <- function(gene_cnv_data) {
    # Calculate frequency of CNV events per gene
    gene_freq <- gene_cnv_data %>%
        group_by(Gene) %>%
        summarize(
            amp_freq = sum(CN > 2.5, na.rm=TRUE) / n(),
            del_freq = sum(CN < 1.5, na.rm=TRUE) / n(),
            total_freq = amp_freq + del_freq
        ) %>%
        arrange(desc(total_freq)) %>%
        head(30)  # Top 30 genes
    
    # Create bubble plot
    p <- ggplot(gene_freq, aes(x=amp_freq, y=del_freq)) +
        geom_point(aes(size=total_freq), alpha=0.6) +
        geom_text(aes(label=Gene), size=3, vjust=2) +
        theme_minimal() +
        labs(x="Amplification Frequency", y="Deletion Frequency",
             title="CNV Frequency Bubble Plot",
             size="Total Frequency")
    
    pdf("cnv_bubble_plot.pdf", width=12, height=10)
    print(p)
    dev.off()
}

# Function to create oncoplot-style heatmap
plotCNVHeatmap <- function(gene_cnv_data, top_n_genes=30) {
    # Get top frequently altered genes
    top_genes <- gene_cnv_data %>%
        group_by(Gene) %>%
        summarize(
            altered = sum(CN > 2.5 | CN < 1.5, na.rm=TRUE) / n()
        ) %>%
        arrange(desc(altered)) %>%
        head(top_n_genes) %>%
        pull(Gene)
    
    # Prepare matrix for heatmap - handle duplicates by taking mean
    cnv_matrix <- gene_cnv_data %>%
        filter(Gene %in% top_genes) %>%
        group_by(Sample, Gene) %>%
        summarize(CN = mean(CN, na.rm=TRUE), .groups='drop') %>%
        spread(Sample, CN) %>%
        as.data.frame()
    
    rownames(cnv_matrix) <- cnv_matrix$Gene
    cnv_matrix$Gene <- NULL
    cnv_matrix <- as.matrix(cnv_matrix)
    
    # Define colors
    col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
    
    # Create heatmap
    pdf("cnv_heatmap.pdf", width=15, height=10)
    Heatmap(cnv_matrix,
            name = "Copy Number",
            col = col_fun,
            show_column_names = FALSE,
            column_title = "Samples",
            row_title = "Genes",
            cluster_columns = TRUE,
            cluster_rows = TRUE)
    dev.off()
}

# Execute the visualizations
plotChromosomeCNV(segments)
plotCNVBubble(gene_cnv)
plotCNVHeatmap(gene_cnv) 