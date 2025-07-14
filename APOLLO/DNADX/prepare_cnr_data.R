#!/usr/bin/env Rscript

# Script to prepare CNR data for DNADX analysis
# This script converts CNR files to a format compatible with DNADX

# Load required libraries
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)

# Set paths
cnr_dir <- "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Annotated_cnr_file_500Kb"
output_dir <- "./output"
dnadx_dir <- "./DNA-based-predictors-of-non-genetic-cancer-phenotypes"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load helper functions from DNADX
source(file.path(dnadx_dir, "Rscript", "helper.R"))

# Function to read CNR file and extract copy number data
read_cnr_file <- function(file_path) {
  # Read CNR file
  cnr_data <- fread(file_path)
  
  # Extract sample name from file name
  sample_name <- gsub(".annotated.cnr", "", basename(file_path))
  
  # Check if the required columns exist
  call_col <- paste0(sample_name, "_29.Corrected_Call")
  cn_col <- paste0(sample_name, "_29.Corrected_Copy_Number")
  
  if (!call_col %in% colnames(cnr_data) || !cn_col %in% colnames(cnr_data)) {
    cat("Warning: Required columns not found in", file_path, "\n")
    # Try to find alternative column names
    call_cols <- grep("Call", colnames(cnr_data), value = TRUE)
    cn_cols <- grep("Copy_Number", colnames(cnr_data), value = TRUE)
    
    if (length(call_cols) > 0 && length(cn_cols) > 0) {
      call_col <- call_cols[1]
      cn_col <- cn_cols[1]
      cat("Using alternative columns:", call_col, "and", cn_col, "\n")
    } else {
      cat("Error: Cannot find suitable columns in", file_path, "\n")
      return(NULL)
    }
  }
  
  # Create a data frame first to ensure all components have the same length
  gr_df <- data.frame(
    seqnames = cnr_data$chromosome,
    start = cnr_data$start,
    end = cnr_data$end,
    log2 = cnr_data$log2,
    call = cnr_data[[call_col]],
    copy_number = cnr_data[[cn_col]],
    gene = cnr_data$gene
  )
  
  # Create a GRanges object for genomic regions
  gr <- GRanges(
    seqnames = gr_df$seqnames,
    ranges = IRanges(start = gr_df$start, end = gr_df$end)
  )
  
  # Add metadata columns
  mcols(gr) <- DataFrame(
    log2 = gr_df$log2,
    call = gr_df$call,
    copy_number = gr_df$copy_number,
    gene = gr_df$gene
  )
  
  # Return the GRanges object with sample name
  return(list(sample = sample_name, data = gr))
}

# Function to convert CNR data to gene-level copy number
cnr_to_gene_level <- function(cnr_list) {
  # Remove any NULL entries from the list
  cnr_list <- cnr_list[!sapply(cnr_list, is.null)]
  
  if (length(cnr_list) == 0) {
    stop("No valid CNR data found")
  }
  
  # Extract gene information from CNR files
  gene_info <- lapply(cnr_list, function(x) {
    if (is.null(x) || is.null(x$data)) return(NULL)
    
    data.frame(
      chromosome = as.character(seqnames(x$data)),
      start = start(x$data),
      end = end(x$data),
      gene = mcols(x$data)$gene,
      log2 = mcols(x$data)$log2,
      call = mcols(x$data)$call,
      copy_number = mcols(x$data)$copy_number,
      sample = x$sample,
      stringsAsFactors = FALSE
    )
  })
  
  # Remove any NULL entries
  gene_info <- gene_info[!sapply(gene_info, is.null)]
  
  if (length(gene_info) == 0) {
    stop("No valid gene information extracted from CNR files")
  }
  
  # Combine all gene information
  all_gene_info <- do.call(rbind, gene_info)
  
  # Handle missing or NA values
  all_gene_info$gene[is.na(all_gene_info$gene)] <- "-"
  all_gene_info$log2[is.na(all_gene_info$log2)] <- 0
  all_gene_info$copy_number[is.na(all_gene_info$copy_number)] <- 2  # Default diploid
  
  # Split gene column (which may contain multiple genes)
  gene_data <- all_gene_info %>%
    separate_rows(gene, sep = ",") %>%
    filter(gene != "-") %>%  # Remove regions with no gene annotation
    filter(gene != "") %>%   # Remove empty gene names
    group_by(gene, sample) %>%
    summarize(
      mean_log2 = mean(log2, na.rm = TRUE),
      mean_copy_number = mean(copy_number, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = gene,
      names_from = sample,
      values_from = c(mean_log2, mean_copy_number)
    )
  
  # Create separate matrices for log2 and copy number
  log2_cols <- grep("mean_log2_", colnames(gene_data), value = TRUE)
  cn_cols <- grep("mean_copy_number_", colnames(gene_data), value = TRUE)
  
  log2_matrix <- as.matrix(gene_data[, log2_cols, drop = FALSE])
  rownames(log2_matrix) <- gene_data$gene
  colnames(log2_matrix) <- gsub("mean_log2_", "", log2_cols)
  
  cn_matrix <- as.matrix(gene_data[, cn_cols, drop = FALSE])
  rownames(cn_matrix) <- gene_data$gene
  colnames(cn_matrix) <- gsub("mean_copy_number_", "", cn_cols)
  
  return(list(log2_matrix = log2_matrix, cn_matrix = cn_matrix))
}

# Main execution
cat("Reading CNR files...\n")
cnr_files <- list.files(cnr_dir, pattern = "*.annotated.cnr", full.names = TRUE)

# Process a subset of files for testing if there are many files
if (length(cnr_files) > 10) {
  cat("Found", length(cnr_files), "CNR files. Processing all files...\n")
}

cnr_data_list <- lapply(cnr_files, function(file) {
  cat("Processing", basename(file), "...\n")
  tryCatch({
    read_cnr_file(file)
  }, error = function(e) {
    cat("Error processing", file, ":", e$message, "\n")
    return(NULL)
  })
})

cat("Converting to gene-level copy number...\n")
gene_level_data <- cnr_to_gene_level(cnr_data_list)

# Save the processed data
cat("Saving processed data...\n")
save(gene_level_data, file = file.path(output_dir, "gene_level_cn_data.rda"))

cat("Data preparation complete!\n") 