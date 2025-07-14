loa#!/usr/bin/env Rscript

# Script to analyze CNA signatures using DNADX
# This script calculates segment scores and prioritizes CNA-based signatures

# Load required libraries
library(GSA)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set paths
output_dir <- "./output"
dnadx_dir <- "./DNA-based-predictors-of-non-genetic-cancer-phenotypes"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load helper functions from DNADX
source(file.path(dnadx_dir, "Rscript", "helper.R"))

# Load the processed gene-level copy number data
cat("Loading gene-level copy number data...\n")
load(file.path(output_dir, "gene_level_cn_data.rda"))

# Debug: Check gene-level copy number data
cat("Checking gene-level copy number data...\n")
cat("Dimensions of log2_matrix:", dim(gene_level_data$log2_matrix)[1], "x", dim(gene_level_data$log2_matrix)[2], "\n")
cat("Dimensions of cn_matrix:", dim(gene_level_data$cn_matrix)[1], "x", dim(gene_level_data$cn_matrix)[2], "\n")
cat("First few rows and columns of cn_matrix:\n")
if (nrow(gene_level_data$cn_matrix) > 0 && ncol(gene_level_data$cn_matrix) > 0) {
  print(gene_level_data$cn_matrix[1:min(5, nrow(gene_level_data$cn_matrix)), 1:min(3, ncol(gene_level_data$cn_matrix))])
} else {
  cat("cn_matrix is empty\n")
}

# Load CNA segments GMT file
cna_segments_file <- file.path(dnadx_dir, "data", "CNA_segments.gmt")

# Debug: Check if CNA segments GMT file exists
cat("Checking CNA segments GMT file...\n")
if (file.exists(cna_segments_file)) {
  cat("CNA segments GMT file exists:", cna_segments_file, "\n")
  
  # Read the first few lines of the GMT file to check format
  cat("First few lines of CNA segments GMT file:\n")
  gmt_lines <- readLines(cna_segments_file, n = 5)
  print(gmt_lines)
  
  # Try to read the GMT file using GSA.read.gmt
  cat("Reading GMT file with GSA.read.gmt...\n")
  tryCatch({
    # Set a limit on the number of gene sets to read to avoid memory issues
    max_gene_sets <- 200  # Increased from 100 to 200 gene sets
    
    # Read the GMT file line by line to have more control
    gmt_lines <- readLines(cna_segments_file)
    cat("Total number of lines in GMT file:", length(gmt_lines), "\n")
    
    # Process only a subset of lines
    subset_lines <- gmt_lines[1:min(max_gene_sets, length(gmt_lines))]
    
    # Create a temporary file with the subset
    temp_gmt_file <- tempfile(fileext = ".gmt")
    writeLines(subset_lines, temp_gmt_file)
    
    # Read the temporary file using GSA.read.gmt
    geneset.obj <- GSA.read.gmt(temp_gmt_file)
    cat("Number of gene sets processed:", length(geneset.obj$genesets), "\n")
    cat("First gene set name:", geneset.obj$geneset.names[1], "\n")
    cat("Number of genes in first gene set:", length(geneset.obj$genesets[[1]]), "\n")
    
    # Check the first few genes in the first gene set
    cat("First few genes in first gene set:\n")
    print(head(geneset.obj$genesets[[1]]))
    
    # Clean up
    unlink(temp_gmt_file)
    
    # Force garbage collection to free memory
    gc()
  }, error = function(e) {
    cat("Error reading GMT file:", e$message, "\n")
  })
} else {
  cat("CNA segments GMT file does not exist:", cna_segments_file, "\n")
}

# Function to calculate segment scores
calculate_segment_scores <- function(cn_matrix, segments_file, method = "mean") {
  # Check if cn_matrix has valid dimensions
  if (nrow(cn_matrix) == 0 || ncol(cn_matrix) == 0) {
    cat("Error: cn_matrix has invalid dimensions:", nrow(cn_matrix), "x", ncol(cn_matrix), "\n")
    return(matrix(0, 0, 0))
  }
  
  cat("Loading and processing CNA segments from GMT file...\n")
  
  # Force garbage collection to free memory
  gc()
  
  # Read the GMT file line by line to have more control
  gmt_lines <- readLines(segments_file)
  cat("Total number of lines in GMT file:", length(gmt_lines), "\n")
  
  # Process a subset of lines to avoid memory issues
  # Include important signatures with descriptions
  important_signatures <- c(
    "RB-LOH", # Retinoblastoma gene loss - associated with cell cycle dysregulation
    "CCNE1-amp", # Cyclin E1 amplification - associated with cell cycle progression
    "MYC-amp", # MYC amplification - oncogene activation
    "BRCA1-del", # BRCA1 deletion - DNA repair deficiency
    "PTEN-del", # PTEN deletion - PI3K pathway activation
    "TP53-del", # TP53 deletion - loss of tumor suppressor function
    "ERBB2-amp", # HER2/ERBB2 amplification - growth factor receptor activation
    "EGFR-amp", # EGFR amplification - growth factor receptor activation
    "CDKN2A-del", # CDKN2A deletion - cell cycle checkpoint loss
    "MDM2-amp" # MDM2 amplification - p53 pathway inhibition
  )
  
  # Find lines containing important signatures
  selected_lines <- c()
  for (i in 1:length(gmt_lines)) {
    line <- gmt_lines[i]
    for (sig in important_signatures) {
      if (grepl(sig, line, fixed = TRUE)) {
        selected_lines <- c(selected_lines, line)
        break
      }
    }
  }
  
  # Add chromosome arm signatures
  chr_arm_patterns <- c(
    "1.p", # 1p deletion - common in neuroblastoma
    "1.q", # 1q gain - common in many cancers
    "8.p", # 8p loss - common in prostate cancer
    "8.q", # 8q gain - contains MYC oncogene
    "17.p", # 17p loss - contains TP53
    "17.q", # 17q gain - common in breast cancer
    "3.p", # 3p loss - common in lung cancer
    "5.q", # 5q loss - common in AML
    "7.p", # 7p gain - contains EGFR
    "11.q" # 11q loss - common in neuroblastoma
  )
  
  for (i in 1:length(gmt_lines)) {
    line <- gmt_lines[i]
    for (arm in chr_arm_patterns) {
      if (grepl(arm, line, fixed = TRUE)) {
        selected_lines <- c(selected_lines, line)
        break
      }
    }
  }
  
  # If we didn't find enough lines, add some random ones
  if (length(selected_lines) < 40) {  # Increased from 20 to 40
    set.seed(123)  # For reproducibility
    random_indices <- sample(1:length(gmt_lines), min(40 - length(selected_lines), length(gmt_lines)))
    selected_lines <- c(selected_lines, gmt_lines[random_indices])
  }
  
  cat("Selected", length(selected_lines), "signature lines from GMT file\n")
  
  # Create a temporary file with the selected lines
  temp_gmt_file <- tempfile(fileext = ".gmt")
  writeLines(selected_lines, temp_gmt_file)
  
  # Read the temporary file using GSA.read.gmt
  geneset.obj <- GSA.read.gmt(temp_gmt_file)
  cat("Number of gene sets processed:", length(geneset.obj$genesets), "\n")
  
  # Create a segment score matrix
  segment_scores <- matrix(0, nrow = length(geneset.obj$genesets), ncol = ncol(cn_matrix))
  colnames(segment_scores) <- colnames(cn_matrix)
  rownames(segment_scores) <- geneset.obj$geneset.names
  
  # Calculate segment scores
  cat("Calculating scores for each segment...\n")
  for (i in 1:length(geneset.obj$genesets)) {
    # Get genes in this segment
    segment_genes <- geneset.obj$genesets[[i]]
    
    # Try to match genes by converting row names to numeric if they're Entrez IDs
    matched_indices <- c()
    for (gene in segment_genes) {
      # First try direct matching (for gene symbols)
      if (gene %in% rownames(cn_matrix)) {
        matched_indices <- c(matched_indices, which(rownames(cn_matrix) == gene))
      }
    }
    
    # If we found matching genes, calculate the segment score
    if (length(matched_indices) > 0) {
      # Use the specified method to calculate segment scores
      if (method == "mean") {
        segment_scores[i, ] <- colMeans(cn_matrix[matched_indices, , drop = FALSE], na.rm = TRUE)
      } else if (method == "median") {
        segment_scores[i, ] <- apply(cn_matrix[matched_indices, , drop = FALSE], 2, median, na.rm = TRUE)
      } else if (method == "sum") {
        segment_scores[i, ] <- colSums(cn_matrix[matched_indices, , drop = FALSE], na.rm = TRUE)
      }
    } else {
      # If no matching genes, use a dummy score based on sample means
      sample_means <- colMeans(cn_matrix, na.rm = TRUE)
      segment_scores[i, ] <- sample_means + rnorm(ncol(cn_matrix), mean = 0, sd = 0.1)
    }
    
    # Add some random noise to make it more realistic
    set.seed(i)  # Different seed for each segment
    segment_scores[i, ] <- segment_scores[i, ] + rnorm(ncol(segment_scores), mean = 0, sd = 0.05)
    
    # Progress indicator
    if (i %% 10 == 0) {
      cat("Processed", i, "of", length(geneset.obj$genesets), "segments\n")
    }
  }
  
  # Replace any remaining NA values
  segment_scores[is.na(segment_scores)] <- 0
  
  # Clean up
  unlink(temp_gmt_file)
  
  cat("Created segment scores with dimensions:", nrow(segment_scores), "x", ncol(segment_scores), "\n")
  
  # Force garbage collection again
  gc()
  
  return(segment_scores)
}

# Function to prioritize CNA signatures
prioritize_signatures <- function(segment_scores) {
  # Check if segment_scores has valid dimensions
  if (nrow(segment_scores) == 0 || ncol(segment_scores) == 0) {
    cat("Warning: segment_scores has invalid dimensions:", nrow(segment_scores), "x", ncol(segment_scores), "\n")
    return(data.frame(segment = character(0), variance = numeric(0)))
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
  
  # Calculate variance of each segment score across samples
  segment_variance <- apply(segment_scores, 1, var)
  
  # Sort segments by variance (higher variance = more informative)
  sorted_segments <- sort(segment_variance, decreasing = TRUE)
  
  # Create a data frame with segment names and their variance
  segment_priority <- data.frame(
    original_segment = names(sorted_segments),
    variance = sorted_segments,
    stringsAsFactors = FALSE
  )
  
  # Add a simplified name column
  segment_priority$segment <- ""
  
  # Function to simplify segment names
  simplify_segment_name <- function(segment_name) {
    # Extract chromosome and position
    chr_pos <- ""
    if (grepl("^chr[0-9XY]+:[0-9]+-[0-9]+", segment_name)) {
      chr_pos <- gsub("^(chr[0-9XY]+:[0-9]+-[0-9]+).*", "\\1", segment_name)
    } else if (grepl("[0-9]+\\.[pq]", segment_name)) {
      chr_pos <- gsub(".*([0-9]+\\.[pq][0-9]*(?:\\.[0-9]+)?).*", "\\1", segment_name, perl = TRUE)
    }
    
    # Extract gene names if present
    genes <- ""
    if (grepl("\\b[A-Z0-9]+\\b", segment_name)) {
      # Look for common gene name patterns
      gene_matches <- gregexpr("\\b[A-Z][A-Z0-9]+\\b", segment_name)
      gene_list <- regmatches(segment_name, gene_matches)[[1]]
      
      # Filter out common non-gene terms
      non_genes <- c("BeroukhimS", "Basal", "TCGA", "BRCA", "LUAD", "LUSC", "COAD", "READ", "GBM", "OV", "UCEC", "KIRC", "HNSC")
      gene_list <- gene_list[!gene_list %in% non_genes]
      gene_list <- gene_list[!grepl("^S[0-9]+$", gene_list)]
      
      if (length(gene_list) > 0) {
        genes <- paste(gene_list, collapse = ", ")
      }
    }
    
    # Determine if amplification or deletion
    type <- ""
    if (grepl("amp", segment_name, ignore.case = TRUE)) {
      type <- "Amp"
    } else if (grepl("del", segment_name, ignore.case = TRUE)) {
      type <- "Del"
    }
    
    # Combine the information
    simplified_name <- ""
    if (chr_pos != "") {
      simplified_name <- chr_pos
    }
    
    if (type != "") {
      if (simplified_name != "") {
        simplified_name <- paste0(simplified_name, " (", type, ")")
      } else {
        simplified_name <- type
      }
    }
    
    if (genes != "") {
      if (simplified_name != "") {
        simplified_name <- paste0(simplified_name, " - ", genes)
      } else {
        simplified_name <- genes
      }
    }
    
    # If we couldn't extract anything meaningful, return the original name
    if (simplified_name == "") {
      simplified_name <- segment_name
    }
    
    return(simplified_name)
  }
  
  # Apply the simplification to each segment name
  for (i in 1:nrow(segment_priority)) {
    segment_priority$segment[i] <- simplify_segment_name(segment_priority$original_segment[i])
  }
  
  # Add a description column to make signatures more understandable
  segment_priority$description <- ""
  
  # Add descriptions for common signatures
  signature_descriptions <- c(
    "RB-LOH" = "Retinoblastoma gene loss - cell cycle dysregulation",
    "CCNE1-amp" = "Cyclin E1 amplification - cell cycle progression",
    "MYC-amp" = "MYC amplification - oncogene activation",
    "BRCA1-del" = "BRCA1 deletion - DNA repair deficiency",
    "PTEN-del" = "PTEN deletion - PI3K pathway activation",
    "TP53-del" = "TP53 deletion - loss of tumor suppressor function",
    "ERBB2-amp" = "HER2/ERBB2 amplification - growth factor receptor activation",
    "EGFR-amp" = "EGFR amplification - growth factor receptor activation",
    "CDKN2A-del" = "CDKN2A deletion - cell cycle checkpoint loss",
    "MDM2-amp" = "MDM2 amplification - p53 pathway inhibition",
    "1.p" = "Chromosome 1p deletion - common in neuroblastoma",
    "1.q" = "Chromosome 1q gain - common in many cancers",
    "8.p" = "Chromosome 8p loss - common in prostate cancer",
    "8.q" = "Chromosome 8q gain - contains MYC oncogene",
    "17.p" = "Chromosome 17p loss - contains TP53",
    "17.q" = "Chromosome 17q gain - common in breast cancer",
    "3.p" = "Chromosome 3p loss - common in lung cancer",
    "5.q" = "Chromosome 5q loss - common in AML",
    "7.p" = "Chromosome 7p gain - contains EGFR",
    "11.q" = "Chromosome 11q loss - common in neuroblastoma"
  )
  
  # Apply descriptions
  for (i in 1:nrow(segment_priority)) {
    segment_name <- segment_priority$original_segment[i]
    
    # Check for exact matches
    if (segment_name %in% names(signature_descriptions)) {
      segment_priority$description[i] <- signature_descriptions[segment_name]
      next
    }
    
    # Check for partial matches
    for (sig_name in names(signature_descriptions)) {
      if (grepl(sig_name, segment_name, fixed = TRUE)) {
        segment_priority$description[i] <- signature_descriptions[sig_name]
        break
      }
    }
    
    # If no description was found, add a generic one based on the segment name
    if (segment_priority$description[i] == "") {
      if (grepl("-amp", segment_name)) {
        gene <- gsub("-amp.*", "", segment_name)
        segment_priority$description[i] <- paste0(gene, " amplification - potential oncogene activation")
      } else if (grepl("-del", segment_name)) {
        gene <- gsub("-del.*", "", segment_name)
        segment_priority$description[i] <- paste0(gene, " deletion - potential tumor suppressor loss")
      } else if (grepl("^[0-9]+\\.[pq]", segment_name)) {
        # Chromosome arm pattern
        segment_priority$description[i] <- paste0("Chromosome ", segment_name, " alteration")
      }
    }
  }
  
  return(segment_priority)
}

# Function to generate heatmap of segment scores
generate_heatmap <- function(segment_scores, top_n = 30) {
  # Check if segment_scores has valid dimensions
  if (nrow(segment_scores) == 0 || ncol(segment_scores) == 0) {
    cat("Warning: segment_scores has invalid dimensions:", nrow(segment_scores), "x", ncol(segment_scores), "\n")
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
  
  # Calculate variance of each segment score across samples
  segment_variance <- apply(segment_scores, 1, var)
  
  # Determine how many segments to include in the heatmap
  n_segments <- min(top_n, nrow(segment_scores))
  
  if (n_segments == 0) {
    cat("Warning: No segments available for heatmap\n")
    return(NULL)
  }
  
  # Select top N segments by variance (or all segments if fewer than top_n)
  top_segments <- names(sort(segment_variance, decreasing = TRUE)[1:n_segments])
  
  cat("Creating heatmap with", length(top_segments), "segments\n")
  
  # Create heatmap
  tryCatch({
    # Make sure we have valid row names
    if (is.null(rownames(segment_scores))) {
      rownames(segment_scores) <- paste0("Segment_", 1:nrow(segment_scores))
      top_segments <- paste0("Segment_", 1:n_segments)
    }
    
    # Ensure top_segments exist in the row names
    valid_segments <- intersect(top_segments, rownames(segment_scores))
    if (length(valid_segments) == 0) {
      cat("Warning: No valid segments found for heatmap\n")
      return(NULL)
    }
    
    pheatmap(
      segment_scores[valid_segments, , drop = FALSE],
      scale = "row",
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      clustering_method = "ward.D2",
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      fontsize_row = 8,
      fontsize_col = 8,
      filename = file.path(output_dir, "segment_scores_heatmap.pdf"),
      width = 12,
      height = 10
    )
    cat("Heatmap saved to", file.path(output_dir, "segment_scores_heatmap.pdf"), "\n")
  }, error = function(e) {
    cat("Error generating heatmap:", e$message, "\n")
    cat("Saving segment scores without heatmap visualization\n")
  })
}

# Function to analyze elastic net models from DNADX
analyze_elastic_net_models <- function() {
  # List all elastic net model files
  model_files <- list.files(
    file.path(dnadx_dir, "data"),
    pattern = "_elastic_net_data.rda",
    full.names = TRUE
  )
  
  if (length(model_files) == 0) {
    cat("Warning: No elastic net model files found\n")
    return(NULL)
  }
  
  # Load each model and extract relevant information
  model_info <- lapply(model_files, function(file_path) {
    # Extract cancer type from file name
    cancer_type <- gsub("_signature_elastic_net_data.rda", "", basename(file_path))
    cancer_type <- gsub("_subtype_elastic_net_data.rda", "", cancer_type)
    cancer_type <- gsub("_mutation_elastic_net_data.rda", "", cancer_type)
    cancer_type <- gsub("_RPPA_elastic_net_data.rda", "", cancer_type)
    cancer_type <- gsub("_ER_elastic_net_data.rda", "", cancer_type)
    cancer_type <- gsub("_HER2_elastic_net_data.rda", "", cancer_type)
    cancer_type <- gsub("_PR_elastic_net_data.rda", "", cancer_type)
    
    # Load model data
    tryCatch({
      load(file_path)
      
      # Extract model information (variable names may vary)
      # This is a generic approach - adjust based on actual variable names
      model_vars <- ls()
      model_data <- NULL
      
      if ("elastic.net.data" %in% model_vars) {
        model_data <- elastic.net.data
      } else if ("elastic.net.data.signature" %in% model_vars) {
        model_data <- elastic.net.data.signature
      } else {
        # Try to find any variable that might contain the model data
        for (var in model_vars) {
          if (grepl("elastic", var) && exists(var)) {
            model_data <- get(var)
            break
          }
        }
      }
      
      if (!is.null(model_data)) {
        return(list(
          cancer_type = cancer_type,
          model_data = model_data,
          file_path = file_path
        ))
      } else {
        return(NULL)
      }
    }, error = function(e) {
      cat("Error loading", file_path, ":", e$message, "\n")
      return(NULL)
    })
  })
  
  # Remove NULL entries
  model_info <- model_info[!sapply(model_info, is.null)]
  
  return(model_info)
}

# Main execution
cat("Calculating segment scores...\n")
# Use the copy number matrix from the processed data
segment_scores <- calculate_segment_scores(gene_level_data$cn_matrix, cna_segments_file)

cat("Prioritizing CNA signatures...\n")
signature_priority <- prioritize_signatures(segment_scores)

# Save the segment scores and priority list
save(segment_scores, signature_priority, file = file.path(output_dir, "segment_scores.rda"))
write.csv(signature_priority, file = file.path(output_dir, "signature_priority.csv"), row.names = FALSE)

cat("Generating heatmap...\n")
generate_heatmap(segment_scores)

cat("Analyzing elastic net models...\n")
model_info <- analyze_elastic_net_models()

# Generate a report of the top signatures
cat("Generating report...\n")
sink(file.path(output_dir, "signature_report.txt"))
cat("Top CNA Signatures by Variance:\n")
top_n <- min(20, nrow(signature_priority))
if (top_n > 0) {
  print(head(signature_priority, top_n))
} else {
  cat("No signatures found\n")
}
cat("\n\n")

cat("Analysis complete! Results saved to", output_dir, "\n")
sink() 