#!/usr/bin/env Rscript

# Script to group samples by patient and create new plots
# This script reads the existing analysis results and creates new visualizations
# with samples grouped by patient

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
library(dendextend)
library(readxl)  # Add readxl library for Excel file reading

# Set paths
input_dir <- "./output"
output_dir <- "./output_v2"
clinical_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Clinical_information/TFx_hg38_1000_500kb.csv"
molecular_type_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Clinical_information/Molecular_type.xlsx"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check if required files exist
segment_scores_file <- file.path(input_dir, "segment_scores.rda")

if (!file.exists(segment_scores_file)) {
  cat("Error: Segment scores file not found:", segment_scores_file, "\n")
  cat("Please run the analyze_signatures.R script first.\n")
  quit(status = 1)
}

# Load the segment scores and priority list
cat("Loading analysis results...\n")
load(segment_scores_file)

# Read clinical information
cat("Reading clinical information...\n")
clinical_data <- read.csv(clinical_file, skip = 1, stringsAsFactors = FALSE)

# Read molecular type information
cat("Reading molecular type information...\n")
molecular_data <- read_excel(molecular_type_file)
molecular_data <- as.data.frame(molecular_data)
colnames(molecular_data)[1] <- "Patient"  # Rename first column to Patient

# Convert molecular type columns to factors
molecular_data$HR_HER2 <- factor(apply(molecular_data[, c("HR+HER2-", "HR+HER2+", "TNBC")], 1, 
                                      function(x) {
                                        if(x[1] == 1) return("HR+HER2-")
                                        if(x[2] == 1) return("HR+HER2+")
                                        if(x[3] == 1) return("TNBC")
                                        return(NA)
                                      }))

# Create a mapping of patient to molecular type
molecular_type_map <- setNames(molecular_data$HR_HER2, molecular_data$Patient)

# Print molecular type mapping for debugging
cat("\nMolecular type mapping:\n")
print(molecular_type_map)

# Function to extract sample IDs from clinical data
extract_sample_ids <- function(clinical_data) {
  # Initialize empty list to store patient-sample mappings
  patient_samples <- list()
  
  # Process each row (patient)
  for (i in 1:nrow(clinical_data)) {
    patient_id <- clinical_data$Patient_Code[i]
    samples <- c()
    
    # Check each progression column
    progression_cols <- c(2, 5, 8, 11, 14, 17) # Column indices for progression sample IDs
    
    for (col in progression_cols) {
      if (col <= ncol(clinical_data) && !is.na(clinical_data[i, col]) && clinical_data[i, col] != "") {
        # Extract sample ID (remove any text in parentheses)
        sample_id <- gsub(" \\(.*\\)$", "", clinical_data[i, col])
        samples <- c(samples, sample_id)
      }
    }
    
    if (length(samples) > 0) {
      patient_samples[[patient_id]] <- samples
    }
  }
  
  return(patient_samples)
}

# Extract patient-sample mappings
patient_samples <- extract_sample_ids(clinical_data)

# Create a data frame with patient-sample mappings
sample_patient_df <- data.frame(
  Sample = unlist(patient_samples),
  Patient = rep(names(patient_samples), sapply(patient_samples, length)),
  stringsAsFactors = FALSE
)

# Print summary of patient-sample mappings
cat("Found", length(patient_samples), "patients with", nrow(sample_patient_df), "samples\n")
for (patient in names(patient_samples)) {
  cat("Patient", patient, "has", length(patient_samples[[patient]]), "samples:", 
      paste(patient_samples[[patient]], collapse=", "), "\n")
}

# Function to reorder samples by patient
reorder_samples_by_patient <- function(segment_scores, patient_samples) {
  # Get current sample names
  current_samples <- colnames(segment_scores)
  
  # Print diagnostic information
  cat("\nDiagnostic information for sample matching:\n")
  cat("Current sample names in segment_scores:\n")
  print(head(current_samples, 10))
  cat("...\n")
  
  # Create a mapping of sample names to standardized format
  sample_mapping <- data.frame(
    original = current_samples,
    standardized = current_samples,
    stringsAsFactors = FALSE
  )
  
  # Try to match samples with different formats
  for (i in 1:nrow(sample_mapping)) {
    orig <- sample_mapping$original[i]
    # Try different transformations to match with clinical data
    # 1. As is
    # 2. Remove leading/trailing whitespace
    # 3. Extract TP pattern
    if (grepl("TP[0-9]+-M[0-9]+", orig)) {
      sample_mapping$standardized[i] <- gsub(".*?(TP[0-9]+-M[0-9]+).*", "\\1", orig)
    } else {
      sample_mapping$standardized[i] <- trimws(orig)
    }
  }
  
  # Print the standardized sample names
  cat("\nStandardized sample names (first 10):\n")
  print(head(sample_mapping, 10))
  cat("...\n")
  
  # Create a mapping for samples with underscores
  underscore_mapping <- list()
  
  # Map TP18-M643_1 to TP18-M643
  underscore_mapping[["TP18-M643_1"]] <- "TP18-M643"
  
  # For ULP_3 and ULP_4, both have samples with TP18-M681 as base
  # We need to handle this special case
  # First, check if TP18-M681 appears twice in the data
  m681_indices <- which(sample_mapping$standardized == "TP18-M681")
  
  if (length(m681_indices) == 2) {
    cat("\nFound two instances of TP18-M681, assigning them to ULP_3 and ULP_4\n")
    # Assign the first one to ULP_4 (TP18-M681_1)
    underscore_mapping[["TP18-M681_1"]] <- sample_mapping$original[m681_indices[1]]
    # Assign the second one to ULP_3 (TP18-M681_2)
    underscore_mapping[["TP18-M681_2"]] <- sample_mapping$original[m681_indices[2]]
  } else if (length(m681_indices) == 1) {
    cat("\nFound only one instance of TP18-M681, cannot uniquely assign to ULP_3 or ULP_4\n")
    # If there's only one, we'll assign it to ULP_4 (arbitrary choice)
    underscore_mapping[["TP18-M681_1"]] <- sample_mapping$original[m681_indices[1]]
  }
  
  # Create a new order based on patient grouping
  new_order <- c()
  patient_annotations <- c()
  
  # Sort patient IDs numerically
  patient_ids <- names(patient_samples)
  patient_numbers <- as.numeric(gsub("ULP_", "", patient_ids))
  sorted_indices <- order(patient_numbers)
  sorted_patient_ids <- patient_ids[sorted_indices]
  
  cat("\nProcessing patients in numerical order:", paste(sorted_patient_ids, collapse=", "), "\n")
  
  # Process each patient in numerical order
  for (patient in sorted_patient_ids) {
    patient_sample_ids <- patient_samples[[patient]]
    
    # Find which samples from this patient are in our data
    matched_indices <- c()
    matched_samples <- c()
    for (sample_id in patient_sample_ids) {
      # First check if this is a special case with underscore
      if (sample_id %in% names(underscore_mapping)) {
        mapped_id <- underscore_mapping[[sample_id]]
        matches <- which(sample_mapping$original == mapped_id)
        if (length(matches) > 0) {
          cat("  Using special mapping for", sample_id, "->", mapped_id, "\n")
        }
      } else {
        # Try to find matches in our standardized sample names
        matches <- which(sample_mapping$standardized == sample_id)
      }
      
      if (length(matches) > 0) {
        matched_indices <- c(matched_indices, matches)
        matched_samples <- c(matched_samples, sample_id)
      }
    }
    
    if (length(matched_indices) > 0) {
      # Sort the matched samples chronologically
      # Extract TP year and M number for sorting
      sample_info <- data.frame(
        index = matched_indices,
        sample_id = matched_samples,
        stringsAsFactors = FALSE
      )
      
      # Extract year and M number for chronological sorting
      sample_info$year <- as.numeric(gsub("TP([0-9]+)-M.*", "\\1", sample_info$sample_id))
      sample_info$m_number <- as.numeric(gsub("TP[0-9]+-M([0-9]+).*", "\\1", sample_info$sample_id))
      
      # Sort by year, then by M number
      sample_info <- sample_info[order(sample_info$year, sample_info$m_number), ]
      
      # Add these samples to our new order in chronological order
      patient_samples_in_data <- sample_mapping$original[sample_info$index]
      new_order <- c(new_order, patient_samples_in_data)
      
      # Add patient annotation
      patient_annotations <- c(patient_annotations, 
                              rep(patient, length(patient_samples_in_data)))
    } else {
      cat("No matches found for patient", patient, "samples:", paste(patient_sample_ids, collapse=", "), "\n")
    }
  }
  
  # Identify unmatched samples but don't include them
  unmatched_samples <- setdiff(current_samples, new_order)
  if (length(unmatched_samples) > 0) {
    cat("\nExcluding unmatched samples (", length(unmatched_samples), "):\n")
    print(unmatched_samples)
  }
  
  # Return the new order and annotations
  return(list(
    new_order = new_order,
    patient_annotations = patient_annotations
  ))
}

# Reorder samples by patient
reordering <- reorder_samples_by_patient(segment_scores, patient_samples)
new_sample_order <- reordering$new_order
patient_annotations <- reordering$patient_annotations

# Create a reordered segment scores matrix
segment_scores_reordered <- segment_scores[, new_sample_order, drop = FALSE]

# Function to create a heatmap with patient annotations
create_patient_annotated_heatmap <- function(segment_scores, patient_annotations, top_n = 30) {
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
  
  # Create annotation data frame for patients and molecular types
  annotation_col <- data.frame(
    Patient = factor(patient_annotations),
    Molecular_Type = factor(sapply(patient_annotations, function(x) molecular_type_map[x]))
  )
  rownames(annotation_col) <- colnames(segment_scores)
  
  # Print annotation data frame for debugging
  cat("\nAnnotation data frame:\n")
  print(head(annotation_col))
  
  # Define colors for patient annotation
  n_patients <- length(unique(patient_annotations))
  patient_colors <- colorRampPalette(brewer.pal(min(9, n_patients), "Set1"))(n_patients)
  names(patient_colors) <- unique(patient_annotations)
  
  # Define colors for molecular type annotation
  molecular_colors <- c(
    "HR+HER2-" = "#FF8C00",  # Dark Orange
    "HR+HER2+" = "#9932CC",  # Dark Orchid (Purple)
    "TNBC" = "#008B8B"       # Dark Cyan (Teal)
  )
  
  # Print molecular colors for debugging
  cat("\nMolecular colors:\n")
  print(molecular_colors)
  
  annotation_colors <- list(
    Patient = patient_colors,
    Molecular_Type = molecular_colors
  )
  
  # Print annotation colors for debugging
  cat("\nAnnotation colors:\n")
  print(annotation_colors)
  
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
    
    # Create a clustering order for rows only (segments)
    # Use a different approach to calculate distance
    segment_data <- segment_scores[valid_segments, , drop = FALSE]
    
    # Calculate dimensions for square cells
    # Adjust width based on number of samples and segments
    n_samples <- ncol(segment_data)
    n_valid_segments <- nrow(segment_data)
    
    # Calculate width to ensure square cells and enough space for legend
    # Base width on number of samples with some extra space for legend
    width <- max(18, n_samples * 0.25 + 6)  # Increased base width to 18 inches
    # Height based on number of segments to maintain square cells
    height <- max(12, n_valid_segments * 0.25 + 2)
    
    # Create the heatmap with custom clustering and square cells
    pheatmap(
      segment_data,
      scale = "row",
      cluster_rows = TRUE,  # Cluster rows (segments)
      cluster_cols = FALSE, # Don't cluster columns (samples) to keep patient groups together
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      fontsize_row = 8,
      fontsize_col = 8,
      cellwidth = 12,  # Set cell width for square cells
      cellheight = 12, # Set cell height for square cells
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      filename = file.path(output_dir, "segment_scores_heatmap_by_patient.pdf"),
      width = width,
      height = height,
      main = "Segment Scores Heatmap (Grouped by Patient)\nZ-scores represent standardized segment scores across samples",
      show_rownames = TRUE,
      show_colnames = TRUE,
      border_color = NA,
      display_numbers = FALSE,
      number_format = "%.2f",
      number_color = "black",
      fontsize = 8,
      fontsize_number = 8,
      silent = FALSE
    )
    
    # Create a second heatmap with a gap between patient groups
    # First, create a list to store the patient groups
    patient_groups <- list()
    unique_patients <- unique(patient_annotations)
    
    # Get the indices for each patient group
    current_index <- 1
    gap_positions <- c()
    
    for (patient in unique_patients) {
      patient_indices <- which(patient_annotations == patient)
      if (length(patient_indices) > 0) {
        patient_groups[[patient]] <- patient_indices
        current_index <- current_index + length(patient_indices)
        if (current_index <= length(patient_annotations)) {
          gap_positions <- c(gap_positions, current_index - 0.5)
        }
      }
    }
    
    # Create a heatmap with gaps between patient groups
    pdf(file.path(output_dir, "segment_scores_heatmap_by_patient_with_gaps.pdf"), width = width, height = height)
    
    # Use the pheatmap function with gaps
    pheatmap(
      segment_data,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      fontsize_row = 8,
      fontsize_col = 8,
      cellwidth = 12,  # Set cell width for square cells
      cellheight = 12, # Set cell height for square cells
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      gaps_col = gap_positions,
      main = "Segment Scores Heatmap (Grouped by Patient with Gaps)\nZ-scores represent standardized segment scores across samples",
      show_rownames = TRUE,
      show_colnames = TRUE,
      border_color = NA,
      display_numbers = FALSE,
      number_format = "%.2f",
      number_color = "black",
      fontsize = 8,
      fontsize_number = 8,
      silent = FALSE
    )
    
    dev.off()
    
    # Create a third heatmap with patient labels
    pdf(file.path(output_dir, "segment_scores_heatmap_by_patient_labeled.pdf"), width = width, height = height + 2)
    
    # Plot the heatmap directly
    pheatmap(
      segment_data,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      fontsize_row = 8,
      fontsize_col = 8,
      cellwidth = 12,  # Set cell width for square cells
      cellheight = 12, # Set cell height for square cells
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      gaps_col = gap_positions,
      main = "Segment Scores Heatmap (Grouped by Patient)\nZ-scores represent standardized segment scores across samples",
      show_rownames = TRUE,
      show_colnames = TRUE,
      border_color = NA,
      display_numbers = FALSE,
      number_format = "%.2f",
      number_color = "black",
      fontsize = 8,
      fontsize_number = 8,
      silent = FALSE
    )
    
    dev.off()
    
    # Create a fourth heatmap with patient labels and molecular type below
    pdf(file.path(output_dir, "segment_scores_heatmap_with_patient_labels.pdf"), width = width, height = height + 6)
    
    # Create a layout with space for the heatmap and patient labels
    layout(matrix(c(1, 2, 3), nrow = 3), heights = c(5, 1, 1))
    
    # Plot the heatmap in the top panel
    par(mar = c(0, 4, 4, 2))
    image(t(scale(t(as.matrix(segment_data)))), 
          col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
          breaks = seq(-3, 3, length.out = 101),
          axes = FALSE,
          main = "Segment Scores Heatmap (Grouped by Patient)\nZ-scores represent standardized segment scores across samples")
    
    # Add grid lines with more visible appearance
    abline(v = seq(0, 1, length.out = ncol(segment_data) + 1), 
           h = seq(0, 1, length.out = nrow(segment_data) + 1),
           col = "gray50", lty = "solid", lwd = 0.5)
    
    # Add a color legend for the heatmap
    legend_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(5)
    legend("right", legend = c("-3", "-1.5", "0", "1.5", "3"), 
           fill = legend_colors, title = "Z-score\n(Row-scaled)\nBlue = Low\nRed = High", 
           cex = 0.8, box.lwd = 1)
    
    # Add patient labels in the middle panel
    par(mar = c(2, 4, 0, 2))
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", 
         xlim = c(0, 1), ylim = c(0, 1))
    
    # Add grid lines for patient labels with more visible appearance
    abline(v = seq(0, 1, length.out = length(patient_groups) + 1),
           col = "gray50", lty = "solid", lwd = 0.5)
    
    # Add patient labels
    patient_positions <- c()
    patient_labels <- c()
    
    for (patient in names(patient_groups)) {
      indices <- patient_groups[[patient]]
      if (length(indices) > 0) {
        center_pos <- mean(indices) / ncol(segment_scores)
        patient_positions <- c(patient_positions, center_pos)
        patient_labels <- c(patient_labels, patient)
      }
    }
    
    text(patient_positions, rep(0.5, length(patient_positions)), 
         patient_labels, cex = 1.2, font = 2)
    
    # Add molecular type labels in the bottom panel
    par(mar = c(2, 4, 0, 2))
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", 
         xlim = c(0, 1), ylim = c(0, 1))
    
    # Add grid lines for molecular type labels with more visible appearance
    abline(v = seq(0, 1, length.out = length(patient_groups) + 1),
           col = "gray50", lty = "solid", lwd = 0.5)
    
    # Add molecular type labels with colors
    molecular_positions <- c()
    molecular_labels <- c()
    molecular_colors_vec <- c()
    
    for (patient in names(patient_groups)) {
      indices <- patient_groups[[patient]]
      if (length(indices) > 0) {
        center_pos <- mean(indices) / ncol(segment_scores)
        molecular_positions <- c(molecular_positions, center_pos)
        mol_type <- molecular_type_map[patient]
        molecular_labels <- c(molecular_labels, as.character(mol_type))
        molecular_colors_vec <- c(molecular_colors_vec, molecular_colors[as.character(mol_type)])
      }
    }
    
    # Draw colored rectangles for molecular types
    rect_height <- 0.4
    rect_y_bottom <- 0.5 - rect_height/2
    rect_y_top <- 0.5 + rect_height/2
    
    for (i in seq_along(molecular_positions)) {
      if (!is.na(molecular_labels[i])) {
        rect_width <- 0.04
        rect_x_left <- molecular_positions[i] - rect_width/2
        rect_x_right <- molecular_positions[i] + rect_width/2
        rect(rect_x_left, rect_y_bottom, rect_x_right, rect_y_top,
             col = molecular_colors_vec[i], border = "black")
      }
    }
    
    # Add legend for molecular types
    legend("right", 
           legend = names(molecular_colors),
           fill = molecular_colors,
           title = "Molecular Type",
           cex = 0.8,
           box.lwd = 1)
    
    dev.off()
    
    cat("Heatmap saved to", file.path(output_dir, "segment_scores_heatmap_by_patient.pdf"), "\n")
    cat("Heatmap with gaps saved to", file.path(output_dir, "segment_scores_heatmap_by_patient_with_gaps.pdf"), "\n")
    cat("Heatmap with patient labels saved to", file.path(output_dir, "segment_scores_heatmap_by_patient_labeled.pdf"), "\n")
    cat("Heatmap with patient labels and molecular type saved to", file.path(output_dir, "segment_scores_heatmap_with_patient_labels.pdf"), "\n")
  }, error = function(e) {
    cat("Error generating heatmap:", e$message, "\n")
    cat("Saving segment scores without heatmap visualization\n")
  })
}

# Function to create a correlation heatmap with patient annotations
create_patient_annotated_correlation_heatmap <- function(segment_scores, patient_annotations, top_n = 40) {
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
  
  # Create annotation data frame for patients
  annotation_col <- data.frame(
    Patient = factor(patient_annotations)
  )
  rownames(annotation_col) <- colnames(segment_scores)
  
  # Define colors for patient annotation
  n_patients <- length(unique(patient_annotations))
  patient_colors <- colorRampPalette(brewer.pal(min(9, n_patients), "Set1"))(n_patients)
  names(patient_colors) <- unique(patient_annotations)
  
  annotation_colors <- list(
    Patient = patient_colors
  )
  
  # Calculate correlation matrix
  tryCatch({
    cor_matrix <- cor(t(segment_scores[valid_segments, , drop = FALSE]), method = "pearson")
    
    # Calculate dimensions for square cells
    n_valid_segments <- length(valid_segments)
    width <- max(18, n_valid_segments * 0.25 + 6)  # Increased base width
    height <- width  # Make it square
    
    # Create heatmap
    pheatmap(
      cor_matrix,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
      fontsize_row = 8,
      fontsize_col = 8,
      cellwidth = 12,  # Set cell width for square cells
      cellheight = 12, # Set cell height for square cells
      filename = file.path(output_dir, "signature_correlation_heatmap_by_patient.pdf"),
      width = width,
      height = height,
      main = "Signature Correlation Heatmap\nPearson correlation between genomic segments",
      legend_labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
      legend_breaks = seq(-1, 1, 0.5),
      legend = TRUE,
      legend_title = "Pearson Correlation\nBlue = Negative\nYellow = Neutral\nRed = Positive"
    )
    cat("Correlation heatmap saved to", file.path(output_dir, "signature_correlation_heatmap_by_patient.pdf"), "\n")
  }, error = function(e) {
    cat("Error creating correlation heatmap:", e$message, "\n")
  })
}

# Function to create a PCA plot with patient annotations
create_patient_annotated_pca_plot <- function(segment_scores, patient_annotations) {
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
    
    # Add patient information
    pca_data$Patient <- factor(patient_annotations)
    pca_data$Sample <- colnames(segment_scores)
    
    # Calculate variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
    
    # Define colors for patients
    n_patients <- length(unique(patient_annotations))
    patient_colors <- colorRampPalette(brewer.pal(min(9, n_patients), "Set1"))(n_patients)
    names(patient_colors) <- unique(patient_annotations)
    
    # Create PCA plot
    p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Patient, label = Sample)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_text(size = 2, vjust = -1, show.legend = FALSE) +
      scale_color_manual(values = patient_colors) +
      labs(
        title = "PCA of Samples Based on CNA Signatures (Grouped by Patient)",
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
        caption = "Values represent principal component scores derived from segment scores"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.key.size = unit(0.8, "cm"),  # Increase legend key size
        legend.text = element_text(size = 10),  # Increase legend text size
        plot.caption = element_text(size = 10, hjust = 0)  # Add caption explanation
      )
    
    # Save the plot with increased width for legend
    ggsave(file.path(output_dir, "pca_plot_by_patient.pdf"), p, width = 16, height = 12)
    cat("PCA plot saved to", file.path(output_dir, "pca_plot_by_patient.pdf"), "\n")
    
    return(p)
  }, error = function(e) {
    cat("Error creating PCA plot:", e$message, "\n")
    return(NULL)
  })
}

# Function to create a sample clustering dendrogram with patient annotations
create_patient_annotated_dendrogram <- function(segment_scores, patient_annotations) {
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
    
    # Convert to dendrogram object
    dend <- as.dendrogram(hc)
    
    # Create a color mapping for patients
    n_patients <- length(unique(patient_annotations))
    patient_colors <- colorRampPalette(brewer.pal(min(9, n_patients), "Set1"))(n_patients)
    names(patient_colors) <- unique(patient_annotations)
    
    # Create a color vector for the labels
    label_colors <- patient_colors[patient_annotations[match(labels(dend), colnames(segment_scores))]]
    
    # Color the labels
    dend <- set(dend, "labels_col", label_colors)
    
    # Create dendrogram plot with increased width
    pdf(file.path(output_dir, "sample_dendrogram_by_patient.pdf"), width = 18, height = 12)
    
    # Plot the dendrogram
    plot(dend, main = "Sample Clustering Based on CNA Signatures (Colored by Patient)", 
         xlab = "", sub = "Euclidean distance based on segment scores")
    
    # Add a legend with increased size
    legend("topright", legend = names(patient_colors), 
           fill = patient_colors, title = "Patient", cex = 1.0, 
           box.lwd = 2, box.col = "black")
    
    dev.off()
    
    cat("Sample dendrogram saved to", file.path(output_dir, "sample_dendrogram_by_patient.pdf"), "\n")
  }, error = function(e) {
    cat("Error creating sample dendrogram:", e$message, "\n")
  })
}

# Main execution
cat("Creating patient-annotated visualizations...\n")

# Create heatmap with patient annotations
create_patient_annotated_heatmap(segment_scores_reordered, patient_annotations)

# Create correlation heatmap
create_patient_annotated_correlation_heatmap(segment_scores_reordered, patient_annotations)

# Create PCA plot with patient annotations
create_patient_annotated_pca_plot(segment_scores_reordered, patient_annotations)

# Create sample dendrogram with patient annotations
create_patient_annotated_dendrogram(segment_scores_reordered, patient_annotations)

# Copy the signature priority file to the new output directory
file.copy(file.path(input_dir, "signature_priority.csv"), 
          file.path(output_dir, "signature_priority.csv"))

cat("Patient-grouped analysis complete! Results saved to", output_dir, "\n") 