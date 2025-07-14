#!/usr/bin/env Rscript

# Survival Analysis for CNA Signatures
# This script performs survival analysis on all signatures from DNADX analysis
# - Treats samples individually for PFS
# - Treats patients individually for OS
# - Calculates hazard ratios
# - Analyzes all 40 signatures but only plots those with p < 0.05
# - Uses distinct colors for high/low groups with clear annotations

# Load required libraries
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(corrplot)  # For correlation plots

# Set paths
output_dir <- "../output"
survival_dir <- "."
sample_info_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Clinical_information/TFx_hg38_1000_500kb.csv"
patient_survival_file <- "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Clinical_information/Filtered_Survival_Data.csv"

# Create output directory if it doesn't exist
dir.create(survival_dir, showWarnings = FALSE, recursive = TRUE)

# Load segment scores and signature priority
cat("Loading segment scores and signature priority...\n")
load(file.path(output_dir, "segment_scores.rda"))
signature_priority <- read.csv(file.path(output_dir, "signature_priority.csv"), stringsAsFactors = FALSE)

# Add diagnostic information about segment scores
cat("Dimensions of segment_scores matrix:", dim(segment_scores)[1], "rows x", dim(segment_scores)[2], "columns\n")
cat("Number of samples in segment_scores:", ncol(segment_scores), "\n")

# Extract sample names from segment scores
sample_names <- colnames(segment_scores)
sample_names <- gsub(".annotated", "", sample_names)
cat("Total number of samples after processing names:", length(sample_names), "\n")

# Function to extract gene annotations from signature names
extract_gene_annotation <- function(signature_name) {
  # First check if it's one of the manually annotated signatures
  if (grepl("chr1_168438419_168880930", signature_name)) {
    return("XCL1, XCL2, DPT")
  } else if (grepl("chr14_56307190_56415610", signature_name)) {
    return("No genes")
  } else if (grepl("chr4_15313739_17632477", signature_name)) {
    return("C1QTNF7, FBXL5, CD38, FGFBP1, TAPT1, LDB2")
  } else if (grepl("chr9_21489625_22474701", signature_name)) {
    return("CDKN2A, CDKN2B")
  }
  
  # If not one of the manually annotated signatures, use the original logic
  if (grepl("____", signature_name)) {
    gene_part <- strsplit(signature_name, "____")[[1]][2]
    return(gene_part)
  } else {
    # Extract chromosome and coordinates for potential annotation lookup
    parts <- strsplit(signature_name, "_")[[1]]
    chr <- parts[1]
    # Try to find start and end positions
    coords <- grep("^[0-9]+$", parts, value = TRUE)
    if (length(coords) >= 2) {
      return(paste0(chr, ":", coords[1], "-", coords[2]))
    }
    return("")
  }
}

# Load sample information
cat("Loading sample information...\n")
sample_info <- fread(sample_info_file, skip = 1, header = TRUE)
cat("Number of rows in sample_info:", nrow(sample_info), "\n")

# Load patient survival data
cat("Loading patient survival data...\n")
patient_survival <- fread(patient_survival_file, header = TRUE)
# Remove the first column if it's an index
if (is.numeric(patient_survival[[1]])) {
  patient_survival <- patient_survival[, -1, with = FALSE]
}
cat("Number of rows in patient_survival:", nrow(patient_survival), "\n")

# Print column names to debug
cat("Sample info columns:", paste(colnames(sample_info), collapse = ", "), "\n")
cat("Patient survival columns:", paste(colnames(patient_survival), collapse = ", "), "\n")

# Extract sample names from segment scores
sample_names <- colnames(segment_scores)
sample_names <- gsub(".annotated", "", sample_names)
cat("Total number of samples in segment scores:", length(sample_names), "\n")

# Function to extract sample progression time from sample name
extract_sample_progression_time <- function(sample_name) {
  # Extract the patient code and progression time from the sample name
  for (i in 1:nrow(sample_info)) {
    patient_code <- sample_info[[1]][i]
    
    # Check each progression column
    for (j in seq(2, ncol(sample_info), 3)) {
      if (j <= ncol(sample_info) && !is.na(sample_info[[j]][i]) && grepl(sample_name, sample_info[[j]][i])) {
        # Find the corresponding progression time in patient_survival
        patient_idx <- which(patient_survival$Patient_Code == patient_code)
        if (length(patient_idx) > 0) {
          # Calculate which progression this is (1st, 2nd, etc.)
          prog_num <- (j - 2) / 3 + 1
          
          # Get the corresponding progression time column
          prog_col <- paste0("Time from Dx to ", switch(prog_num, "1st", "2nd", "3rd", "4th", "5th", "6th"), " Progression (months)")
          
          if (prog_col %in% colnames(patient_survival)) {
            prog_time <- patient_survival[[prog_col]][patient_idx]
            return(list(
              patient_code = patient_code,
              progression_time = as.numeric(prog_time),
              os_time = as.numeric(patient_survival$`OS (months)`[patient_idx])
            ))
          }
        }
        break
      }
    }
  }
  return(list(patient_code = NA, progression_time = NA, os_time = NA))
}

# Map sample names to progression times
cat("Mapping sample names to progression times...\n")
sample_progression <- lapply(sample_names, extract_sample_progression_time)
names(sample_progression) <- sample_names

# Count samples with valid progression data
valid_progression <- sum(!sapply(sample_progression, function(x) is.na(x$progression_time)))
cat("Number of samples with valid progression data:", valid_progression, "\n")

# Print samples with missing progression data
missing_progression <- sample_names[sapply(sample_progression, function(x) is.na(x$progression_time))]
if (length(missing_progression) > 0) {
  cat("Samples with missing progression data:", paste(missing_progression, collapse=", "), "\n")
}

# Create a data frame with sample names and progression times (for PFS analysis)
sample_pfs_data <- data.frame(
  Sample = sample_names,
  Patient_Code = sapply(sample_progression, function(x) x$patient_code),
  PFS = sapply(sample_progression, function(x) x$progression_time)
)

# Count samples with valid PFS data
valid_pfs <- sum(!is.na(sample_pfs_data$PFS))
cat("Number of samples with valid PFS data:", valid_pfs, "\n")

# Identify and print samples with missing PFS data
missing_pfs_samples <- sample_pfs_data$Sample[is.na(sample_pfs_data$PFS)]
cat("Number of samples with missing PFS data:", length(missing_pfs_samples), "\n")
cat("Samples with missing PFS data:", paste(missing_pfs_samples, collapse=", "), "\n")

# Create a data frame with patient codes and OS times (for OS analysis)
# For patients with multiple samples, we'll use the patient-level OS data
patient_os_data <- patient_survival %>%
  select(Patient_Code, `OS (months)`) %>%
  rename(OS = `OS (months)`)

# Define color schemes for high/low groups
high_color <- "#FF0000"  # Red for high
low_color <- "#0000FF"   # Blue for low

# Function to create an enhanced KM plot with better annotations
create_km_plot <- function(fit, title, xlab = "Time (months)", ylab = "Survival Probability", 
                          hr_text = "", n_high = 0, n_low = 0, signature_name = "", analysis_type = "PFS") {
  # Extract gene annotation
  gene_annotation <- extract_gene_annotation(signature_name)
  
  # Create subtitle with gene annotation if available
  subtitle <- hr_text
  if (gene_annotation != "") {
    subtitle <- paste0(subtitle, "\nGenes: ", gene_annotation)
  }
  
  # Use the full signature name for the title with analysis type
  plot_title <- paste0(signature_name, " (", analysis_type, ")")
  
  # Extract data from survfit object
  surv_data <- data.frame(
    time = fit$time,
    surv = fit$surv,
    upper = fit$upper,
    lower = fit$lower,
    strata = rep(names(fit$strata), fit$strata)
  )
  
  # Fix strata names to ensure they match exactly with "High" and "Low"
  surv_data$strata <- gsub("Group=", "", surv_data$strata)
  
  # Create plot with explicit color mapping
  p <- ggplot(surv_data, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_step(size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    labs(
      title = plot_title,  # Use the full signature name with analysis type
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      caption = paste0("High group (n=", n_high, "), Low group (n=", n_low, ")")
    ) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size = 11, face = "bold"),  # Even smaller title font
      plot.subtitle = element_text(size = 9, face = "italic"),  # Smaller subtitle font
      plot.caption = element_text(size = 8, hjust = 0),  # Smaller caption font
      legend.position = "top",
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.margin = margin(6, 6, 6, 6),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 9),  # Smaller axis title font
      axis.text = element_text(size = 8)  # Smaller axis text font
    ) +
    # Explicitly set colors for High and Low groups
    scale_color_manual(values = c("High" = high_color, "Low" = low_color)) +
    scale_fill_manual(values = c("High" = high_color, "Low" = low_color)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = seq(0, max(surv_data$time, na.rm = TRUE) + 50, 50))
  
  # Return the plot directly without the risk table
  return(p)
}

# Function to perform PFS analysis for a signature (using individual samples)
perform_pfs_analysis <- function(signature_idx, signature_name, signature_scores) {
  cat("Performing PFS analysis for signature:", signature_name, "\n")
  
  # Create a data frame with sample names and signature scores
  survival_data <- data.frame(
    Sample = colnames(signature_scores),
    Score = as.numeric(signature_scores[signature_idx, ])
  )
  
  # Print diagnostic information about initial sample count
  cat("  Initial sample count:", nrow(survival_data), "\n")
  
  # Merge with sample progression information
  survival_data <- merge(survival_data, sample_pfs_data, by = "Sample", all.x = TRUE)
  
  # Print diagnostic information after merging
  cat("  Sample count after merging with PFS data:", nrow(survival_data), "\n")
  cat("  Samples with missing PFS data:", sum(is.na(survival_data$PFS)), "\n")
  
  # Remove samples without progression data
  survival_data <- survival_data[!is.na(survival_data$PFS), ]
  
  # Print diagnostic information after filtering
  cat("  Final sample count for analysis:", nrow(survival_data), "\n")
  
  # Check if we have enough data
  if (nrow(survival_data) < 5) {
    cat("Warning: Not enough data for PFS analysis of signature", signature_name, "\n")
    return(list(
      signature = signature_name,
      pval = NA,
      hr = NA,
      hr_lower = NA,
      hr_upper = NA,
      survival_data = survival_data
    ))
  }
  
  # Determine high/low groups based on median score
  median_score <- median(survival_data$Score, na.rm = TRUE)
  survival_data$Group <- factor(ifelse(survival_data$Score >= median_score, "High", "Low"), levels = c("High", "Low"))
  
  # Count samples in each group
  n_high <- sum(survival_data$Group == "High")
  n_low <- sum(survival_data$Group == "Low")
  
  # Print diagnostic information about group sizes
  cat("  High group:", n_high, "samples, Low group:", n_low, "samples\n")
  
  # Check if we have both groups
  if (length(unique(survival_data$Group)) < 2) {
    cat("Warning: Only one group for PFS analysis of signature", signature_name, "\n")
    return(list(
      signature = signature_name,
      pval = NA,
      hr = NA,
      hr_lower = NA,
      hr_upper = NA,
      survival_data = survival_data
    ))
  }
  
  # Create survival object
  surv_obj <- Surv(time = survival_data$PFS, event = rep(1, nrow(survival_data)))
  
  # Fit survival model
  fit <- survfit(surv_obj ~ Group, data = survival_data)
  
  # Log-rank test
  log_rank <- survdiff(surv_obj ~ Group, data = survival_data)
  pval <- 1 - pchisq(log_rank$chisq, df = 1)
  
  # Calculate hazard ratio using Cox proportional hazards model
  cox_model <- tryCatch({
    coxph(surv_obj ~ Group, data = survival_data)
  }, error = function(e) {
    cat("Error in Cox model for PFS analysis of", signature_name, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(cox_model)) {
    cox_summary <- summary(cox_model)
    hr <- exp(cox_summary$coefficients[1, 1])
    hr_lower <- exp(cox_summary$coefficients[1, 1] - 1.96 * cox_summary$coefficients[1, 3])
    hr_upper <- exp(cox_summary$coefficients[1, 1] + 1.96 * cox_summary$coefficients[1, 3])
    
    hr_text <- sprintf("HR (Low vs High) = %.2f (95%% CI: %.2f-%.2f), p = %.3f", hr, hr_lower, hr_upper, pval)
  } else {
    hr <- NA
    hr_lower <- NA
    hr_upper <- NA
    hr_text <- paste0("p = ", format(pval, digits = 3), ", HR not available")
  }
  
  # Create KM plot if p-value < 0.05
  if (!is.na(pval) && pval < 0.05) {
    tryCatch({
      pfs_plot <- create_km_plot(
        fit,
        title = paste0("Progression-Free Survival by ", signature_name),
        hr_text = hr_text,
        n_high = n_high,
        n_low = n_low,
        signature_name = signature_name,
        analysis_type = "PFS"
      )
      
      pdf(file.path(survival_dir, paste0(signature_name, "_PFS.pdf")), width = 8, height = 6)
      print(pfs_plot)
      dev.off()
    }, error = function(e) {
      cat("Error creating PFS plot for", signature_name, ":", e$message, "\n")
    })
  }
  
  # Return results
  return(list(
    signature = signature_name,
    pval = pval,
    hr = hr,
    hr_lower = hr_lower,
    hr_upper = hr_upper,
    n_high = n_high,
    n_low = n_low,
    survival_data = survival_data
  ))
}

# Function to perform OS analysis for a signature (using patient-level data)
perform_os_analysis <- function(signature_idx, signature_name, signature_scores) {
  cat("Performing OS analysis for signature:", signature_name, "\n")
  
  # Create a data frame with sample names and signature scores
  sample_data <- data.frame(
    Sample = colnames(signature_scores),
    Score = as.numeric(signature_scores[signature_idx, ])
  )
  
  # Add patient codes
  sample_data <- merge(sample_data, sample_pfs_data[, c("Sample", "Patient_Code")], by = "Sample", all.x = TRUE)
  
  # Remove samples without patient codes
  sample_data <- sample_data[!is.na(sample_data$Patient_Code), ]
  
  # For patients with multiple samples, use the maximum score
  patient_data <- sample_data %>%
    group_by(Patient_Code) %>%
    summarize(Score = max(Score, na.rm = TRUE))
  
  # Merge with patient OS data
  survival_data <- merge(patient_data, patient_os_data, by = "Patient_Code", all.x = TRUE)
  
  # Remove patients without OS data
  survival_data <- survival_data[!is.na(survival_data$OS), ]
  
  # Check if we have enough data
  if (nrow(survival_data) < 5) {
    cat("Warning: Not enough data for OS analysis of signature", signature_name, "\n")
    return(list(
      signature = signature_name,
      pval = NA,
      hr = NA,
      hr_lower = NA,
      hr_upper = NA,
      survival_data = survival_data
    ))
  }
  
  # Determine high/low groups based on median score
  median_score <- median(survival_data$Score, na.rm = TRUE)
  survival_data$Group <- factor(ifelse(survival_data$Score >= median_score, "High", "Low"), levels = c("High", "Low"))
  
  # Count patients in each group
  n_high <- sum(survival_data$Group == "High")
  n_low <- sum(survival_data$Group == "Low")
  
  # Check if we have both groups
  if (length(unique(survival_data$Group)) < 2) {
    cat("Warning: Only one group for OS analysis of signature", signature_name, "\n")
    return(list(
      signature = signature_name,
      pval = NA,
      hr = NA,
      hr_lower = NA,
      hr_upper = NA,
      survival_data = survival_data
    ))
  }
  
  # Create survival object
  surv_obj <- Surv(time = survival_data$OS, event = rep(1, nrow(survival_data)))
  
  # Fit survival model
  fit <- survfit(surv_obj ~ Group, data = survival_data)
  
  # Log-rank test
  log_rank <- survdiff(surv_obj ~ Group, data = survival_data)
  pval <- 1 - pchisq(log_rank$chisq, df = 1)
  
  # Calculate hazard ratio using Cox proportional hazards model
  cox_model <- tryCatch({
    coxph(surv_obj ~ Group, data = survival_data)
  }, error = function(e) {
    cat("Error in Cox model for OS analysis of", signature_name, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(cox_model)) {
    cox_summary <- summary(cox_model)
    hr <- exp(cox_summary$coefficients[1, 1])
    hr_lower <- exp(cox_summary$coefficients[1, 1] - 1.96 * cox_summary$coefficients[1, 3])
    hr_upper <- exp(cox_summary$coefficients[1, 1] + 1.96 * cox_summary$coefficients[1, 3])
    
    hr_text <- sprintf("HR (Low vs High) = %.2f (95%% CI: %.2f-%.2f), p = %.3f", hr, hr_lower, hr_upper, pval)
  } else {
    hr <- NA
    hr_lower <- NA
    hr_upper <- NA
    hr_text <- paste0("p = ", format(pval, digits = 3), ", HR not available")
  }
  
  # Create KM plot if p-value < 0.05
  if (!is.na(pval) && pval < 0.05) {
    tryCatch({
      os_plot <- create_km_plot(
        fit,
        title = paste0("Overall Survival by ", signature_name),
        hr_text = hr_text,
        n_high = n_high,
        n_low = n_low,
        signature_name = signature_name,
        analysis_type = "OS"
      )
      
      pdf(file.path(survival_dir, paste0(signature_name, "_OS.pdf")), width = 8, height = 6)
      print(os_plot)
      dev.off()
    }, error = function(e) {
      cat("Error creating OS plot for", signature_name, ":", e$message, "\n")
    })
  }
  
  # Return results
  return(list(
    signature = signature_name,
    pval = pval,
    hr = hr,
    hr_lower = hr_lower,
    hr_upper = hr_upper,
    n_high = n_high,
    n_low = n_low,
    survival_data = survival_data
  ))
}

# Perform survival analysis for all signatures
cat("Performing survival analysis for all signatures...\n")
pfs_results <- list()
os_results <- list()

# Use all signatures (not just top 10)
for (i in 1:nrow(signature_priority)) {
  signature_idx <- which(rownames(segment_scores) == signature_priority$original_segment[i])
  if (length(signature_idx) > 0) {
    signature_name <- gsub("[^a-zA-Z0-9]", "_", signature_priority$segment[i])
    
    # For the first signature, print detailed diagnostic information
    if (i == 1) {
      cat("\nDetailed analysis for first signature:", signature_name, "\n")
      
      # Create a data frame with sample names and signature scores
      survival_data <- data.frame(
        Sample = colnames(segment_scores),
        Score = as.numeric(segment_scores[signature_idx, ])
      )
      
      cat("Initial sample count:", nrow(survival_data), "\n")
      
      # Merge with sample progression information
      survival_data <- merge(survival_data, sample_pfs_data, by = "Sample", all.x = TRUE)
      
      cat("Sample count after merging with PFS data:", nrow(survival_data), "\n")
      cat("Samples with missing PFS data:", sum(is.na(survival_data$PFS)), "\n")
      
      # List samples with missing PFS data
      missing_pfs_samples <- survival_data$Sample[is.na(survival_data$PFS)]
      if (length(missing_pfs_samples) > 0) {
        cat("Samples missing PFS data:", paste(missing_pfs_samples, collapse=", "), "\n")
      }
      
      # Remove samples without progression data
      survival_data <- survival_data[!is.na(survival_data$PFS), ]
      
      cat("Final sample count for analysis:", nrow(survival_data), "\n")
      
      # Determine high/low groups based on median score
      median_score <- median(survival_data$Score, na.rm = TRUE)
      survival_data$Group <- factor(ifelse(survival_data$Score >= median_score, "High", "Low"), levels = c("High", "Low"))
      
      # Count samples in each group
      n_high <- sum(survival_data$Group == "High")
      n_low <- sum(survival_data$Group == "Low")
      
      cat("High group:", n_high, "samples, Low group:", n_low, "samples\n")
      cat("Median score used for stratification:", median_score, "\n")
      
      # Continue with regular analysis
    }
    
    # Perform PFS analysis (using individual samples)
    tryCatch({
      pfs_results[[i]] <- perform_pfs_analysis(signature_idx, signature_name, segment_scores)
    }, error = function(e) {
      cat("Error in PFS analysis for", signature_name, ":", e$message, "\n")
    })
    
    # Perform OS analysis (using patient-level data)
    tryCatch({
      os_results[[i]] <- perform_os_analysis(signature_idx, signature_name, segment_scores)
    }, error = function(e) {
      cat("Error in OS analysis for", signature_name, ":", e$message, "\n")
    })
  } else {
    cat("Warning: Signature not found in segment scores:", signature_priority$original_segment[i], "\n")
  }
}

# Filter out NULL results
pfs_results <- pfs_results[!sapply(pfs_results, is.null)]
os_results <- os_results[!sapply(os_results, is.null)]

# Create summary tables
if (length(pfs_results) > 0) {
  pfs_summary <- data.frame(
    Signature = sapply(pfs_results, function(x) x$signature),
    Gene_Annotation = sapply(pfs_results, function(x) extract_gene_annotation(x$signature)),
    PFS_pvalue = sapply(pfs_results, function(x) x$pval),
    PFS_HR = sapply(pfs_results, function(x) x$hr),
    PFS_HR_lower = sapply(pfs_results, function(x) x$hr_lower),
    PFS_HR_upper = sapply(pfs_results, function(x) x$hr_upper),
    High_n = sapply(pfs_results, function(x) ifelse(is.null(x$n_high), NA, x$n_high)),
    Low_n = sapply(pfs_results, function(x) ifelse(is.null(x$n_low), NA, x$n_low))
  )
  
  # Sort by p-value
  pfs_summary <- pfs_summary[order(pfs_summary$PFS_pvalue), ]
  
  # Save summary table
  write.csv(pfs_summary, file.path(survival_dir, "pfs_summary.csv"), row.names = FALSE)
  
  # Print significant results
  cat("\nSignificant PFS results (p < 0.05):\n")
  print(pfs_summary[pfs_summary$PFS_pvalue < 0.05, ])
}

if (length(os_results) > 0) {
  os_summary <- data.frame(
    Signature = sapply(os_results, function(x) x$signature),
    Gene_Annotation = sapply(os_results, function(x) extract_gene_annotation(x$signature)),
    OS_pvalue = sapply(os_results, function(x) x$pval),
    OS_HR = sapply(os_results, function(x) x$hr),
    OS_HR_lower = sapply(os_results, function(x) x$hr_lower),
    OS_HR_upper = sapply(os_results, function(x) x$hr_upper),
    High_n = sapply(os_results, function(x) ifelse(is.null(x$n_high), NA, x$n_high)),
    Low_n = sapply(os_results, function(x) ifelse(is.null(x$n_low), NA, x$n_low))
  )
  
  # Sort by p-value
  os_summary <- os_summary[order(os_summary$OS_pvalue), ]
  
  # Save summary table
  write.csv(os_summary, file.path(survival_dir, "os_summary.csv"), row.names = FALSE)
  
  # Print significant results
  cat("\nSignificant OS results (p < 0.05):\n")
  print(os_summary[os_summary$OS_pvalue < 0.05, ])
}

# Create a combined summary table
if (length(pfs_results) > 0 && length(os_results) > 0) {
  combined_summary <- merge(pfs_summary, os_summary, by = c("Signature", "Gene_Annotation"), all = TRUE, suffixes = c("_PFS", "_OS"))
  write.csv(combined_summary, file.path(survival_dir, "combined_summary.csv"), row.names = FALSE)
  
  # Create a combined PDF with all significant results
  significant_signatures <- combined_summary[(!is.na(combined_summary$PFS_pvalue) & combined_summary$PFS_pvalue < 0.05) | 
                                            (!is.na(combined_summary$OS_pvalue) & combined_summary$OS_pvalue < 0.05), ]
  
  if (nrow(significant_signatures) > 0) {
    cat("\nCreating combined plot for all significant signatures...\n")
    
    # Create individual plots for each significant signature
    all_plots <- list()
    plot_count <- 0
    
    # Store signature indices for correlation analysis
    sig_indices <- list()
    
    for (i in 1:nrow(significant_signatures)) {
      sig_name <- significant_signatures$Signature[i]
      
      # Find the signature index in segment_scores
      signature_idx <- which(rownames(segment_scores) == signature_priority$original_segment[which(gsub("[^a-zA-Z0-9]", "_", signature_priority$segment) == sig_name)])
      if (length(signature_idx) > 0) {
        sig_indices[[sig_name]] <- signature_idx
      }
      
      # PFS plot if significant
      if (!is.na(significant_signatures$PFS_pvalue[i]) && significant_signatures$PFS_pvalue[i] < 0.05) {
        pfs_idx <- which(sapply(pfs_results, function(x) x$signature) == sig_name)
        
        if (length(pfs_idx) > 0) {
          result <- pfs_results[[pfs_idx]]
          
          surv_data <- result$survival_data
          surv_obj <- Surv(time = surv_data$PFS, event = rep(1, nrow(surv_data)))
          fit <- survfit(surv_obj ~ Group, data = surv_data)
          
          hr_text <- sprintf("HR (Low vs High) = %.2f (95%% CI: %.2f-%.2f), p = %.3f", 
                            result$hr, result$hr_lower, result$hr_upper, result$pval)
          
          tryCatch({
            plot_count <- plot_count + 1
            all_plots[[plot_count]] <- create_km_plot(
              fit,
              title = paste0("Progression-Free Survival by ", sig_name),
              hr_text = hr_text,
              n_high = result$n_high,
              n_low = result$n_low,
              signature_name = sig_name,
              analysis_type = "PFS"
            )
          }, error = function(e) {
            cat("Error creating combined PFS plot for", sig_name, ":", e$message, "\n")
          })
        }
      }
      
      # OS plot if significant
      if (!is.na(significant_signatures$OS_pvalue[i]) && significant_signatures$OS_pvalue[i] < 0.05) {
        os_idx <- which(sapply(os_results, function(x) x$signature) == sig_name)
        
        if (length(os_idx) > 0) {
          result <- os_results[[os_idx]]
          
          surv_data <- result$survival_data
          surv_obj <- Surv(time = surv_data$OS, event = rep(1, nrow(surv_data)))
          fit <- survfit(surv_obj ~ Group, data = surv_data)
          
          hr_text <- sprintf("HR (Low vs High) = %.2f (95%% CI: %.2f-%.2f), p = %.3f", 
                            result$hr, result$hr_lower, result$hr_upper, result$pval)
          
          tryCatch({
            plot_count <- plot_count + 1
            all_plots[[plot_count]] <- create_km_plot(
              fit,
              title = paste0("Overall Survival by ", sig_name),
              hr_text = hr_text,
              n_high = result$n_high,
              n_low = result$n_low,
              signature_name = sig_name,
              analysis_type = "OS"
            )
          }, error = function(e) {
            cat("Error creating combined OS plot for", sig_name, ":", e$message, "\n")
          })
        }
      }
    }
    
    # Create a multi-panel figure with all plots
    if (length(all_plots) > 0) {
      # Calculate grid dimensions
      n_plots <- length(all_plots)
      n_cols <- min(2, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      # Create the combined plot
      pdf(file.path(survival_dir, "all_significant_signatures.pdf"), width = 12, height = 8)
      
      # Use grid.arrange to create a multi-panel figure
      do.call(gridExtra::grid.arrange, c(all_plots, ncol = n_cols))
      
      # Also save individual plots to the combined PDF
      for (p in all_plots) {
        print(p)
      }
      
      dev.off()
      
      # Create a separate summary figure with all plots
      pdf(file.path(survival_dir, "summary_figure.pdf"), width = 12, height = 10)
      do.call(gridExtra::grid.arrange, c(all_plots, ncol = n_cols))
      dev.off()
    }
    
    # Create correlation plot for significant signatures
    if (length(sig_indices) > 0) {
      cat("\nCalculating correlations between significant signatures...\n")
      
      # Extract scores for significant signatures
      sig_scores <- data.frame(matrix(ncol = length(sig_indices), nrow = ncol(segment_scores)))
      colnames(sig_scores) <- names(sig_indices)
      
      for (i in 1:length(sig_indices)) {
        sig_name <- names(sig_indices)[i]
        sig_scores[, i] <- as.numeric(segment_scores[sig_indices[[sig_name]], ])
      }
      
      # Add sample names
      sig_scores$Sample <- colnames(segment_scores)
      
      # Merge with PFS data to filter out samples without progression data
      sig_scores <- merge(sig_scores, sample_pfs_data[, c("Sample", "PFS")], by = "Sample", all.x = TRUE)
      sig_scores <- sig_scores[!is.na(sig_scores$PFS), ]
      
      # Calculate correlation matrix
      cor_matrix <- cor(sig_scores[, names(sig_indices)], method = "pearson")
      
      # Create correlation plot
      pdf(file.path(survival_dir, "signature_correlations.pdf"), width = 8, height = 7)
      
      # Create a more informative title for each signature
      sig_labels <- sapply(names(sig_indices), function(name) {
        gene_annotation <- extract_gene_annotation(name)
        if (gene_annotation != "") {
          return(paste0(name, "\n(", gene_annotation, ")"))
        } else {
          return(name)
        }
      })
      
      # Plot correlation matrix
      corrplot(cor_matrix, method = "circle", type = "upper", 
               tl.col = "black", tl.srt = 45, tl.cex = 0.8,
               col = colorRampPalette(c("#0000FF", "white", "#FF0000"))(100),
               title = "Pearson Correlation Between Significant Signatures",
               mar = c(0, 0, 2, 0))
      
      # Also create a scatter plot matrix for more detailed view
      pairs(sig_scores[, names(sig_indices)], 
            main = "Scatter Plot Matrix of Significant Signatures",
            pch = 19, col = rgb(0, 0, 1, 0.5),
            labels = sig_labels)
      
      dev.off()
      
      # Create a heatmap of the signature scores
      sig_scores_matrix <- as.matrix(sig_scores[, names(sig_indices)])
      rownames(sig_scores_matrix) <- sig_scores$Sample
      
      # Scale the data for better visualization
      sig_scores_scaled <- scale(sig_scores_matrix)
      
      # Create heatmap
      pdf(file.path(survival_dir, "signature_heatmap.pdf"), width = 10, height = 8)
      heatmap(sig_scores_scaled, 
              main = "Heatmap of Significant Signature Scores",
              col = colorRampPalette(c("blue", "white", "red"))(100),
              scale = "none",  # Already scaled
              margins = c(10, 5),
              cexRow = 0.5,    # Smaller text for row labels
              cexCol = 0.8)    # Slightly smaller text for column labels
      dev.off()
    }
  }
}

cat("Survival analysis complete! Results are available in the survival directory.\n") 