# Script to check for a specific signature in signature_score.rda

signature_to_check <- "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP"
results_file <- "results/signature_score.rda"

if (file.exists(results_file)) {
  load(results_file) # This should load 'signature_score'
  
  if (exists("signature_score")) {
    cat(paste("Checking for signature:", signature_to_check, "\n"))
    
    if (signature_to_check %in% rownames(signature_score)) {
      cat("Signature FOUND in signature_score.rda.\n")
      score_values <- signature_score[signature_to_check, ]
      num_total_values <- length(score_values)
      num_na_values <- sum(is.na(score_values))
      cat(paste("  Total values:", num_total_values, "\n"))
      cat(paste("  Number of NA values:", num_na_values, "\n"))
      cat(paste("  Number of non-NA values:", num_total_values - num_na_values, "\n"))
      if (num_total_values - num_na_values > 0) {
        cat("  First few non-NA values:\n")
        print(head(score_values[!is.na(score_values)], 5))
      } else {
        cat("  All values are NA for this signature.\n")
      }
    } else {
      cat("Signature NOT FOUND in signature_score.rda.\n")
      cat("Available signatures (first 10):
")
      print(head(rownames(signature_score),10))
    }
  } else {
    cat("'signature_score' object not found in", results_file, "\n")
  }
} else {
  cat("Results file not found:", results_file, "\nPlease ensure signature_score_and_segment_score_calculation.R has run successfully.\n")
} 