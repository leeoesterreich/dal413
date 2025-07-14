# Script to load and display Elastic Net model parameters (alpha and lambda)

# Load the saved model results
results_file <- "results/elastic_net_models.rda"
if (file.exists(results_file)) {
  load(results_file) # This should load 'model_results'
  
  if (exists("model_results") && length(model_results) > 0) {
    cat("--- Elastic Net Model Parameters (Best Tune) ---\n")
    
    for (pheno_name in names(model_results)) {
      cat(paste("\nPhenotype:", pheno_name, "\n"))
      result_item <- model_results[[pheno_name]]
      
      # Check if this phenotype was skipped due to an error during modeling
      if (!is.null(result_item$error)) {
        cat(paste("  Error during modeling:", result_item$error, "\n"))
      } else if (is.null(result_item$model) || is.null(result_item$model$bestTune)) {
        cat("  Model or bestTune not found for this phenotype.\n")
      } else {
        best_alpha <- result_item$model$bestTune$alpha
        best_lambda <- result_item$model$bestTune$lambda
        cat(paste("  Best Alpha:", best_alpha, "\n"))
        cat(paste("  Best Lambda:", best_lambda, "\n"))
      }
    }
  } else {
    cat("No model results found in", results_file, "or the list is empty.\n")
  }
} else {
  cat("Results file not found:", results_file, "\nPlease run the Elastic_Net_modeling.R script first.\n")
} 