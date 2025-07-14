library(sampling)
library(caret)
library(glmnet)
library(ROCR)
library(doMC)
registerDoMC(cores = 8)
set.seed(1992)
source("Rscript/helper.R")

# Load the data
load('results/signature_score.rda')
CNdata <- read.table("Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt", header=TRUE, sep="\t", check.names=FALSE)
CN_score <- as.matrix(CNdata[,-1])
rownames(CN_score) <- CNdata[,1]

# Initialize results list
model_results <- list()

# Run elastic net for each signature
for(pheno in rownames(signature_score)) {
    score <- unlist(signature_score[pheno,])
    
    # Handle NAs in score for quantile calculation and binarization
    quant_val <- quantile(score, 0.67, na.rm = TRUE)
    
    score_bi <- rep(NA, length(score)) # Initialize with NAs
    non_na_indices <- !is.na(score)
    
    if (is.na(quant_val)) {
        # This happens if all scores are NA or too few non-NA for quantile
        warning(paste("Quantile could not be computed for phenotype:", pheno, "due to NAs or insufficient data. Skipping this phenotype."))
        model_results[[pheno]] <- list(error = "Could not compute quantile due to NAs/insufficient data")
        next # Skip to the next phenotype
    }
    
    score_bi[non_na_indices] <- ifelse(score[non_na_indices] >= quant_val, 1, 0)
    
    # Identify samples with non-NA binarized scores for modeling
    valid_samples_idx <- which(!is.na(score_bi))
    
    if(length(valid_samples_idx) < 10) { # Arbitrary threshold: need enough samples for splitting and modeling
        warning(paste("Too few valid (non-NA) samples for phenotype:", pheno, "after binarization. Skipping this phenotype."))
        model_results[[pheno]] <- list(error = "Too few valid samples after binarization")
        next
    }

    # Subset CN_score and score_bi to only valid samples
    current_CN_score <- CN_score[, valid_samples_idx, drop = FALSE]
    current_score_bi <- score_bi[valid_samples_idx]

    # Ensure CN_score has rows (genes) and current_CN_score has columns (samples)
    # The original script had CN_score[train,] which implies samples are rows in CN_score
    # Let's adjust based on the Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt format (genes x samples)
    # So, current_CN_score should be samples x genes for caret.
    # Transpose current_CN_score before splitting.
    # However, caret usually expects predictors as samples x features.
    # Let's re-verify: CNdata is genes x samples. CN_score is genes x samples.
    # trainX <- CN_score[train,] implies train is indices for *columns* if CN_score is genes x samples.
    # This seems to be a mismatch with how balancedstratification might work or how caret expects data.

    # The original script: trainX <- CN_score[train,]
    # If CN_score is genes x samples, this selects *columns* (samples).
    # Let's assume CN_score is genes x samples as loaded.
    # For caret, trainX should be samples x features (genes).
    # So, we need to transpose the subset of CN_score corresponding to the current samples.

    # Filter samples for balancing_variables if it exists and is used
    # The original script used balancedstratification(balancing_variables,strata,pik,comment = F)
    # balancing_variables is not loaded in the current version. We'll do a simple random split.
    # If balancing_variables is needed, it must be loaded and sample-matched like other data.

    set.seed(1992) # for reproducibility
    train_fraction <- 0.7
    num_valid_samples <- length(current_score_bi)
    train_size <- floor(train_fraction * num_valid_samples)
    
    shuffled_indices <- sample(num_valid_samples)
    train_indices_in_valid <- shuffled_indices[1:train_size]
    test_indices_in_valid <- shuffled_indices[(train_size + 1):num_valid_samples]

    # trainX should be samples x genes
    trainX <- t(current_CN_score[, train_indices_in_valid, drop = FALSE])
    testX <- t(current_CN_score[, test_indices_in_valid, drop = FALSE])
    trainY <- current_score_bi[train_indices_in_valid]
    testY <- current_score_bi[test_indices_in_valid]

    if(nrow(trainX) < 2 || nrow(testX) < 2 || length(unique(trainY)) < 2) {
        warning(paste("Not enough data or variance in outcome for phenotype:", pheno, "after splitting. Skipping."))
        model_results[[pheno]] <- list(error = "Insufficient data/variance after splitting")
        next
    }

    glmnet_obj <- caret_wrap(trainX,trainY,testX,testY,bi = T)
    
    # Get predictions
    pred_train <- predict(glmnet_obj,newdata = trainX,type = 'prob')
    pred_test <- predict(glmnet_obj,newdata = testX,type = 'prob')
    pred_train_obj <- prediction(pred_train$pos,labels = trainY)
    pred_test_obj <- prediction(pred_test$pos,labels = testY)
    
    # Calculate AUC
    auc_train <- signif(performance(pred_train_obj,measure = 'auc')@y.values[[1]][1],2)
    auc_test <- signif(performance(pred_test_obj,measure = 'auc')@y.values[[1]][1],2)
    
    # Get feature coefficients
    beta <- as.matrix(coef(glmnet_obj$finalModel,glmnet_obj$bestTune$lambda))
    beta <- beta[-1,1]
    
    # Store results
    model_results[[pheno]] <- list(
        model = glmnet_obj,
        train_auc = auc_train,
        test_auc = auc_test,
        coefficients = beta
    )
}

# Save all results
save(model_results, file='results/elastic_net_models.rda')






