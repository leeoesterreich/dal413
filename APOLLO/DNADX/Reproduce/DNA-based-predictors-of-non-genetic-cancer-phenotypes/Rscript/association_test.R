source("Rscript/helper.R")

# Test for associations between gene signatures and DNA CNA
# Given a signature_score matrix and 
# gene-level CNA matrix CN_score: continuous CNA score for each gene, 
# CN_gain: binary matrix where 1 is copy number gain for a gene and 0 is no gain
# CN_loss: binary matrix where 1 is copy number loss for a gene and 0 is no loss
# Fisher's exact test and spearman correlation test/linear model control for subtypes
# choose a gene signature to test, for example RB-LOH signature

# Load the required data
load('results/signature_score.rda')
# load('results/segment_score.rda') # Not directly used for CN_matrix for sigCNTest

# Read CNA data for binary matrices and CN_matrix (CN_score in sigCNTest)
CNdata <- read.table("Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt", header=TRUE, sep="\t", check.names=FALSE)
CN_matrix <- as.matrix(CNdata[,-1])
rownames(CN_matrix) <- CNdata[,1]

# --- Debugging: Print structure of signature_score and CN_matrix ---
print("--- signature_score (RNA-based) --- ")
print(paste("Dimensions:", paste(dim(signature_score), collapse=" x ")))
print("First 5 column names (samples):")
print(colnames(signature_score)[1:min(5, ncol(signature_score))])
print("Number of columns (samples):")
print(ncol(signature_score))

print("--- CN_matrix (CNA-based) --- ")
print(paste("Dimensions:", paste(dim(CN_matrix), collapse=" x ")))
print("First 5 column names (samples):")
print(colnames(CN_matrix)[1:min(5, ncol(CN_matrix))])
print("Number of columns (samples):")
print(ncol(CN_matrix))

# --- Match samples between signature_score and CN_matrix ---
common_samples <- intersect(colnames(signature_score), colnames(CN_matrix))
if (length(common_samples) == 0) {
  stop("No common samples found between signature_score and CN_matrix. Cannot proceed.")
}
if (length(common_samples) < ncol(signature_score) || length(common_samples) < ncol(CN_matrix)){
    warning("Mismatch in samples between signature_score and CN_matrix. Subsetting to common samples.")
}

print(paste("Number of common samples:", length(common_samples)))
signature_score_matched <- signature_score[, common_samples, drop=FALSE]
CN_matrix_matched <- CN_matrix[, common_samples, drop=FALSE]

print("--- signature_score_matched (RNA-based, after matching) --- ")
print(paste("Dimensions:", paste(dim(signature_score_matched), collapse=" x ")))
print("Number of columns (samples):")
print(ncol(signature_score_matched))

print("--- CN_matrix_matched (CNA-based, after matching) --- ")
print(paste("Dimensions:", paste(dim(CN_matrix_matched), collapse=" x ")))
print("Number of columns (samples):")
print(ncol(CN_matrix_matched))
# ---------------------------------------------------------------------

# Create binary matrices for gains and losses from the matched CN_matrix
CN_gain <- CN_matrix_matched > 0.3
CN_loss <- CN_matrix_matched < -0.3
# Ensure they remain matrices with correct dimensions and names
CN_gain <- matrix(as.numeric(CN_gain), nrow=nrow(CN_matrix_matched), dimnames=dimnames(CN_matrix_matched))
CN_loss <- matrix(as.numeric(CN_loss), nrow=nrow(CN_matrix_matched), dimnames=dimnames(CN_matrix_matched))

# Test associations for each signature using matched data
results <- list()
# Use signature_score_matched for the loop
for(pheno in rownames(signature_score_matched)) {
    score <- unlist(signature_score_matched[pheno,])
    # Pass CN_matrix_matched to sigCNTest, which expects the full gene x sample matrix
    p_value <- sigCNTest(score, CN_matrix_matched, CN_gain, CN_loss)
    
    # Benjamini-Hochberg correct p values
    p_value_adj <- apply(p_value, 2, function(x){return(p.adjust(x, method = "BH"))})
    # log10 transform adjusted p values
    log_p <- -log(p_value_adj,10)
    
    results[[pheno]] <- list(
        raw_pvalues = p_value,
        adjusted_pvalues = p_value_adj,
        log_pvalues = log_p
    )
}

# Save results
save(results, file='results/association_test_results.rda')
