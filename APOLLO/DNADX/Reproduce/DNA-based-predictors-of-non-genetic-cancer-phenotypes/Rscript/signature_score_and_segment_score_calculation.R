# This script is for calculating gene signature scores based on RNA and segment scores based on gene-level DNA CNA

# Install and load required packages
if (!require("readxl")) {
    install.packages("readxl", repos="http://cran.r-project.org")
    library(readxl)
}
if (!require("GSA")) {
    install.packages("GSA", repos="http://cran.r-project.org")
    library(GSA)
}

# Source helper functions
source("Rscript/helper.R")

# Load input data
edata <- read.table("HiSeqV2", header=TRUE, sep="\t", check.names=FALSE)
CNdata <- read.table("Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt", header=TRUE, sep="\t", check.names=FALSE)

# Convert data frames to matrices with proper row names
edata_matrix <- as.matrix(edata[,-1])  # Remove first column (gene names)
rownames(edata_matrix) <- as.character(edata[[1]])  # Use first column as row names

CNdata_matrix <- as.matrix(CNdata[,-1])  # Remove first column (gene names)
rownames(CNdata_matrix) <- as.character(CNdata[[1]])  # Use first column as row names

# --- Debug: Check GMT loading and initial signature calculation ---
cat("Attempting to load GMT file: data/gene_signatures_20170111.gmt\n")
temp_geneset_obj <- tryCatch(GSA.read.gmt("data/gene_signatures_20170111.gmt"), error = function(e) { cat("Error loading GMT file:", e$message, "\n"); NULL })

if (!is.null(temp_geneset_obj)) {
    cat("GMT file loaded. Number of gene sets found:", length(temp_geneset_obj$genesets), "\n")
    cat("First 10 gene set names from GMT:
")
    print(head(temp_geneset_obj$geneset.names, 10))
    # Check if your specific signature is in the loaded names
    target_gsea_sig <- "GSEA_Median_GP17_Basal_signaling_r0_958_SMID_BREAST_CANCER_BASAL_UP"
    if (target_gsea_sig %in% temp_geneset_obj$geneset.names) {
        cat(paste("Target signature '", target_gsea_sig, "' FOUND in GMT loaded names.\n"))
    } else {
        cat(paste("Target signature '", target_gsea_sig, "' NOT FOUND in GMT loaded names. Checking for similar names...\n"))
        # Very simple partial match check (case-insensitive)
        matches <- temp_geneset_obj$geneset.names[grepl("BASAL_SIGNALING", toupper(temp_geneset_obj$geneset.names))]
        if(length(matches) > 0) {
            cat("Potential matches found in GMT for 'BASAL_SIGNALING':\n")
            print(matches)
        }
    }
} else {
    cat("Failed to load or parse GMT file. calc_signatures will likely fail or produce empty results for GMT signatures.\n")
}
# -------------------------------------------------------------------

# Calculate signature scores
cat("Calculating signature scores from GMT...\n")
signature_score <- calc_signatures(edata_matrix, "data/gene_signatures_20170111.gmt", method = "median")

# --- Debug: Check signature_score immediately after calc_signatures ---
cat("\n--- signature_score immediately after calc_signatures (from GMT) ---\n")
if (exists("signature_score") && !is.null(signature_score) && !is.null(dim(signature_score))) {
    cat(paste("Dimensions:", paste(dim(signature_score), collapse=" x "), "\n"))
    cat("First 10 rownames (signatures) from GMT processing:
")
    print(head(rownames(signature_score), 10))
    # Check for your specific signature again
    if (target_gsea_sig %in% rownames(signature_score)) {
        cat(paste("Target signature '", target_gsea_sig, "' FOUND in initial signature_score object.\n"))
    } else {
        cat(paste("Target signature '", target_gsea_sig, "' NOT FOUND in initial signature_score object.\n"))
    }
} else {
    cat("signature_score object is NULL, empty, or not a matrix/data.frame after calc_signatures.\n")
}
cat("---------------------------------------------------------------------\n\n")
# ---------------------------------------------------------------------

# find NA signatures
NAsig <- rownames(signature_score[is.na(signature_score[,1]),])
# remove NA signatures
i <- match(NAsig,rownames(signature_score))
if (length(i) > 0) {
    signature_score <- signature_score[-i,]
}

# CD103_Ratio
CD103_pos <- signature_score[rownames(signature_score)=="CD103_Positive_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_neg <- signature_score[rownames(signature_score)=="CD103_Negative_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_ratio <- CD103_pos - CD103_neg # log2 scale division 
signature_score <- rbind(signature_score,CD103_ratio)
rownames(signature_score)[nrow(signature_score)] <- "CD103_Ratio_Cancer.Cell.2014_PMID.25446897"

# differentiation score
print("Reading differentiation score training data...")
diff_centroid <- read.table("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data/special_gene_signature_training_sets/UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035.txt", sep = '\t', header = T, row.names = 1, check.names = F)
print("Structure of diff_centroid:")
print(str(diff_centroid))
print("Head of diff_centroid:")
print(head(diff_centroid))
print("Structure of edata_matrix:")
print(str(edata_matrix))
print("Head of edata_matrix:")
print(head(edata_matrix))

diff_score <- assignDiffScore.dwd(diff_centroid,edata_matrix)
signature_score <- rbind(signature_score,diff_score)
rownames(signature_score)[nrow(signature_score)] <- "UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035"

# oncotype DX score
oncotype <- GHI_RS(edata_matrix)
signature_score <- rbind(signature_score,oncotype)
rownames(signature_score)[nrow(signature_score)] <- "GHI_RS_Model_NJEM.2004_PMID.15591335"

# Save results in the results directory
save(signature_score,file = 'results/signature_score.rda')

# Calculate segment scores
segment_score <- calc_segments(CNdata_matrix,'data/CNA_segments.gmt',method = 'mean')

# Save segment scores in the results directory
save(segment_score,file = 'results/segment_score.rda')





