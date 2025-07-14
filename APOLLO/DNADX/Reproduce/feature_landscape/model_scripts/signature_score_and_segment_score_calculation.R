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
edata <- read_excel("../41467_2019_13588_MOESM2_ESM.xlsx", sheet = 1)
CNdata <- read_excel("../41467_2019_13588_MOESM3_ESM.xlsx", sheet = 1)

# Print column names to check structure
print("edata column names:")
print(colnames(edata))
print("\nCNdata column names:")
print(colnames(CNdata))

# Convert data frames to matrices with proper row names
edata_matrix <- as.matrix(edata[,-1])  # Remove first column (gene names)
rownames(edata_matrix) <- as.character(edata[[1]])  # Use first column as row names

CNdata_matrix <- as.matrix(CNdata[,-1])  # Remove first column (gene names)
rownames(CNdata_matrix) <- as.character(CNdata[[1]])  # Use first column as row names

# Given a gene expression data matrix (gene X sample): edata
# rows are genes in Entrez ID and columns are samples
# run calc_signatures
signature_score <- calc_signatures(edata_matrix,"data/gene_signatures_20170111.gmt",method = "median")

# find NA signatures
NAsig <- rownames(signature_score[is.na(signature_score[,1]),])
# remove NA signatures
i <- match(NAsig,rownames(signature_score))
signature_score <- signature_score[-i,]

# CD103_Ratio
CD103_pos <- signature_score[rownames(signature_score)=="CD103_Positive_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_neg <- signature_score[rownames(signature_score)=="CD103_Negative_Median_Cancer.Cell.2014_PMID.25446897"]
CD103_ratio <- CD103_pos - CD103_neg # log2 scale division 
signature_score <- rbind(signature_score,CD103_ratio)
rownames(signature_score)[nrow(signature_score)] <- "CD103_Ratio_Cancer.Cell.2014_PMID.25446897"

# differentiation score
diff_centroid <- read.table("data/special_gene_signature_training_sets/UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035.txt",sep = '\t',header = T,row.names = 1,check.names = F)
diff_score <- assignDiffScore.dwd(diff_centroid,edata_matrix)
signature_score <- rbind(signature_score,diff_score)
rownames(signature_score)[nrow(signature_score)] <- "UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035"

# oncotype DX score
oncotype <- GHI_RS(edata_matrix)
signature_score <- rbind(signature_score,diff_score)
rownames(signature_score)[nrow(signature_score)] <- "GHI_RS_Model_NJEM.2004_PMID.15591335"

# Save results in the results directory
save(signature_score,file = 'results/signature_score.rda')

# For special signatures calculated as correlation to predetermined gene centroids
# All traing sets files are included in the ~/data/special_gene_signature_training_sets folder
# For calculation of such special signature, DWD was used to merge current edata with 
# training set, then DWD prediction tool was used to compute pearson/spearman correlation/euclidean distance
# for each sample

# Given a gene-level CNA score matrix (gene X sample): CNdata
# rows are genes in Entrez ID and columns are samples
segment_score <- calc_segments(CNdata_matrix,'data/CNA_segments.gmt',method = 'mean')

# Save segment scores in the results directory
save(segment_score,file = 'results/segment_score.rda')





