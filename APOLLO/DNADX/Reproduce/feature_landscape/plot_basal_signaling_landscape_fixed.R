library(readxl)

# Source helper functions
source("model_scripts/helper.R")

# Read the Excel file containing the elastic net weights
excel_path <- "41467_2019_13588_MOESM4_ESM_Elastic_Net_gene_signatures.xlsx"
df <- read_excel(excel_path)

# Find the Basal signaling signature row
signature_name <- "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP"
signature_row <- which(df[[1]] == signature_name)

print(paste("Found signature at row:", signature_row))

# Extract the weights for this signature (skip first column which contains signature names)
weights <- as.numeric(df[signature_row, 2:ncol(df)])
names(weights) <- colnames(df)[2:ncol(df)]

# Remove NA weights
weights <- weights[!is.na(weights)]

print(paste("Number of features with weights:", length(weights)))
print(paste("Number of non-zero weights:", sum(weights != 0)))

# Load segment annotation
load("segment_anno.rda")  # This should contain segment_anno and vertical_lines

# Filter weights to only include segments that exist in segment_anno and have valid coordinates
valid_segments <- intersect(names(weights), rownames(segment_anno))
print(paste("Number of valid segments found in annotation:", length(valid_segments)))

# Check for invalid coordinates in segment_anno
invalid_coords <- c()
for (seg in valid_segments) {
  start_pos <- segment_anno[seg, 2]
  end_pos <- segment_anno[seg, 3]
  if (is.na(start_pos) || is.na(end_pos) || !is.finite(start_pos) || !is.finite(end_pos)) {
    invalid_coords <- c(invalid_coords, seg)
  }
}

if (length(invalid_coords) > 0) {
  print(paste("Removing", length(invalid_coords), "segments with invalid coordinates"))
  valid_segments <- setdiff(valid_segments, invalid_coords)
}

# Filter weights to only valid segments
weights_filtered <- weights[valid_segments]
print(paste("Number of segments after filtering:", length(weights_filtered)))
print(paste("Number of non-zero weights after filtering:", sum(weights_filtered != 0)))

# Load vertical lines for chromosome boundaries if not in segment_anno.rda
if (!exists("vertical_lines")) {
  # Define chromosome boundaries (these would typically be loaded from data)
  vertical_lines <- c(1, 2229, 4329, 5870, 7265, 8537, 9573, 10649, 11602, 12478, 
                     13253, 14026, 14633, 15287, 15952, 16649, 17315, 17952, 
                     18609, 19201, 19822, 20401, 21030, 21653, 22276, 24776)
}

# Plot the feature landscape
plot_seg_ss(weights_filtered, paste(signature_name, "_fixed", sep=""))

print("Feature landscape plot saved!") 