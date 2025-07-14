library(readxl)

# Source helper functions
source("model_scripts/helper.R")

# Load segment annotation
load("segment_anno.rda")

# Create proper vertical_lines for chromosomes 1-22, X
vertical_lines <- c(0, 2.4725, 4.9020, 6.8970, 8.8097, 10.6183, 12.3273, 13.9155, 15.3782, 16.7809,
                   18.1346, 19.4791, 20.8026, 21.9440, 23.0077, 24.0111, 24.8994, 25.6871, 26.4483,
                   27.0864, 27.7108, 28.1802, 28.6771, 30.2265)

print("=== SETUP ===")
print(paste("Using vertical_lines with", length(vertical_lines), "boundaries"))
print(paste("Loaded segment_anno with", nrow(segment_anno), "segments"))

# Read the corrected Excel file
excel_path <- "41467_2019_13588_MOESM4_ESM_Elastic_Net_gene_signatures.xlsx"
df <- read_excel(excel_path)

# Find the Basal signaling signature row
signature_name <- "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP"
signature_row <- which(df[[1]] == signature_name)
print(paste("Found signature at row:", signature_row))

# Extract weights
weights <- as.numeric(df[signature_row, 2:ncol(df)])
segment_names <- colnames(df)[2:ncol(df)]
names(weights) <- segment_names

# Filter to matching segments
matching_segments <- intersect(segment_names, rownames(segment_anno))
weights_filtered <- weights[matching_segments]
weights_filtered <- weights_filtered[!is.na(weights_filtered)]

print(paste("Number of valid weights:", length(weights_filtered)))
print(paste("Number of non-zero weights:", sum(weights_filtered != 0)))

# Enhanced plotting function matching original paper style
plot_basal_final <- function(beta, main) {
  print("Creating final plot matching original paper style...")
  
  total_gene <- 30.2265
  min_y <- min(beta, na.rm = TRUE)
  max_y <- max(beta, na.rm = TRUE)
  
  # Separate whole arm from segment features
  wholearm_indices <- grep('wholearm', names(beta))
  if (length(wholearm_indices) > 0) {
    beta_wholearm <- beta[wholearm_indices]
    beta_segments <- beta[-wholearm_indices]
  } else {
    beta_wholearm <- numeric(0)
    beta_segments <- beta
  }
  
  print(paste("Found", length(beta_wholearm), "whole-arm features"))
  print(paste("Found", length(beta_segments), "segment features"))
  
  # Create the plot
  png(filename = paste(main, '.png', sep = ''), width = 28, height = 5, res = 150, units = 'in')
  par(cex.axis = 2, cex.lab = 2.5, mai = c(0.6, 1.5, 0.6, 0.5))
  
  # Initialize empty plot
  plot(0, 0, type = 'n', xlim = c(0, total_gene), ylim = c(min_y * 1.1, max_y * 1.1), 
       xaxs = "i", xaxt = 'n', yaxt = 'n', ylab = '', xlab = "", bty = 'n')
  
  # Add whole-arm background shading first
  if (length(beta_wholearm) > 0) {
    for (i in 1:length(beta_wholearm)) {
      wholearm_name <- names(beta_wholearm)[i]
      wholearm_anno <- segment_anno[wholearm_name, ]
      
      # Determine if it's p or q arm and positive or negative
      arm_start <- wholearm_anno[2]
      arm_end <- wholearm_anno[3]
      weight_val <- beta_wholearm[i]
      
      if (weight_val > 0) {
        # Pink background for positive (gains)
        rect(arm_start, min_y * 1.1, arm_end, max_y * 1.1, 
             col = rgb(1, 0.8, 0.8, alpha = 0.7), border = NA)
      } else if (weight_val < 0) {
        # Light blue background for negative (losses)
        rect(arm_start, min_y * 1.1, arm_end, max_y * 1.1, 
             col = rgb(0.8, 0.8, 1, alpha = 0.7), border = NA)
      }
    }
  }
  
  # Create coordinate vectors for all segments
  all_x <- c()
  all_y <- c()
  all_colors <- c()
  
  # Process segment features
  if (length(beta_segments) > 0) {
    for (i in 1:length(beta_segments)) {
      seg_name <- names(beta_segments)[i]
      if (seg_name %in% rownames(segment_anno)) {
        seg_anno <- segment_anno[seg_name, ]
        seg_mid <- (seg_anno[2] + seg_anno[3]) / 2
        weight_val <- beta_segments[i]
        
        all_x <- c(all_x, seg_mid)
        all_y <- c(all_y, weight_val)
        
        if (weight_val > 0) {
          all_colors <- c(all_colors, 'red')
        } else {
          all_colors <- c(all_colors, 'blue')
        }
      }
    }
  }
  
  # Process whole-arm features (plot as thick lines)
  if (length(beta_wholearm) > 0) {
    for (i in 1:length(beta_wholearm)) {
      wholearm_name <- names(beta_wholearm)[i]
      wholearm_anno <- segment_anno[wholearm_name, ]
      wholearm_mid <- (wholearm_anno[2] + wholearm_anno[3]) / 2
      weight_val <- beta_wholearm[i]
      
      all_x <- c(all_x, wholearm_mid)
      all_y <- c(all_y, weight_val)
      
      if (weight_val > 0) {
        all_colors <- c(all_colors, 'red')
      } else {
        all_colors <- c(all_colors, 'blue')
      }
    }
  }
  
  # Plot all features as histogram-style bars
  if (length(all_x) > 0) {
    # Create histogram-style plot
    for (i in 1:length(all_x)) {
      lines(c(all_x[i], all_x[i]), c(0, all_y[i]), 
            col = all_colors[i], lwd = 3)
    }
  }
  
  # Add chromosome boundaries
  abline(v = vertical_lines, lwd = 1, col = 'gray60')
  abline(h = 0, lwd = 1.5, col = 'black')
  
  # Add y-axis
  y_ticks <- pretty(c(min_y, max_y), n = 5)
  axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = 2)
  
  # Calculate text positions for chromosome labels
  chrom_labels <- c(1:22, 'X')
  text_positions <- c()
  for (i in 1:(length(vertical_lines)-1)) {
    mid_pos <- (vertical_lines[i] + vertical_lines[i+1]) / 2
    text_positions <- c(text_positions, mid_pos)
  }
  
  # Ensure we have the right number of positions
  if (length(text_positions) > length(chrom_labels)) {
    text_positions <- text_positions[1:length(chrom_labels)]
  }
  
  # Add chromosome labels
  mtext(chrom_labels, side = 1, at = text_positions, line = 1.5, cex = 2.5)
  mtext('Weight', side = 2, line = 4, cex = 2.5)
  
  # Add title
  title(main = paste("Basal signaling signature", "17", sep = ""), cex.main = 2.5)
  
  dev.off()
  print(paste("Enhanced plot saved as:", paste(main, '.png', sep = '')))
}

# Plot with enhanced function
if (length(weights_filtered) > 0 && sum(weights_filtered != 0) > 0) {
  plot_basal_final(weights_filtered, "Basal_signaling_final")
} else {
  print("No valid weights to plot")
} 