library(maftools)

# Path to GISTIC output directories
ilc_gistic_folder <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb_Broad/GISTIC/gistic_run"
idc_gistic_folder <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Tendo/GISTIC"

# Read the GISTIC results for both ILC and IDC
ilc_gistic <- readGistic(
  gisticDir = ilc_gistic_folder,
  isTCGA = FALSE
)

idc_gistic <- readGistic(
  gisticDir = idc_gistic_folder,
  isTCGA = FALSE
)

# Print summary of loaded data
cat("ILC GISTIC Results:\n")
print(ilc_gistic)
cat("\nIDC GISTIC Results:\n")
print(idc_gistic)

# Create output directory for plots
dir.create("gistic_plots", showWarnings = FALSE)

# Generate comparison plots for amplifications and deletions
pdf("gistic_plots/ilc_idc_comparison_amplifications.pdf", width = 12, height = 8)
# Set larger text size for the plot
par(cex = 1.3)  # Increase overall text size
coGisticChromPlot(gistic1 = idc_gistic, gistic2 = ilc_gistic, 
                  g1Name = "IDC", g2Name = "ILC", 
                  type = 'Amp',
                  markBands = "all",
                  ref.build = "hg38")
dev.off()

pdf("gistic_plots/ilc_idc_comparison_deletions.pdf", width = 12, height = 8)
# Set larger text size for the plot
par(cex = 1.3)  # Increase overall text size
coGisticChromPlot(gistic1 = idc_gistic, gistic2 = ilc_gistic, 
                  g1Name = "IDC", g2Name = "ILC", 
                  type = 'Del',
                  markBands = "all",
                  ref.build = "hg38")
dev.off()

# Reset par to default
par(cex = 1)

# Generate individual plots for ILC with adjusted parameters
pdf("gistic_plots/ilc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = ilc_gistic, 
                markBands = "all", 
                ref.build = "hg38",
                amp.thresh = 0.1,  # Lower amplitude threshold to catch more subtle amplifications
                del.thresh = -0.1,  # Adjusted deletion threshold for consistency
                max.size = 2)  # Adjust point size for better visibility
dev.off()

pdf("gistic_plots/ilc_bubble_plot.pdf", width = 12, height = 8)
gisticBubblePlot(gistic = ilc_gistic, 
                 markBands = "all",
                 amp.thresh = 0.1,  # Lower amplitude threshold
                 del.thresh = -0.1)  # Adjusted deletion threshold
dev.off()

pdf("gistic_plots/ilc_oncoplot_top20.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = ilc_gistic, 
               top = 20,
               amp.thresh = 0.1,  # Lower amplitude threshold
               del.thresh = -0.1)  # Adjusted deletion threshold
dev.off()

# Generate individual plots for IDC with adjusted parameters
pdf("gistic_plots/idc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = idc_gistic, 
                markBands = "all",
                ref.build = "hg38",
                amp.thresh = 0.1,  # Lower amplitude threshold
                del.thresh = -0.1,  # Adjusted deletion threshold
                max.size = 2)  # Adjust point size for better visibility
dev.off()

pdf("gistic_plots/idc_bubble_plot.pdf", width = 12, height = 8)
gisticBubblePlot(gistic = idc_gistic, 
                 markBands = "all",
                 amp.thresh = 0.1,  # Lower amplitude threshold
                 del.thresh = -0.1)  # Adjusted deletion threshold
dev.off()

pdf("gistic_plots/idc_oncoplot_top20.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = idc_gistic, 
               top = 20,
               amp.thresh = 0.1,  # Lower amplitude threshold
               del.thresh = -0.1)  # Adjusted deletion threshold
dev.off()

# Save summary statistics
sink("gistic_plots/gistic_comparison_summary.txt")
cat("GISTIC Analysis Comparison Summary\n")
cat("================================\n\n")
cat("ILC GISTIC Results:\n")
cat("------------------\n")
print(ilc_gistic)
cat("\nIDC GISTIC Results:\n")
cat("------------------\n")
print(idc_gistic)
sink() 