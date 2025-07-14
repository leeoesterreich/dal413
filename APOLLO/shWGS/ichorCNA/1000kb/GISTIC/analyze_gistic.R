library(maftools)

# Path to GISTIC output directories
ilc_gistic_folder <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/GISTIC/gistic_run"
idc_gistic_folder <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/GISTIC/gistic_run_NST/"

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
dir.create("ILC_IDC_comparison", showWarnings = FALSE)

# Generate comparison plots for amplifications and deletions
pdf("ILC_IDC_comparison/ilc_idc_comparison_amplifications.pdf", width = 12, height = 8)
# Set larger text size for the plot
par(cex = 1.3)  # Increase overall text size
coGisticChromPlot(gistic1 = idc_gistic, gistic2 = ilc_gistic, 
                  g1Name = "IDC", g2Name = "ILC", 
                  type = 'Amp',
                  markBands = "all")
dev.off()

pdf("ILC_IDC_comparison/ilc_idc_comparison_deletions.pdf", width = 12, height = 8)
# Set larger text size for the plot
par(cex = 1.3)  # Increase overall text size
coGisticChromPlot(gistic1 = idc_gistic, gistic2 = ilc_gistic, 
                  g1Name = "IDC", g2Name = "ILC", 
                  type = 'Del',
                  markBands = "all")
dev.off()

# Reset par to default
par(cex = 1)

# Generate individual plots for ILC
pdf("ILC_IDC_comparison/ilc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = ilc_gistic, markBands = "all")
dev.off()

pdf("ILC_IDC_comparison/ilc_bubble_plot.pdf", width = 12, height = 8)
tryCatch({
  gisticBubblePlot(gistic = ilc_gistic, markBands = "all")
}, error = function(e) {
  cat("Could not generate ILC bubble plot: ", conditionMessage(e), "\n")
})
dev.off()

pdf("ILC_IDC_comparison/ilc_oncoplot_top20.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = ilc_gistic, top = 20)
dev.off()

# Generate individual plots for IDC
pdf("ILC_IDC_comparison/idc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = idc_gistic, markBands = "all")
dev.off()

pdf("ILC_IDC_comparison/idc_bubble_plot.pdf", width = 12, height = 8)
tryCatch({
  gisticBubblePlot(gistic = idc_gistic, markBands = "all")
}, error = function(e) {
  cat("Could not generate IDC bubble plot: ", conditionMessage(e), "\n")
})
dev.off()

pdf("ILC_IDC_comparison/idc_oncoplot_top20.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = idc_gistic, top = 20)
dev.off()

# Save summary statistics
sink("ILC_IDC_comparison/gistic_comparison_summary.txt")
cat("GISTIC Analysis Comparison Summary\n")
cat("================================\n\n")
cat("ILC GISTIC Results:\n")
cat("------------------\n")
print(ilc_gistic)
cat("\nIDC GISTIC Results:\n")
cat("------------------\n")
print(idc_gistic)
sink() 