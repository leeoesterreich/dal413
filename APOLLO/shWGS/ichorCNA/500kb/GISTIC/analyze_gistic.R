library(maftools)

# Path to your GISTIC output directory
gistic_res_folder <- "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/500kb/GISTIC/GISTIC/gistic_run"

# Read the GISTIC results
gistic_results <- readGistic(
  gisticDir = gistic_res_folder,
  isTCGA = FALSE
)

# Print summary of loaded data
print(gistic_results)

# Create output directory for plots
dir.create("gistic_plots", showWarnings = FALSE)

# Generate and save various visualizations
pdf("gistic_plots/chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = gistic_results, markBands = "all")
dev.off()

pdf("gistic_plots/bubble_plot.pdf", width = 12, height = 8)
gisticBubblePlot(gistic = gistic_results, markBands = "all")
dev.off()

# Generate separate plots for amplifications and deletions
pdf("gistic_plots/oncoplot_top20.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = gistic_results, top = 20)
dev.off()

# Save summary statistics
sink("gistic_plots/gistic_summary.txt")
cat("GISTIC Analysis Summary\n")
cat("======================\n\n")
print(gistic_results)
sink() 