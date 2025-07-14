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
dir.create("../gistic_plots", showWarnings = FALSE)

# Generate comparison plots for amplifications and deletions
pdf("../gistic_plots/ilc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = ilc_gistic, markBands = "all")
dev.off()

pdf("../gistic_plots/idc_chromosome_plot.pdf", width = 12, height = 8)
gisticChromPlot(gistic = idc_gistic, markBands = "all")
dev.off()

# Generate bubble plots
pdf("../gistic_plots/ilc_bubble_plot.pdf", width = 12, height = 8)
gisticBubblePlot(gistic = ilc_gistic)
dev.off()

pdf("../gistic_plots/idc_bubble_plot.pdf", width = 12, height = 8)
gisticBubblePlot(gistic = idc_gistic)
dev.off()

# Generate oncoplot-style visualization
pdf("../gistic_plots/ilc_oncoplot.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = ilc_gistic)
dev.off()

pdf("../gistic_plots/idc_oncoplot.pdf", width = 12, height = 8)
gisticOncoPlot(gistic = idc_gistic)
dev.off()

# Compare significant alterations between ILC and IDC
ilc_sig <- ilc_gistic@gis.scores
idc_sig <- idc_gistic@gis.scores

# Write comparison summary
sink("../gistic_plots/gistic_comparison_summary.txt")
cat("Significant alterations in ILC:\n")
print(ilc_sig[ilc_sig$q.value < 0.25,])
cat("\nSignificant alterations in IDC:\n")
print(idc_sig[idc_sig$q.value < 0.25,])
sink() 