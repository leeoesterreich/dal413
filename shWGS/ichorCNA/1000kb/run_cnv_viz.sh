#!/bin/bash
#SBATCH --job-name=cnv_viz
#SBATCH --output=cnv_viz_%j.out
#SBATCH --error=cnv_viz_%j.err
#SBATCH --time=2:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dal413@pitt.edu

# Load required modules
module load gcc/12.2.0
module load r/4.4.0

# Install required R packages if not already installed
R --quiet --no-save << EOF
if (!require("ggplot2")) install.packages("ggplot2", repos="http://cran.us.r-project.org")
if (!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org")
if (!require("tidyr")) install.packages("tidyr", repos="http://cran.us.r-project.org")
if (!require("ComplexHeatmap")) {
    if (!require("BiocManager")) install.packages("BiocManager", repos="http://cran.us.r-project.org")
    BiocManager::install("ComplexHeatmap")
}
if (!require("circlize")) install.packages("circlize", repos="http://cran.us.r-project.org")
if (!require("RColorBrewer")) install.packages("RColorBrewer", repos="http://cran.us.r-project.org")
EOF

# Run the R script
Rscript visualize_cnv.R 