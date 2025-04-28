# Set user library path
.libPaths(c(paste0(Sys.getenv("HOME"), "/R/library"), .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="http://cran.us.r-project.org", lib=paste0(Sys.getenv("HOME"), "/R/library"))
}

# Create directory if it doesn't exist
dir.create(paste0(Sys.getenv("HOME"), "/R/library"), recursive=TRUE, showWarnings=FALSE)

# Install required packages
BiocManager::install("maftools", lib=paste0(Sys.getenv("HOME"), "/R/library"))
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", lib=paste0(Sys.getenv("HOME"), "/R/library"))
install.packages("ggplot2", repos="http://cran.us.r-project.org", lib=paste0(Sys.getenv("HOME"), "/R/library")) 