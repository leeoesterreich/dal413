#!/bin/bash
#SBATCH --job-name=seg_cnv_analysis
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
module load python/3.7.0

echo "Starting CNV comparison analysis for 11q13.3"

# --- Define region of interest ---
# TODO: Please replace these with the correct coordinates for 11q13.3
export CHR="11"
export START_POS="69100000" # Placeholder
export END_POS="69500000"   # Placeholder

# --- File Paths ---
CLINICAL_FILE="brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt"
MASTER_SEG_FILE="brca_tcga.seg"
IDC_SAMPLES="idc_samples.txt"
ILC_SAMPLES="ilc_samples.txt"
IDC_SEG="idc.seg"
ILC_SEG="ilc.seg"
ANALYSIS_SCRIPT="analyze_region.py"
PLOT_FILE="11q13.3_IDC_vs_ILC_boxplot.png"

# --- Step 1: Create sample lists for IDC and ILC ---
echo "Creating sample lists..."
# Get header line to find sample ID column index
header=$(head -n 5 "$CLINICAL_FILE" | tail -n 1)
sample_col=$(echo "$header" | tr '\t' '\n' | grep -nx "SAMPLE_ID" | cut -d: -f1)
oncotree_col=$(echo "$header" | tr '\t' '\n' | grep -nx "ONCOTREE_CODE" | cut -d: -f1)

# Extract sample IDs based on Oncotree code
awk -F'\\t' -v smp_col="$sample_col" -v onc_col="$oncotree_col" 'NR > 5 && $onc_col == "IDC"' "$CLINICAL_FILE" | cut -f$smp_col > "$IDC_SAMPLES"
awk -F'\\t' -v smp_col="$sample_col" -v onc_col="$oncotree_col" 'NR > 5 && $onc_col == "ILC"' "$CLINICAL_FILE" | cut -f$smp_col > "$ILC_SAMPLES"

echo "Found $(wc -l < "$IDC_SAMPLES") IDC samples and $(wc -l < "$ILC_SAMPLES") ILC samples."

# --- Step 2: Filter the master segmentation file ---
echo "Filtering segmentation file for IDC and ILC samples..."
# Add header to new seg files
head -n 1 "$MASTER_SEG_FILE" > "$IDC_SEG"
head -n 1 "$MASTER_SEG_FILE" > "$ILC_SEG"

# Use grep to efficiently filter based on sample lists
grep -F -f "$IDC_SAMPLES" "$MASTER_SEG_FILE" >> "$IDC_SEG"
grep -F -f "$ILC_SAMPLES" "$MASTER_SEG_FILE" >> "$ILC_SEG"

# --- Step 3: Run the analysis using the standalone Python script ---
echo "Running Python analysis script..."
python analyze_region.py \
    --chr "$CHR" \
    --start "$START_POS" \
    --end "$END_POS" \
    --idc_seg_file "$IDC_SEG" \
    --ilc_seg_file "$ILC_SEG" \
    --plot_file "$PLOT_FILE"

echo "Analysis complete."

# --- Clean up intermediate files ---
rm "$IDC_SAMPLES" "$ILC_SAMPLES" "$IDC_SEG" "$ILC_SEG" 