import pandas as pd
import os
import sys

# Read CPT6 gene lists
def read_gene_list(file_path):
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"Successfully read {len(genes)} genes from {file_path}")
        return genes
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        sys.exit(1)

# Read gene lists
gain_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/gain.csv')
shallow_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/shallowdeletion.csv')
homo_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/homodeletion.csv')

# Combine all genes
all_genes = set(gain_genes + shallow_genes + homo_genes)
print(f"Total unique genes to filter: {len(all_genes)}")

# Read TCGA IDC gistic file
try:
    print("Reading TCGA IDC gistic file...")
    tcga_idc_df = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic.txt', sep='\t')
    print(f"Successfully read TCGA IDC data with {len(tcga_idc_df)} rows and {len(tcga_idc_df.columns)} columns")
except Exception as e:
    print(f"Error reading TCGA IDC gistic file: {str(e)}")
    sys.exit(1)

# Filter for genes in CPT6 lists
print("Filtering genes...")
filtered_df = tcga_idc_df[tcga_idc_df['Gene Symbol'].isin(all_genes)]

# Save filtered file
output_file = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic_filtered.txt'
try:
    print(f"Saving filtered file to {output_file}...")
    filtered_df.to_csv(output_file, sep='\t', index=False)
    print(f"Successfully saved filtered file with {len(filtered_df)} rows")
except Exception as e:
    print(f"Error saving filtered file: {str(e)}")
    sys.exit(1)

# Print summary of genes found/not found
found_genes = set(filtered_df['Gene Symbol'])
not_found = all_genes - found_genes
print(f"\nGenes not found in TCGA IDC data: {len(not_found)}")
for gene in not_found:
    print(f"- {gene}")

print("\nScript completed successfully!") 