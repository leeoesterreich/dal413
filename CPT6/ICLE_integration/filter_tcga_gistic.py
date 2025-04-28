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

def filter_gistic_file(input_file, output_file, genes_to_filter):
    try:
        print(f"\nProcessing {input_file}...")
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully read data with {len(df)} rows and {len(df.columns)} columns")
        
        filtered_df = df[df['Gene Symbol'].isin(genes_to_filter)]
        print(f"Filtering genes...")
        
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"Successfully saved filtered file to {output_file} with {len(filtered_df)} rows")
        
        found_genes = set(filtered_df['Gene Symbol'])
        not_found = genes_to_filter - found_genes
        print(f"Genes not found: {len(not_found)}")
        for gene in not_found:
            print(f"- {gene}")
            
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        sys.exit(1)

# Read gene lists
gain_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/gain.csv')
shallow_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/shallowdeletion.csv')
homo_genes = read_gene_list('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/homodeletion.csv')

# Combine all genes
all_genes = set(gain_genes + shallow_genes + homo_genes)
print(f"\nTotal unique genes to filter: {len(all_genes)}")

# Process TCGA ILC data
filter_gistic_file(
    '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic.txt',
    '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic_filtered.txt',
    all_genes
)

# Process TCGA IDC data
filter_gistic_file(
    '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic.txt',
    '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic_filtered.txt',
    all_genes
)

print("\nScript completed successfully!") 