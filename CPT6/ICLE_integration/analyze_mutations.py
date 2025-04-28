#!/usr/bin/env python3
# Script to analyze mutations across cell lines and convert mouse genes to human homologs

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# File paths
human_cell_lines_maf = "/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf"
human_ilc_maf = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/TCGA_ILC_BRCA_MAF.maf"
mouse_cpt6_maf = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf"
oncogene_list_file = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/Oncogene_CPT6.csv"
output_dir = "results"

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Check if input files exist
def check_input_files():
    missing_files = []
    for file_path in [human_cell_lines_maf, human_ilc_maf, mouse_cpt6_maf, oncogene_list_file]:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print("Error: The following input files are missing:")
        for file_path in missing_files:
            print("  - {}".format(file_path))
        sys.exit(1)

# Function to convert mouse genes to human homologs using a simple approach
def convert_mouse_to_human_genes(mouse_maf):
    print("Converting mouse genes to human homologs...")
    
    try:
        # Read the mouse MAF file
        mouse_df = pd.read_csv(mouse_maf, sep='\t', comment='#')
        
        # Get mouse gene symbols
        mouse_genes = mouse_df['Hugo_Symbol'].unique().tolist()
        
        # Simple conversion: capitalize first letter and keep the rest
        # This is a simplified approach since we don't have gseapy
        mouse_to_human = {}
        for gene in mouse_genes:
            try:
                # Most mouse genes are named the same as human genes but with different capitalization
                # Mouse: first letter capitalized, rest lowercase (e.g., Trp53)
                # Human: all uppercase (e.g., TP53)
                if gene.startswith('LOC') or gene.startswith('Gm'):
                    # Skip LOC and Gm genes which are often mouse-specific
                    mouse_to_human[gene] = gene
                else:
                    # For known important genes, provide manual mapping
                    if gene.lower() == 'trp53':
                        mouse_to_human[gene] = 'TP53'
                    elif gene.lower() == 'erbb2':
                        mouse_to_human[gene] = 'ERBB2'
                    elif gene.lower() == 'kras':
                        mouse_to_human[gene] = 'KRAS'
                    elif gene.lower() == 'braf':
                        mouse_to_human[gene] = 'BRAF'
                    elif gene.lower() == 'egfr':
                        mouse_to_human[gene] = 'EGFR'
                    else:
                        # Default: uppercase the gene name
                        mouse_to_human[gene] = gene.upper()
            except:
                mouse_to_human[gene] = gene  # Keep original if any error occurs
        
        # Create a new column with human homologs
        mouse_df['Human_Homolog'] = mouse_df['Hugo_Symbol'].map(mouse_to_human)
        
        # Replace mouse genes with human homologs
        mouse_df['Original_Hugo_Symbol'] = mouse_df['Hugo_Symbol']
        mouse_df['Hugo_Symbol'] = mouse_df['Human_Homolog']
        
        # Save the converted MAF file
        converted_maf = os.path.join(output_dir, "CPT6_human_homologs.maf")
        mouse_df.to_csv(converted_maf, sep='\t', index=False)
        
        print("Converted MAF file saved to {}".format(converted_maf))
        return converted_maf
    except Exception as e:
        print("Error converting mouse genes to human homologs: {}".format(str(e)))
        sys.exit(1)

# Function to merge MAF files
def merge_maf_files(human_cell_lines_maf, human_ilc_maf, converted_mouse_maf):
    print("Merging MAF files...")
    
    try:
        # Read the MAF files
        human_cell_lines_df = pd.read_csv(human_cell_lines_maf, sep='\t', comment='#')
        human_ilc_df = pd.read_csv(human_ilc_maf, sep='\t', comment='#')
        mouse_df = pd.read_csv(converted_mouse_maf, sep='\t', comment='#')
        
        # Add a source column to identify the origin of each mutation
        human_cell_lines_df['Source'] = 'ICLE_Cell_Lines'
        human_ilc_df['Source'] = 'TCGA_ILC'
        mouse_df['Source'] = 'CPT6'
        
        # Ensure all dataframes have the same columns
        common_columns = set(human_cell_lines_df.columns) & set(human_ilc_df.columns) & set(mouse_df.columns)
        
        if not common_columns:
            print("Error: No common columns found between MAF files.")
            print("Human cell lines columns:", human_cell_lines_df.columns.tolist())
            print("Human ILC columns:", human_ilc_df.columns.tolist())
            print("Mouse CPT6 columns:", mouse_df.columns.tolist())
            sys.exit(1)
        
        # Merge the dataframes
        merged_df = pd.concat([
            human_cell_lines_df[list(common_columns)],
            human_ilc_df[list(common_columns)],
            mouse_df[list(common_columns)]
        ], ignore_index=True)
        
        # Save the merged MAF file
        merged_maf = os.path.join(output_dir, "merged_maf.maf")
        merged_df.to_csv(merged_maf, sep='\t', index=False)
        
        print("Merged MAF file saved to {}".format(merged_maf))
        return merged_maf
    except Exception as e:
        print("Error merging MAF files: {}".format(str(e)))
        sys.exit(1)

# Function to analyze mutation frequencies of oncogenes
def analyze_oncogene_mutations(merged_maf, oncogene_list_file):
    print("Analyzing oncogene mutations...")
    
    try:
        # Read the oncogene list
        oncogenes = pd.read_csv(oncogene_list_file, header=None)[0].tolist()
        
        # Read the merged MAF file
        maf_df = pd.read_csv(merged_maf, sep='\t')
        
        # Filter for oncogenes
        oncogene_maf = maf_df[maf_df['Hugo_Symbol'].isin(oncogenes)]
        
        if oncogene_maf.empty:
            print("Warning: No mutations found for the specified oncogenes.")
            return
        
        # Count mutations per gene per source
        mutation_counts = oncogene_maf.groupby(['Hugo_Symbol', 'Source']).size().reset_index(name='Mutation_Count')
        
        # Pivot the data for better visualization
        pivot_counts = mutation_counts.pivot(index='Hugo_Symbol', columns='Source', values='Mutation_Count').fillna(0)
        
        # Save the mutation counts
        pivot_counts.to_csv(os.path.join(output_dir, "oncogene_mutation_counts.csv"))
        
        # Create a heatmap of mutation frequencies
        plt.figure(figsize=(12, len(oncogenes) * 0.4))
        sns.heatmap(pivot_counts, annot=True, cmap='YlOrRd', fmt='g')
        plt.title('Oncogene Mutation Frequencies Across Samples')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "oncogene_mutation_heatmap.png"), dpi=300)
        plt.close()
        
        print("Oncogene mutation analysis saved to {}".format(output_dir))
    except Exception as e:
        print("Error analyzing oncogene mutations: {}".format(str(e)))

# Function to analyze tumor mutation burden
def analyze_tumor_mutation_burden(human_cell_lines_maf, human_ilc_maf, mouse_cpt6_maf, converted_mouse_maf):
    print("Analyzing tumor mutation burden...")
    
    try:
        # Function to calculate TMB for a MAF file
        def calculate_tmb(maf_file, source_name):
            try:
                maf_df = pd.read_csv(maf_file, sep='\t', comment='#')
                # Group by sample and count mutations
                if 'Tumor_Sample_Barcode' not in maf_df.columns:
                    print("Warning: 'Tumor_Sample_Barcode' column not found in {}".format(maf_file))
                    return pd.DataFrame(columns=['Tumor_Sample_Barcode', 'Mutation_Count', 'Source'])
                
                sample_counts = maf_df.groupby('Tumor_Sample_Barcode').size().reset_index(name='Mutation_Count')
                sample_counts['Source'] = source_name
                return sample_counts
            except Exception as e:
                print("Error calculating TMB for {}: {}".format(source_name, str(e)))
                return pd.DataFrame(columns=['Tumor_Sample_Barcode', 'Mutation_Count', 'Source'])
        
        # Calculate TMB for each source
        human_cell_lines_tmb = calculate_tmb(human_cell_lines_maf, 'ICLE_Cell_Lines')
        human_ilc_tmb = calculate_tmb(human_ilc_maf, 'TCGA_ILC')
        mouse_tmb = calculate_tmb(mouse_cpt6_maf, 'CPT6_Mouse')
        converted_mouse_tmb = calculate_tmb(converted_mouse_maf, 'CPT6_Human_Homologs')
        
        # Combine TMB data
        combined_tmb = pd.concat([human_cell_lines_tmb, human_ilc_tmb, mouse_tmb, converted_mouse_tmb])
        
        if combined_tmb.empty:
            print("Warning: No TMB data available.")
            return
        
        # Save TMB data
        combined_tmb.to_csv(os.path.join(output_dir, "tumor_mutation_burden.csv"), index=False)
        
        # Create a boxplot of TMB by source
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='Source', y='Mutation_Count', data=combined_tmb)
        plt.title('Tumor Mutation Burden Across Sample Sources')
        plt.ylabel('Mutation Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "tmb_boxplot.png"), dpi=300)
        plt.close()
        
        print("Tumor mutation burden analysis saved to {}".format(output_dir))
    except Exception as e:
        print("Error analyzing tumor mutation burden: {}".format(str(e)))

# Main function to run the analysis
def main():
    try:
        # Check if input files exist
        check_input_files()
        
        # Convert mouse genes to human homologs
        converted_mouse_maf = convert_mouse_to_human_genes(mouse_cpt6_maf)
        
        # Merge MAF files
        merged_maf = merge_maf_files(human_cell_lines_maf, human_ilc_maf, converted_mouse_maf)
        
        # Analyze oncogene mutations
        analyze_oncogene_mutations(merged_maf, oncogene_list_file)
        
        # Analyze tumor mutation burden
        analyze_tumor_mutation_burden(human_cell_lines_maf, human_ilc_maf, mouse_cpt6_maf, converted_mouse_maf)
        
        print("Analysis complete. Results saved to the 'results' directory.")
    except Exception as e:
        print("Error in main function: {}".format(str(e)))
        sys.exit(1)

if __name__ == "__main__":
    main() 