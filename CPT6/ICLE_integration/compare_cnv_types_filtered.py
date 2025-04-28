import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse
from matplotlib import cm
from scipy.stats import fisher_exact

# Set up argument parser
parser = argparse.ArgumentParser(description='Compare CNV types between TCGA ILC, TCGA IDC, and ICLE datasets for CPT6 genes')
parser.add_argument('--icle_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/data_CNA_gistic.txt", 
                    help='Path to ICLE copy number alteration file')
parser.add_argument('--tcga_ilc_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic_filtered.txt", 
                    help='Path to filtered TCGA ILC copy number alteration file')
parser.add_argument('--tcga_idc_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic_filtered.txt", 
                    help='Path to filtered TCGA IDC copy number alteration file')
parser.add_argument('--output_dir', type=str, default="results_combined", 
                    help='Directory to save results')
args = parser.parse_args()

# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Define CNV types and colors
cnv_colors = {
    'Shallow_Deletion': '#ADD8E6',  # Light blue
    'Homozygous_Deletion': '#0000FF',  # Dark blue
    'Gain': '#FF0000',  # Red
    'No_CNV': '#F0F0F0'  # Light Gray
}

def process_cnv_data(cna_path):
    """Process CNV data from gistic file."""
    df = pd.read_csv(cna_path, sep='\t')
    samples = [col for col in df.columns if col not in ['Hugo_Symbol', 'Gene Symbol', 'Entrez_Gene_Id']]
    sample_count = len(samples)
    print(f"Found {sample_count} samples in {cna_path}")
    return df, samples, sample_count

def calculate_cnv_frequencies(df, samples, sample_count):
    """Calculate CNV frequencies for each gene."""
    frequencies = {}
    gene_col = 'Hugo_Symbol' if 'Hugo_Symbol' in df.columns else 'Gene Symbol'
    
    for _, row in df.iterrows():
        gene = row[gene_col]
        cnv_counts = {
            'Shallow_Deletion': 0,
            'Homozygous_Deletion': 0,
            'Gain': 0,
            'No_CNV': 0
        }
        for sample in samples:
            value = row[sample]
            if value == -1:
                cnv_counts['Shallow_Deletion'] += 1
            elif value == -2:
                cnv_counts['Homozygous_Deletion'] += 1
            elif value == 1:
                cnv_counts['Gain'] += 1
            else:
                cnv_counts['No_CNV'] += 1
        frequencies[gene] = {k: v/sample_count for k, v in cnv_counts.items()}
    return frequencies

def get_gene_categories(icle_freq):
    """Determine gene categories based on ICLE frequencies."""
    categories = {
        'Gain': [],
        'Shallow_Deletion': [],
        'Homozygous_Deletion': []
    }
    
    for gene, freqs in icle_freq.items():
        max_freq = max(freqs['Gain'], freqs['Shallow_Deletion'], freqs['Homozygous_Deletion'])
        if max_freq > 0:
            if freqs['Gain'] == max_freq:
                categories['Gain'].append(gene)
            elif freqs['Shallow_Deletion'] == max_freq:
                categories['Shallow_Deletion'].append(gene)
            elif freqs['Homozygous_Deletion'] == max_freq:
                categories['Homozygous_Deletion'].append(gene)
    
    return categories

def plot_cnv_comparison(icle_freq, tcga_ilc_freq, tcga_idc_freq, output_dir):
    """Create CNV comparison plot."""
    # Get common genes
    common_genes = set(icle_freq.keys()) & set(tcga_ilc_freq.keys()) & set(tcga_idc_freq.keys())
    print(f"Found {len(common_genes)} common genes")
    
    # Get gene categories
    categories = get_gene_categories(icle_freq)
    print("\nGene Categories:")
    print(f"Gain genes: {', '.join(categories['Gain'])}")
    print(f"Shallow deletion genes: {', '.join(categories['Shallow_Deletion'])}")
    print(f"Homozygous deletion genes: {', '.join(categories['Homozygous_Deletion'])}")
    
    # Sort genes by frequency in ICLE
    gene_freqs = [(gene, icle_freq[gene]['Gain'] + icle_freq[gene]['Shallow_Deletion'] + icle_freq[gene]['Homozygous_Deletion']) 
                  for gene in common_genes]
    gene_freqs.sort(key=lambda x: x[1], reverse=True)
    sorted_genes = [gene for gene, _ in gene_freqs]
    
    # Create figure
    plt.figure(figsize=(15, 8))
    
    # Plot bars for each dataset
    x = np.arange(len(sorted_genes))
    width = 0.25
    
    # Plot ICLE frequencies
    icle_gains = [icle_freq[gene]['Gain'] for gene in sorted_genes]
    icle_shallow = [icle_freq[gene]['Shallow_Deletion'] for gene in sorted_genes]
    icle_homo = [icle_freq[gene]['Homozygous_Deletion'] for gene in sorted_genes]
    
    # Plot TCGA ILC frequencies
    tcga_ilc_gains = [tcga_ilc_freq[gene]['Gain'] for gene in sorted_genes]
    tcga_ilc_shallow = [tcga_ilc_freq[gene]['Shallow_Deletion'] for gene in sorted_genes]
    tcga_ilc_homo = [tcga_ilc_freq[gene]['Homozygous_Deletion'] for gene in sorted_genes]
    
    # Plot TCGA IDC frequencies
    tcga_idc_gains = [tcga_idc_freq[gene]['Gain'] for gene in sorted_genes]
    tcga_idc_shallow = [tcga_idc_freq[gene]['Shallow_Deletion'] for gene in sorted_genes]
    tcga_idc_homo = [tcga_idc_freq[gene]['Homozygous_Deletion'] for gene in sorted_genes]
    
    # Plot stacked bars
    plt.bar(x - width, icle_gains, width, label='ICLE Gain', color=cnv_colors['Gain'])
    plt.bar(x - width, icle_shallow, width, bottom=icle_gains, label='ICLE Shallow Deletion', color=cnv_colors['Shallow_Deletion'])
    plt.bar(x - width, icle_homo, width, bottom=np.array(icle_gains) + np.array(icle_shallow), 
            label='ICLE Homozygous Deletion', color=cnv_colors['Homozygous_Deletion'])
    
    plt.bar(x, tcga_ilc_gains, width, label='TCGA ILC Gain', color=cnv_colors['Gain'])
    plt.bar(x, tcga_ilc_shallow, width, bottom=tcga_ilc_gains, label='TCGA ILC Shallow Deletion', color=cnv_colors['Shallow_Deletion'])
    plt.bar(x, tcga_ilc_homo, width, bottom=np.array(tcga_ilc_gains) + np.array(tcga_ilc_shallow), 
            label='TCGA ILC Homozygous Deletion', color=cnv_colors['Homozygous_Deletion'])
    
    plt.bar(x + width, tcga_idc_gains, width, label='TCGA IDC Gain', color=cnv_colors['Gain'])
    plt.bar(x + width, tcga_idc_shallow, width, bottom=tcga_idc_gains, label='TCGA IDC Shallow Deletion', color=cnv_colors['Shallow_Deletion'])
    plt.bar(x + width, tcga_idc_homo, width, bottom=np.array(tcga_idc_gains) + np.array(tcga_idc_shallow), 
            label='TCGA IDC Homozygous Deletion', color=cnv_colors['Homozygous_Deletion'])
    
    # Customize plot
    plt.xlabel('Genes')
    plt.ylabel('Frequency')
    plt.title('CNV Comparison between ICLE, TCGA ILC, and TCGA IDC (CPT6 Genes)')
    plt.xticks(x, sorted_genes, rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add sample size annotations
    plt.text(0.02, 0.98, f'ICLE (n={len(icle_samples)})', transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.text(0.02, 0.92, f'TCGA ILC (n={len(tcga_ilc_samples)})', transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.text(0.02, 0.86, f'TCGA IDC (n={len(tcga_idc_samples)})', transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add gene category annotations
    y_pos = 0.80
    for category, genes in categories.items():
        plt.text(0.02, y_pos, f'{category} genes:', transform=plt.gca().transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        y_pos -= 0.04
        for gene in genes:
            plt.text(0.02, y_pos, f'  {gene}', transform=plt.gca().transAxes, 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            y_pos -= 0.04
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig(os.path.join(output_dir, 'cnv_comparison_filtered.pdf'), bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(output_dir, 'cnv_comparison_filtered.png'), bbox_inches='tight', dpi=300)
    plt.close()

def main():
    # Process CNV data
    icle_df, icle_samples, icle_sample_count = process_cnv_data(args.icle_cna)
    tcga_ilc_df, tcga_ilc_samples, tcga_ilc_sample_count = process_cnv_data(args.tcga_ilc_cna)
    tcga_idc_df, tcga_idc_samples, tcga_idc_sample_count = process_cnv_data(args.tcga_idc_cna)
    
    # Calculate frequencies
    icle_freq = calculate_cnv_frequencies(icle_df, icle_samples, icle_sample_count)
    tcga_ilc_freq = calculate_cnv_frequencies(tcga_ilc_df, tcga_ilc_samples, tcga_ilc_sample_count)
    tcga_idc_freq = calculate_cnv_frequencies(tcga_idc_df, tcga_idc_samples, tcga_idc_sample_count)
    
    # Create comparison plot
    plot_cnv_comparison(icle_freq, tcga_ilc_freq, tcga_idc_freq, args.output_dir)
    
    print("\nScript completed successfully!")

if __name__ == "__main__":
    main() 