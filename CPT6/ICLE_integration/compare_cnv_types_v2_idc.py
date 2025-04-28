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
parser = argparse.ArgumentParser(description='Compare CNV types between TCGA ILC, TCGA IDC, and ICLE datasets')
parser.add_argument('--icle_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/data_CNA_gistic.txt", 
                    help='Path to ICLE copy number alteration file')
parser.add_argument('--tcga_ilc_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic.txt", 
                    help='Path to TCGA ILC copy number alteration file')
parser.add_argument('--tcga_idc_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/IDC_data_CNA_gistic.txt", 
                    help='Path to TCGA IDC copy number alteration file')
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

# Define comparison bar colors using viridis colormap
viridis = plt.cm.viridis
comparison_colors = {
    'ICLE': viridis(0.2),
    'TCGA_ILC': viridis(0.5),
    'TCGA_IDC': viridis(0.8)
}

def process_icle_cnv(icle_cna_path):
    """Process ICLE CNV data."""
    icle_df = pd.read_csv(icle_cna_path, sep='\t')
    icle_samples = [col for col in icle_df.columns if col not in ['Hugo_Symbol', 'Entrez_Gene_Id']]
    icle_sample_count = len(icle_samples)
    print(f"Found {icle_sample_count} ICLE samples")
    return icle_df, icle_samples, icle_sample_count

def process_tcga_cnv(tcga_cna_path):
    """Process TCGA CNV data."""
    tcga_df = pd.read_csv(tcga_cna_path, sep='\t')
    tcga_samples = [col for col in tcga_df.columns if col not in ['Hugo_Symbol', 'Entrez_Gene_Id']]
    tcga_sample_count = len(tcga_samples)
    print(f"Found {tcga_sample_count} TCGA samples")
    return tcga_df, tcga_samples, tcga_sample_count

def calculate_cnv_frequencies(df, samples, sample_count):
    """Calculate CNV frequencies for each gene."""
    frequencies = {}
    for _, row in df.iterrows():
        gene = row['Hugo_Symbol']
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

def calculate_commonality_score(freq_dict1, freq_dict2):
    """Calculate commonality score between two frequency dictionaries."""
    common_genes = set(freq_dict1.keys()) & set(freq_dict2.keys())
    if not common_genes:
        return 0
    scores = []
    for gene in common_genes:
        score = 0
        for cnv_type in ['Shallow_Deletion', 'Homozygous_Deletion', 'Gain']:
            score += min(freq_dict1[gene][cnv_type], freq_dict2[gene][cnv_type])
        scores.append(score)
    return np.mean(scores)

def plot_cnv_comparison(icle_freq, tcga_ilc_freq, tcga_idc_freq, output_dir):
    """Create CNV comparison plot."""
    # Get common genes
    common_genes = set(icle_freq.keys()) & set(tcga_ilc_freq.keys()) & set(tcga_idc_freq.keys())
    print(f"Found {len(common_genes)} common genes")
    
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
    plt.title('CNV Comparison between ICLE, TCGA ILC, and TCGA IDC')
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
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig(os.path.join(output_dir, 'cnv_comparison_all_datasets.pdf'), bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(output_dir, 'cnv_comparison_all_datasets.png'), bbox_inches='tight', dpi=300)
    plt.close()

def main():
    # Process CNV data
    icle_df, icle_samples, icle_sample_count = process_icle_cnv(args.icle_cna)
    tcga_ilc_df, tcga_ilc_samples, tcga_ilc_sample_count = process_tcga_cnv(args.tcga_ilc_cna)
    tcga_idc_df, tcga_idc_samples, tcga_idc_sample_count = process_tcga_cnv(args.tcga_idc_cna)
    
    # Calculate frequencies
    icle_freq = calculate_cnv_frequencies(icle_df, icle_samples, icle_sample_count)
    tcga_ilc_freq = calculate_cnv_frequencies(tcga_ilc_df, tcga_ilc_samples, tcga_ilc_sample_count)
    tcga_idc_freq = calculate_cnv_frequencies(tcga_idc_df, tcga_idc_samples, tcga_idc_sample_count)
    
    # Create comparison plot
    plot_cnv_comparison(icle_freq, tcga_ilc_freq, tcga_idc_freq, args.output_dir)
    
    # Calculate and print commonality scores
    icle_tcga_ilc_score = calculate_commonality_score(icle_freq, tcga_ilc_freq)
    icle_tcga_idc_score = calculate_commonality_score(icle_freq, tcga_idc_freq)
    tcga_ilc_tcga_idc_score = calculate_commonality_score(tcga_ilc_freq, tcga_idc_freq)
    
    print(f"\nCommonality Scores:")
    print(f"ICLE vs TCGA ILC: {icle_tcga_ilc_score:.3f}")
    print(f"ICLE vs TCGA IDC: {icle_tcga_idc_score:.3f}")
    print(f"TCGA ILC vs TCGA IDC: {tcga_ilc_tcga_idc_score:.3f}")

if __name__ == "__main__":
    main() 