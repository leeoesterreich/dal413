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
from collections import Counter

# Set up argument parser
parser = argparse.ArgumentParser(description='Compare CNV types between TCGA, CPT6, and ICLE datasets')
parser.add_argument('--cpt6_homodeletion', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/homodeletion.csv", 
                    help='Path to CPT6 homodeletion gene list')
parser.add_argument('--cpt6_shallowdeletion', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/shallowdeletion.csv", 
                    help='Path to CPT6 shallow deletion gene list')
parser.add_argument('--cpt6_gain', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/CNV_summary/gain.csv", 
                    help='Path to CPT6 gain gene list')
parser.add_argument('--icle_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/data_CNA_gistic.txt", 
                    help='Path to ICLE copy number alteration file')
parser.add_argument('--tcga_cna', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic.txt", 
                    help='Path to TCGA ILC copy number alteration file')
parser.add_argument('--output_dir', type=str, default="results", 
                    help='Directory to save results')
args = parser.parse_args()

# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Define CNV types and colors
cnv_colors = {
    'Shallow_Deletion': '#ADD8E6',  # Light blue
    'Homozygous_Deletion': '#0000FF',  # Dark blue
    'Any_Deletion': '#6495ED',      # Medium blue
    'Gain': '#FF0000',  # Red
    'No_CNV': '#F0F0F0'  # Light Gray
}

# Define comparison bar colors
comparison_colors = {
    'ICLE': '#ADD8E6',  # Light blue (same as Shallow_Deletion)
    'TCGA': '#0000FF'   # Dark blue (same as Homozygous_Deletion)
}

# Define the specific genes to focus on
TARGET_GENES = ['CDKN2A', 'CDKN2B', 'CDH1', 'PTEN']

def read_gene_list(file_path):
    """Read gene list from a CSV file"""
    try:
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        return genes
    except Exception as e:
        print(f"Error reading gene list from {file_path}: {str(e)}")
        return []

def process_cpt6_cnv(homodeletion_path, shallowdeletion_path, gain_path):
    """Process CPT6 CNV data from gene lists"""
    print(f"Processing CPT6 CNV data from gene lists")
    
    # Read gene lists
    homodeletion_genes = read_gene_list(homodeletion_path)
    shallowdeletion_genes = read_gene_list(shallowdeletion_path)
    gain_genes = read_gene_list(gain_path)
    
    # Create CNV dictionary for CPT6
    cpt6_cnv = {}
    for gene in homodeletion_genes:
        cpt6_cnv[gene] = 'Homozygous_Deletion'
    for gene in shallowdeletion_genes:
        cpt6_cnv[gene] = 'Shallow_Deletion'
    for gene in gain_genes:
        cpt6_cnv[gene] = 'Gain'
    
    return homodeletion_genes, shallowdeletion_genes, gain_genes, cpt6_cnv

def process_tcga_cnv(tcga_cna_path, selected_genes):
    """Process TCGA CNV data for selected genes"""
    print(f"Processing TCGA CNV data from {tcga_cna_path}")
    
    try:
        # Read TCGA CNA file
        tcga_df = pd.read_csv(tcga_cna_path, sep='\t')
        
        # Get sample names
        tcga_samples = [col for col in tcga_df.columns if col != 'Gene Symbol']
        tcga_sample_count = len(tcga_samples)
        print(f"Found {tcga_sample_count} TCGA samples")
        
        # Initialize CNV frequency dictionary
        tcga_cnv_freq = {}
        
        # Process each selected gene
        for gene in selected_genes:
            if gene in tcga_df['Gene Symbol'].values:
                gene_row = tcga_df[tcga_df['Gene Symbol'] == gene]
                
                # Count different CNV types for this gene
                homdel_count = sum(1 for sample in tcga_samples if gene_row[sample].values[0] == -2)
                shalldel_count = sum(1 for sample in tcga_samples if gene_row[sample].values[0] == -1)
                gain_count = sum(1 for sample in tcga_samples if gene_row[sample].values[0] >= 1)
                
                # Calculate frequencies
                tcga_cnv_freq[gene] = {
                    'Homozygous_Deletion': homdel_count / tcga_sample_count * 100,
                    'Shallow_Deletion': shalldel_count / tcga_sample_count * 100,
                    'Gain': gain_count / tcga_sample_count * 100,
                    'Any_Deletion': (homdel_count + shalldel_count) / tcga_sample_count * 100
                }
            else:
                # If gene not found, initialize with zeros
                tcga_cnv_freq[gene] = {
                    'Homozygous_Deletion': 0,
                    'Shallow_Deletion': 0,
                    'Gain': 0,
                    'Any_Deletion': 0
                }
        
        return tcga_cnv_freq, tcga_samples
    
    except Exception as e:
        print(f"Error processing TCGA CNV data: {str(e)}")
        return {}, []

def process_icle_cnv(icle_cna_path, selected_genes):
    """Process ICLE CNV data for selected genes"""
    print(f"Processing ICLE CNV data from {icle_cna_path}")
    
    try:
        # Read ICLE CNA file
        icle_df = pd.read_csv(icle_cna_path, sep='\t')
        
        # Get sample names
        icle_samples = [col for col in icle_df.columns if col not in ['Hugo_Symbol', 'Entrez_Gene_Id']]
        
        # Initialize CNV dictionary
        icle_cnv = {}
        
        # Process selected genes
        for gene in selected_genes:
            if gene in icle_df['Hugo_Symbol'].values:
                gene_row = icle_df[icle_df['Hugo_Symbol'] == gene]
                icle_cnv[gene] = {}
                
                for sample in icle_samples:
                    cnv_value = gene_row[sample].values[0]
                    
                    if cnv_value == -2:
                        icle_cnv[gene][sample] = 'Homozygous_Deletion'
                    elif cnv_value == -1:
                        icle_cnv[gene][sample] = 'Shallow_Deletion'
                    elif cnv_value >= 1:
                        icle_cnv[gene][sample] = 'Gain'
                    else:
                        icle_cnv[gene][sample] = 'No_CNV'
            else:
                # If gene not found, initialize with No_CNV for all samples
                icle_cnv[gene] = {sample: 'No_CNV' for sample in icle_samples}
        
        # Calculate frequency of CNV types for each gene
        gene_cnv_freq = {}
        for gene in selected_genes:
            if gene in icle_cnv:
                homdel_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
                shalldel_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
                gain_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Gain')
                
                gene_cnv_freq[gene] = {
                    'Homozygous_Deletion': homdel_count / len(icle_samples) * 100,
                    'Shallow_Deletion': shalldel_count / len(icle_samples) * 100,
                    'Gain': gain_count / len(icle_samples) * 100,
                    'Any_Deletion': (homdel_count + shalldel_count) / len(icle_samples) * 100
                }
        
        return icle_cnv, icle_samples, gene_cnv_freq
    
    except Exception as e:
        print(f"Error processing ICLE CNV data: {str(e)}")
        return {}, [], {}

def calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, tcga_samples, target_genes):
    """Calculate CNV frequencies and perform Fisher's exact test for target genes"""
    print("Calculating CNV frequencies and performing statistical tests")
    
    gene_stats = {}
    icle_total = len(icle_samples)
    tcga_total = len(tcga_samples)
    
    # Calculate statistics for both homozygous and any deletion
    for gene in target_genes:
        # Calculate for homozygous deletion
        homdel_icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
        homdel_icle_freq = (homdel_icle_count / icle_total) * 100
        homdel_tcga_freq = tcga_cnv_freq.get(gene, {}).get('Homozygous_Deletion', 0)
        
        # Calculate for any deletion (homozygous + shallow)
        anydel_icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] in ['Homozygous_Deletion', 'Shallow_Deletion'])
        anydel_icle_freq = (anydel_icle_count / icle_total) * 100
        anydel_tcga_freq = tcga_cnv_freq.get(gene, {}).get('Any_Deletion', 0)
        
        # Perform Fisher's exact test for homozygous deletion
        homdel_pvalue = 1
        homdel_odds_ratio = 1
        if homdel_icle_count > 0 or homdel_tcga_freq > 0:
            homdel_tcga_count = int(np.round(homdel_tcga_freq * tcga_total / 100))
            homdel_table = [[homdel_icle_count, icle_total - homdel_icle_count],
                    [homdel_tcga_count, tcga_total - homdel_tcga_count]]
            homdel_odds_ratio, homdel_pvalue = fisher_exact(homdel_table)
        
        # Perform Fisher's exact test for any deletion
        anydel_pvalue = 1
        anydel_odds_ratio = 1
        if anydel_icle_count > 0 or anydel_tcga_freq > 0:
            anydel_tcga_count = int(np.round(anydel_tcga_freq * tcga_total / 100))
            anydel_table = [[anydel_icle_count, icle_total - anydel_icle_count],
                    [anydel_tcga_count, tcga_total - anydel_tcga_count]]
            anydel_odds_ratio, anydel_pvalue = fisher_exact(anydel_table)
        
        # Store statistics for both types
        gene_stats[f"{gene}_homdel"] = {
            'gene': gene,
            'icle_freq': homdel_icle_freq,
            'tcga_freq': homdel_tcga_freq,
            'pvalue': homdel_pvalue,
            'odds_ratio': homdel_odds_ratio,
            'cnv_type': 'Homozygous_Deletion',
            'display_name': f"{gene} (HomDel)"
        }
        
        gene_stats[f"{gene}_anydel"] = {
            'gene': gene,
            'icle_freq': anydel_icle_freq,
            'tcga_freq': anydel_tcga_freq,
            'pvalue': anydel_pvalue,
            'odds_ratio': anydel_odds_ratio,
            'cnv_type': 'Any_Deletion',
            'display_name': f"{gene} (Any Del)"
        }
    
    return gene_stats

def create_combined_deletion_plot(icle_cnv, icle_samples, tcga_samples, gene_stats, target_genes):
    """Create plot comparing homozygous and any deletion for target genes"""
    print("Creating combined deletion comparison plot")
    
    # Create a list of entries to display
    display_entries = []
    for gene in target_genes:
        display_entries.append(f"{gene}_homdel")
        display_entries.append(f"{gene}_anydel")
    
    # Sort entries by gene name
    sorted_entries = sorted(display_entries, key=lambda x: (target_genes.index(x.split('_')[0]), x.split('_')[1]))
    
    # Create figure
    fig_height = 8
    fig_width = 12
    fig, ax = plt.figure(figsize=(fig_width, fig_height), dpi=100), plt.axes()
    
    # Set up positions
    y_pos = np.arange(len(sorted_entries))
    bar_height = 0.35
    
    # Plot ICLE frequencies
    icle_bars = ax.barh(y_pos + bar_height/2,
                        [gene_stats.get(entry, {}).get('icle_freq', 0) for entry in sorted_entries],
                        height=bar_height, color=comparison_colors['ICLE'], alpha=0.8,
                        label=f'ICLE (n={len(icle_samples)})')
    
    # Add ICLE frequency annotations
    for i, entry in enumerate(sorted_entries):
        if entry in gene_stats:
            freq = gene_stats[entry]['icle_freq']
            if freq > 0:
                # Calculate actual count
                count = int(round(freq * len(icle_samples) / 100))
                ax.text(freq + 0.5, y_pos[i] + bar_height/2,
                        f'{count}/{len(icle_samples)}', va='center', ha='left', fontsize=8)
    
    # Plot TCGA frequencies
    tcga_bars = ax.barh(y_pos - bar_height/2,
                        [gene_stats.get(entry, {}).get('tcga_freq', 0) for entry in sorted_entries],
                        height=bar_height, color=comparison_colors['TCGA'], alpha=0.8,
                        label=f'TCGA (n={len(tcga_samples)})')
    
    # Add TCGA frequency annotations (as percentages)
    for i, entry in enumerate(sorted_entries):
        if entry in gene_stats:
            freq = gene_stats[entry]['tcga_freq']
            if freq > 0:
                ax.text(freq + 0.5, y_pos[i] - bar_height/2,
                        f'{freq:.1f}%', va='center', ha='left', fontsize=8)
    
    # Add significance stars
    for i, entry in enumerate(sorted_entries):
        if entry in gene_stats:
            stats = gene_stats[entry]
            if stats['pvalue'] < 0.05 and (stats['odds_ratio'] > 2 or stats['odds_ratio'] < 0.5):
                stars = '*' * sum([stats['pvalue'] < cutoff for cutoff in [0.05, 0.01, 0.001]])
                ax.text(max(stats['icle_freq'], stats['tcga_freq']) + 5,
                        i, stars, va='center', ha='left', fontsize=10)
    
    # Customize plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels([gene_stats.get(entry, {}).get('display_name', entry) for entry in sorted_entries], fontsize=10)
    
    # Set axis limits
    ax.set_xlim(0, 100)
    ax.set_xlabel('Alteration Frequency (%)', fontsize=12)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, axis='x')
    
    # Add legend with Fisher's exact test information
    legend_text = ('* p<0.05 & fold change >2\n(Fisher\'s exact test)')
    legend = ax.legend(loc='upper right', bbox_to_anchor=(1, 1), 
                      frameon=True, fontsize=10)
    legend.set_title(legend_text, prop={'size': 10})
    
    # Add title
    plt.title('Comparison of Homozygous and Any Deletion Frequencies for Key Genes', fontsize=14, pad=20)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig(os.path.join(args.output_dir, "combined_deletion_comparison.pdf"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.savefig(os.path.join(args.output_dir, "combined_deletion_comparison.png"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.close()

def print_combined_deletion_summary(icle_cnv, icle_samples, gene_stats, target_genes):
    """Print summary of combined deletion analysis"""
    print("\nCombined Deletion Analysis Summary:")
    
    for gene in target_genes:
        print(f"\n{gene}:")
        
        # Homozygous deletion stats
        homdel_entry = f"{gene}_homdel"
        if homdel_entry in gene_stats:
            homdel_stats = gene_stats[homdel_entry]
            homdel_icle_count = int(round(homdel_stats['icle_freq'] * len(icle_samples) / 100))
            print(f"  Homozygous Deletion: {homdel_icle_count}/{len(icle_samples)} ICLE samples ({homdel_stats['icle_freq']:.1f}%)")
            print(f"  Homozygous Deletion: {homdel_stats['tcga_freq']:.1f}% TCGA samples")
            print(f"  p-value: {homdel_stats['pvalue']:.4f}, Odds Ratio: {homdel_stats['odds_ratio']:.2f}")
        
        # Any deletion stats
        anydel_entry = f"{gene}_anydel"
        if anydel_entry in gene_stats:
            anydel_stats = gene_stats[anydel_entry]
            anydel_icle_count = int(round(anydel_stats['icle_freq'] * len(icle_samples) / 100))
            print(f"  Any Deletion: {anydel_icle_count}/{len(icle_samples)} ICLE samples ({anydel_stats['icle_freq']:.1f}%)")
            print(f"  Any Deletion: {anydel_stats['tcga_freq']:.1f}% TCGA samples")
            print(f"  p-value: {anydel_stats['pvalue']:.4f}, Odds Ratio: {anydel_stats['odds_ratio']:.2f}")
        
        # List ICLE samples with deletions
        print("\n  ICLE samples with deletions:")
        for sample in icle_samples:
            cnv_type = icle_cnv[gene][sample]
            if cnv_type in ['Homozygous_Deletion', 'Shallow_Deletion']:
                print(f"    {sample}: {cnv_type}")

def main():
    try:
        # Process CPT6 CNV data
        cpt6_homodeletion_genes, cpt6_shallowdeletion_genes, cpt6_gain_genes, cpt6_cnv = process_cpt6_cnv(
            args.cpt6_homodeletion, args.cpt6_shallowdeletion, args.cpt6_gain
        )
        
        # Process TCGA CNV data for target genes
        tcga_cnv_freq, tcga_samples = process_tcga_cnv(args.tcga_cna, TARGET_GENES)
        
        # Process ICLE CNV data for target genes
        icle_cnv, icle_samples, gene_cnv_freq = process_icle_cnv(args.icle_cna, TARGET_GENES)
        
        if not icle_samples:
            print("Error: No ICLE samples found.")
            sys.exit(1)
        
        # Calculate CNV frequencies and perform statistical tests
        gene_stats = calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, tcga_samples, TARGET_GENES)
        
        # Create combined deletion comparison plot
        create_combined_deletion_plot(icle_cnv, icle_samples, tcga_samples, gene_stats, TARGET_GENES)
        
        # Print combined deletion summary
        print_combined_deletion_summary(icle_cnv, icle_samples, gene_stats, TARGET_GENES)
        
        print(f"\nAnalysis complete. Results saved to the '{args.output_dir}' directory.")
    except Exception as e:
        print(f"Error in main function: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 