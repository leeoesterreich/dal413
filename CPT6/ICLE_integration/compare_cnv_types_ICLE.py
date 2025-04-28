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
    'Gain': '#FF0000',  # Red
    'No_CNV': '#F0F0F0'  # Light Gray
}

# Define comparison bar colors
comparison_colors = {
    'ICLE': '#ADD8E6',  # Light blue (same as Shallow_Deletion)
    'TCGA': '#0000FF'   # Dark blue (same as Homozygous_Deletion)
}

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
        tcga_sample_count = 179
        print(f"Using fixed TCGA sample size of {tcga_sample_count}")
        
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

def calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, tcga_samples, selected_genes, cpt6_cnv):
    """Calculate CNV frequencies and perform Fisher's exact test"""
    print("Calculating CNV frequencies and performing statistical tests")
    
    gene_stats = {}
    icle_total = len(icle_samples)
    tcga_total = 179
    
    # Calculate combined frequencies for sorting
    combined_freq = {}
    for gene in selected_genes:
        cpt6_type = cpt6_cnv.get(gene, 'No_CNV')
        if cpt6_type == 'No_CNV':
            continue
            
        # For homozygous deletion genes
        if cpt6_type == 'Homozygous_Deletion':
            icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
            icle_freq = (icle_count / icle_total) * 100
            tcga_freq = tcga_cnv_freq.get(gene, {}).get('Homozygous_Deletion', 0)
        # For shallow deletion genes
        elif cpt6_type == 'Shallow_Deletion':
            icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
            icle_freq = (icle_count / icle_total) * 100
            tcga_freq = tcga_cnv_freq.get(gene, {}).get('Shallow_Deletion', 0)
        # For gain genes
        elif cpt6_type == 'Gain':
            icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Gain')
            icle_freq = (icle_count / icle_total) * 100
            tcga_freq = tcga_cnv_freq.get(gene, {}).get('Gain', 0)
            
        # Store combined frequency for sorting
        combined_freq[gene] = {
            'type': cpt6_type,
            'combined_freq': icle_freq + tcga_freq
        }
    
    # Get top 10 genes for each CNV type (except homozygous deletion which we keep all)
    top_homdel_genes = [gene for gene, data in combined_freq.items() if data['type'] == 'Homozygous_Deletion']
    
    top_shalldel_genes = sorted(
        [gene for gene, data in combined_freq.items() if data['type'] == 'Shallow_Deletion'],
        key=lambda x: combined_freq[x]['combined_freq'],
        reverse=True
    )[:10]
    
    top_gain_genes = sorted(
        [gene for gene, data in combined_freq.items() if data['type'] == 'Gain'],
        key=lambda x: combined_freq[x]['combined_freq'],
        reverse=True
    )[:10]
    
    # Combine selected genes
    filtered_genes = top_homdel_genes + top_shalldel_genes + top_gain_genes
    
    # Calculate statistics for selected genes
    for gene in filtered_genes:
        # Get CPT6 CNV type
        cpt6_type = cpt6_cnv.get(gene, 'No_CNV')
        
        # Calculate ICLE frequency for the same CNV type
        if gene in icle_cnv:
            # For homozygous deletion genes
            if cpt6_type == 'Homozygous_Deletion':
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Homozygous_Deletion', 0)
            # For shallow deletion genes
            elif cpt6_type == 'Shallow_Deletion':
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Shallow_Deletion', 0)
            # For gain genes
            elif cpt6_type == 'Gain':
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Gain')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Gain', 0)
            
            # Perform Fisher's exact test
            pvalue = 1
            odds_ratio = 1
            if icle_count > 0 or tcga_freq > 0:
                tcga_count = int(np.round(tcga_freq * tcga_total / 100))
                table = [[icle_count, icle_total - icle_count],
                        [tcga_count, tcga_total - tcga_count]]
                odds_ratio, pvalue = fisher_exact(table)
            
            gene_stats[gene] = {
                'icle_freq': icle_freq,
                'tcga_freq': tcga_freq,
                'pvalue': pvalue,
                'odds_ratio': odds_ratio,
                'cnv_type': cpt6_type
            }
    
    return gene_stats, filtered_genes

def create_cnv_comparison_plot(cpt6_cnv, icle_cnv, icle_samples, tcga_samples, gene_stats, selected_genes):
    """Create CNV comparison plot"""
    print("Creating CNV comparison plot")
    
    # Sort genes by CNV type and frequency
    sorted_genes = sorted(
        [gene for gene in selected_genes if gene in gene_stats],
        key=lambda x: (
            0 if cpt6_cnv.get(x) == 'Homozygous_Deletion' else 
            (1 if cpt6_cnv.get(x) == 'Shallow_Deletion' else 2),
            -gene_stats.get(x, {}).get('tcga_freq', 0)
        )
    )
    
    # Sort ICLE samples by CNV count
    cnv_counts = {}
    for sample in icle_samples:
        count = sum(1 for gene in sorted_genes if icle_cnv[gene][sample] != 'No_CNV')
        cnv_counts[sample] = count
    
    sorted_icle_samples = sorted(icle_samples, key=lambda x: cnv_counts[x])
    
    # Create figure with gridspec layout
    fig_height = min(max(len(sorted_genes) * 0.3 + 4, 12), 24)
    fig_width = min(max(len(sorted_icle_samples) * 0.25 + 14, 24), 30)
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=100)
    
    # Create a gridspec to handle the layout - adjust width ratios for extended x-axis
    gs = plt.GridSpec(2, 20, height_ratios=[2, 15], width_ratios=[0.6]*10 + [0.4]*10, figure=fig, hspace=0.1)
    
    # Create plots
    ax_top = fig.add_subplot(gs[0, :10])  # Mutation counts
    ax_main = fig.add_subplot(gs[1, :10])  # Main oncoprint
    ax_comp = fig.add_subplot(gs[1, 10:])  # Comparison plot - wider now
    
    # Define box width for consistency
    box_width = 0.7
    
    # Calculate positions for CPT6 and ICLE cell lines
    cpt6_position = -1
    gap = 1  # Gap between CPT6 and ICLE
    icle_positions = np.arange(gap, gap + len(sorted_icle_samples))
    
    # Plot ICLE cell line CNV counts
    ax_top.bar(icle_positions, [cnv_counts[s] for s in sorted_icle_samples], 
              color='lightgray', alpha=0.5, width=box_width)
    
    # Add count labels
    for i, count in enumerate([cnv_counts[s] for s in sorted_icle_samples]):
        ax_top.text(icle_positions[i], count, str(count), ha='center', va='bottom')
    
    # Add "CNV Counts" annotation
    ax_top.text(-0.5, np.mean([cnv_counts[s] for s in sorted_icle_samples]), 
                "CNV\nCounts", ha='right', va='center', fontsize=10, fontweight='bold')
    
    # Add legend for CNV types
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.7)
                      for cnv_type, color in cnv_colors.items() if cnv_type != 'No_CNV']
    legend_labels = [k for k in cnv_colors.keys() if k != 'No_CNV']
    
    # Position the CNV type legend
    ax_top.legend(legend_elements, legend_labels,
                 loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=3, fontsize=8)
    
    # Set the same x-limits for top and main plots
    ax_top.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    
    # Plot settings
    bar_height = 0.6
    y_pos = np.arange(len(sorted_genes))
    
    # Add extra padding to y-axis
    y_padding = 0.5
    
    # Plot CPT6 CNVs
    for i, gene in enumerate(sorted_genes):
        if gene in cpt6_cnv:
            cnv_type = cpt6_cnv[gene]
            # Center the box around the cpt6_position
            box_start = cpt6_position - box_width/2
            ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                          box_width, bar_height,
                                          facecolor=cnv_colors.get(cnv_type, cnv_colors['No_CNV']),
                                          alpha=0.7))
    
    # Plot ICLE CNVs with background boxes
    for i, gene in enumerate(sorted_genes):
        # First plot gray background boxes for all samples
        for j, pos in enumerate(icle_positions):
            # Center the box around the position
            box_start = pos - box_width/2
            ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                          box_width, bar_height,
                                          facecolor=cnv_colors['No_CNV'],
                                          alpha=0.3))
        
        # Then plot CNV boxes
        for j, sample in enumerate(sorted_icle_samples):
            if gene in icle_cnv and sample in icle_cnv[gene]:
                cnv_type = icle_cnv[gene][sample]
                if cnv_type != 'No_CNV':
                    # Center the box around the position
                    box_start = icle_positions[j] - box_width/2
                    ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                                 box_width, bar_height,
                                                 facecolor=cnv_colors.get(cnv_type, cnv_colors['No_CNV']),
                                                 alpha=0.7))
    
    # Customize main plot
    ax_main.set_yticks(y_pos)
    ax_main.set_yticklabels(sorted_genes, fontsize=9)
    
    # Set y-axis limits with padding
    ax_main.set_ylim(-y_padding, len(sorted_genes) - 1 + y_padding)
    
    # Set x-axis limits and ticks
    ax_main.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    
    # Add all x-ticks at the exact positions
    all_tick_positions = [cpt6_position] + list(icle_positions)
    all_tick_labels = ['CPT6'] + [s.split('_')[0] for s in sorted_icle_samples]
    
    ax_main.set_xticks(all_tick_positions)
    ax_main.set_xticklabels(all_tick_labels, rotation=45, ha='right', va='top', fontsize=8)
    
    # Add x-axis label for ICLE cell lines
    ax_main.text(gap + len(sorted_icle_samples)/2, -len(sorted_genes)*0.15, 'ICLE Cell Lines', 
                ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Plot comparison bars
    bar_height = 0.35
    
    # Calculate maximum frequency for scaling
    max_freq = max(
        max(gene_stats.get(gene, {}).get('icle_freq', 0) for gene in sorted_genes),
        max(gene_stats.get(gene, {}).get('tcga_freq', 0) for gene in sorted_genes)
    )
    
    # Plot ICLE frequencies
    icle_bars = ax_comp.barh(y_pos + bar_height/2,
                            [gene_stats.get(gene, {}).get('icle_freq', 0) for gene in sorted_genes],
                            height=bar_height, color=comparison_colors['ICLE'], alpha=0.8,
                            label=f'ICLE (n={len(icle_samples)})')
    
    # Add ICLE frequency annotations
    for i, gene in enumerate(sorted_genes):
        if gene in gene_stats:
            freq = gene_stats[gene]['icle_freq']
            if freq > 0:
                # Calculate actual count
                count = int(round(freq * len(icle_samples) / 100))
                ax_comp.text(freq + 0.5, y_pos[i] + bar_height/2,
                            f'{count}/{len(icle_samples)}', va='center', ha='left', fontsize=7)
    
    # Plot TCGA frequencies
    tcga_bars = ax_comp.barh(y_pos - bar_height/2,
                            [gene_stats.get(gene, {}).get('tcga_freq', 0) for gene in sorted_genes],
                            height=bar_height, color=comparison_colors['TCGA'], alpha=0.8,
                            label=f'TCGA (n=179)')
    
    # Add TCGA frequency annotations (as percentages)
    for i, gene in enumerate(sorted_genes):
        if gene in gene_stats:
            freq = gene_stats[gene]['tcga_freq']
            if freq > 0:
                ax_comp.text(freq + 0.5, y_pos[i] - bar_height/2,
                            f'{freq:.1f}%', va='center', ha='left', fontsize=7)
    
    # Add significance stars - move more to the right
    for i, gene in enumerate(sorted_genes):
        if gene in gene_stats:
            stats = gene_stats[gene]
            if stats['pvalue'] < 0.05 and (stats['odds_ratio'] > 2 or stats['odds_ratio'] < 0.5):
                stars = '*' * sum([stats['pvalue'] < cutoff for cutoff in [0.05, 0.01, 0.001]])
                ax_comp.text(max(stats['icle_freq'], stats['tcga_freq']) + 5,  # Moved more to the right
                            i, stars, va='center', ha='left', fontsize=8)
    
    # Customize comparison plot
    max_freq_extended = 100  # Extend to 100%
    ax_comp.set_xlim(0, max_freq_extended)
    ax_comp.set_ylim(ax_main.get_ylim())
    ax_comp.set_xlabel('Alteration Frequency (%)', fontsize=10)
    ax_comp.grid(True, linestyle='--', alpha=0.3, axis='x')
    ax_comp.set_yticks([])
    ax_comp.spines['top'].set_visible(False)
    ax_comp.spines['right'].set_visible(False)
    
    # Add x-ticks at regular intervals
    ax_comp.set_xticks([0, 20, 40, 60, 80, 100])
    
    # Add legend with Fisher's exact test information
    legend_text = ('* p<0.05 & fold change >2\n(Fisher\'s exact test)')
    legend = ax_comp.legend(loc='upper right', bbox_to_anchor=(1, 1), 
                          frameon=True, fontsize=8)
    legend.set_title(legend_text, prop={'size': 8})
    
    # Add title indicating gene selection
    plt.suptitle('CNV Comparison of CPT6 Genes (Top 10 Common CNVs)', fontsize=14, y=0.98)
    
    # Adjust layout with more space
    plt.subplots_adjust(right=0.85, top=0.95, bottom=0.15)  # Added bottom margin
    
    # Save plot
    plt.savefig(os.path.join(args.output_dir, "cnv_comparison_top10.pdf"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.savefig(os.path.join(args.output_dir, "cnv_comparison_top10.png"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.close()

def print_cnv_summary(icle_cnv, icle_samples, selected_genes, cnv_counts, cpt6_cnv):
    """Print summary of CNVs"""
    print("\nCNV Summary:")
    
    print("\nCNVs per cell line:")
    for sample in icle_samples:
        print(f"{sample.split('_')[0]}: {cnv_counts[sample]} CNVs")
    
    print("\nSelected Genes CNV Status:")
    for gene in selected_genes:
        cpt6_type = cpt6_cnv.get(gene, 'No_CNV')
        if cpt6_type != 'No_CNV':
            cnv_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == cpt6_type)
            print(f"{gene} ({cpt6_type}): {cnv_samples}/{len(icle_samples)} samples with same CNV in ICLE")
    
    print("\nICLE CNVs by cell line:")
    for gene in selected_genes:
        if gene not in cpt6_cnv:
            continue
            
        cnv_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] != 'No_CNV')
        print(f"\n{gene} ({cpt6_cnv[gene]}): {cnv_samples}/{len(icle_samples)} samples")
        for sample in icle_samples:
            cnv_type = icle_cnv[gene][sample]
            if cnv_type != 'No_CNV':
                print(f"  {sample}: {cnv_type}")

def main():
    try:
        # Process CPT6 CNV data
        cpt6_homodeletion_genes, cpt6_shallowdeletion_genes, cpt6_gain_genes, cpt6_cnv = process_cpt6_cnv(
            args.cpt6_homodeletion, args.cpt6_shallowdeletion, args.cpt6_gain
        )
        
        # Use all gene types from CPT6
        selected_genes = cpt6_homodeletion_genes + cpt6_shallowdeletion_genes + cpt6_gain_genes
        
        # Remove duplicates while preserving order
        unique_selected_genes = []
        for gene in selected_genes:
            if gene not in unique_selected_genes:
                unique_selected_genes.append(gene)
        
        # Process TCGA CNV data for selected genes
        tcga_cnv_freq, tcga_samples = process_tcga_cnv(args.tcga_cna, unique_selected_genes)
        
        # Process ICLE CNV data for selected genes
        icle_cnv, icle_samples, gene_cnv_freq = process_icle_cnv(args.icle_cna, unique_selected_genes)
        
        if not icle_samples:
            print("Error: No ICLE samples found.")
            sys.exit(1)
        
        # Calculate CNV frequencies and perform statistical tests
        gene_stats, filtered_genes = calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, tcga_samples, 
                                             unique_selected_genes, cpt6_cnv)
        
        # Calculate CNV counts for each cell line
        cnv_counts = {}
        for sample in icle_samples:
            count = sum(1 for gene in filtered_genes if icle_cnv[gene][sample] != 'No_CNV')
            cnv_counts[sample] = count
        
        # Create CNV comparison plot
        create_cnv_comparison_plot(cpt6_cnv, icle_cnv, icle_samples, tcga_samples, gene_stats, filtered_genes)
        
        # Print CNV summary
        print_cnv_summary(icle_cnv, icle_samples, filtered_genes, cnv_counts, cpt6_cnv)
        
        print(f"\nAnalysis complete. Results saved to the '{args.output_dir}' directory.")
    except Exception as e:
        print(f"Error in main function: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 