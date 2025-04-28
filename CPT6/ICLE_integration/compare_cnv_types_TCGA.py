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

def process_tcga_cnv(tcga_cna_path):
    """Process TCGA CNV data and select top genes"""
    print(f"Processing TCGA CNV data from {tcga_cna_path}")
    
    try:
        # Read TCGA CNA file
        tcga_df = pd.read_csv(tcga_cna_path, sep='\t')
        
        # Get sample names
        tcga_samples = [col for col in tcga_df.columns if col != 'Gene Symbol']
        tcga_sample_count = 179
        print(f"Using fixed TCGA sample size of {tcga_sample_count}")
        
        # Initialize CNV frequency dictionaries
        tcga_cnv_freq = {}
        
        # Process each gene
        for _, row in tcga_df.iterrows():
            gene = row['Gene Symbol']
            
            # Count different CNV types for this gene
            homdel_count = sum(1 for sample in tcga_samples if row[sample] == -2)
            shalldel_count = sum(1 for sample in tcga_samples if row[sample] == -1)
            gain_count = sum(1 for sample in tcga_samples if row[sample] >= 1)
            
            # Calculate frequencies
            if gene not in tcga_cnv_freq:
                tcga_cnv_freq[gene] = {
                    'Homozygous_Deletion': 0,
                    'Shallow_Deletion': 0,
                    'Gain': 0
                }
            
            # Update frequencies
            tcga_cnv_freq[gene]['Homozygous_Deletion'] = homdel_count / tcga_sample_count * 100
            tcga_cnv_freq[gene]['Shallow_Deletion'] = shalldel_count / tcga_sample_count * 100
            tcga_cnv_freq[gene]['Gain'] = gain_count / tcga_sample_count * 100
        
        # Select top 10 homozygously deleted genes
        homdel_genes = sorted(
            [(gene, data['Homozygous_Deletion']) for gene, data in tcga_cnv_freq.items() if data['Homozygous_Deletion'] > 0],
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        # Select top 10 shallow deleted genes
        shalldel_genes = sorted(
            [(gene, data['Shallow_Deletion']) for gene, data in tcga_cnv_freq.items() if data['Shallow_Deletion'] > 0],
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        # Select top 10 gained genes
        gained_genes = sorted(
            [(gene, data['Gain']) for gene, data in tcga_cnv_freq.items() if data['Gain'] > 0],
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        # Extract gene names
        top_homdel_genes = [gene for gene, _ in homdel_genes]
        top_shalldel_genes = [gene for gene, _ in shalldel_genes]
        top_gained_genes = [gene for gene, _ in gained_genes]
        
        # Combine selected genes
        selected_genes = top_homdel_genes + top_shalldel_genes + top_gained_genes
        
        return tcga_cnv_freq, selected_genes, top_homdel_genes + top_shalldel_genes, top_gained_genes
    
    except Exception as e:
        print(f"Error processing TCGA CNV data: {str(e)}")
        return {}, [], [], []

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

def calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, selected_genes, top_deleted_genes, top_gained_genes):
    """Calculate CNV frequencies and perform Fisher's exact test"""
    print("Calculating CNV frequencies and performing statistical tests")
    
    gene_stats = {}
    icle_total = len(icle_samples)
    tcga_total = 179  # TCGA sample count
    
    for gene in selected_genes:
        # Calculate ICLE frequencies
        if gene in icle_cnv:
            # For homozygous deletion genes, compare with TCGA homozygous deletions
            if gene in top_deleted_genes[:4]:  # First 4 are homozygous deletions
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Homozygous_Deletion', 0)
                cnv_type = 'Homozygous_Deletion'
            # For shallow deletion genes, compare with TCGA shallow deletions
            elif gene in top_deleted_genes[4:]:  # Rest are shallow deletions
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Shallow_Deletion', 0)
                cnv_type = 'Shallow_Deletion'
            # For gain genes, count only gains
            elif gene in top_gained_genes:
                icle_count = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Gain')
                icle_freq = (icle_count / icle_total) * 100
                tcga_freq = tcga_cnv_freq.get(gene, {}).get('Gain', 0)
                cnv_type = 'Gain'
            else:
                continue
            
            # Perform Fisher's exact test
            pvalue = 1
            odds_ratio = 1
            if icle_count > 0 and tcga_freq > 0:
                tcga_count = int(np.round(tcga_freq * tcga_total / 100))
                table = [[icle_count, icle_total - icle_count],
                        [tcga_count, tcga_total - tcga_count]]
                odds_ratio, pvalue = fisher_exact(table)
            
            gene_stats[gene] = {
                'icle_freq': icle_freq,
                'tcga_freq': tcga_freq,
                'pvalue': pvalue,
                'odds_ratio': odds_ratio,
                'cnv_type': cnv_type
            }
    
    return gene_stats

def calculate_commonality_score(gene, icle_cnv, icle_samples, cpt6_cnv):
    """Calculate a commonality score between CPT6 and ICLE for a gene"""
    if gene not in cpt6_cnv or gene not in icle_cnv:
        return 0
    
    # Get CPT6 CNV type
    cpt6_type = cpt6_cnv[gene]
    
    # Count matching CNV types in ICLE
    matching_count = 0
    for sample in icle_samples:
        icle_type = icle_cnv[gene][sample]
        
        # For shallow deletions, only match shallow deletions
        if cpt6_type == 'Shallow_Deletion' and icle_type == 'Shallow_Deletion':
            matching_count += 1
        # For homozygous deletions, match both homozygous and shallow deletions
        elif cpt6_type == 'Homozygous_Deletion' and icle_type in ['Homozygous_Deletion', 'Shallow_Deletion']:
            matching_count += 1
        # For gains, only exact matches count
        elif cpt6_type == 'Gain' and icle_type == 'Gain':
            matching_count += 1
    
    # Calculate percentage of ICLE samples matching CPT6
    return (matching_count / len(icle_samples)) * 100

def create_cnv_comparison_plot(cpt6_cnv, icle_cnv, icle_samples, gene_stats, selected_genes, top_deleted_genes, top_gained_genes):
    """Create CNV comparison plot"""
    print("Creating CNV comparison plot")
    
    # Calculate commonality scores for each gene
    commonality_scores = {}
    for gene in selected_genes:
        commonality_scores[gene] = calculate_commonality_score(gene, icle_cnv, icle_samples, cpt6_cnv)
    
    # Separate genes into groups
    homodeletion_genes = [gene for gene in selected_genes if gene in top_deleted_genes[:4]]
    shallow_deletion_genes = [gene for gene in selected_genes if gene in top_deleted_genes[4:]]
    gain_genes = [gene for gene in selected_genes if gene in top_gained_genes]
    
    # Sort each group by commonality score in ascending order
    sorted_homodeletion = sorted(homodeletion_genes, key=lambda x: commonality_scores[x])
    sorted_shallow_deletion = sorted(shallow_deletion_genes, key=lambda x: commonality_scores[x])
    sorted_gain = sorted(gain_genes, key=lambda x: commonality_scores[x])
    
    # Combine the sorted groups while maintaining the group order
    sorted_genes = sorted_homodeletion + sorted_shallow_deletion + sorted_gain
    
    # Sort ICLE samples by CNV count
    cnv_counts = {}
    for sample in icle_samples:
        count = sum(1 for gene in selected_genes if icle_cnv[gene][sample] != 'No_CNV')
        cnv_counts[sample] = count
    
    sorted_icle_samples = sorted(icle_samples, key=lambda x: cnv_counts[x])
    
    # Create figure with gridspec layout
    fig_height = min(max(len(sorted_genes) * 0.25 + 4, 12), 24)
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
    ax_main.set_yticklabels(sorted_genes, fontsize=8)
    
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
    
    # Get viridis colors for comparison bars
    # viridis_colors = plt.cm.viridis([0.2, 0.8])
    
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
                            f'{count}/{len(icle_samples)}', va='center', ha='left', fontsize=6)
    
    # Plot TCGA frequencies
    tcga_bars = ax_comp.barh(y_pos - bar_height/2,
                            [gene_stats.get(gene, {}).get('tcga_freq', 0) for gene in sorted_genes],
                            height=bar_height, color=comparison_colors['TCGA'], alpha=0.8,
                            label=f'TCGA (n=179)')
    
    # Add TCGA frequency annotations
    for i, gene in enumerate(sorted_genes):
        if gene in gene_stats:
            freq = gene_stats[gene]['tcga_freq']
            if freq > 0:
                ax_comp.text(freq + 0.5, y_pos[i] - bar_height/2,
                            f'{freq:.1f}%', va='center', ha='left', fontsize=6)
    
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
    legend_text = ('* p<0.05 & fold change >2\n(Fisher\'s exact test)\nHomozygous Deletions: CPT6 vs TCGA\nShallow Deletions: CPT6 vs TCGA & ICLE\nGains: CPT6 vs TCGA')
    legend = ax_comp.legend(loc='upper right', bbox_to_anchor=(1, 1), 
                          frameon=True, fontsize=8)
    legend.set_title(legend_text, prop={'size': 8})
    
    # Add title indicating gene selection
    plt.suptitle('CNV Comparison: Genes Grouped by Type and Ranked by Commonality (Ascending)\nHomozygous Deletions: CPT6 vs TCGA, Shallow Deletions: CPT6 vs TCGA & ICLE', fontsize=14, y=0.98)
    
    # Adjust layout with more space
    plt.subplots_adjust(right=0.85, top=0.95, bottom=0.15)  # Added bottom margin
    
    # Save plot
    plt.savefig(os.path.join(args.output_dir, "cnv_comparison_specific_types.pdf"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.savefig(os.path.join(args.output_dir, "cnv_comparison_specific_types.png"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.close()

def print_cnv_summary(icle_cnv, icle_samples, selected_genes, cnv_counts, cpt6_homodeletion_genes, common_deleted_genes, common_gained_genes, cpt6_cnv):
    """Print summary of CNVs"""
    print("\nCNV Summary:")
    
    print("\nCNVs per cell line:")
    for sample in icle_samples:
        print(f"{sample.split('_')[0]}: {cnv_counts[sample]} CNVs")
    
    # Calculate commonality scores
    commonality_scores = {}
    for gene in selected_genes:
        if gene in cpt6_cnv:
            cpt6_type = cpt6_cnv[gene]
            matching_count = 0
            for sample in icle_samples:
                icle_type = icle_cnv[gene][sample]
                # For shallow deletions, only match shallow deletions
                if cpt6_type == 'Shallow_Deletion' and icle_type == 'Shallow_Deletion':
                    matching_count += 1
                # For homozygous deletions, match both homozygous and shallow deletions
                elif cpt6_type == 'Homozygous_Deletion' and icle_type in ['Homozygous_Deletion', 'Shallow_Deletion']:
                    matching_count += 1
                # For gains, only exact matches count
                elif cpt6_type == 'Gain' and icle_type == 'Gain':
                    matching_count += 1
            commonality_scores[gene] = (matching_count / len(icle_samples)) * 100
        else:
            commonality_scores[gene] = 0
    
    # Sort genes within each group by commonality score in ascending order
    sorted_homodeletion = sorted(cpt6_homodeletion_genes, key=lambda x: commonality_scores[x])
    sorted_deletion = sorted(common_deleted_genes, key=lambda x: commonality_scores[x])
    sorted_gain = sorted(common_gained_genes, key=lambda x: commonality_scores[x])
    
    print("\nCPT6 Homozygously Deleted Genes (sorted by commonality with ICLE in ascending order, compared with TCGA homozygous deletions):")
    for gene in sorted_homodeletion:
        homo_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Homozygous_Deletion')
        shallow_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
        print(f"{gene}: {homo_samples}/{len(icle_samples)} homozygous deletions, {shallow_samples}/{len(icle_samples)} shallow deletions in ICLE (Commonality: {commonality_scores[gene]:.1f}%)")
    
    print("\nCommon Deleted Genes (TCGA & CPT6) (sorted by commonality with ICLE in ascending order, comparing shallow deletions):")
    for gene in sorted_deletion:
        shallow_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Shallow_Deletion')
        print(f"{gene}: {shallow_samples}/{len(icle_samples)} samples with shallow deletion in ICLE (Commonality: {commonality_scores[gene]:.1f}%)")
    
    print("\nCommon Gained Genes (TCGA & CPT6) (sorted by commonality with ICLE in ascending order):")
    for gene in sorted_gain:
        gain_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] == 'Gain')
        print(f"{gene}: {gain_samples}/{len(icle_samples)} samples with gain in ICLE (Commonality: {commonality_scores[gene]:.1f}%)")
    
    print("\nICLE CNVs by cell line (grouped by type and sorted by commonality with CPT6 in ascending order):")
    # Combine the sorted groups while maintaining the group order
    all_sorted_genes = sorted_homodeletion + sorted_deletion + sorted_gain
    
    for gene in all_sorted_genes:
        cnv_samples = sum(1 for sample in icle_samples if icle_cnv[gene][sample] != 'No_CNV')
        gene_type = "Homozygous Deletion" if gene in cpt6_homodeletion_genes else ("Shallow Deletion" if gene in common_deleted_genes else "Gain")
        print(f"\n{gene} ({gene_type}): {cnv_samples}/{len(icle_samples)} samples (Commonality with CPT6: {commonality_scores[gene]:.1f}%)")
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
        
        # Exclude FAS and SMARCA2 from gain genes
        cpt6_gain_genes = [gene for gene in cpt6_gain_genes if gene not in ['FAS', 'SMARCA2']]
        if 'FAS' in cpt6_cnv:
            del cpt6_cnv['FAS']
        if 'SMARCA2' in cpt6_cnv:
            del cpt6_cnv['SMARCA2']
        
        # Process TCGA CNV data
        tcga_cnv_freq, _, tcga_top_deleted_genes, tcga_top_gained_genes = process_tcga_cnv(args.tcga_cna)
        
        # Create a list of all CPT6 genes
        all_cpt6_genes = cpt6_homodeletion_genes + cpt6_shallowdeletion_genes + cpt6_gain_genes
        
        # Find CPT6 genes that are also in TCGA data
        cpt6_tcga_genes = []
        for gene in all_cpt6_genes:
            if gene in tcga_cnv_freq:
                cpt6_tcga_genes.append(gene)
        
        # Calculate TCGA frequencies for CPT6 genes
        cpt6_tcga_freq = {}
        for gene in cpt6_tcga_genes:
            if gene in cpt6_homodeletion_genes:
                cpt6_tcga_freq[gene] = {
                    'freq': tcga_cnv_freq[gene]['Homozygous_Deletion'],
                    'type': 'Homozygous_Deletion'
                }
            elif gene in cpt6_shallowdeletion_genes:
                cpt6_tcga_freq[gene] = {
                    'freq': tcga_cnv_freq[gene]['Shallow_Deletion'],
                    'type': 'Shallow_Deletion'
                }
            elif gene in cpt6_gain_genes:
                cpt6_tcga_freq[gene] = {
                    'freq': tcga_cnv_freq[gene]['Gain'],
                    'type': 'Gain'
                }
        
        # Sort CPT6 genes by TCGA frequency
        sorted_deleted_genes = sorted(
            [(gene, cpt6_tcga_freq[gene]['freq']) for gene in cpt6_tcga_genes 
             if gene in cpt6_shallowdeletion_genes and cpt6_tcga_freq[gene]['type'] == 'Shallow_Deletion'],
            key=lambda x: x[1],
            reverse=True
        )
        
        sorted_gained_genes = sorted(
            [(gene, cpt6_tcga_freq[gene]['freq']) for gene in cpt6_tcga_genes 
             if gene in cpt6_gain_genes and cpt6_tcga_freq[gene]['type'] == 'Gain'],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Select top 10 deleted and gained genes
        top_deleted_genes = [gene for gene, _ in sorted_deleted_genes[:10]]
        top_gained_genes = [gene for gene, _ in sorted_gained_genes[:10]]
        
        # Exclude TLX1 from gained genes
        top_gained_genes = [gene for gene in top_gained_genes if gene != 'TLX1']
        
        # Combine all selected genes - always include all homodeletion genes
        all_genes = cpt6_homodeletion_genes + top_deleted_genes + top_gained_genes
        selected_genes = []
        for gene in all_genes:
            if gene not in selected_genes:
                selected_genes.append(gene)
        
        # Process ICLE CNV data for selected genes
        icle_cnv, icle_samples, gene_cnv_freq = process_icle_cnv(args.icle_cna, selected_genes)
        
        if not icle_samples:
            print("Error: No ICLE samples found.")
            sys.exit(1)
        
        # Calculate CNV frequencies and perform statistical tests
        gene_stats = calculate_cnv_frequencies(icle_cnv, icle_samples, tcga_cnv_freq, selected_genes, 
                                             cpt6_homodeletion_genes + top_deleted_genes, top_gained_genes)
        
        # Calculate CNV counts for each cell line
        cnv_counts = {}
        for sample in icle_samples:
            count = sum(1 for gene in selected_genes if icle_cnv[gene][sample] != 'No_CNV')
            cnv_counts[sample] = count
        
        # Create CNV comparison plot
        create_cnv_comparison_plot(cpt6_cnv, icle_cnv, icle_samples, gene_stats, selected_genes, 
                                 cpt6_homodeletion_genes + top_deleted_genes, top_gained_genes)
        
        # Print CNV summary
        print_cnv_summary(icle_cnv, icle_samples, selected_genes, cnv_counts, 
                        cpt6_homodeletion_genes, top_deleted_genes, top_gained_genes, cpt6_cnv)
        
        print(f"\nAnalysis complete. Results saved to the '{args.output_dir}' directory.")
    except Exception as e:
        print(f"Error in main function: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 