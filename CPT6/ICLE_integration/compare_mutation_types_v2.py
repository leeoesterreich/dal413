import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import sys
import argparse
from matplotlib.colors import ListedColormap
from matplotlib import cm
from scipy.stats import fisher_exact

# Set up argument parser
parser = argparse.ArgumentParser(description='Compare mutation types between CPT6 and ICLE cell lines')
parser.add_argument('--cpt6_maf', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/CPT6.maf", 
                    help='Path to CPT6 MAF file')
parser.add_argument('--cpt6_ref', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/agg_sorted.csv", 
                    help='Path to CPT6 mutation reference file')
parser.add_argument('--icle_maf', type=str, default="/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf", 
                    help='Path to ICLE MAF file')
parser.add_argument('--icle_cna', type=str, default="/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_CNA_gistic.txt", 
                    help='Path to ICLE copy number alteration file')
parser.add_argument('--tcga_freq', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/TCGA2015_Metabirc2012_16_MSK2018_19_20_ILC.txt",
                    help='Path to TCGA frequency data')
parser.add_argument('--oncogene_list', type=str, default="/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/Oncogene_CPT6.csv", 
                    help='Path to oncogene list file')
parser.add_argument('--output_dir', type=str, default="results", 
                    help='Directory to save results')
args = parser.parse_args()

# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Map CPT6 consequences to cBioPortal variant classifications
consequence_map = {
    'frameshift_variant': 'Frame_Shift_Del',
    'frameshift_insertion': 'Frame_Shift_Ins',
    'frameshift_deletion': 'Frame_Shift_Del',
    'inframe_deletion': 'In_Frame_Del',
    'inframe_insertion': 'In_Frame_Ins',
    'missense_variant': 'Missense_Mutation',
    'splice_acceptor_variant': 'Splice_Site',
    'splice_donor_variant': 'Splice_Site',
    'stop_gained': 'Nonsense_Mutation',
    'stop_lost': 'Nonstop_Mutation',
    'start_lost': 'Translation_Start_Site',
    'protein_altering_variant': 'Missense_Mutation'
}

# Define mutation types and create viridis color map
mutation_types = [
    'Frame_Shift_Del', 
    'Frame_Shift_Ins', 
    'In_Frame_Del', 
    'In_Frame_Ins', 
    'Missense_Mutation', 
    'Splice_Site', 
    'Nonsense_Mutation', 
    'Nonstop_Mutation', 
    'Translation_Start_Site', 
    'Silent', 
    'Other',
    'Homozygous_Deletion',
    'High_Amplification',
    'No_Mutation'
]

# Generate viridis colors for mutation types
viridis = cm.get_cmap('viridis', len(mutation_types) - 1)  # -1 to exclude No_Mutation
color_map = {mutation_types[i]: viridis(i / (len(mutation_types) - 2)) for i in range(len(mutation_types) - 1)}
color_map['No_Mutation'] = '#F0F0F0'  # Light Gray for unmutated

def read_oncogene_list(file_path):
    """Read the oncogene list from a CSV file"""
    try:
        return pd.read_csv(file_path, header=None)[0].tolist()
    except Exception as e:
        print(f"Error reading oncogene list: {str(e)}")
        return []

def convert_mouse_to_human_genes(gene):
    """Convert mouse gene to human homolog"""
    # Simple conversion for known genes
    if gene.lower() == 'trp53':
        return 'TP53'
    elif gene.lower() in ['erbb2', 'kras', 'braf', 'egfr']:
        return gene.upper()
    # Default: uppercase the gene name
    return gene.upper()

def get_mutated_genes(cpt6_mutations):
    """Get list of genes that have mutations in CPT6"""
    return [gene for gene, mutations in cpt6_mutations.items() 
            if mutations and mutations[0] != 'No_Mutation']

def process_cpt6_mutations_from_ref(ref_file_path, oncogenes):
    """Process CPT6 mutations from reference file"""
    print(f"Processing CPT6 mutations from reference file {ref_file_path}")
    cpt6_mutations = {gene: [] for gene in oncogenes}
    
    try:
        # Read the CPT6 reference file - force reload
        cpt6_ref_df = pd.read_csv(ref_file_path)
        
        # Process each mutation
        for _, row in cpt6_ref_df.iterrows():
            gene_symbol = row['Gene Symbol']
            if gene_symbol in oncogenes:
                consequence = row['Consequence']
                variant_type = consequence_map.get(consequence, 'Other')
                cpt6_mutations[gene_symbol].append(variant_type)
        
        # For genes with no mutations, add a placeholder
        for gene in oncogenes:
            if not cpt6_mutations[gene]:
                cpt6_mutations[gene].append('No_Mutation')
        
        # Filter to keep only genes with mutations
        mutated_genes = get_mutated_genes(cpt6_mutations)
        return {gene: cpt6_mutations[gene] for gene in mutated_genes}
    except Exception as e:
        print(f"Error processing CPT6 mutations from reference: {str(e)}")
        return {}

def process_cpt6_mutations(cpt6_maf_path, oncogenes):
    """Process CPT6 mutations from MAF file"""
    print(f"Processing CPT6 mutations from {cpt6_maf_path}")
    cpt6_mutations = {gene: [] for gene in oncogenes}
    
    try:
        # Read the CPT6 MAF file
        cpt6_df = pd.read_csv(cpt6_maf_path, sep='\t', comment='#')
        
        # Convert mouse genes to human homologs
        cpt6_df['Human_Gene'] = cpt6_df['Hugo_Symbol'].apply(convert_mouse_to_human_genes)
        
        # Filter for oncogenes
        for _, row in cpt6_df.iterrows():
            human_gene = row['Human_Gene']
            if human_gene in oncogenes:
                # Map the consequence to cBioPortal variant classification
                consequence = row.get('Consequence', '').split('&')[0]  # Take first consequence if multiple
                variant_type = consequence_map.get(consequence, 'Other')
                cpt6_mutations[human_gene].append(variant_type)
        
        # For genes with no mutations, add a placeholder
        for gene in oncogenes:
            if not cpt6_mutations[gene]:
                cpt6_mutations[gene].append('No_Mutation')
                
        return cpt6_mutations
    except Exception as e:
        print(f"Error processing CPT6 mutations: {str(e)}")
        # Return empty mutations with placeholders
        return {gene: ['No_Mutation'] for gene in oncogenes}

def process_icle_mutations(icle_maf_path, icle_cna_path, oncogenes):
    """Process ICLE mutations and copy number alterations"""
    print(f"Processing ICLE mutations from {icle_maf_path}")
    print(f"Processing ICLE copy number alterations from {icle_cna_path}")
    
    try:
        # Read the ICLE MAF file
        icle_df = pd.read_csv(icle_maf_path, sep='\t')
        
        # Read the ICLE CNA file
        cna_df = pd.read_csv(icle_cna_path, sep='\t')
        
        # Get unique samples from MAF file
        icle_samples = sorted(set(icle_df['Tumor_Sample_Barcode']))
        
        # Initialize mutations dictionary
        icle_mutations = {gene: {sample: [] for sample in icle_samples} for gene in oncogenes}
        
        # Fill in actual mutations from MAF file
        for gene in oncogenes:
            gene_muts = icle_df[icle_df['Hugo_Symbol'] == gene]
            for _, row in gene_muts.iterrows():
                sample = row['Tumor_Sample_Barcode']
                mut_type = row['Variant_Classification']
                icle_mutations[gene][sample].append(mut_type)
        
        # Process copy number alterations
        if 'Hugo_Symbol' in cna_df.columns:
            # Filter for oncogenes
            cna_oncogenes = cna_df[cna_df['Hugo_Symbol'].isin(oncogenes)]
            
            # Map sample names from CNA to MAF format
            cna_samples = [col for col in cna_df.columns if col not in ['Hugo_Symbol', 'Entrez_Gene_Id']]
            sample_map = {}
            for cna_sample in cna_samples:
                for maf_sample in icle_samples:
                    if cna_sample in maf_sample:
                        sample_map[cna_sample] = maf_sample
            
            # Add copy number alterations to mutations
            for _, row in cna_oncogenes.iterrows():
                gene = row['Hugo_Symbol']
                for cna_sample in cna_samples:
                    if cna_sample in sample_map:
                        maf_sample = sample_map[cna_sample]
                        cna_value = row[cna_sample]
                        
                        # Only include homozygous deletions (-2) and high-level amplifications (2)
                        if cna_value == -2:
                            icle_mutations[gene][maf_sample].append('Homozygous_Deletion')
                        elif cna_value == 2:
                            icle_mutations[gene][maf_sample].append('High_Amplification')
        
        return icle_mutations, icle_samples
    except Exception as e:
        print(f"Error processing ICLE mutations and CNAs: {str(e)}")
        return {}, []

def calculate_mutation_counts(icle_mutations, icle_samples):
    """Calculate mutation counts for each gene and cell line"""
    gene_freq = {}
    mutation_counts = {}
    mutations_per_cell_line = {}
    
    # Calculate mutations per cell line
    for sample in icle_samples:
        count = 0
        for gene in icle_mutations:
            if icle_mutations[gene][sample]:  # If there are any mutations for this gene in this sample
                count += 1  # Count each gene only once per sample
        mutations_per_cell_line[sample] = count
    
    # Calculate gene frequencies
    for gene in icle_mutations:
        mut_samples = sum(1 for sample in icle_samples if icle_mutations[gene][sample])
        gene_freq[gene] = (mut_samples / len(icle_samples)) * 100
        mutation_counts[gene] = mut_samples
    
    return gene_freq, mutation_counts, mutations_per_cell_line

def calculate_gene_frequencies(icle_mutations, icle_samples, tcga_data):
    """Calculate gene frequencies and perform Fisher's exact test"""
    gene_stats = {}
    icle_total = len(icle_samples)
    tcga_total = 643  # TCGA sample count
    
    # Convert TCGA frequencies to float
    tcga_data['Freq'] = tcga_data['Freq'].str.rstrip('%').astype(float)
    
    for gene in icle_mutations:
        # Calculate ICLE frequency
        mut_samples = sum(1 for sample in icle_samples if icle_mutations[gene][sample])
        icle_freq = (mut_samples / icle_total) * 100
        
        # Get TCGA frequency if available
        tcga_freq = 0
        pvalue = 1
        odds_ratio = 1
        if gene in tcga_data['Gene'].values:
            tcga_freq = tcga_data[tcga_data['Gene'] == gene]['Freq'].iloc[0]
            
            # Only perform Fisher's exact test if ICLE has mutations
            if mut_samples > 0:
                tcga_count = int(np.round(tcga_freq * tcga_total / 100))
                table = [[mut_samples, icle_total - mut_samples],
                        [tcga_count, tcga_total - tcga_count]]
                odds_ratio, pvalue = fisher_exact(table)
        
        gene_stats[gene] = {
            'icle_freq': icle_freq,
            'tcga_freq': tcga_freq,
            'pvalue': pvalue,
            'odds_ratio': odds_ratio
        }
    
    return gene_stats

def create_oncoprint_plot(cpt6_mutations, icle_mutations, icle_samples, gene_freq, mutation_counts, mutations_per_cell_line):
    """Create oncoprint-style plot comparing CPT6 and ICLE mutations"""
    print("Creating oncoprint comparison plot")
    
    # Read TCGA frequency data
    tcga_data = pd.read_csv(args.tcga_freq, sep='\t')
    gene_stats = calculate_gene_frequencies(icle_mutations, icle_samples, tcga_data)
    
    # Sort genes by ICLE frequency
    sorted_genes = sorted(cpt6_mutations.keys(), 
                        key=lambda x: gene_stats[x]['icle_freq'],
                        reverse=True)
    
    # Sort ICLE samples by mutation count (ascending)
    sorted_icle_samples = sorted(icle_samples, key=lambda x: mutations_per_cell_line[x])
    
    # Create figure with gridspec layout
    fig_height = min(max(len(sorted_genes) * 0.25 + 4, 12), 24)  # Adjusted scaling
    fig_width = min(max(len(sorted_icle_samples) * 0.25 + 14, 24), 30)  # Increased width for legend
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=100)
    
    # Create a gridspec to handle the layout with more space for legend
    gs = plt.GridSpec(2, 20, height_ratios=[2, 15], width_ratios=[0.6]*12 + [0.4]*8, figure=fig, hspace=0.1)
    
    # Create plots
    ax_top = fig.add_subplot(gs[0, :12])  # Mutation counts
    ax_main = fig.add_subplot(gs[1, :12])  # Main oncoprint
    ax_comp = fig.add_subplot(gs[1, 12:])  # Comparison plot
    
    # Define box width for consistency
    box_width = 0.7
    
    # Calculate positions for CPT6 and ICLE cell lines
    cpt6_position = -1
    gap = 1  # Gap between CPT6 and ICLE
    icle_positions = np.arange(gap, gap + len(sorted_icle_samples))
    
    # Get mutation counts for sorted samples
    mutation_counts_sorted = [mutations_per_cell_line[s] for s in sorted_icle_samples]
    
    # Plot bars
    ax_top.bar(icle_positions, mutation_counts_sorted, 
              color='lightgray', alpha=0.5, width=box_width)
    
    # Add count labels
    for i, count in enumerate(mutation_counts_sorted):
        ax_top.text(icle_positions[i], count + 0.5, str(count), 
                   ha='center', va='bottom', fontsize=8)
    
    # Customize top plot
    ax_top.set_ylim(0, max(mutation_counts_sorted) * 1.2)  # Add 20% padding for labels
    ax_top.set_yticks([])  # Remove y-ticks
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_visible(False)  # Hide left spine
    
    # Add "Mutation Gene Count" annotation to the left of the top plot
    ax_top.text(-1.5, np.mean(mutation_counts_sorted), 
                "Mutation\nGene Count", ha='right', va='center', fontsize=10, fontweight='bold')
    
    # Set the same x-limits for top and main plots
    ax_top.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    ax_top.set_xticks([])  # Remove x-ticks from top plot
    
    # Plot settings
    bar_height = 0.6  # Reduced height to prevent overlap
    y_pos = np.arange(len(sorted_genes))
    
    # Add extra padding to y-axis to show full boxes
    y_padding = 0.5
    
    # Plot CPT6 mutations
    for i, gene in enumerate(sorted_genes):
        mutations = cpt6_mutations[gene]
        unique_muts = set(mutations)
        if 'No_Mutation' in unique_muts and len(unique_muts) > 1:
            unique_muts.remove('No_Mutation')
        
        mut_box_width = box_width / len(unique_muts) if unique_muts else box_width
        for j, mut_type in enumerate(unique_muts):
            if mut_type != 'No_Mutation':
                # Center the box around the cpt6_position
                box_start = cpt6_position - box_width/2 + j * mut_box_width
                ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                              mut_box_width, bar_height,
                                              facecolor=color_map.get(mut_type, color_map['Other']),
                                              alpha=0.7))
    
    # Plot ICLE mutations with background boxes
    for i, gene in enumerate(sorted_genes):
        # First plot gray background boxes for all samples
        for j, pos in enumerate(icle_positions):
            # Center the box around the position
            box_start = pos - box_width/2
            ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                          box_width, bar_height,
                                          facecolor=color_map['No_Mutation'],
                                          alpha=0.3))
        
        # Then plot mutation boxes
        for j, sample in enumerate(sorted_icle_samples):
            mutations = icle_mutations[gene][sample]
            if mutations:
                unique_muts = set(mutations)
                mut_box_width = box_width / len(unique_muts)
                for k, mut_type in enumerate(unique_muts):
                    # Center the box around the position
                    box_start = icle_positions[j] - box_width/2 + k * mut_box_width
                    ax_main.add_patch(plt.Rectangle((box_start, y_pos[i] - bar_height/2),
                                                 mut_box_width, bar_height,
                                                 facecolor=color_map.get(mut_type, color_map['Other']),
                                                 alpha=0.7))
    
    # Customize main plot
    ax_main.set_yticks(y_pos)
    ax_main.set_yticklabels(sorted_genes, fontsize=8)  # Added back gene labels with smaller font
    
    # Set y-axis limits with padding to show full boxes
    ax_main.set_ylim(-y_padding, len(sorted_genes) - 1 + y_padding)
    
    # Set x-axis limits and ticks
    ax_main.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    
    # Add all x-ticks at the exact positions (CPT6 and ICLE cell lines)
    all_tick_positions = [cpt6_position] + list(icle_positions)
    all_tick_labels = ['CPT6'] + [s.split('_')[0] for s in sorted_icle_samples]
    
    ax_main.set_xticks(all_tick_positions)
    ax_main.set_xticklabels(all_tick_labels, rotation=45, ha='right', va='top', fontsize=8)
    
    # Add x-axis label for ICLE cell lines
    ax_main.text(gap + len(sorted_icle_samples)/2, -len(sorted_genes)*0.08, 'ICLE Cell Lines', 
                ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Plot comparison bars
    bar_height = 0.35
    
    # Calculate maximum frequency for scaling
    max_freq = max(max(gene_stats[gene]['icle_freq'] for gene in sorted_genes),
                  max(gene_stats[gene]['tcga_freq'] for gene in sorted_genes))
    
    # Get viridis colors for comparison bars
    viridis_colors = plt.cm.viridis([0.2, 0.8])  # Get two distinct colors from viridis
    
    # Plot ICLE frequencies
    icle_bars = ax_comp.barh(y_pos + bar_height/2,
                            [gene_stats[gene]['icle_freq'] for gene in sorted_genes],
                            height=bar_height, color=viridis_colors[0], alpha=0.8,
                            label=f'ICLE (n=17)')
    
    # Add ICLE frequency annotations
    for i, gene in enumerate(sorted_genes):
        freq = gene_stats[gene]['icle_freq']
        mut_samples = sum(1 for sample in icle_samples if icle_mutations[gene][sample])
        if freq > 0:  # Only annotate non-zero frequencies
            ax_comp.text(freq + 0.5, y_pos[i] + bar_height/2,
                        f'{mut_samples}/17', va='center', ha='left', fontsize=6)

    # Plot TCGA frequencies
    tcga_bars = ax_comp.barh(y_pos - bar_height/2,
                            [gene_stats[gene]['tcga_freq'] for gene in sorted_genes],
                            height=bar_height, color=viridis_colors[1], alpha=0.8,
                            label=f'TCGA (n=643)')
    
    # Add TCGA frequency annotations
    for i, gene in enumerate(sorted_genes):
        freq = gene_stats[gene]['tcga_freq']
        if freq > 0:  # Only annotate non-zero frequencies
            ax_comp.text(freq + 0.5, y_pos[i] - bar_height/2,
                        f'{freq:.1f}%', va='center', ha='left', fontsize=6)

    # Add significance stars with more spacing
    for i, gene in enumerate(sorted_genes):
        stats = gene_stats[gene]
        if stats['icle_freq'] > 0 and stats['pvalue'] < 0.05 and (stats['odds_ratio'] > 2 or stats['odds_ratio'] < 0.5):
            stars = '*' * sum([stats['pvalue'] < cutoff for cutoff in [0.05, 0.01, 0.001]])
            ax_comp.text(max(stats['icle_freq'], stats['tcga_freq']) + 4,
                        i, stars, va='center', ha='left', fontsize=8)
    
    # Customize comparison plot
    max_freq_extended = 35  # Reduced from 40 to 35 to make x-axis shorter
    ax_comp.set_xlim(0, max_freq_extended)  # Extended x-axis with fixed value
    ax_comp.set_ylim(ax_main.get_ylim())
    ax_comp.set_xlabel('Alteration Frequency (%)', fontsize=10)
    ax_comp.grid(True, linestyle='--', alpha=0.3, axis='x')
    ax_comp.set_yticks([])
    ax_comp.spines['top'].set_visible(False)
    ax_comp.spines['right'].set_visible(False)
    
    # Add x-ticks at regular intervals
    ax_comp.set_xticks([0, 5, 10, 15, 20, 25, 30, 35])  # Updated ticks to match 35% limit with 5% intervals
    
    # Add legend with Fisher's exact test information - remove first two lines
    legend_text = ('* p<0.05 & fold change >2\n(Fisher\'s exact test)')
    legend = ax_comp.legend(loc='upper right', bbox_to_anchor=(1, 1), 
                          frameon=True, fontsize=8)
    legend.set_title(legend_text, prop={'size': 8})
    
    # Adjust layout with more space
    plt.subplots_adjust(right=0.85)  # Add more space on the right for legend
    
    # Save plot with reduced DPI and extra space for legend
    plt.savefig(os.path.join(args.output_dir, "oncoprint_comparison_v2.pdf"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.savefig(os.path.join(args.output_dir, "oncoprint_comparison_v2.png"), 
                bbox_inches='tight', dpi=150, pad_inches=0.5)
    plt.close()

def print_summary(cpt6_mutations, icle_mutations, icle_samples, mutations_per_cell_line, sorted_genes):
    """Print summary of mutations"""
    print("\nMutation Summary:")
    print("\nMutations per cell line:")
    for sample in icle_samples:
        print(f"{sample.split('_')[0]}: {mutations_per_cell_line[sample]} mutations")

    print("\nCPT6 mutations:")
    for gene in sorted_genes:
        muts = cpt6_mutations[gene]
        print(f"{gene}: {', '.join(set(muts))}")

    print("\nICLE mutations by cell line:")
    for gene in sorted_genes:
        mut_samples = sum(1 for sample in icle_samples if icle_mutations[gene][sample])
        print(f"\n{gene}: {mut_samples}/{len(icle_samples)} samples")
        for sample in icle_samples:
            muts = icle_mutations[gene][sample]
            if muts:
                print(f"  {sample}: {', '.join(set(muts))}")

def main():
    try:
        # Read oncogene list
        oncogenes = read_oncogene_list(args.oncogene_list)
        if not oncogenes:
            print("Error: No oncogenes found in the list.")
            sys.exit(1)
        
        # Process CPT6 mutations from reference file and get only mutated genes
        cpt6_mutations = process_cpt6_mutations_from_ref(args.cpt6_ref, oncogenes)
        mutated_genes = list(cpt6_mutations.keys())  # Get only genes with mutations
        
        # Process ICLE mutations and copy number alterations for mutated genes only
        icle_mutations, icle_samples = process_icle_mutations(args.icle_maf, args.icle_cna, mutated_genes)
        if not icle_samples:
            print("Error: No ICLE samples found.")
            sys.exit(1)
        
        # Calculate mutation counts
        gene_freq, mutation_counts, mutations_per_cell_line = calculate_mutation_counts(icle_mutations, icle_samples)
        
        # Sort genes by mutation frequency
        sorted_genes = sorted(mutated_genes, key=lambda x: gene_freq.get(x, 0), reverse=True)
        
        # Create oncoprint plot
        create_oncoprint_plot(cpt6_mutations, icle_mutations, icle_samples, gene_freq, mutation_counts, mutations_per_cell_line)
        
        # Print summary
        print_summary(cpt6_mutations, icle_mutations, icle_samples, mutations_per_cell_line, sorted_genes)
        
        print(f"\nAnalysis complete. Results saved to the '{args.output_dir}' directory.")
    except Exception as e:
        print(f"Error in main function: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 