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
                
        return cpt6_mutations
    except Exception as e:
        print(f"Error processing CPT6 mutations from reference: {str(e)}")
        # Return empty mutations with placeholders
        return {gene: ['No_Mutation'] for gene in oncogenes}

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
        count = sum(1 for gene in icle_mutations if icle_mutations[gene][sample])
        mutations_per_cell_line[sample] = count
    
    # Calculate gene frequencies
    for gene in icle_mutations:
        mut_samples = sum(1 for sample in icle_samples if icle_mutations[gene][sample])
        gene_freq[gene] = (mut_samples / len(icle_samples)) * 100
        mutation_counts[gene] = mut_samples
    
    return gene_freq, mutation_counts, mutations_per_cell_line

def create_oncoprint_plot(cpt6_mutations, icle_mutations, icle_samples, gene_freq, mutation_counts, mutations_per_cell_line):
    """Create oncoprint-style plot comparing CPT6 and ICLE mutations"""
    print("Creating oncoprint comparison plot")
    
    # Sort genes by mutation frequency
    sorted_genes = sorted(cpt6_mutations.keys(), key=lambda x: gene_freq.get(x, 0), reverse=True)
    
    # Sort ICLE samples by mutation count (ascending)
    sorted_icle_samples = sorted(icle_samples, key=lambda x: mutations_per_cell_line[x])
    
    # Create figure with gridspec layout
    fig_height = max(len(sorted_genes) * 0.4 + 6, 14)  # Increased minimum height
    fig_width = max(len(sorted_icle_samples) * 0.3 + 12, 24)  # Minimum width of 24 inches
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create a gridspec to handle the layout
    gs = plt.GridSpec(3, 12, height_ratios=[0.5, 1, 6], figure=fig)
    
    # Add title at the top
    ax_title = fig.add_subplot(gs[0, :])
    ax_title.text(0.5, 0.5, 'Comparison of Mutation Types between CPT6 and ICLE Cell Lines', 
                 ha='center', va='center', fontsize=16, fontweight='bold')
    ax_title.axis('off')
    
    # Create top plot for cell line mutation counts - ensure it aligns with main plot
    ax_top = fig.add_subplot(gs[1, :10])  # Use the same width as the main plot
    ax_main = fig.add_subplot(gs[2, :10])  # Main plot
    ax_freq = fig.add_subplot(gs[2, 10:])  # Frequency plot
    
    # Define box width for consistency
    box_width = 0.8
    
    # Calculate positions for CPT6 and ICLE cell lines
    cpt6_position = -1
    gap = 1  # Gap between CPT6 and ICLE
    icle_positions = np.arange(gap, gap + len(sorted_icle_samples))
    
    # Plot only ICLE cell line mutation counts - align exactly with columns
    ax_top.bar(icle_positions, [mutations_per_cell_line[s] for s in sorted_icle_samples], 
              color='lightgray', alpha=0.5, width=box_width)
    
    # Add count labels for ICLE only
    for i, count in enumerate([mutations_per_cell_line[s] for s in sorted_icle_samples]):
        ax_top.text(icle_positions[i], count, str(count), ha='center', va='bottom')
    
    # Add "Mutation Counts" annotation to the left of the top plot - moved a bit to the right
    ax_top.text(-1.0, np.mean([mutations_per_cell_line[s] for s in sorted_icle_samples]), 
                "Mutation\nCounts", ha='right', va='center', fontsize=10, fontweight='bold')
    
    # Set the same x-limits for top and main plots to ensure alignment
    ax_top.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    
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
    
    # Plot mutation frequency bars and text
    max_count = len(sorted_icle_samples)
    for i, gene in enumerate(sorted_genes):
        count = mutation_counts[gene]
        ax_freq.barh(y_pos[i], count, height=bar_height, color='lightgray', alpha=0.5)
        # Position text to the right of the bar
        ax_freq.text(count + 0.2, y_pos[i], f'{count}/{len(sorted_icle_samples)}',
                    va='center', ha='left', fontsize=10)
    
    # Customize main plot
    ax_main.set_yticks(y_pos)
    ax_main.set_yticklabels(sorted_genes)
    
    # Set y-axis limits with padding to show full boxes
    ax_main.set_ylim(-y_padding, len(sorted_genes) - 1 + y_padding)
    
    # Set x-axis limits and ticks
    ax_main.set_xlim(-1.5, len(sorted_icle_samples) + gap + 0.5)
    
    # Add all x-ticks at the exact positions (CPT6 and ICLE cell lines)
    all_tick_positions = [cpt6_position] + list(icle_positions)
    all_tick_labels = ['CPT6'] + [s.split('_')[0] for s in sorted_icle_samples]
    
    ax_main.set_xticks(all_tick_positions)
    ax_main.set_xticklabels(all_tick_labels, rotation=45, ha='right', va='top')
    
    # Add x-axis label for ICLE cell lines - moved further down to prevent overlap
    ax_main.text(gap + len(sorted_icle_samples)/2, -len(sorted_genes)*0.08, 'ICLE Cell Lines', 
                ha='center', va='top', fontsize=12, fontweight='bold')
    
    # Customize frequency plot
    ax_freq.set_xlim(-1, max_count + 2)  # Extended to make room for text
    ax_freq.set_ylim(ax_main.get_ylim())
    ax_freq.set_yticks([])
    ax_freq.set_xticks([])
    ax_freq.spines['left'].set_visible(False)
    ax_freq.spines['right'].set_visible(False)
    ax_freq.spines['top'].set_visible(False)
    ax_freq.spines['bottom'].set_visible(False)
    
    # Add legend at the bottom of the plot with further reduced space
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.7)
                      for mut_type, color in color_map.items() if mut_type != 'No_Mutation']
    ax_main.legend(legend_elements, [k for k in color_map.keys() if k != 'No_Mutation'], 
                  loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4)
    
    # Remove unnecessary spines and grid
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.grid(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot with extra space at bottom for legend
    plt.savefig(os.path.join(args.output_dir, "oncoprint_comparison.pdf"), bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(args.output_dir, "oncoprint_comparison.png"), bbox_inches='tight', dpi=300)
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
        
        # Process CPT6 mutations from reference file
        cpt6_mutations = process_cpt6_mutations_from_ref(args.cpt6_ref, oncogenes)
        
        # Process ICLE mutations and copy number alterations
        icle_mutations, icle_samples = process_icle_mutations(args.icle_maf, args.icle_cna, oncogenes)
        if not icle_samples:
            print("Error: No ICLE samples found.")
            sys.exit(1)
        
        # Calculate mutation counts
        gene_freq, mutation_counts, mutations_per_cell_line = calculate_mutation_counts(icle_mutations, icle_samples)
        
        # Sort genes by mutation frequency
        sorted_genes = sorted(oncogenes, key=lambda x: gene_freq.get(x, 0), reverse=True)
        
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