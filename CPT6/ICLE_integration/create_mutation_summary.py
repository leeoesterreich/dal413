import pandas as pd
import os

# Read the MAF file
maf_file = '/bgfs/alee/LO_LAB/General/cbioportal/ICLE_2023_v3/data_mutation.maf'
df = pd.read_csv(maf_file, sep='\t')

# Get total number of unique samples
total_samples = len(df['Tumor_Sample_Barcode'].unique())

# Calculate mutation frequencies per gene
gene_mutations = {}
for gene in df['Hugo_Symbol'].unique():
    gene_df = df[df['Hugo_Symbol'] == gene]
    
    # Count unique samples with mutations in this gene
    mutated_samples = len(gene_df['Tumor_Sample_Barcode'].unique())
    total_mutations = len(gene_df)
    frequency = (mutated_samples / total_samples) * 100
    
    gene_mutations[gene] = {
        'total_mutations': total_mutations,
        'mutated_samples': mutated_samples,
        'frequency': frequency
    }

# Sort genes by frequency
sorted_genes = sorted(gene_mutations.items(), key=lambda x: x[1]['frequency'], reverse=True)

# Create output directory if it doesn't exist
os.makedirs('results', exist_ok=True)

# Write summary to file
with open('results/ICLE_mutation_frequency_summary.txt', 'w') as f:
    # Write header
    f.write("Gene\tMutSig(Q-value)\t# Mut\t#\tProfiled Samples\tFreq\tIs Cancer Gene (source: OncoKB)\n")
    
    # Write data for each gene
    for gene, data in sorted_genes:
        # Format: Gene, MutSig(Q-value), # Mut, # Mutated Samples, Total Samples, Freq%, Is Cancer Gene
        f.write(f"{gene}\t\t{data['total_mutations']}\t{data['mutated_samples']}\t{total_samples}\t{data['frequency']:.1f}%\t\n")