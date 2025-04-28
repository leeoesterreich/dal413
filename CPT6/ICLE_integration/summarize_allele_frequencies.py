import pandas as pd
import os

# Read the merged MAF file
maf_file = 'results/merged_maf.maf'
df = pd.read_csv(maf_file, sep='\t', comment='#')

# Count unique samples
total_samples = len(df['Tumor_Sample_Barcode'].unique())

# Group by gene and calculate frequencies
gene_summary = df.groupby('Hugo_Symbol').agg({
    'Tumor_Sample_Barcode': lambda x: len(set(x)),  # Number of unique samples with mutation
    'Variant_Classification': 'count'  # Total number of mutations
}).reset_index()

# Rename columns
gene_summary.columns = ['Gene', '# Mut', 'Profiled Samples']

# Calculate frequency
gene_summary['Freq'] = (gene_summary['Profiled Samples'] / total_samples * 100).round(1).astype(str) + '%'

# Sort by number of affected samples
gene_summary = gene_summary.sort_values('Profiled Samples', ascending=False)

# Add empty columns for MutSig and Is Cancer Gene
gene_summary.insert(1, 'MutSig(Q-value)', '')
gene_summary.insert(5, 'Is Cancer Gene (source: OncoKB)', '')

# Save to file
output_file = 'results/ICLE_mutation_frequency_summary.txt'
gene_summary.to_csv(output_file, sep='\t', index=False)

print(f"Analysis complete. Results saved to {output_file}")
print(f"Total number of samples analyzed: {total_samples}") 