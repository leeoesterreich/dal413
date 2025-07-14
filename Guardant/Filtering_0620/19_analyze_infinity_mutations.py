import pandas as pd
import numpy as np
import sys
from io import StringIO
import datetime

# Create a string buffer to capture all output
output_buffer = StringIO()
original_stdout = sys.stdout
sys.stdout = output_buffer

def analyze_gene_mutations(idc_df, ilc_df, gene):
    print(f"\n{gene} Mutation Types:")
    
    # IDC Analysis
    idc_data = idc_df[idc_df['Gene'] == gene]
    idc_total = len(idc_data)
    print(f"  IDC ({idc_total} total mutations):")
    if idc_total > 0:
        type_counts = idc_data['Type'].value_counts().sort_index()
        for mut_type, count in type_counts.items():
            print(f"    - {mut_type}: {count} mutations")
    
    # ILC Analysis
    ilc_data = ilc_df[ilc_df['Gene'] == gene]
    ilc_total = len(ilc_data)
    print(f"  ILC ({ilc_total} total mutations):")
    if ilc_total > 0:
        type_counts = ilc_data['Type'].value_counts().sort_index()
        for mut_type, count in type_counts.items():
            print(f"    - {mut_type}: {count} mutations")
    
    # Combined Analysis
    combined_data = pd.concat([idc_data, ilc_data])
    combined_total = len(combined_data)
    print(f"  Combined ({combined_total} total mutations):")
    if combined_total > 0:
        type_counts = combined_data['Type'].value_counts().sort_index()
        for mut_type, count in type_counts.items():
            print(f"    - {mut_type}: {count} mutations")

def analyze_mutations(idc_df, ilc_df):
    # Add timestamp
    print(f"Analysis performed on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nReading Infinity data files...")
    
    # Overall statistics
    print("\nOverall Statistics:")
    print("-" * 20)
    
    # IDC
    idc_unique_mutations = idc_df[['Gene', 'Alteration']].drop_duplicates()
    idc_unique_genes = idc_df['Gene'].dropna().unique()
    print(f"\nIDC Infinity Cohort:")
    print(f"  Unique Mutations Identified: {len(idc_unique_mutations)}")
    print(f"  Unique Genes Mutated: {len(idc_unique_genes)}")
    
    # ILC
    ilc_unique_mutations = ilc_df[['Gene', 'Alteration']].drop_duplicates()
    ilc_unique_genes = ilc_df['Gene'].dropna().unique()
    print(f"\nILC Infinity Cohort:")
    print(f"  Unique Mutations Identified: {len(ilc_unique_mutations)}")
    print(f"  Unique Genes Mutated: {len(ilc_unique_genes)}")
    
    # Combined
    combined_df = pd.concat([idc_df, ilc_df])
    combined_unique_mutations = combined_df[['Gene', 'Alteration']].drop_duplicates()
    combined_unique_genes = combined_df['Gene'].dropna().unique()
    print(f"\nCombined Cohorts:")
    print(f"  Total Unique Mutations Identified: {len(combined_unique_mutations)}")
    print(f"  Total Unique Genes Mutated: {len(combined_unique_genes)}")
    
    # Analyze key genes
    print("\nDetailed Gene Analysis:")
    print("-" * 20)
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    for gene in key_genes:
        analyze_gene_mutations(idc_df, ilc_df, gene)

# Read the Infinity data files
idc_infinity = pd.read_csv('IDC_genomic_infinity.csv')
ilc_infinity = pd.read_csv('ILC_genomic_infinity.csv')

# Perform analysis
analyze_mutations(idc_infinity, ilc_infinity)

# Save summary to file
print("\nGenerating summary report...")

# Create summary DataFrame
summary_data = {
    'Cohort': ['IDC', 'ILC', 'Combined'],
    'Unique Mutations Identified': [
        len(idc_infinity[['Gene', 'Alteration']].drop_duplicates()),
        len(ilc_infinity[['Gene', 'Alteration']].drop_duplicates()),
        len(pd.concat([idc_infinity, ilc_infinity])[['Gene', 'Alteration']].drop_duplicates())
    ],
    'Unique Genes Mutated': [
        len(idc_infinity['Gene'].unique()),
        len(ilc_infinity['Gene'].unique()),
        len(pd.concat([idc_infinity, ilc_infinity])['Gene'].unique())
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('infinity_mutation_summary.csv', index=False)
print("Summary saved to: infinity_mutation_summary.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('infinity_mutation_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: infinity_mutation_analysis_report.txt") 