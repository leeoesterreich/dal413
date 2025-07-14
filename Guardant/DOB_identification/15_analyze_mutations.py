import pandas as pd
import numpy as np
import sys
from io import StringIO
import datetime
import os

# Create output directory if it doesn't exist
os.makedirs('Vus_filtering', exist_ok=True)

# Create a string buffer to capture all output
output_buffer = StringIO()
original_stdout = sys.stdout
sys.stdout = output_buffer

def analyze_mutations(df, cohort_name):
    print(f"\n{cohort_name} Cohort Analysis (Filtered for Vus=FALSE):")
    print("-" * (len(cohort_name) + 35))
    
    # Count unique Gene-Alteration combinations
    gene_alt_combinations = df[['Gene', 'Alteration']].drop_duplicates()
    print(f"\nTotal unique Gene-Alteration combinations: {len(gene_alt_combinations)}")
    
    # Count unique genes with mutations
    unique_genes = df['Gene'].dropna().unique()
    print(f"Total unique genes with mutations: {len(unique_genes)}")
    
    # Function to check if a gene is in the string (for fusion genes)
    def contains_gene(gene_string, target_gene):
        if pd.isna(gene_string):
            return False
        return target_gene in str(gene_string)
    
    # Count specific genes
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    for gene in key_genes:
        gene_count = df[df['Gene'].apply(lambda x: contains_gene(x, gene))]['Effective Patient ID'].nunique()
        total_mutations = len(df[df['Gene'].apply(lambda x: contains_gene(x, gene))])
        print(f"\n{gene}:")
        print(f"- Unique patients with mutations: {gene_count}")
        print(f"- Total mutations (including repeats): {total_mutations}")
        
        # Get unique alterations for this gene
        gene_alterations = df[df['Gene'].apply(lambda x: contains_gene(x, gene))][['Gene', 'Alteration']].drop_duplicates()
        print("- Unique alterations:")
        for _, row in gene_alterations.iterrows():
            print(f"  * {row['Gene']}: {row['Alteration']}")

# Add timestamp
print(f"Analysis performed on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Read the data files
print("\nReading data files...")
try:
    # Read files with the first row as column names
    idc_genomic = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/IDC_genomic.csv', 
                             encoding='latin-1', 
                             low_memory=False)
    ilc_genomic = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/ILC_genomic.csv', 
                             encoding='latin-1', 
                             low_memory=False)
except FileNotFoundError as e:
    print(f"Error: An input file was not found. {e}")
    sys.exit(1)

# Print column names for debugging
print("\nAvailable columns in IDC genomic data:")
for col in idc_genomic.columns:
    print(f"'{col}'")

print("\nAvailable columns in ILC genomic data:")
for col in ilc_genomic.columns:
    print(f"'{col}'")

# Filter for Vus = FALSE
idc_genomic = idc_genomic[idc_genomic['Vus'] == "FALSE"]
ilc_genomic = ilc_genomic[ilc_genomic['Vus'] == "FALSE"]

# Analyze both cohorts
analyze_mutations(idc_genomic, "IDC")
analyze_mutations(ilc_genomic, "ILC")

# Save detailed mutation data
print("\nGenerating detailed mutation report...")

def generate_mutation_details(df, cohort):
    # Get mutation counts by gene
    gene_counts = df['Gene'].value_counts().reset_index()
    gene_counts.columns = ['Gene', 'Total_Mutations']
    gene_counts['Cohort'] = cohort
    gene_counts['Unique_Patients'] = df.groupby('Gene')['Effective Patient ID'].nunique().reindex(gene_counts['Gene']).values
    return gene_counts

idc_details = generate_mutation_details(idc_genomic, "IDC")
ilc_details = generate_mutation_details(ilc_genomic, "ILC")

# Combine and save details
all_details = pd.concat([idc_details, ilc_details])
all_details = all_details.sort_values(['Cohort', 'Total_Mutations'], ascending=[True, False])
all_details.to_csv('mutation_details.csv', index=False)
print("Detailed mutation report saved to: Vus_filtering/mutation_details.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('mutation_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: Vus_filtering/mutation_analysis_report.txt") 