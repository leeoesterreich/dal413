import pandas as pd
import numpy as np
from collections import defaultdict
from io import StringIO
import sys

# Create string buffers to capture output
infinity_summary_buffer = StringIO()
infinity_detail_buffer = StringIO()

# Define known BRCA fusions
BRCA_FUSIONS = {
    'BRCA1-B4GALNT2', 
    'BRCA1-IGR', 
    'BRCA1-ZNF563', 
    'BRCA2-STARD13'
}

# Define known PALB2 fusions
PALB2_FUSIONS = {
    'PALB2-IGR',
    'PALB2-LONP2',
    'PALB2-SHISA9'
}

def analyze_test_counts(df, cohort_name, output_buffer):
    # Group by patient ID
    test_counts = df.groupby('Effective Patient ID')['Total Test'].first()
    
    # Basic statistics
    total_patients = len(test_counts)
    total_tests = test_counts.sum()
    avg_tests = total_tests / total_patients if total_patients > 0 else 0
    
    # Patients with multiple tests
    multiple_tests = test_counts[test_counts > 1]
    multiple_test_median = multiple_tests.median() if len(multiple_tests) > 0 else 0
    
    print(f"\n{cohort_name} Infinity Test Count Analysis:", file=output_buffer)
    print(f"Total unique patients: {total_patients}", file=output_buffer)
    print(f"Total tests: {total_tests}", file=output_buffer)
    print(f"Average tests per patient: {avg_tests:.2f}", file=output_buffer)
    print(f"Patients with multiple tests: {len(multiple_tests)} ({(len(multiple_tests)/total_patients*100):.1f}%)", file=output_buffer)
    if len(multiple_tests) > 0:
        print(f"For patients with multiple tests - median number of tests: {multiple_test_median}", file=output_buffer)
    
    return total_patients, total_tests, len(multiple_tests)

def analyze_mutations(df, cohort_name, summary_buffer, detail_buffer):
    # Count unique Gene-Alteration combinations
    gene_alt_combinations = df[['Gene', 'Alteration']].drop_duplicates()
    unique_gene_alt_count = len(gene_alt_combinations)
    
    # Count unique genes with mutations
    unique_genes = df['Gene'].dropna().unique()
    unique_gene_count = len(unique_genes)
    
    def contains_gene(gene_string, target_gene):
        if pd.isna(gene_string):
            return False
        gene_str = str(gene_string)
        return target_gene == gene_str or gene_str.startswith(f"{target_gene}-") or gene_str.startswith(f"{target_gene}_")
    
    def is_fusion(gene_string, target_gene):
        if pd.isna(gene_string):
            return False
        gene_str = str(gene_string)
        return gene_str.startswith(f"{target_gene}-") or gene_str.startswith(f"{target_gene}_")
    
    # Get unique patients with any mutation in key genes
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    patients_with_mutations = set()
    gene_patient_counts = {}
    gene_to_patients = {}
    patient_to_genes = defaultdict(set)
    
    for gene in key_genes:
        gene_mask = df['Gene'].apply(lambda x: contains_gene(x, gene))
        gene_data = df[gene_mask]
        if len(gene_data) > 0:
            gene_patients = set(gene_data['Effective Patient ID'].unique())
            gene_patient_counts[gene] = len(gene_patients)
            gene_to_patients[gene] = gene_patients
            patients_with_mutations.update(gene_patients)
            for patient in gene_patients:
                patient_to_genes[patient].add(gene)
    
    # Analyze mutation counts per patient
    mutation_count_distribution = defaultdict(int)
    for patient, genes in patient_to_genes.items():
        mutation_count_distribution[len(genes)] += 1
    
    # Print summary
    print(f"\n{cohort_name} Infinity Patient Summary:", file=summary_buffer)
    print(f"Total unique patients with any mutation in key genes: {len(patients_with_mutations)}", file=summary_buffer)
    print("\nNumber of genes with mutations per patient:", file=summary_buffer)
    for num_genes in sorted(mutation_count_distribution.keys()):
        count = mutation_count_distribution[num_genes]
        percentage = (count / len(patients_with_mutations)) * 100
        print(f"- {count} patients ({percentage:.1f}%) have mutations in {num_genes} gene{'s' if num_genes != 1 else ''}", file=summary_buffer)
    
    print("\nPatients per gene:", file=summary_buffer)
    for gene in key_genes:
        if gene in gene_patient_counts:
            percentage = (gene_patient_counts[gene] / len(patients_with_mutations)) * 100
            print(f"- {gene}: {gene_patient_counts[gene]} patients ({percentage:.1f}% of patients with mutations)", file=summary_buffer)
    
    # Analyze co-mutations
    print("\nCommon co-occurring mutations:", file=summary_buffer)
    for i, gene1 in enumerate(key_genes):
        if gene1 not in gene_to_patients:
            continue
        for gene2 in key_genes[i+1:]:
            if gene2 not in gene_to_patients:
                continue
            common_patients = gene_to_patients[gene1].intersection(gene_to_patients[gene2])
            if len(common_patients) > 0:
                percentage = (len(common_patients) / len(patients_with_mutations)) * 100
                print(f"- {gene1} + {gene2}: {len(common_patients)} patients ({percentage:.1f}% of patients with mutations)", file=summary_buffer)
    
    # Analyze mutation types
    print(f"\n{cohort_name} Mutation Type Analysis:", file=detail_buffer)
    for gene in key_genes:
        gene_mask = df['Gene'].apply(lambda x: contains_gene(x, gene))
        gene_data = df[gene_mask]
        
        if len(gene_data) == 0:
            continue
            
        print(f"\n{gene} Analysis:", file=detail_buffer)
        
        # Count mutation types
        mutation_types = defaultdict(int)
        for _, row in gene_data.iterrows():
            mut_type = row['Type'] if pd.notna(row['Type']) else 'unknown'
            mutation_types[mut_type] += 1
        
        print("Mutation types:", file=detail_buffer)
        for mut_type, count in sorted(mutation_types.items()):
            print(f"- {mut_type}: {count}", file=detail_buffer)
        
        # List all unique variants
        unique_variants = gene_data[['Alteration', 'Type']].drop_duplicates()
        print(f"\nUnique variants ({len(unique_variants)}):", file=detail_buffer)
        for _, variant in unique_variants.iterrows():
            print(f"- {variant['Alteration']} ({variant['Type']})", file=detail_buffer)

def main():
    # Read the Infinity data
    idc_infinity = pd.read_csv('IDC_genomic_infinity_corrected.csv')
    ilc_infinity = pd.read_csv('ILC_genomic_infinity_corrected.csv')
    
    # Analyze IDC Infinity data
    analyze_test_counts(idc_infinity, "IDC", infinity_summary_buffer)
    analyze_mutations(idc_infinity, "IDC", infinity_summary_buffer, infinity_detail_buffer)
    
    # Analyze ILC Infinity data
    analyze_test_counts(ilc_infinity, "ILC", infinity_summary_buffer)
    analyze_mutations(ilc_infinity, "ILC", infinity_summary_buffer, infinity_detail_buffer)
    
    # Save the summaries to files
    with open('infinity_analysis_summary_corrected.txt', 'w') as f:
        f.write(infinity_summary_buffer.getvalue())
    
    with open('infinity_detailed_mutation_analysis_corrected.txt', 'w') as f:
        f.write(infinity_detail_buffer.getvalue())
    
    print("Analysis complete. Results saved to infinity_analysis_summary_corrected.txt and infinity_detailed_mutation_analysis_corrected.txt")

if __name__ == "__main__":
    main() 