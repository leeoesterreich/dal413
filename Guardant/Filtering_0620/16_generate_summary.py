import pandas as pd
import numpy as np
import datetime
from io import StringIO
import sys
from collections import defaultdict

print("Script starting...")  # Debug print

# Create string buffers to capture output
summary_buffer = StringIO()
detail_buffer = StringIO()
original_stdout = sys.stdout

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
    # Group by patient ID and get their Total Test values
    gh_id_col = 'GH_ID'
    test_counts = df.groupby(gh_id_col)['Total Test'].first()
    
    # Basic statistics
    total_patients = len(test_counts)
    total_tests = test_counts.sum()
    avg_tests = total_tests / total_patients if total_patients > 0 else 0
    
    # Patients with multiple tests
    multiple_tests = test_counts[test_counts > 1]
    multiple_test_median = multiple_tests.median() if len(multiple_tests) > 0 else 0
    
    print(f"\n{cohort_name} Test Count Analysis:", file=output_buffer)
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
        # Check if it's a fusion pattern (gene-something or gene_something)
        return gene_str.startswith(f"{target_gene}-") or gene_str.startswith(f"{target_gene}_")
    
    # Get unique patients with any mutation in key genes
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    patients_with_mutations = set()
    gene_patient_counts = {}
    gene_to_patients = {}  # Dictionary to store patients for each gene
    patient_to_genes = defaultdict(set)  # Dictionary to store genes for each patient
    
    for gene in key_genes:
        gene_mask = df['Gene'].apply(lambda x: contains_gene(x, gene))
        gene_data = df[gene_mask]
        if len(gene_data) > 0:
            gene_patients = set(gene_data['GH_ID'].unique())
            gene_patient_counts[gene] = len(gene_patients)
            gene_to_patients[gene] = gene_patients
            patients_with_mutations.update(gene_patients)
            # Record which genes each patient has mutations in
            for patient in gene_patients:
                patient_to_genes[patient].add(gene)
    
    # Analyze how many genes each patient has mutations in
    mutation_count_distribution = defaultdict(int)
    for patient, genes in patient_to_genes.items():
        mutation_count_distribution[len(genes)] += 1
    
    # Print summary of patients with any mutation
    print(f"\n{cohort_name} Patient Summary:", file=summary_buffer)
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
    
    # Analyze common co-mutations
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
    
    # Function to analyze mutations for a specific gene
    def analyze_gene_mutations(gene):
        # Get all rows for this gene
        gene_mask = df['Gene'].apply(lambda x: contains_gene(x, gene))
        gene_data = df[gene_mask]
        
        if len(gene_data) == 0:
            return None
        
        # Count unique patients
        unique_patients = len(gene_data['GH_ID'].unique())
        
        # First identify fusions by looking at Gene column patterns
        fusions = []
        fusion_patients = set()
        
        for _, row in gene_data.iterrows():
            gene_name = str(row['Gene'])
            # Check if it's a fusion pattern in the Gene column
            if is_fusion(gene_name, gene):
                fusions.append({
                    'name': gene_name,
                    'patient_id': row['GH_ID']
                })
                fusion_patients.add(row['GH_ID'])
        
        # Group mutations by type and count patients
        mutation_counts = defaultdict(lambda: {'count': 0, 'patients': set()})
        
        for _, row in gene_data.iterrows():
            gene_name = str(row['Gene'])
            mutation_type = row['Type'] if pd.notna(row['Type']) else 'unknown'
            alteration = row['Alteration'] if pd.notna(row['Alteration']) else 'unknown'
            patient_id = row['GH_ID']
            
            # Skip if it's a fusion as we've already counted it
            if is_fusion(gene_name, gene):
                continue
                
            key = f"{alteration} ({mutation_type})"
            mutation_counts[key]['count'] += 1
            mutation_counts[key]['patients'].add(patient_id)
        
        # Convert to list of tuples and sort by count
        mutation_list = [(k, len(v['patients']), v['count']) for k, v in mutation_counts.items()]
        mutation_list.sort(key=lambda x: x[1], reverse=True)
        
        # Count mutation types
        mutation_types = defaultdict(set)
        for _, row in gene_data.iterrows():
            gene_name = str(row['Gene'])
            if is_fusion(gene_name, gene):
                mutation_types['fusion'].add(gene_name)
            else:
                mut_type = row['Type'] if pd.notna(row['Type']) else 'unknown'
                mutation_types[mut_type].add(f"{row['Gene']}_{row['Alteration']}")
        
        return {
            'unique_patients': unique_patients,
            'mutations': mutation_list,
            'mutation_types': {k: len(v) for k, v in mutation_types.items()},
            'fusions': fusions,
            'fusion_patients': len(fusion_patients)
        }
    
    # Analyze key genes
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    gene_analyses = {}
    
    print(f"\n{cohort_name} Mutation Analysis:", file=summary_buffer)
    print(f"Total unique gene-alteration combinations: {unique_gene_alt_count}", file=summary_buffer)
    print(f"Total unique genes with mutations: {unique_gene_count}", file=summary_buffer)
    
    print(f"\n{cohort_name} Mutation Analysis:", file=detail_buffer)
    print(f"Total unique gene-alteration combinations: {unique_gene_alt_count}", file=detail_buffer)
    print(f"Total unique genes with mutations: {unique_gene_count}", file=detail_buffer)
    
    # Print detailed mutation analysis for each key gene
    for gene in key_genes:
        analysis = analyze_gene_mutations(gene)
        if analysis:
            gene_analyses[gene] = analysis
            
            # Summary output (vital stats only)
            print(f"\n{gene} Summary:", file=summary_buffer)
            print(f"Total patients: {analysis['unique_patients']}", file=summary_buffer)
            if analysis['fusions']:
                print(f"Fusion mutations: {analysis['fusion_patients']} patients", file=summary_buffer)
                for fusion in analysis['fusions']:
                    print(f"- {fusion['name']}", file=summary_buffer)
            print("Mutation types:", file=summary_buffer)
            for mut_type, count in analysis['mutation_types'].items():
                print(f"- {mut_type}: {count} variants", file=summary_buffer)
            
            # Detailed output
            print(f"\n{gene} Analysis:", file=detail_buffer)
            print(f"Total patients with {gene} mutations: {analysis['unique_patients']}", file=detail_buffer)
            
            if analysis['fusions']:
                print(f"\nFusion mutations ({analysis['fusion_patients']} patients):", file=detail_buffer)
                for fusion in analysis['fusions']:
                    print(f"- {fusion['name']}", file=detail_buffer)
            
            print("\nMutation type counts:", file=detail_buffer)
            for mut_type, count in analysis['mutation_types'].items():
                print(f"- {mut_type}: {count} different variants", file=detail_buffer)
            
            print("\nDetailed mutations (by patient count):", file=detail_buffer)
            for mutation, patient_count, total_count in analysis['mutations']:
                print(f"- {mutation}: {patient_count} patients ({total_count} total occurrences)", file=detail_buffer)
    
    return gene_analyses, unique_gene_alt_count, unique_gene_count

def main():
    # Read the data
    idc_df = pd.read_csv('IDC_genomic_with_histology.csv')
    ilc_df = pd.read_csv('ILC_genomic_with_histology.csv')
    
    # Analyze test counts
    print("\nAnalyzing test counts...", file=summary_buffer)
    idc_patients, idc_tests, idc_multiple = analyze_test_counts(idc_df, "IDC", summary_buffer)
    ilc_patients, ilc_tests, ilc_multiple = analyze_test_counts(ilc_df, "ILC", summary_buffer)
    
    # Print combined test statistics
    print("\nCombined Test Statistics:", file=summary_buffer)
    print(f"Total unique patients: {idc_patients + ilc_patients}", file=summary_buffer)
    print(f"Total tests: {idc_tests + ilc_tests}", file=summary_buffer)
    print(f"Total patients with multiple tests: {idc_multiple + ilc_multiple}", file=summary_buffer)
    
    # Analyze mutations
    print("\nAnalyzing mutations...", file=summary_buffer)
    idc_analysis, idc_combinations, idc_genes = analyze_mutations(idc_df, "IDC", summary_buffer, detail_buffer)
    ilc_analysis, ilc_combinations, ilc_genes = analyze_mutations(ilc_df, "ILC", summary_buffer, detail_buffer)
    
    # Print combined mutation statistics
    print("\nCombined Mutation Statistics:", file=summary_buffer)
    print(f"Total unique gene-alteration combinations: {idc_combinations + ilc_combinations}", file=summary_buffer)
    all_genes = set(idc_df['Gene'].dropna()).union(set(ilc_df['Gene'].dropna()))
    print(f"Total unique genes: {len(all_genes)}", file=summary_buffer)
    
    # Add mutation type summary to detailed analysis
    print("\n" + "="*50, file=detail_buffer)
    print("UNIQUE VARIANT COUNTS BY MUTATION TYPE", file=detail_buffer)
    print("="*50, file=detail_buffer)
    
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    
    total_variants = 0
    
    for gene in key_genes:
        print(f"\n{gene}:", file=detail_buffer)
        
        # Count unique variants in IDC
        idc_variants = defaultdict(set)
        idc_mask = idc_df['Gene'].apply(lambda x: str(x).startswith(gene))
        idc_data = idc_df[idc_mask]
        for _, row in idc_data.iterrows():
            if pd.notna(row['Type']) and pd.notna(row['Alteration']):
                idc_variants[row['Type']].add(row['Alteration'])
        
        # Count unique variants in ILC
        ilc_variants = defaultdict(set)
        ilc_mask = ilc_df['Gene'].apply(lambda x: str(x).startswith(gene))
        ilc_data = ilc_df[ilc_mask]
        for _, row in ilc_data.iterrows():
            if pd.notna(row['Type']) and pd.notna(row['Alteration']):
                ilc_variants[row['Type']].add(row['Alteration'])
        
        # Combine variants
        all_types = set(idc_variants.keys()).union(set(ilc_variants.keys()))
        gene_total = 0
        
        print("  IDC:", file=detail_buffer)
        idc_total = 0
        for mut_type in sorted(all_types):
            if len(idc_variants[mut_type]) > 0:
                print(f"    - {mut_type}: {len(idc_variants[mut_type])} unique variants", file=detail_buffer)
                idc_total += len(idc_variants[mut_type])
        print(f"    Total: {idc_total} unique variants", file=detail_buffer)
        
        print("  ILC:", file=detail_buffer)
        ilc_total = 0
        for mut_type in sorted(all_types):
            if len(ilc_variants[mut_type]) > 0:
                print(f"    - {mut_type}: {len(ilc_variants[mut_type])} unique variants", file=detail_buffer)
                ilc_total += len(ilc_variants[mut_type])
        print(f"    Total: {ilc_total} unique variants", file=detail_buffer)
        
        # Count combined unique variants
        combined_variants = defaultdict(set)
        for mut_type in all_types:
            combined_variants[mut_type] = idc_variants[mut_type].union(ilc_variants[mut_type])
        
        print("  Combined:", file=detail_buffer)
        combined_total = 0
        for mut_type in sorted(all_types):
            if len(combined_variants[mut_type]) > 0:
                print(f"    - {mut_type}: {len(combined_variants[mut_type])} unique variants", file=detail_buffer)
                combined_total += len(combined_variants[mut_type])
        print(f"    Total: {combined_total} unique variants", file=detail_buffer)
        
        total_variants += combined_total
    
    print(f"\nTotal unique variants across all genes: {total_variants}", file=detail_buffer)
    
    print("\n" + "="*50, file=detail_buffer)
    print("DETAILED MUTATION COUNTS", file=detail_buffer)
    print("="*50, file=detail_buffer)
    
    for gene in key_genes:
        print(f"\n{gene} Mutation Types:", file=detail_buffer)
        
        # Get mutation types for IDC
        idc_types = defaultdict(int)
        idc_mask = idc_df['Gene'].apply(lambda x: str(x).startswith(gene))
        idc_data = idc_df[idc_mask]
        for _, row in idc_data.iterrows():
            if pd.notna(row['Type']):
                idc_types[row['Type']] += 1
        
        # Get mutation types for ILC
        ilc_types = defaultdict(int)
        ilc_mask = ilc_df['Gene'].apply(lambda x: str(x).startswith(gene))
        ilc_data = ilc_df[ilc_mask]
        for _, row in ilc_data.iterrows():
            if pd.notna(row['Type']):
                ilc_types[row['Type']] += 1
        
        # Combine and get all unique types
        all_types = set(list(idc_types.keys()) + list(ilc_types.keys()))
        
        print(f"  IDC ({len(idc_data)} total mutations):", file=detail_buffer)
        for type_name in sorted(all_types):
            if idc_types[type_name] > 0:
                print(f"    - {type_name}: {idc_types[type_name]} mutations", file=detail_buffer)
        
        print(f"  ILC ({len(ilc_data)} total mutations):", file=detail_buffer)
        for type_name in sorted(all_types):
            if ilc_types[type_name] > 0:
                print(f"    - {type_name}: {ilc_types[type_name]} mutations", file=detail_buffer)
        
        print(f"  Combined ({len(idc_data) + len(ilc_data)} total mutations):", file=detail_buffer)
        for type_name in sorted(all_types):
            total = idc_types[type_name] + ilc_types[type_name]
            if total > 0:
                print(f"    - {type_name}: {total} mutations", file=detail_buffer)
    
    # Write to files
    with open('combined_analysis_summary.txt', 'w') as f:
        f.write(summary_buffer.getvalue())
    
    with open('detailed_mutation_analysis.txt', 'w') as f:
        f.write(detail_buffer.getvalue())

if __name__ == "__main__":
    main() 