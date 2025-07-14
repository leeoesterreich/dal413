import pandas as pd
import os

def analyze_infinity_proportions(idc_path, ilc_path, output_path):
    """
    Calculates patient proportions for key gene mutations in the Infinity cohort,
    relative to the total number of patients in that cohort.

    Args:
        idc_path (str): Path to the IDC Infinity genomic data CSV.
        ilc_path (str): Path to the ILC Infinity genomic data CSV.
        output_path (str): Path to save the output summary CSV.
    """
    try:
        idc_df = pd.read_csv(idc_path, low_memory=False)
        ilc_df = pd.read_csv(ilc_path, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: A required input file was not found. Please run script 17 first. Details: {e.filename}")
        return

    combined_df = pd.concat([idc_df, ilc_df], ignore_index=True)
    all_results = []

    # Define gene groups in the main function scope
    gene_groups = {
        'Clinically Relevant Gene Mutations': ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2'],
        'PIK3CA': ['PIK3CA'],
        'ESR1': ['ESR1'],
        'ERBB2': ['ERBB2'],
        'BRCA 1/2': ['BRCA1', 'BRCA2'],
        'PALB2': ['PALB2']
    }

    def calculate_stats_for_cohort(df, cohort_name):
        """Helper function to perform calculation for a given cohort."""
        if df.empty: return
        total_patients_in_cohort = df['GH_ID'].nunique()
        if total_patients_in_cohort == 0: return

        print(f"\\n--- {cohort_name} Infinity Cohort Analysis (Total Patients = {total_patients_in_cohort}) ---")

        def get_patient_count_for_genes(gene_list):
            """Gets the unique patient count for a list of genes."""
            mask = df['Gene'].apply(lambda x: any(str(x).startswith(g) for g in gene_list))
            return df[mask]['GH_ID'].nunique()

        for name, genes in gene_groups.items():
            patient_count = get_patient_count_for_genes(genes)
            percentage = (patient_count / total_patients_in_cohort) * 100
            
            all_results.append({
                'Cohort': cohort_name,
                'Metric': name,
                '# Patients with Mutation': patient_count,
                '% of Total Patients in Cohort': f"{percentage:.1f}%"
            })
            print(f"- {name}: {patient_count} patients ({percentage:.1f}%)")

    # Run analysis
    calculate_stats_for_cohort(idc_df, "IDC")
    calculate_stats_for_cohort(ilc_df, "ILC")
    calculate_stats_for_cohort(combined_df, "Total")
    
    # Save results to CSV
    if all_results:
        results_df = pd.DataFrame(all_results)
        # Pivot the table to get the desired format
        pivot_df = results_df.pivot_table(
            index='Metric', 
            columns='Cohort', 
            values=['# Patients with Mutation', '% of Total Patients in Cohort'],
            aggfunc='first'
        )
        
        # Reorder and format columns to match the desired output
        pivot_df = pivot_df[[
            ('# Patients with Mutation', 'IDC'), ('% of Total Patients in Cohort', 'IDC'),
            ('# Patients with Mutation', 'ILC'), ('% of Total Patients in Cohort', 'ILC'),
            ('# Patients with Mutation', 'Total'), ('% of Total Patients in Cohort', 'Total')
        ]]
        
        # Reorder rows to match the desired output
        pivot_df = pivot_df.reindex(gene_groups.keys())
        
        pivot_df.to_csv(output_path)
        print(f"\\nResults successfully saved to '{output_path}'")

if __name__ == "__main__":
    IDC_FILE = "final_genomic_cohorts/IDC_genomic_infinity.csv"
    ILC_FILE = "final_genomic_cohorts/ILC_genomic_infinity.csv"
    OUTPUT_FILE = "final_genomic_cohorts/infinity_gene_proportions_summary.csv"
    analyze_infinity_proportions(IDC_FILE, ILC_FILE, OUTPUT_FILE) 