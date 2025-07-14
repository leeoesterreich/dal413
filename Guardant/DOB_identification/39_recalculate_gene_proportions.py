import pandas as pd
import os

def recalculate_proportions(idc_path, ilc_path, output_path):
    """
    Calculates the number and percentage of patients with mutations in key genes,
    relative to the TOTAL number of patients in each cohort, and saves to a CSV.

    Args:
        idc_path (str): Path to the IDC genomic data CSV.
        ilc_path (str): Path to the ILC genomic data CSV.
        output_path (str): Path to save the output summary CSV.
    """
    try:
        idc_df = pd.read_csv(idc_path, low_memory=False)
        ilc_df = pd.read_csv(ilc_path, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found: {e.filename}")
        return

    combined_df = pd.concat([idc_df, ilc_df], ignore_index=True)

    # Define the list of clinically relevant genes
    key_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']
    all_results = []

    def calculate_stats_for_cohort(df, cohort_name):
        """Helper function to perform calculation for a given cohort."""
        if df.empty:
            print(f"\nCohort '{cohort_name}' is empty or could not be processed.")
            return

        total_patients_in_cohort = df['GH_ID'].nunique()
        if total_patients_in_cohort == 0:
            print(f"\nCohort '{cohort_name}' has no patients.")
            return
            
        print(f"\n--- {cohort_name} Cohort Analysis (Total Patients = {total_patients_in_cohort}) ---")

        # Function to check if a gene from the dataframe matches any in our list
        def contains_gene(gene_string, target_gene):
            if pd.isna(gene_string):
                return False
            # Handles exact matches (e.g., 'BRCA1') and fusions (e.g., 'BRCA1-...')
            return str(gene_string).startswith(target_gene)

        # Calculate stats for each gene
        for gene in key_genes:
            gene_mask = df['Gene'].apply(lambda x: contains_gene(x, gene))
            patients_with_gene_mutation = df[gene_mask]['GH_ID'].nunique()
            
            # Calculate percentage relative to the total cohort size
            percentage = (patients_with_gene_mutation / total_patients_in_cohort) * 100
            
            # Store results for CSV
            all_results.append({
                'Cohort': cohort_name,
                'Total Patients in Cohort': total_patients_in_cohort,
                'Gene': gene,
                '# Patients with Mutation': patients_with_gene_mutation,
                '% of Total Patients in Cohort': f"{percentage:.1f}%"
            })
            print(f"- {gene}: {patients_with_gene_mutation} patients ({percentage:.1f}%)")

    # Run the analysis for each cohort
    calculate_stats_for_cohort(idc_df, "IDC")
    calculate_stats_for_cohort(ilc_df, "ILC")
    calculate_stats_for_cohort(combined_df, "Combined Total")
    
    # Save the results to a CSV file
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_path, index=False)
        print(f"\nResults successfully saved to '{output_path}'")
    else:
        print("\nNo results were generated to save.")


if __name__ == "__main__":
    IDC_FILE = "final_genomic_cohorts/IDC_genomic.csv"
    ILC_FILE = "final_genomic_cohorts/ILC_genomic.csv"
    OUTPUT_FILE = "final_genomic_cohorts/gene_proportions_by_cohort.csv"
    recalculate_proportions(IDC_FILE, ILC_FILE, OUTPUT_FILE) 