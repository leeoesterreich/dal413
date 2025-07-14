import pandas as pd
import os

def analyze_clinical_relevance(idc_path, ilc_path):
    """
    Calculates the number and percentage of patients with clinically relevant mutations
    for IDC, ILC, and a combined total cohort.

    Args:
        idc_path (str): Path to the IDC genomic data CSV.
        ilc_path (str): Path to the ILC genomic data CSV.
    """
    try:
        idc_df = pd.read_csv(idc_path, low_memory=False)
        ilc_df = pd.read_csv(ilc_path, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found: {e.filename}")
        return

    combined_df = pd.concat([idc_df, ilc_df], ignore_index=True)

    # Define the list of clinically relevant genes
    relevant_genes = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']

    def calculate_stats_for_cohort(df, cohort_name):
        """Helper function to perform calculation for a given cohort dataframe."""
        if df.empty:
            print(f"\\nCohort '{cohort_name}' is empty or could not be processed.")
            return

        total_patients_in_cohort = df['GH_ID'].nunique()
        if total_patients_in_cohort == 0:
            print(f"\\nCohort '{cohort_name}' has no patients.")
            return
            
        # Function to check if a gene from the dataframe matches any in our list
        # This handles exact matches (e.g., 'BRCA1') and fusions (e.g., 'BRCA1-...')
        def is_relevant(gene_string):
            if pd.isna(gene_string):
                return False
            return any(str(gene_string).startswith(g) for g in relevant_genes)

        # Filter the dataframe to only rows with clinically relevant mutations
        relevant_mutations_df = df[df['Gene'].apply(is_relevant)]

        # Find the number of unique patients who have these mutations
        patients_with_relevant_mutation = relevant_mutations_df['GH_ID'].nunique()

        # Calculate the percentage
        percentage = (patients_with_relevant_mutation / total_patients_in_cohort) * 100

        # --- Print Results for the Cohort ---
        print(f"\\n--- {cohort_name} Cohort Analysis ---")
        print(f"Total Patients: {total_patients_in_cohort}")
        print(f"Patients with Clinically Relevant Gene Mutations: {patients_with_relevant_mutation}")
        print(f"Percentage: {percentage:.1f}%")

    # Run the analysis for each cohort
    calculate_stats_for_cohort(idc_df, "IDC")
    calculate_stats_for_cohort(ilc_df, "ILC")
    calculate_stats_for_cohort(combined_df, "Combined Total")


if __name__ == "__main__":
    IDC_FILE = "final_genomic_cohorts/IDC_genomic.csv"
    ILC_FILE = "final_genomic_cohorts/ILC_genomic.csv"
    analyze_clinical_relevance(IDC_FILE, ILC_FILE) 