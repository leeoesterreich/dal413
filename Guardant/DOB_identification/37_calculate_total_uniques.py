import pandas as pd
import os

def calculate_total_uniques(idc_path, ilc_path):
    """
    Calculates the total unique mutations and genes from the ILC and IDC cohort files.

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

    # Combine the two cohorts into a single dataframe
    combined_df = pd.concat([idc_df, ilc_df], ignore_index=True)
    print(f"Total rows in combined data: {len(combined_df)}")
    print(f"Total unique patients in combined data: {combined_df['GH_ID'].nunique()}")

    # 1. Calculate Unique Mutations Identified (unique Gene + Alteration pairs)
    unique_mutations_count = combined_df[['Gene', 'Alteration']].drop_duplicates().shape[0]

    # 2. Calculate Unique Genes Mutated
    unique_genes_count = combined_df['Gene'].nunique()

    # --- Print Summary ---
    print("\\n--- Total Cohort Analysis ---")
    print(f"Unique Mutations Identified: {unique_mutations_count}")
    print(f"Unique Genes Mutated:        {unique_genes_count}")

if __name__ == "__main__":
    IDC_FILE = "final_genomic_cohorts/IDC_genomic.csv"
    ILC_FILE = "final_genomic_cohorts/ILC_genomic.csv"
    calculate_total_uniques(IDC_FILE, ILC_FILE) 