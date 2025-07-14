import pandas as pd
import os

def investigate_high_tmb(idc_file):
    """
    Investigates and prints details for IDC cases with high TMB scores.

    Args:
        idc_file (str): Path to the IDC genomic infinity CSV file.
    """
    try:
        idc_df = pd.read_csv(idc_file, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error loading file: {e}")
        return

    tmb_col = 'TMB Score'
    patient_col = 'Effective Patient ID'

    # Ensure TMB score is numeric
    idc_df[tmb_col] = pd.to_numeric(idc_df[tmb_col], errors='coerce')
    idc_df.dropna(subset=[tmb_col], inplace=True)

    # Filter for high TMB scores
    high_tmb_threshold = 100
    high_tmb_df = idc_df[idc_df[tmb_col] > high_tmb_threshold]

    if high_tmb_df.empty:
        print(f"No IDC cases found with a TMB score greater than {high_tmb_threshold}.")
        return

    print(f"Found {len(high_tmb_df)} tests from {high_tmb_df[patient_col].nunique()} unique patients with TMB Score > {high_tmb_threshold}:\n")
    print("Effective Patient ID | TMB Score")
    print("-------------------- | -----------")
    
    # Sort by TMB score for clarity
    sorted_high_tmb = high_tmb_df.sort_values(by=tmb_col, ascending=False)
    
    for index, row in sorted_high_tmb.iterrows():
        print(f"{row[patient_col]:<20} | {row[tmb_col]:.2f}")


if __name__ == '__main__':
    base_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts'
    idc_infinity_file = os.path.join(base_path, 'IDC_genomic_infinity.csv')
    investigate_high_tmb(idc_infinity_file) 