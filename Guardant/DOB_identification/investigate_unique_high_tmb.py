import pandas as pd
import os

def investigate_unique_high_tmb(idc_file):
    """
    Finds unique IDC patients with a TMB score greater than 100,
    handling duplicates by taking the first test record for each patient.

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

    # Ensure TMB score is numeric and drop rows without a valid score
    idc_df[tmb_col] = pd.to_numeric(idc_df[tmb_col], errors='coerce')
    valid_tmb_df = idc_df.dropna(subset=[tmb_col])

    # --- This is the key step from the reference script ---
    # De-duplicate to ensure one record per unique patient (keeping the first)
    unique_patient_df = valid_tmb_df.drop_duplicates(subset=[patient_col], keep='first')

    # Now, filter for high TMB scores from the unique patient list
    high_tmb_threshold = 100
    high_tmb_unique_patients_df = unique_patient_df[unique_patient_df[tmb_col] > high_tmb_threshold]

    if high_tmb_unique_patients_df.empty:
        print(f"No unique IDC patients found with a TMB score greater than {high_tmb_threshold}.")
        return

    print(f"Found {len(high_tmb_unique_patients_df)} unique patients with TMB Score > {high_tmb_threshold}:\n")
    print("Effective Patient ID | TMB Score")
    print("-------------------- | -----------")
    
    # Sort by TMB score for clarity
    sorted_high_tmb = high_tmb_unique_patients_df.sort_values(by=tmb_col, ascending=False)
    
    for index, row in sorted_high_tmb.iterrows():
        print(f"{row[patient_col]:<20} | {row[tmb_col]:.2f}")


if __name__ == '__main__':
    base_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts'
    idc_infinity_file = os.path.join(base_path, 'IDC_genomic_infinity.csv')
    investigate_unique_high_tmb(idc_infinity_file) 