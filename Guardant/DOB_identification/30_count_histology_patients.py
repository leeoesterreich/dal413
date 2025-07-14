import pandas as pd
import sys

def count_histology_patients(ilc_file, idc_file):
    """
    Counts the unique patients in the ILC and IDC cohorts and calculates
    their overlap and total unique count.

    Args:
        ilc_file (str): Path to the ILC cohort data.
        idc_file (str): Path to the IDC cohort data.
    """
    try:
        df_ilc = pd.read_csv(ilc_file, low_memory=False)
        df_idc = pd.read_csv(idc_file, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # Get unique patient ID sets from each cohort
    ilc_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_ids = set(df_idc['Effective Patient ID'].unique())

    # Calculate counts
    num_ilc = len(ilc_ids)
    num_idc = len(idc_ids)
    num_overlap = len(ilc_ids.intersection(idc_ids))
    total_unique = len(ilc_ids.union(idc_ids))

    # --- Print Report ---
    print("--- Histology Cohort Patient Count ---")
    print(f"Unique patients in ILC cohort: {num_ilc}")
    print(f"Unique patients in IDC cohort: {num_idc}")
    print(f"Patients present in BOTH cohorts (Overlap): {num_overlap}")
    print("-" * 38)
    print(f"Total unique patients across both cohorts: {total_unique}")
    print(f"(Calculation: {num_ilc} + {num_idc} - {num_overlap} = {total_unique})")

if __name__ == "__main__":
    ILC_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_FILE = 'idc_cohort_with_patient_id.csv'
    
    count_histology_patients(ILC_FILE, IDC_FILE) 