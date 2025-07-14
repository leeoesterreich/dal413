import pandas as pd
import sys

def verify_patient_subsets(serial_file, ilc_file, idc_file):
    """
    Verifies if all patients from ILC and IDC cohorts are present in the main
    serial genomic data file.

    Args:
        serial_file (str): Path to the main serial genomic data CSV.
        ilc_file (str): Path to the ILC cohort CSV.
        idc_file (str): Path to the IDC cohort CSV.
    """
    try:
        df_serial = pd.read_csv(serial_file)
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # Get unique patient IDs from each file
    serial_patient_ids = set(df_serial['Effective Patient ID'].unique())
    ilc_patient_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_patient_ids = set(df_idc['Effective Patient ID'].unique())

    print("--- Verification Report ---")

    # Check for ILC cohort
    ilc_missing = ilc_patient_ids - serial_patient_ids
    if not ilc_missing:
        print("✅ All patients from the ILC cohort are present in the serial test data.")
    else:
        print(f"❌ Found {len(ilc_missing)} patients from the ILC cohort NOT present in the serial test data.")
        print("   Missing ILC Patient IDs:", sorted(list(ilc_missing)))

    # Check for IDC cohort
    idc_missing = idc_patient_ids - serial_patient_ids
    if not idc_missing:
        print("✅ All patients from the IDC cohort are present in the serial test data.")
    else:
        print(f"❌ Found {len(idc_missing)} patients from the IDC cohort NOT present in the serial test data.")
        print("   Missing IDC Patient IDs:", sorted(list(idc_missing)))

    print("\nVerification complete.")


if __name__ == "__main__":
    SERIAL_GENOMIC_FILE = 'serial_genomic_with_patient_id.csv'
    ILC_COHORT_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_COHORT_FILE = 'idc_cohort_with_patient_id.csv'
    
    verify_patient_subsets(SERIAL_GENOMIC_FILE, ILC_COHORT_FILE, IDC_COHORT_FILE) 