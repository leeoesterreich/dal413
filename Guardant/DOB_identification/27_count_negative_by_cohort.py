import pandas as pd
import sys

def count_negative_by_cohort(serial_file, ilc_file, idc_file):
    """
    Counts the number of patients in ILC and IDC cohorts who have no
    alterations detected across all their tests.

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

    # Get patient ID sets for each cohort
    ilc_patient_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_patient_ids = set(df_idc['Effective Patient ID'].unique())

    # Find patients who have at least one 'True' detection
    patients_with_true = set(df_serial[df_serial['Aleration Detected?'] == True]['Effective Patient ID'].unique())

    # Find patients in each cohort who are NOT in the 'patients_with_true' set
    ilc_all_negative = ilc_patient_ids - patients_with_true
    idc_all_negative = idc_patient_ids - patients_with_true

    print("--- Report on Patients Without Alteration Detections ---")
    print(f"Number of ILC patients with no alterations detected: {len(ilc_all_negative)}")
    if len(ilc_all_negative) > 0:
        print("  ILC Patient IDs with no alterations:", sorted(list(ilc_all_negative)))

    print(f"\nNumber of IDC patients with no alterations detected: {len(idc_all_negative)}")
    if len(idc_all_negative) > 0:
        print("  IDC Patient IDs with no alterations:", sorted(list(idc_all_negative)))

    print("\nAnalysis complete.")


if __name__ == "__main__":
    SERIAL_GENOMIC_FILE = 'serial_genomic_with_patient_id.csv'
    ILC_COHORT_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_COHORT_FILE = 'idc_cohort_with_patient_id.csv'
    
    count_negative_by_cohort(SERIAL_GENOMIC_FILE, ILC_COHORT_FILE, IDC_COHORT_FILE) 