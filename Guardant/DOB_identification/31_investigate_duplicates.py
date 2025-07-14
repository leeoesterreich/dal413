import pandas as pd
import sys

def investigate_duplicates(ilc_file, idc_file):
    """
    Investigates duplicate patient entries within the ILC and IDC cohort files.

    Args:
        ilc_file (str): Path to the ILC cohort data.
        idc_file (str): Path to the IDC cohort data.
    """
    print("--- Duplicate Patient ID Investigation ---")

    def analyze_file(file_path, cohort_name):
        try:
            df = pd.read_csv(file_path, low_memory=False)
        except FileNotFoundError:
            print(f"\nError: {cohort_name} file not found at {file_path}")
            return

        total_rows = len(df)
        unique_patients = df['Effective Patient ID'].nunique()

        print(f"\n--- Analysis for {cohort_name} Cohort ---")
        print(f"Total number of rows: {total_rows}")
        print(f"Number of unique patients: {unique_patients}")

        if total_rows > unique_patients:
            duplicates = df[df.duplicated(subset=['Effective Patient ID'], keep=False)]
            duplicate_counts = duplicates['Effective Patient ID'].value_counts()
            
            print(f"Discrepancy found: {total_rows - unique_patients} extra row(s) corresponding to {len(duplicate_counts)} patient(s).")
            print("Patients with multiple records:")
            for patient_id, count in duplicate_counts.items():
                print(f"  - Patient ID: {patient_id} appears {count} times.")
        else:
            print("No duplicate patient entries found in this cohort.")

    analyze_file(ilc_file, "ILC")
    analyze_file(idc_file, "IDC")


if __name__ == "__main__":
    ILC_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_FILE = 'idc_cohort_with_patient_id.csv'
    
    investigate_duplicates(ILC_FILE, IDC_FILE) 