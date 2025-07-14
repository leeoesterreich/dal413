import pandas as pd
import sys

def investigate_patient(patient_id, serial_genomic_file, mapping_file):
    """
    Loads all test data for a specific patient and prints their
    serial testing information.
    """
    print(f"--- Investigating Patient: {patient_id} ---")
    try:
        serial_genomic_df = pd.read_csv(serial_genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except Exception as e:
        print(f"Error reading files: {e}", file=sys.stderr)
        return

    full_df = pd.merge(serial_genomic_df, mapping_df, on='GH_ID', how='left')
    patient_df = full_df[full_df['Effective Patient ID'] == patient_id]

    if patient_df.empty:
        print(f"No records found for patient {patient_id}.")
        return
        
    print(f"Found {len(patient_df)} total records for this patient.")
    
    unique_serial_tests = patient_df['Serial Test'].unique()
    unique_report_dates = patient_df['Final Report Date'].unique()
    
    print(f"Unique 'Serial Test' values found: {unique_serial_tests}")
    print(f"Unique 'Final Report Date' values found: {unique_report_dates}")
    print("-" * 30)

if __name__ == '__main__':
    serial_file = '../Guardant Project - Serial Genomic Data.csv'
    mapping_file = '../final_genomic_cohorts/complete_patient_mapping.csv'
    investigate_patient('P869', serial_file, mapping_file) 