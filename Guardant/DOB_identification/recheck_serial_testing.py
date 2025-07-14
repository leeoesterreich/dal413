import pandas as pd
import sys

def recheck_serial_testing(multi_hit_file, serial_genomic_file, mapping_file):
    """
    Re-analyzes patients from a given list to provide a comprehensive view of their
    serial testing history from the complete genomic data.

    Args:
        multi_hit_file (str): Path to the CSV of patients with multiple non-VUS ERBB2 hits.
        serial_genomic_file (str): Path to the complete serial genomic data CSV.
        mapping_file (str): Path to the patient mapping CSV file.
    """
    try:
        multi_hit_df = pd.read_csv(multi_hit_file)
        serial_genomic_df = pd.read_csv(serial_genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading files: {e}", file=sys.stderr)
        return

    # 1. Get the list of 12 patients to analyze
    patient_ids = multi_hit_df['Effective Patient ID'].unique()
    print(f"--- Re-analyzing Serial Testing for {len(patient_ids)} Patients ---")
    print("Patient IDs being analyzed:", patient_ids.tolist())
    print("-" * 30 + "\n")

    # 2. Prepare the full dataset
    full_merged_df = pd.merge(serial_genomic_df, mapping_df, on='GH_ID', how='left')

    # 3. Filter for only the tests belonging to our 12 patients
    patient_all_tests_df = full_merged_df[full_merged_df['Effective Patient ID'].isin(patient_ids)].copy()
    
    # 4. Apply the sanity check filter requested
    initial_rows = len(patient_all_tests_df)
    patient_all_tests_df = patient_all_tests_df[patient_all_tests_df['Serial Test'] <= patient_all_tests_df['Total Test']]
    final_rows = len(patient_all_tests_df)
    if initial_rows > final_rows:
        print(f"Data Quality Note: Removed {initial_rows - final_rows} records where 'Serial Test' > 'Total Test'.\n")

    # 5. Analyze for Serial Testing
    
    # Definition A: Patients with tests on multiple report dates
    print("--- Analysis 1: Serial Testing by Multiple Test Dates ---")
    patient_all_tests_df['Final Report Date'] = pd.to_datetime(patient_all_tests_df['Final Report Date'], errors='coerce')
    serial_by_date = patient_all_tests_df.groupby('Effective Patient ID')['Final Report Date'].nunique()
    patients_serial_by_date = serial_by_date[serial_by_date > 1]
    
    if not patients_serial_by_date.empty:
        print(f"Found {len(patients_serial_by_date)} patients with tests on multiple dates:")
        print(patients_serial_by_date)
    else:
        print("No patients found with tests on multiple dates.")
    print("-" * 30 + "\n")

    # Definition B: Patients with 'Serial Test' column value > 1
    print("--- Analysis 2: Serial Testing by 'Serial Test' Column Value ---")
    serial_by_column = patient_all_tests_df.groupby('Effective Patient ID')['Serial Test'].max()
    patients_serial_by_column = serial_by_column[serial_by_column > 1]
    
    if not patients_serial_by_column.empty:
        print(f"Found {len(patients_serial_by_column)} patients with a 'Serial Test' value > 1:")
        print(patients_serial_by_column)
    else:
        print("No patients found with a 'Serial Test' value > 1.")
    print("-" * 30 + "\n")


if __name__ == '__main__':
    multi_hit_patient_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/multi_hit_non_vus_erbb2_patients.csv'
    serial_data_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/Guardant Project - Serial Genomic Data.csv'
    patient_mapping_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/complete_patient_mapping.csv'
    
    recheck_serial_testing(multi_hit_patient_file, serial_data_file, patient_mapping_file) 