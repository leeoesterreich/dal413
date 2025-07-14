import pandas as pd
import sys

def analyze_test_patterns(serial_data_file, mapping_file):
    """
    Analyzes the serial testing data to understand the relationship between
    serial test number, total tests, and when positive alterations are detected.

    Args:
        serial_data_file (str): Path to the serial genomic data.
        mapping_file (str): Path to the complete patient mapping file.
    """
    try:
        # --- 1. Load Data ---
        print("Step 1: Loading data files...")
        patient_map_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        serial_df = pd.read_csv(serial_data_file)
        print("...data loading complete.")

        # --- 2. Data Preparation ---
        print("Step 2: Preparing and merging data...")
        # Clean column names
        serial_df.columns = serial_df.columns.str.strip()
        # Ensure correct data types
        serial_df['Serial Test'] = pd.to_numeric(serial_df['Serial Test'], errors='coerce')
        serial_df['Total Test'] = pd.to_numeric(serial_df['Total Test'], errors='coerce')
        # Explicitly handle boolean conversion for the detection column
        serial_df['Alteration_Detected'] = serial_df['Aleration Detected?'].apply(lambda x: True if x == 'TRUE' else (False if x == 'FALSE' else None))

        # Merge with mapping to get Effective Patient ID
        merged_df = pd.merge(serial_df, patient_map_df[['GH_ID', 'Effective Patient ID']], on='GH_ID', how='left')
        merged_df.dropna(subset=['Effective Patient ID', 'Serial Test', 'Total Test'], inplace=True)
        print("...data prepared.")

        # --- Debug: Print data for a specific patient ---
        target_gh_id = 'A0192917' # A patient we know exists
        patient_specific_data = merged_df[merged_df['GH_ID'] == target_gh_id]
        print(f"\n--- Data for GH_ID: {target_gh_id} ---")
        print(patient_specific_data[['GH_ID', 'Effective Patient ID', 'Serial Test', 'Total Test', 'Aleration Detected?']].to_string())
        print("-------------------------------------\n")

        # --- 3. Analyze Patterns ---
        print("\nStep 3: Analyzing test patterns for patients with multiple tests...")
        
        # Get a list of patients with more than one test
        patient_test_counts = merged_df.groupby('Effective Patient ID')['Serial Test'].nunique()
        print("\n--- Debug: Patient Test Counts ---")
        print(patient_test_counts.head(10))
        print(f"Max test count found for any patient: {patient_test_counts.max()}")
        print("----------------------------------\n")
        multi_test_patients = patient_test_counts[patient_test_counts > 1].index

        if len(multi_test_patients) == 0:
            print("No patients with multiple tests were found.")
            sys.exit(0)

        print(f"Found {len(multi_test_patients)} patients with more than one test. Analyzing a sample...")

        # Filter for these patients
        multi_test_df = merged_df[merged_df['Effective Patient ID'].isin(multi_test_patients)]

        # For each patient, find out which of their tests had a positive result
        # A test is 'positive' if any of its gene rows has Alteration_Detected == True
        test_positivity = multi_test_df.groupby(['Effective Patient ID', 'Serial Test'])['Alteration_Detected'].any().reset_index()
        test_positivity.rename(columns={'Alteration_Detected': 'Has_Positive_Alteration'}, inplace=True)

        # Merge this back with the total test info
        # We only need one row per test, so let's get the Total Test value
        test_summary = multi_test_df[['Effective Patient ID', 'Serial Test', 'Total Test']].drop_duplicates()
        final_summary = pd.merge(test_positivity, test_summary, on=['Effective Patient ID', 'Serial Test'])

        # --- 4. Report Findings ---
        print("\n--- Analysis Summary ---")
        print("This table shows, for each test taken by a multi-test patient, whether it contained a positive result.")
        
        # Displaying a sample of 15 patients (or all if fewer than 15)
        sample_patients = pd.Series(final_summary['Effective Patient ID'].unique()).sample(min(15, len(multi_test_patients)))
        
        print(final_summary[final_summary['Effective Patient ID'].isin(sample_patients)].to_string())

        print("\n--- Overall Statistics ---")
        # Check for your hypothesis: Do positive alterations only appear on the last test?
        # A 'last test' is where Serial Test == Total Test
        last_tests = final_summary[final_summary['Serial Test'] == final_summary['Total Test']]
        not_last_tests = final_summary[final_summary['Serial Test'] != final_summary['Total Test']]

        positive_on_last_only = last_tests[last_tests['Has_Positive_Alteration'] == True]
        positive_on_not_last = not_last_tests[not_last_tests['Has_Positive_Alteration'] == True]

        print(f"\nNumber of 'last tests' that were positive: {len(positive_on_last_only)}")
        print(f"Number of 'non-last tests' that were positive: {len(positive_on_not_last)}")
        
        if len(positive_on_not_last) > 0:
            print("\nConclusion: The hypothesis that alterations are ONLY found on the last test appears to be FALSE.")
            print("There are instances of positive results on tests that are not the patient's final test.")
        else:
            print("\nConclusion: The hypothesis that alterations are ONLY found on the last test appears to be TRUE for this dataset.")
        
        print("\nScript finished.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    SERIAL_DATA = 'Guardant Project - Serial Genomic Data.csv'
    MAPPING_FILE = 'final_genomic_cohorts/complete_patient_mapping.csv'
    analyze_test_patterns(SERIAL_DATA, MAPPING_FILE) 