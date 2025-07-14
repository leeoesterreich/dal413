import pandas as pd
import os
import sys

def get_first_positive_tests(serial_data_file, mapping_file, output_dir):
    """
    Filters genomic data to find the complete results for each patient's first 
    test that contains at least one positive alteration.

    Args:
        serial_data_file (str): Path to the serial genomic data.
        mapping_file (str): Path to the complete patient mapping file.
        output_dir (str): Directory to save the final cohort files.
    """
    try:
        # --- 1. Load Data ---
        print("Step 1: Loading data files...")
        patient_map_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        serial_df = pd.read_csv(serial_data_file)
        print("...data loading complete.")

        # --- 2. Data Preparation and Merging ---
        print("Step 2: Preparing and merging data...")
        serial_df.columns = serial_df.columns.str.strip()
        merged_df = pd.merge(serial_df, patient_map_df, on='GH_ID', how='left')
        merged_df.dropna(subset=['Effective Patient ID'], inplace=True)
        
        # Ensure correct data types for key columns
        merged_df['Serial Test'] = pd.to_numeric(merged_df['Serial Test'], errors='coerce')
        # The 'Aleration Detected?' column is already boolean, so we can use it directly.
        # We will rename it for clarity and to handle the '?' character.
        merged_df.rename(columns={'Aleration Detected?': 'Alteration_Detected'}, inplace=True)

        merged_df.dropna(subset=['Serial Test', 'Alteration_Detected'], inplace=True)
        merged_df['Serial Test'] = merged_df['Serial Test'].astype(int)
        print("...data merged and prepared.")

        # --- 3. Identify First Positive Test for Each Patient ---
        print("Step 3: Identifying the first positive test for each patient...")
        
        # Identify which tests had any positive result
        # A 'test' is a unique combo of patient and serial number
        test_positivity = merged_df.groupby(['Effective Patient ID', 'Serial Test'])['Alteration_Detected'].any().reset_index()
        positive_tests = test_positivity[test_positivity['Alteration_Detected'] == True]

        if positive_tests.empty:
            print("No positive tests found in the entire dataset.", file=sys.stderr)
            sys.exit(0)

        # Find the minimum (first) positive serial test number for each patient
        first_positive_test_map = positive_tests.loc[positive_tests.groupby('Effective Patient ID')['Serial Test'].idxmin()]
        first_positive_test_map = first_positive_test_map.rename(columns={'Serial Test': 'First_Positive_Serial_Test'})
        print("...first positive tests identified.")

        # --- 4. Filter Main DataFrame for Complete First Positive Tests ---
        print("Step 4: Filtering for all rows from the first positive tests...")
        
        # Merge this map back to the main dataframe
        final_df = pd.merge(
            merged_df, 
            first_positive_test_map[['Effective Patient ID', 'First_Positive_Serial_Test']],
            on='Effective Patient ID',
            how='inner' # Keep only patients who had a positive test
        )

        # Keep only the rows that match the first positive test number
        final_df = final_df[final_df['Serial Test'] == final_df['First_Positive_Serial_Test']].copy()
        
        # Clean up helper columns
        final_df.drop(columns=['First_Positive_Serial_Test'], inplace=True)
        print("...filtering complete.")

        # --- 5. Split by Histology and Save ---
        print("Step 5: Splitting by histology and saving final files...")
        os.makedirs(output_dir, exist_ok=True)

        ilc_cohort_df = final_df[final_df['Histology'] == 'ILC']
        idc_cohort_df = final_df[final_df['Histology'] == 'IDC']

        ilc_output_path = os.path.join(output_dir, 'ILC_first_positive.csv')
        idc_output_path = os.path.join(output_dir, 'IDC_first_positive.csv')

        ilc_cohort_df.to_csv(ilc_output_path, index=False)
        print(f"-> ILC cohort saved to {ilc_output_path} ({len(ilc_cohort_df)} rows)")
        
        idc_cohort_df.to_csv(idc_output_path, index=False)
        print(f"-> IDC cohort saved to {idc_output_path} ({len(idc_cohort_df)} rows)")

        print("\nScript finished successfully!")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    SERIAL_DATA = 'Guardant Project - Serial Genomic Data.csv'
    MAPPING_FILE = 'final_genomic_cohorts/complete_patient_mapping.csv'
    OUTPUT_DIR = 'serial_filtering'
    
    get_first_positive_tests(SERIAL_DATA, MAPPING_FILE, OUTPUT_DIR) 