import pandas as pd
import os
import sys

def filter_by_serial_test(serial_data_file, mapping_file, output_dir):
    """
    Filters genomic data to find the first positive test for each patient based on
    the lowest 'Serial test #' instead of the collection date.

    Args:
        serial_data_file (str): Path to the serial genomic data.
        mapping_file (str): Path to the complete patient mapping file.
        output_dir (str): Directory to save the final cohort files.
    """
    try:
        # --- 1. Load Data ---
        print("Step 1: Loading data files...")
        try:
            # Handle potential BOM in the mapping file
            patient_map_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        except FileNotFoundError:
            print(f"Error: Mapping file not found at {mapping_file}", file=sys.stderr)
            sys.exit(1)

        try:
            serial_df = pd.read_csv(serial_data_file)
        except FileNotFoundError:
            print(f"Error: Serial data file not found at {serial_data_file}", file=sys.stderr)
            sys.exit(1)
        print("...data loading complete.")

        # --- 2. Data Preparation and Merging ---
        print("Step 2: Preparing and merging data...")
        # Clean column names for easier access
        serial_df.columns = serial_df.columns.str.strip().str.replace(' ', '_').str.replace('#', '_num')
        
        # Merge to get Effective Patient ID and Histology for each test
        merged_df = pd.merge(serial_df, patient_map_df, on='GH_ID', how='left')

        # Drop rows where there's no patient mapping
        merged_df.dropna(subset=['Effective Patient ID'], inplace=True)
        print("...data merged.")
        print("Merged DF columns:", merged_df.columns) # DEBUG

        # --- 3. Filter for First Positive Test ---
        print("Step 3: Filtering for first positive test based on serial number...")
        # Convert test number columns to numeric, coercing errors
        merged_df['Serial_Test'] = pd.to_numeric(merged_df['Serial_Test'], errors='coerce')
        merged_df.dropna(subset=['Serial_Test'], inplace=True) # Drop if conversion failed
        merged_df['Serial_Test'] = merged_df['Serial_Test'].astype(int)

        # Filter for positive results
        positive_tests_df = merged_df[merged_df['Aleration_Detected?'] != 'Not Detected'].copy()

        if positive_tests_df.empty:
            print("Warning: No positive tests found after filtering. Output files will be empty.", file=sys.stderr)
        else:
            # Find the minimum serial test number for each patient
            positive_tests_df['min_serial_test'] = positive_tests_df.groupby('Effective Patient ID')['Serial_Test'].transform('min')
            
            # Filter the dataframe to keep only rows that match the minimum serial number
            first_positive_df = positive_tests_df[positive_tests_df['Serial_Test'] == positive_tests_df['min_serial_test']]
            print("...filtering complete.")

        # --- 4. Split by Histology and Save ---
        print("Step 4: Splitting by histology and saving final files...")
        os.makedirs(output_dir, exist_ok=True)

        # Handle case where no positive tests were found
        final_df = first_positive_df if not positive_tests_df.empty else pd.DataFrame(columns=merged_df.columns)

        # Split into ILC and IDC cohorts
        ilc_cohort_df = final_df[final_df['Histology'] == 'ILC']
        idc_cohort_df = final_df[final_df['Histology'] == 'IDC']

        # Define output paths
        ilc_output_path = os.path.join(output_dir, 'ILC_genomic_filtered.csv')
        idc_output_path = os.path.join(output_dir, 'IDC_genomic_filtered.csv')

        # Save the final files
        ilc_cohort_df.to_csv(ilc_output_path, index=False)
        print(f"-> ILC cohort saved to {ilc_output_path} ({len(ilc_cohort_df)} rows)")
        
        idc_cohort_df.to_csv(idc_output_path, index=False)
        print(f"-> IDC cohort saved to {idc_output_path} ({len(idc_cohort_df)} rows)")

        print("\nScript finished successfully!")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Define file paths
    # Using relative paths for portability
    SERIAL_DATA = 'Guardant Project - Serial Genomic Data.csv'
    # The mapping file is now expected to be in a specific subfolder
    MAPPING_FILE = 'final_genomic_cohorts/complete_patient_mapping.csv'
    OUTPUT_DIR = 'serial_filtering'

    # Run the main function
    filter_by_serial_test(SERIAL_DATA, MAPPING_FILE, OUTPUT_DIR) 