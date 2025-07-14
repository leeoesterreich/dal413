import pandas as pd
import os

def process_final_cohorts(serial_data_path, complete_mapping_path, output_dir):
    """
    Filters serial genomic data for the first positive test by date for each patient,
    using the complete patient mapping file for patient identification and histology information.

    Args:
        serial_data_path (str): Path to the serial genomic data CSV.
        complete_mapping_path (str): Path to the complete patient mapping CSV.
        output_dir (str): Directory to save the final cohort files.
    """
    try:
        print("Step 1: Loading data files...")
        # Load data with specific encoding for the main file
        df_serial = pd.read_csv(serial_data_path, encoding='utf-8-sig')
        df_mapping = pd.read_csv(complete_mapping_path)
        print("...data loading complete.")
    except FileNotFoundError as e:
        print(f"Error: A required file was not found: {e.filename}")
        return

    # --- Data Preparation ---
    print("Step 2: Preparing and merging initial data...")
    # Merge serial data with the complete mapping to get Effective Patient ID and Histology
    df = pd.merge(df_serial, df_mapping[['GH_ID', 'Effective Patient ID', 'Histology']], 
                 on='GH_ID', how='inner')

    # Correct column typo and format date
    if 'Aleration Detected?' in df.columns:
        df.rename(columns={'Aleration Detected?': 'Alteration Detected?'}, inplace=True)
    
    # Convert dates to datetime, keeping the original format
    df['Sample Received Date'] = pd.to_datetime(df['Sample Received Date'], format='%m/%d/%Y')
    df['Sample Received Date Original'] = df['Sample Received Date'].dt.strftime('%m/%d/%Y')
    print("...data preparation complete.")

    # --- Calculate Negative Tests ---
    print("Step 3: Calculating count of tests with only 'False' results...")
    negative_tests_df = df[df['Alteration Detected?'] == False]
    # A unique test is a combination of patient and date
    num_negative_tests = negative_tests_df.drop_duplicates(subset=['Effective Patient ID', 'Sample Received Date']).shape[0]
    print(f"-> Found {num_negative_tests} unique tests where 'Alteration Detected' was False.")

    # --- Find First Positive Test by Date ---
    print("Step 4: Identifying first positive test date for each patient...")
    # Filter for positive tests
    positive_tests = df[df['Alteration Detected?'] == True].copy()
    if positive_tests.empty:
        print("No positive tests found. Cannot proceed.")
        return
        
    # Find the first positive test by date for each patient
    first_positive_tests = positive_tests.loc[positive_tests.groupby('Effective Patient ID')['Sample Received Date'].idxmin()]
    first_positive_tests = first_positive_tests[['Effective Patient ID', 'Sample Received Date']].rename(columns={'Sample Received Date': 'First Positive Date'})
    print("...first positive tests identified.")
    
    # --- Filter and Split ---
    print("Step 5: Filtering for first positive results...")
    # Merge to get the first positive date for each patient
    merged_df = pd.merge(df, first_positive_tests, on='Effective Patient ID')
    
    # Filter for rows where the test date matches the first positive date
    result_df = merged_df[merged_df['Sample Received Date'] == merged_df['First Positive Date']].copy()
    
    # Restore the original date format
    result_df['Sample Received Date'] = result_df['Sample Received Date Original']
    result_df.drop(columns=['First Positive Date', 'Sample Received Date Original'], inplace=True, errors='ignore')
    
    # Split by histology
    df_ilc = result_df[result_df['Histology'] == 'ILC']
    df_idc = result_df[result_df['Histology'] == 'IDC']
    print("...filtering and splitting complete.")

    # --- Save Final Files ---
    print("Step 6: Saving final cohort files...")
    ilc_output_path = os.path.join(output_dir, 'ILC_genomic.csv')
    idc_output_path = os.path.join(output_dir, 'IDC_genomic.csv')
    
    df_ilc.to_csv(ilc_output_path, index=False)
    df_idc.to_csv(idc_output_path, index=False)
    print(f"-> ILC genomic data saved to '{ilc_output_path}' ({len(df_ilc)} rows)")
    print(f"-> IDC genomic data saved to '{idc_output_path}' ({len(df_idc)} rows)")
    
    print("\nScript finished successfully!")


if __name__ == "__main__":
    # Define file paths
    SERIAL_DATA = "Guardant Project - Serial Genomic Data.csv"
    COMPLETE_MAPPING = "complete_patient_mapping.csv"
    OUTPUT_DIR = "final_genomic_cohorts"

    process_final_cohorts(SERIAL_DATA, COMPLETE_MAPPING, OUTPUT_DIR) 