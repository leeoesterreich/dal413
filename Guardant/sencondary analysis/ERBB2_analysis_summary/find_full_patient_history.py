import pandas as pd
import sys

def find_full_patient_history(serial_genomic_file, mapping_file, output_csv_file):
    """
    Finds the complete testing history for a predefined list of patients, adds
    histology information, and saves it to a new CSV file.

    Args:
        serial_genomic_file (str): Path to the complete serial genomic data CSV.
        mapping_file (str): Path to the patient mapping CSV file.
        output_csv_file (str): Path for the output CSV file.
    """
    
    # Adding P869 back into the list to see their full history.
    target_patient_ids = ['P125', 'P2', 'P47', 'P806', 'P811', 'P95', 'P869']

    print("Searching for the complete testing history for the following patients:")
    print(target_patient_ids)
    
    try:
        serial_genomic_df = pd.read_csv(serial_genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except Exception as e:
        print(f"An error occurred while reading files: {e}", file=sys.stderr)
        return

    # Merge to link Effective Patient ID
    full_merged_df = pd.merge(serial_genomic_df, mapping_df, on='GH_ID', how='left')

    # Filter for our target patients
    patient_history_df = full_merged_df[full_merged_df['Effective Patient ID'].isin(target_patient_ids)].copy()

    if patient_history_df.empty:
        print("\nCould not find any testing history for the specified patients.")
        return

    # Sort the data
    patient_history_df.sort_values(by=['Effective Patient ID', 'Final Report Date', 'Serial Test'], inplace=True)
    
    # Define and reorder columns
    final_df = patient_history_df.copy()
    
    # Explicitly select and order the columns for the final output
    final_cols = ['Histology', 'Effective Patient ID', 'GH_ID', 'Final Report Date', 'Serial Test', 'Total Test', 'Gene', 'Alteration', 'Vus']
    # Add other columns from the original dataframe that might be useful, avoiding duplicates
    other_cols = [col for col in patient_history_df.columns if col not in final_cols]
    final_cols.extend(other_cols)
    
    # Ensure all selected columns exist
    final_df = final_df[[col for col in final_cols if col in final_df.columns]]

    try:
        final_df.to_csv(output_csv_file, index=False)
        print(f"\nSuccessfully found {len(patient_history_df)} total records for {len(target_patient_ids)} patients.")
        print(f"The complete testing history has been saved to: {output_csv_file}")
    except Exception as e:
        print(f"\nAn error occurred while saving the CSV file: {e}", file=sys.stderr)


if __name__ == '__main__':
    base_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification'
    serial_data_file = f'{base_path}/Guardant Project - Serial Genomic Data.csv'
    patient_mapping_file = f'{base_path}/final_genomic_cohorts/complete_patient_mapping.csv'
    output_file = f'{base_path}/ERBB2_analysis_summary/serial_patients_full_test_history_v3.csv'
    
    find_full_patient_history(serial_data_file, patient_mapping_file, output_file) 