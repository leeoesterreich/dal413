import pandas as pd
import sys

def get_complete_history_by_gh_id(serial_genomic_file, mapping_file, output_csv_file):
    """
    Finds the complete testing history for patients by first identifying all their
    GH_IDs and then filtering the main data file. This is a more robust method
    that avoids potentially problematic merges.

    Args:
        serial_genomic_file (str): Absolute path to the serial genomic data.
        mapping_file (str): Absolute path to the patient mapping data.
        output_csv_file (str): Absolute path for the final output CSV.
    """
    
    target_patient_ids = ['P125', 'P2', 'P47', 'P806', 'P811', 'P95', 'P869']
    print("--- Generating Final Report (Robust Method) ---")
    print("Target Patients:", target_patient_ids)

    try:
        serial_df = pd.read_csv(serial_genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except Exception as e:
        print(f"Error reading source files: {e}", file=sys.stderr)
        return

    # 1. Get all GH_IDs for our target patients from the mapping file
    target_gh_ids = mapping_df[mapping_df['Effective Patient ID'].isin(target_patient_ids)]['GH_ID'].unique()
    
    if len(target_gh_ids) == 0:
        print("Could not find any GH_IDs for the target patients.")
        return
    
    print(f"Found {len(target_gh_ids)} unique GH_IDs corresponding to the {len(target_patient_ids)} target patients.")

    # 2. Filter the main serial data file using this list of GH_IDs
    final_df = serial_df[serial_df['GH_ID'].isin(target_gh_ids)].copy()

    # 3. Add the 'Effective Patient ID' and 'Histology' back for context
    # We use a map for this to avoid a full merge.
    patient_id_map = mapping_df.set_index('GH_ID')['Effective Patient ID']
    histology_map = mapping_df.set_index('GH_ID')['Histology']
    final_df['Effective Patient ID'] = final_df['GH_ID'].map(patient_id_map)
    final_df['Histology'] = final_df['GH_ID'].map(histology_map)

    # Sort the data for clarity
    final_df.sort_values(by=['Effective Patient ID', 'Serial Test', 'Final Report Date'], inplace=True)
    
    # Reorder columns
    preferred_cols = ['Histology', 'Effective Patient ID', 'GH_ID', 'Serial Test', 'Total Test', 'Final Report Date', 'Gene', 'Alteration', 'Vus']
    other_cols = [col for col in final_df.columns if col not in preferred_cols]
    final_df = final_df[preferred_cols + other_cols]

    try:
        final_df.to_csv(output_csv_file, index=False)
        print(f"\nSuccessfully generated final report with {len(final_df)} records.")
        print(f"File saved to: {output_csv_file}")
    except Exception as e:
        print(f"Error saving final report: {e}", file=sys.stderr)


if __name__ == '__main__':
    base = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification'
    serial_file = f'{base}/Guardant Project - Serial Genomic Data.csv'
    mapping_file = f'{base}/final_genomic_cohorts/complete_patient_mapping.csv'
    output_file = f'{base}/ERBB2_analysis_summary/DEFINITIVELY_FINAL_PATIENT_HISTORY.csv'

    get_complete_history_by_gh_id(serial_file, mapping_file, output_file) 