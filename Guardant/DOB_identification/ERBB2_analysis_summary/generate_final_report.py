import pandas as pd
import sys

def generate_final_report(serial_genomic_file, mapping_file, output_csv_file):
    """
    Consolidates all previous steps to generate a final, complete report
    for a specific list of 7 patients.

    Args:
        serial_genomic_file (str): Absolute path to the serial genomic data.
        mapping_file (str): Absolute path to the patient mapping data.
        output_csv_file (str): Absolute path for the final output CSV.
    """
    
    # The final list of 7 patients we are interested in.
    target_patient_ids = ['P125', 'P2', 'P47', 'P806', 'P811', 'P95', 'P869']
    print("--- Generating Final Consolidated Report ---")
    print("Target Patients:", target_patient_ids)

    try:
        serial_df = pd.read_csv(serial_genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except Exception as e:
        print(f"Error reading source files: {e}", file=sys.stderr)
        return

    # Merge dataframes to get the Effective Patient ID
    full_df = pd.merge(serial_df, mapping_df, on='GH_ID', how='left')
    
    # Filter for the 7 target patients
    final_df = full_df[full_df['Effective Patient ID'].isin(target_patient_ids)].copy()

    if final_df.empty:
        print("Could not find any records for the target patients.")
        return

    # Sort the data for clarity
    final_df.sort_values(by=['Effective Patient ID', 'Serial Test', 'Final Report Date'], inplace=True)

    # Reorder columns to bring important identifiers to the front
    all_cols = final_df.columns.tolist()
    # Define the desired order, the rest will follow
    preferred_cols = ['Histology', 'Effective Patient ID', 'GH_ID', 'Serial Test', 'Total Test', 'Final Report Date', 'Gene', 'Alteration', 'Vus']
    # Create the final column order
    other_cols = [col for col in all_cols if col not in preferred_cols]
    final_col_order = preferred_cols + other_cols
    final_df = final_df[final_col_order]

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
    output_file = f'{base}/ERBB2_analysis_summary/FINAL_COMPLETE_PATIENT_HISTORY.csv'

    generate_final_report(serial_file, mapping_file, output_file) 