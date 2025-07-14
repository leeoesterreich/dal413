import pandas as pd
import sys

def attach_patient_ids(ref_file, cohort_files, output_prefixes):
    """
    Attaches 'Effective Patient ID' to cohort files based on date of birth.

    Args:
        ref_file (str): Path to the reference CSV with patient IDs and DOBs.
        cohort_files (list): A list of paths to the cohort CSV files.
        output_prefixes (list): A list of prefixes for the output files.
    """
    try:
        df_ref = pd.read_csv(ref_file)
    except FileNotFoundError:
        print(f"Error: Reference file not found at {ref_file}")
        sys.exit(1)

    # Prepare the reference dataframe
    # Ensure 'Patient DOB' is string type for consistent matching
    df_ref['Patient DOB'] = df_ref['Patient DOB'].astype(str)
    # Create a mapping from Patient DOB to the first encountered Effective Patient ID
    dob_to_id_map = df_ref.drop_duplicates(subset=['Patient DOB']).set_index('Patient DOB')['Effective Patient ID'].to_dict()
    
    for i, cohort_file in enumerate(cohort_files):
        try:
            df_cohort = pd.read_csv(cohort_file, low_memory=False)
            print(f"Processing {cohort_file}...")

            # Clean up column names by stripping leading/trailing whitespace
            df_cohort.columns = df_cohort.columns.str.strip()

            if 'Date of Birth' not in df_cohort.columns:
                print(f"Warning: 'Date of Birth' column not found in {cohort_file}. Skipping.")
                continue

            # Convert cohort DOB to the same format as the reference DOB (M/D/YYYY)
            # The format %-m/%-d/%Y handles dates without leading zeros (e.g., 7/1/2023)
            df_cohort['Formatted DOB'] = pd.to_datetime(df_cohort['Date of Birth'], format='%Y%m%d', errors='coerce').dt.strftime('%-m/%-d/%Y')
            
            # Map the Effective Patient ID
            df_cohort['Effective Patient ID'] = df_cohort['Formatted DOB'].map(dob_to_id_map)

            # Reorder columns to make 'Effective Patient ID' the first column
            cols = ['Effective Patient ID'] + [col for col in df_cohort if col != 'Effective Patient ID']
            df_cohort = df_cohort[cols].drop(columns=['Formatted DOB']) # Drop the temporary formatted DOB column

            # Save the updated dataframe
            output_file = f"{output_prefixes[i]}_with_patient_id.csv"
            df_cohort.to_csv(output_file, index=False)
            
            # Report on matching
            matched_count = df_cohort['Effective Patient ID'].notna().sum()
            total_count = len(df_cohort)
            print(f"Finished processing {cohort_file}:")
            print(f"  - Matched {matched_count} out of {total_count} records.")
            print(f"  - Output saved to {output_file}\n")

        except FileNotFoundError:
            print(f"Warning: Cohort file not found at {cohort_file}. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing {cohort_file}: {e}")

if __name__ == "__main__":
    REFERENCE_CSV = 'serial_genomic_with_patient_id.csv'
    COHORT_CSVS = [
        'ilc_cohort.csv',
        'idc_cohort.csv'
    ]
    OUTPUT_PREFIXES = ['ilc_cohort', 'idc_cohort']
    
    attach_patient_ids(REFERENCE_CSV, COHORT_CSVS, OUTPUT_PREFIXES) 