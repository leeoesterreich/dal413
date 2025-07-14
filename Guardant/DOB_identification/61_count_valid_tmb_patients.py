import pandas as pd
import os

def count_valid_tmb_patients(file_paths):
    """
    Counts the number of unique patients with a valid (numeric) TMB score in each file.

    Args:
        file_paths (list): A list of paths to the CSV files.
    """
    print("--- Patients with Valid TMB Scores ---")
    
    for file_path in file_paths:
        count = 0
        if not os.path.exists(file_path):
            print(f"Warning: File not found at {file_path}")
        else:
            try:
                df = pd.read_csv(file_path, low_memory=False)

                if 'TMB Score' in df.columns and 'Effective Patient ID' in df.columns:
                    # Convert 'TMB Score' to numeric, turning non-numeric values into NaN
                    df['TMB Score'] = pd.to_numeric(df['TMB Score'], errors='coerce')
                    
                    # Drop rows where 'TMB Score' could not be converted (is NaN)
                    valid_tmb_df = df.dropna(subset=['TMB Score'])
                    
                    # Count the number of unique patients in the cleaned dataframe
                    count = valid_tmb_df['Effective Patient ID'].nunique()
                else:
                    print(f"Warning: Required columns not found in {file_path}")

            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")
        
        print(f"File: {os.path.basename(file_path)}, Unique patients with a valid TMB score: {count}")


if __name__ == "__main__":
    file_paths = [
        'final_genomic_cohorts/ILC_genomic_infinity.csv',
        'final_genomic_cohorts/IDC_genomic_infinity.csv'
    ]
    count_valid_tmb_patients(file_paths) 