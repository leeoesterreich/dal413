import pandas as pd
import os

def calculate_median_tmb_per_file(file_paths):
    """
    Calculates the median TMB score for all unique patients in each specified CSV file.

    Args:
        file_paths (list): A list of paths to the CSV files.
    """
    
    print("--- Median TMB Score per File ---")
    for file_path in file_paths:
        if not os.path.exists(file_path):
            print(f"Warning: File not found at {file_path}")
            continue
            
        try:
            df = pd.read_csv(file_path, low_memory=False)

            if 'TMB Score' not in df.columns or 'Effective Patient ID' not in df.columns:
                print(f"Warning: 'TMB Score' or 'Effective Patient ID' not found in {file_path}")
                continue

            # Convert 'TMB Score' to numeric, coercing errors will turn non-numeric values into NaN
            df['TMB Score'] = pd.to_numeric(df['TMB Score'], errors='coerce')
            
            # Drop rows where 'TMB Score' is NaN to not affect median calculation
            df.dropna(subset=['TMB Score'], inplace=True)

            # Get the TMB score for each unique patient. If a patient has multiple entries, this takes the first one.
            unique_patient_tmb = df.drop_duplicates(subset=['Effective Patient ID'])
            
            # Calculate the median TMB score for the unique patients in the current file
            median_tmb = unique_patient_tmb['TMB Score'].median()
            
            print(f"File: {file_path}, Median TMB Score for all unique patients: {median_tmb:.2f}")

        except Exception as e:
            print(f"An error occurred while processing {file_path}: {e}")

if __name__ == "__main__":
    file_paths = [
        'final_genomic_cohorts/ILC_genomic_infinity.csv',
        'final_genomic_cohorts/IDC_genomic_infinity.csv'
    ]
    calculate_median_tmb_per_file(file_paths) 