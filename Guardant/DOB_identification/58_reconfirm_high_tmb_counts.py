import pandas as pd
import os

def reconfirm_high_tmb_counts(file_paths):
    """
    Counts the number of unique patients with a TMB score > 10, separately for each file and combined.

    Args:
        file_paths (list): A list of paths to the CSV files.
    """
    total_unique_patients = set()
    
    print("--- Patients with TMB Score > 10 ---")
    
    for file_path in file_paths:
        count = 0
        if not os.path.exists(file_path):
            print(f"Warning: File not found at {file_path}")
        else:
            try:
                df = pd.read_csv(file_path, low_memory=False)

                if 'TMB Score' in df.columns and 'Effective Patient ID' in df.columns:
                    df['TMB Score'] = pd.to_numeric(df['TMB Score'], errors='coerce')
                    high_tmb_df = df.dropna(subset=['TMB Score'])
                    high_tmb_df = high_tmb_df[high_tmb_df['TMB Score'] > 10]
                    
                    count = high_tmb_df['Effective Patient ID'].nunique()
                    total_unique_patients.update(high_tmb_df['Effective Patient ID'].unique())
                else:
                    print(f"Warning: Required columns not found in {file_path}")

            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")
        
        print(f"File: {os.path.basename(file_path)}, Unique patients with TMB > 10: {count}")

    print(f"\nTotal unique patients across all files with TMB > 10: {len(total_unique_patients)}")


if __name__ == "__main__":
    file_paths = [
        'final_genomic_cohorts/ILC_genomic_infinity.csv',
        'final_genomic_cohorts/IDC_genomic_infinity.csv'
    ]
    reconfirm_high_tmb_counts(file_paths) 