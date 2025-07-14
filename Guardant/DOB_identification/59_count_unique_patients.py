import pandas as pd
import os

def count_total_unique_patients(file_paths):
    """
    Counts the total number of unique patients for each file and across all files.

    Args:
        file_paths (list): A list of paths to the CSV files.
    """
    total_unique_patients_set = set()
    
    print("--- Unique Patient Counts ---")
    
    for file_path in file_paths:
        count = 0
        if not os.path.exists(file_path):
            print(f"Warning: File not found at {file_path}")
        else:
            try:
                df = pd.read_csv(file_path, low_memory=False)

                if 'Effective Patient ID' in df.columns:
                    # Get unique patient IDs from the current file
                    unique_ids_in_file = df['Effective Patient ID'].unique()
                    count = len(unique_ids_in_file)
                    
                    # Add these unique IDs to the overall set
                    total_unique_patients_set.update(unique_ids_in_file)
                else:
                    print(f"Warning: 'Effective Patient ID' column not found in {file_path}")

            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")
        
        print(f"File: {os.path.basename(file_path)}, Total unique patients: {count}")

    print(f"\nTotal unique patients across all files: {len(total_unique_patients_set)}")


if __name__ == "__main__":
    file_paths = [
        'Vus_filtering/ILC_count_infinity.csv',
        'Vus_filtering/IDC_count_infinity.csv'
    ]
    count_total_unique_patients(file_paths) 