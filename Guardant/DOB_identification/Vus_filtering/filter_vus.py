import pandas as pd
import numpy as np

def filter_vus(input_file, output_file):
    """
    Filters a CSV file to include rows where 'Vus' is FALSE or the value is missing (NaN).
    This is more robust to data entry inconsistencies.

    Args:
        input_file (str): The path to the input CSV file.
        output_file (str): The path to the output CSV file.
    """
    try:
        # Read 'Vus' column as string to avoid automatic type conversion
        df = pd.read_csv(input_file, dtype={'Vus': str}, low_memory=False)

        if 'Vus' in df.columns:
            # Normalize the 'Vus' column: fill NaNs with empty string, strip whitespace, convert to upper case
            vus_normalized = df['Vus'].fillna('').str.strip().str.upper()

            # Keep rows where 'Vus' is 'FALSE' or where it was originally empty/NaN.
            # Clinically relevant mutations often have a blank Vus field.
            filtered_df = df[ (vus_normalized == 'FALSE') | (vus_normalized == '') ]

            filtered_df.to_csv(output_file, index=False)
            print(f"Filtered data from '{input_file}' and saved to '{output_file}'. Kept {len(filtered_df)} rows.")
        else:
            print(f"Error: 'Vus' column not found in '{input_file}'")

    except FileNotFoundError:
        print(f"Error: File not found at '{input_file}'")
    except Exception as e:
        print(f"An error occurred while processing '{input_file}': {e}")

if __name__ == "__main__":
    files_to_process = {
        "IDC_genomic.csv": "IDC_genomic_filtered.csv",
        "ILC_genomic.csv": "ILC_genomic_filtered.csv"
    }

    for input_csv, output_csv in files_to_process.items():
        filter_vus(input_csv, output_csv) 