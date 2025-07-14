import pandas as pd
import numpy as np
import sys

def count_total_unique_in_file(filename):
    """
    Reads a CSV file and counts the total number of unique items across all columns.
    """
    try:
        df = pd.read_csv(filename)
        print(f"--- Analysis for: {filename} ---")

        # Get all values from the dataframe, flatten them into a 1D array, and remove NaNs
        all_values = df.values.ravel()
        all_values_no_nan = all_values[~pd.isna(all_values)]

        # Find the number of unique items in the flattened array
        total_unique_count = len(pd.unique(all_values_no_nan))

        print(f"Total number of unique items across all columns: {total_unique_count}")
        print("-" * (22 + len(filename)) + "\n")

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while processing {filename}: {e}", file=sys.stderr)

if __name__ == "__main__":
    files_to_analyze = [
        "ILC_count.csv",
        "ILC_count_infinity.csv"
    ]
    for f in files_to_analyze:
        count_total_unique_in_file(f) 