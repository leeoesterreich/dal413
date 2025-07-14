import pandas as pd
import sys

def inspect_column(filename, column_name):
    """Prints unique values and their counts for a specific column in a CSV."""
    try:
        # Using low_memory=False to prevent dtype warnings with mixed types
        df = pd.read_csv(filename, low_memory=False)
        if column_name in df.columns:
            print(f"--- Analysis for {filename} ---")
            print(f"Unique values in '{column_name}' column:")
            print(df[column_name].unique())
            print(f"\nValue counts in '{column_name}' column (including empty values):")
            print(df[column_name].value_counts(dropna=False))
            print("-" * 40 + "\n")
        else:
            print(f"Column '{column_name}' not found in '{filename}'.")
    except Exception as e:
        print(f"Error processing {filename}: {e}", file=sys.stderr)

if __name__ == "__main__":
    print("Inspecting 'Vus' column in original data files...")
    inspect_column('IDC_genomic.csv', 'Vus')
    inspect_column('ILC_genomic.csv', 'Vus') 