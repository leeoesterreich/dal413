import pandas as pd

def analyze_combined_data(files):
    """
    Combines multiple CSV files and analyzes the combined data.

    Args:
        files (list): A list of paths to the CSV files.
    """
    try:
        dataframes = [pd.read_csv(f) for f in files]
        combined_df = pd.concat(dataframes, ignore_index=True)

        if 'Gene' in combined_df.columns and 'Alteration' in combined_df.columns:
            unique_genes = combined_df['Gene'].nunique()
            unique_combinations = combined_df[['Gene', 'Alteration']].drop_duplicates().shape[0]

            print("Analysis for combined IDC and ILC data:")
            print(f"  - Number of unique items in 'Gene' column: {unique_genes}")
            print(f"  - Number of unique combinations of 'Gene' and 'Alteration': {unique_combinations}")
        else:
            print("Error: Required columns ('Gene', 'Alteration') not found in the combined data.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    files_to_combine = [
        "IDC_genomic_filtered.csv",
        "ILC_genomic_filtered.csv"
    ]
    analyze_combined_data(files_to_combine) 