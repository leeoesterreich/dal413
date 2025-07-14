import pandas as pd

def analyze_dataframe(df, filename):
    """
    Analyzes a DataFrame to find unique genes and unique Gene/Alteration combinations.

    Args:
        df (pd.DataFrame): The DataFrame to analyze.
        filename (str): The name of the file being analyzed, for printing.
    """
    if 'Gene' in df.columns and 'Alteration' in df.columns:
        unique_genes = df['Gene'].nunique()
        unique_combinations = df[['Gene', 'Alteration']].drop_duplicates().shape[0]

        print(f"Analysis for '{filename}':")
        print(f"  - Number of unique items in 'Gene' column: {unique_genes}")
        print(f"  - Number of unique combinations of 'Gene' and 'Alteration': {unique_combinations}")
        print("-" * 30)
    else:
        print(f"Error: Required columns ('Gene', 'Alteration') not found in '{filename}'")

if __name__ == "__main__":
    filtered_files = [
        "IDC_genomic_filtered.csv",
        "ILC_genomic_filtered.csv"
    ]

    for file in filtered_files:
        try:
            df = pd.read_csv(file)
            analyze_dataframe(df, file)
        except FileNotFoundError:
            print(f"Error: File not found at '{file}'")
        except Exception as e:
            print(f"An error occurred while processing '{file}': {e}") 