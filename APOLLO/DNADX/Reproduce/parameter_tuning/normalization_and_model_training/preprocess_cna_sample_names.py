import pandas as pd
import os

def main():
    base_path = "." 
    cna_input_filename = "cna_segment_score_mean_no_norm.pkl"
    cna_output_filename = "cna_segment_score_mean_no_norm_renamed.pkl"
    
    cna_input_path = os.path.join(base_path, "training_data", cna_input_filename)
    cna_output_path = os.path.join(base_path, "training_data", cna_output_filename)

    print(f"Loading CNA data from: {cna_input_path}")
    try:
        cna_df = pd.read_pickle(cna_input_path)
    except FileNotFoundError:
        print(f"Error: Input CNA file not found at {cna_input_path}")
        return

    print(f"Original CNA data shape: {cna_df.shape}")
    print("Original first 5 CNA sample names:", cna_df.columns.tolist()[:5])

    # Rename columns: replace '.' with '-'
    renamed_columns = [col.replace('.', '-') for col in cna_df.columns]
    cna_df.columns = renamed_columns

    print("Renamed first 5 CNA sample names:", cna_df.columns.tolist()[:5])

    print(f"Saving renamed CNA data to: {cna_output_path}")
    try:
        cna_df.to_pickle(cna_output_path)
        print("Successfully saved renamed CNA data.")
    except Exception as e:
        print(f"Error saving renamed CNA data: {e}")

if __name__ == "__main__":
    main() 