import pandas as pd
import sys

def extract_first_true_tests(input_file, output_file):
    """
    Extracts the first test with a 'True' alteration for each patient using
    the definitive patient IDs.

    Args:
        input_file (str): Path to the updated serial genomic data.
        output_file (str): Path to save the new first-true-test data.
    """
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found.", file=sys.stderr)
        return

    # Filter for rows where 'Aleration Detected?' is True
    df_true = df[df['Aleration Detected?'] == True].copy()

    if df_true.empty:
        print("No records with 'Aleration Detected? == True' found. Output file will be empty.")
        pd.DataFrame(columns=df.columns).to_csv(output_file, index=False)
        return

    # Find the index of the first 'True' test for each patient
    first_true_indices = df_true.loc[df_true.groupby('Effective Patient ID')['Serial Test'].idxmin()].index
    
    first_true_tests_df = df.loc[first_true_indices]

    first_true_tests_df.to_csv(output_file, index=False)
    
    print(f"Successfully extracted {len(first_true_tests_df)} records of the first 'True' test for each patient.")
    print(f"Data saved to {output_file}")

if __name__ == "__main__":
    INPUT_CSV = 'serial_genomic_with_patient_id.csv'
    OUTPUT_CSV = 'first_true_test_per_patient.csv'
    extract_first_true_tests(INPUT_CSV, OUTPUT_CSV) 