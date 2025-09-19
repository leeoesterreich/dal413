import pandas as pd
import sys

def verify_full_history(history_file):
    """
    Reads a patient history file and verifies that the complete
    serial testing history is present for each patient.

    Args:
        history_file (str): Path to the patient testing history CSV file.
    """
    print(f"--- Verifying Test History in: {history_file} ---\n")
    try:
        history_df = pd.read_csv(history_file)
    except FileNotFoundError:
        print(f"Error: The file '{history_file}' was not found in the current directory.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        return

    # Get the list of patients in the file
    patient_ids = history_df['Effective Patient ID'].unique()

    print("Found the following patients in the file:", patient_ids.tolist())
    print("\nVerifying the 'Serial Test' numbers present for each patient...")
    print("-" * 30)

    # Group by patient and list their unique serial test numbers
    # The .sort_values() makes the output cleaner (e.g., [1, 2, 3] instead of [2, 1, 3])
    serial_test_summary = history_df.groupby('Effective Patient ID')['Serial Test'].unique().apply(lambda x: sorted(x))

    if not serial_test_summary.empty:
        print(serial_test_summary.to_string())
    else:
        print("No data to summarize.")
    
    print("-" * 30)
    print("\nVerification complete. The list above shows all test numbers included for each patient.")


if __name__ == '__main__':
    # Using an absolute path to be certain the input file is found.
    input_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/ERBB2_analysis_summary/serial_patients_full_test_history_v3.csv'
    verify_full_history(input_file)