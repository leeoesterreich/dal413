import pandas as pd
import argparse
import os

def assign_patient_ids(input_file, output_file, name_column):
    """
    Assigns a unique numeric Patient ID to each unique name in a specified column.

    Args:
        input_file (str): Path to the input CSV file. #Replace with your path to csv
        output_file (str): Path to save the output CSV file. #Replace your desired output path
        name_column (str): The name of the column containing patient names. #Replace with the exact column name, case sensitive
    """
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    if name_column not in df.columns:
        print(f"Error: Column '{name_column}' not found in the input file.")
        print(f"Available columns are: {list(df.columns)}")
        return

    # Create a unique numeric ID for each unique value in the name column
    df['Patient ID'] = pd.factorize(df[name_column])[0] + 1
    df['Patient ID'] = 'P' + df['Patient ID'].astype(str)

    # Move 'Patient ID' to the first column
    cols = ['Patient ID'] + [col for col in df if col != 'Patient ID']
    df = df[cols]

    # Save the updated dataframe
    try:
        df.to_csv(output_file, index=False)
        print(f"Successfully processed the file.")
        print(f"Output with Patient IDs saved to: '{output_file}'")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Assign a unique numeric Patient ID to each unique name in a CSV file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--input_file',
        type=str,
        required=True,
        help="Path to the input CSV file."
    )
    parser.add_argument(
        '--name_column',
        type=str,
        required=True,
        help="The name of the column containing the names to create IDs from."
    )
    parser.add_argument(
        '--output_file',
        type=str,
        help=(
            "Path to the output CSV file.\n"
            "If not provided, it defaults to '<input_filename>_with_ids.csv'."
        )
    )

    args = parser.parse_args()

    # If no output file is specified, create a default name
    if not args.output_file:
        base, ext = os.path.splitext(args.input_file)
        args.output_file = f"{base}_with_ids{ext}"

    assign_patient_ids(args.input_file, args.output_file, args.name_column) 