import pandas as pd
import sys

def update_serial_file(master_file, raw_serial_file, output_file):
    """
    Applies the definitive Effective Patient IDs to the raw serial data file.

    Args:
        master_file (str): Path to the master patient list (e.g., 'patient_master_list.csv').
        raw_serial_file (str): Path to the original, raw serial genomic data.
        output_file (str): Path to save the updated serial data.
    """
    try:
        print("Reading master patient list and raw serial data...")
        df_master = pd.read_csv(master_file)
        df_raw_serial = pd.read_csv(raw_serial_file)
        print("...done.")
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}", file=sys.stderr)
        sys.exit(1)

    # Create the mapping from GH_ID to the new Effective Patient ID
    id_map = df_master.set_index('GH_ID')['Effective Patient ID'].to_dict()

    print("Mapping new Effective Patient IDs to serial data...")
    df_raw_serial['Effective Patient ID'] = df_raw_serial['GH_ID'].map(id_map)

    # Report on any patients in the serial file that couldn't be mapped
    unmapped_patients = df_raw_serial[df_raw_serial['Effective Patient ID'].isna()]['GH_ID'].nunique()
    if unmapped_patients > 0:
        print(f"Warning: Found {unmapped_patients} GH_IDs in the serial file that were not in the master patient list.")
    
    # Reorder columns to bring the new ID to the front
    cols = ['Effective Patient ID'] + [col for col in df_raw_serial.columns if col != 'Effective Patient ID']
    df_final = df_raw_serial[cols]

    # Save the updated file
    df_final.to_csv(output_file, index=False)
    print(f"Successfully updated serial data. Output saved to {output_file}")


if __name__ == "__main__":
    MASTER_FILE = 'patient_master_list.csv'
    RAW_SERIAL_FILE = 'Guardant Project - Genomic Data_DOB_identification.csv'
    OUTPUT_FILE = 'serial_genomic_with_patient_id.csv' # Overwriting the old one

    update_serial_file(MASTER_FILE, RAW_SERIAL_FILE, OUTPUT_FILE) 