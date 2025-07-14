import pandas as pd
import sys
import os

def analyze_longitudinal_erbb2_mutations(genomic_file, mapping_file, output_file):
    """
    Analyzes serial genomic data to find patients with multiple distinct ERBB2 mutations over time,
    and writes the report to a file.

    Args:
        genomic_file (str): Path to the serial genomic data CSV file.
        mapping_file (str): Path to the patient mapping CSV file.
        output_file (str): Path to the output report file.
    """
    try:
        genomic_df = pd.read_csv(genomic_file, low_memory=False)
        mapping_df = pd.read_csv(mapping_file, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the files: {e}", file=sys.stderr)
        return

    # Merge the two dataframes on GH_ID
    merged_df = pd.merge(genomic_df, mapping_df, on='GH_ID', how='left')

    # Ensure necessary columns are present
    required_cols = ['Effective Patient ID', 'Gene', 'Alteration', 'Final Report Date', 'Histology', 'Vus']
    if not all(col in merged_df.columns for col in required_cols):
        print(f"Error: Missing one or more required columns. Found: {merged_df.columns.tolist()}", file=sys.stderr)
        return

    # Filter for ERBB2 mutations
    erbb2_df = merged_df[merged_df['Gene'] == 'ERBB2'].copy()
    if erbb2_df.empty:
        print("No ERBB2 mutations found in the data.")
        # We can still write this to the report file
        with open(output_file, 'w') as f:
            f.write("No ERBB2 mutations found in the data.\n")
        return
    
    # Convert date column to datetime objects
    erbb2_df['Final Report Date'] = pd.to_datetime(erbb2_df['Final Report Date'], errors='coerce')
    erbb2_df.dropna(subset=['Final Report Date'], inplace=True)

    # Identify distinct alterations for each patient
    patient_mutation_uniqueness = erbb2_df.groupby('Effective Patient ID')['Alteration'].nunique()
    
    # Find patients with more than one distinct mutation
    patients_with_multiple_mutations = patient_mutation_uniqueness[patient_mutation_uniqueness > 1]

    with open(output_file, 'w') as f:
        if patients_with_multiple_mutations.empty:
            f.write("No patients found with multiple distinct ERBB2 mutations over time.\n")
        else:
            num_patients = len(patients_with_multiple_mutations)
            f.write(f"Found {num_patients} unique patient(s) with multiple distinct ERBB2 mutations over time.\n")
            f.write("-" * 30 + "\n")
            f.write("Patient ID and number of unique ERBB2 mutations:\n")
            f.write(patients_with_multiple_mutations.to_string() + "\n")
            f.write("-" * 30 + "\n")

            patient_ids = patients_with_multiple_mutations.index
            detailed_info = erbb2_df[erbb2_df['Effective Patient ID'].isin(patient_ids)].copy()
            detailed_info.sort_values(by=['Effective Patient ID', 'Final Report Date'], inplace=True)
            
            display_cols = ['Effective Patient ID', 'GH_ID', 'Final Report Date', 'Gene', 'Alteration', 'c dot DNA', 'Alteration Name', 'Variant Allele Fraction', 'Histology', 'Vus']
            existing_display_cols = [col for col in display_cols if col in detailed_info.columns]

            f.write("Details of ERBB2 mutations for these patients:\n")
            f.write(detailed_info[existing_display_cols].to_string() + "\n")
    
    print(f"Report successfully generated at: {output_file}")


if __name__ == '__main__':
    genomic_data_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/Guardant Project - Serial Genomic Data.csv'
    patient_mapping_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/complete_patient_mapping.csv'
    output_report_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/erbb2_longitudinal_report_with_details.txt'
    
    analyze_longitudinal_erbb2_mutations(genomic_data_file, patient_mapping_file, output_report_file) 