import pandas as pd
import sys
import os

def analyze_non_vus_erbb2_hits(genomic_file, mapping_file, output_csv_file):
    """
    Analyzes serial genomic data to find patients with multiple non-VUS ERBB2 mutations,
    and saves the filtered data to a new CSV file.

    Args:
        genomic_file (str): Path to the serial genomic data CSV file.
        mapping_file (str): Path to the patient mapping CSV file.
        output_csv_file (str): Path for the output CSV file.
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

    # Merge the dataframes
    merged_df = pd.merge(genomic_df, mapping_df, on='GH_ID', how='left')

    # Filter for ERBB2 mutations
    erbb2_df = merged_df[merged_df['Gene'] == 'ERBB2'].copy()
    if erbb2_df.empty:
        print("No ERBB2 mutations found in the data.")
        return

    # Filter for non-VUS mutations (Vus == FALSE)
    if 'Vus' not in erbb2_df.columns:
        print("Error: 'Vus' column not found in the data.", file=sys.stderr)
        return
    
    non_vus_df = erbb2_df[erbb2_df['Vus'] == 'FALSE'].copy()

    if non_vus_df.empty:
        print("No non-VUS ERBB2 mutations found.")
        return

    # Find patients with multiple non-VUS hits
    mutation_counts = non_vus_df.groupby('Effective Patient ID').size()
    patients_with_multiple_hits = mutation_counts[mutation_counts > 1].index

    if patients_with_multiple_hits.empty:
        print("No patients found with multiple non-VUS ERBB2 mutations.")
        return

    # Filter the non-VUS dataframe to include only these patients
    final_df = non_vus_df[non_vus_df['Effective Patient ID'].isin(patients_with_multiple_hits)].copy()

    # Add a column with the total number of hits for each patient
    final_df['Multiple Hit Count'] = final_df.groupby('Effective Patient ID')['Effective Patient ID'].transform('size')

    # Sort the data for better readability
    final_df.sort_values(by=['Effective Patient ID', 'Final Report Date'], inplace=True)
    
    # Save the final dataframe to a CSV file
    try:
        final_df.to_csv(output_csv_file, index=False)
        print(f"Analysis complete. Report for patients with multiple non-VUS ERBB2 hits saved to:")
        print(output_csv_file)
        
        # Also print a summary to the console
        summary = final_df.groupby('Effective Patient ID')['Multiple Hit Count'].first()
        print("\nSummary of patients with multiple hits:")
        print(summary)

    except Exception as e:
        print(f"An error occurred while saving the CSV file: {e}", file=sys.stderr)


if __name__ == '__main__':
    genomic_data_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/Guardant Project - Serial Genomic Data.csv'
    patient_mapping_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/complete_patient_mapping.csv'
    output_csv = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/multi_hit_non_vus_erbb2_patients.csv'
    
    analyze_non_vus_erbb2_hits(genomic_data_file, patient_mapping_file, output_csv) 