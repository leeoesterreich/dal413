import pandas as pd
import sys

def analyze_patient_subsets(multi_hit_file, mapping_file):
    """
    Analyzes a subset of patients for their histology (ILC/IDC) and
    determines which of them have had serial testing.

    Args:
        multi_hit_file (str): Path to the CSV file of patients with multiple hits.
        mapping_file (str): Path to the patient mapping CSV file.
    """
    try:
        multi_hit_df = pd.read_csv(multi_hit_file)
        mapping_df = pd.read_csv(mapping_file)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the files: {e}", file=sys.stderr)
        return

    # Get the unique list of patients from the multi-hit file
    patient_ids = multi_hit_df['Effective Patient ID'].unique()
    print(f"Analyzing {len(patient_ids)} unique patients from the report.\n")

    # --- Part 1: Analyze Histology ---
    
    # Filter the mapping file for our specific patients
    patient_histology_df = mapping_df[mapping_df['Effective Patient ID'].isin(patient_ids)]
    
    # We only need one histology entry per patient, so we drop duplicates
    unique_histology_df = patient_histology_df.drop_duplicates(subset=['Effective Patient ID'])
    
    if 'Histology' in unique_histology_df.columns:
        histology_counts = unique_histology_df['Histology'].value_counts()
        print("--- Histology Summary ---")
        print("Number of patients by histology type:")
        print(histology_counts)
        print("-" * 25 + "\n")
    else:
        print("Warning: 'Histology' column not found in the mapping file.")


    # --- Part 2: Identify Patients with Serial Testing ---
    
    # A patient has serial testing if they have more than one unique report date.
    if 'Final Report Date' in multi_hit_df.columns:
        # Convert to datetime to ensure proper comparison
        multi_hit_df['Final Report Date'] = pd.to_datetime(multi_hit_df['Final Report Date'], errors='coerce')
        
        # Count unique dates for each patient
        serial_test_counts = multi_hit_df.groupby('Effective Patient ID')['Final Report Date'].nunique()
        
        # Filter for patients with more than one unique date
        patients_with_serial_testing = serial_test_counts[serial_test_counts > 1]
        
        print("--- Serial Testing Summary ---")
        if not patients_with_serial_testing.empty:
            print(f"Found {len(patients_with_serial_testing)} patients with serial testing (tests on multiple dates).")
            print("Patients and their number of unique test dates:")
            print(patients_with_serial_testing)
        else:
            print("No patients with serial testing (tests on multiple dates) were found in this subset.")
        print("-" * 25 + "\n")

    else:
        print("Warning: 'Final Report Date' column not found in the multi-hit file.")


if __name__ == '__main__':
    multi_hit_patient_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/multi_hit_non_vus_erbb2_patients.csv'
    patient_mapping_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts/complete_patient_mapping.csv'
    
    analyze_patient_subsets(multi_hit_patient_file, patient_mapping_file) 