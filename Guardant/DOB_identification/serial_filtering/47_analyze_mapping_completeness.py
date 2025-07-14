import pandas as pd
import sys

def analyze_mapping_completeness(serial_data_file, mapping_file, master_patient_list):
    """
    Investigates if the patient mapping is incomplete by checking if patients
    in the cohort have other GH_IDs in the serial data that are not in the mapping.

    Args:
        serial_data_file (str): Path to the full serial genomic data.
        mapping_file (str): Path to the current complete patient mapping file.
        master_patient_list (str): Path to the master list linking all GH_IDs to Effective PIDs.
    """
    try:
        # --- 1. Load Data ---
        print("Step 1: Loading data files...")
        serial_df = pd.read_csv(serial_data_file)
        mapping_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
        master_list_df = pd.read_csv(master_patient_list)
        print("...data loading complete.")

        # --- 2. Data Preparation ---
        print("Step 2: Preparing data...")
        serial_df.columns = serial_df.columns.str.strip()
        
        # Get the set of GH_IDs that are in our current, histology-confirmed mapping
        mapped_gh_ids = set(mapping_df['GH_ID'])
        
        # Get the set of all patients (Effective Patient IDs) we care about
        cohort_patient_ids = set(mapping_df['Effective Patient ID'])
        print(f"Found {len(mapped_gh_ids)} GH_IDs for {len(cohort_patient_ids)} unique patients in the current mapping.")

        # From the master list, find ALL GH_IDs that belong to our cohort patients
        all_patient_gh_ids_df = master_list_df[master_list_df['Effective Patient ID'].isin(cohort_patient_ids)]
        all_patient_gh_ids = set(all_patient_gh_ids_df['GH_ID'])
        print(f"The master list shows that these patients have a total of {len(all_patient_gh_ids)} unique GH_IDs associated with them.")
        
        # Now, narrow this down to the GH_IDs that actually have data in the serial file
        gh_ids_in_serial_data = set(serial_df['GH_ID'])
        relevant_gh_ids_from_serial = all_patient_gh_ids.intersection(gh_ids_in_serial_data)
        print(f"Of those, {len(relevant_gh_ids_from_serial)} GH_IDs are present in the serial data file.")

        # --- 3. Find Missing GH_IDs ---
        print("\nStep 3: Identifying missed GH_IDs...")
        
        # The missing GH_IDs are those present in the serial data for our patients, but not in our mapping
        missing_gh_ids = relevant_gh_ids_from_serial - mapped_gh_ids

        # --- 4. Report Findings ---
        if not missing_gh_ids:
            print("\n--- Conclusion ---")
            print("Your hypothesis appears to be incorrect.")
            print("All GH_IDs for patients in our cohort that are present in the serial data are already in the mapping file.")
            print("The issue lies elsewhere.")
        else:
            print("\n--- Conclusion ---")
            print("Your hypothesis is CORRECT.")
            print(f"Found {len(missing_gh_ids)} GH_IDs that are in the serial data and belong to patients in our cohorts, but were NOT in the mapping file.")
            print("This confirms the mapping is incomplete and we are likely missing patient data.")
            
            # Show which patients these missing GH_IDs belong to
            missing_gh_patient_map = all_patient_gh_ids_df[all_patient_gh_ids_df['GH_ID'].isin(missing_gh_ids)]
            
            print("\nSample of missed GH_IDs and the patients they belong to:")
            print(missing_gh_patient_map.head(10).to_string())

            # Specifically check for the patient you asked about
            if 'A0229828' in missing_gh_ids:
                print("\n-> CONFIRMED: Patient GH_ID A0229828 was one of the missing IDs.")
                patient_info = missing_gh_patient_map[missing_gh_patient_map['GH_ID'] == 'A0229828']
                print("   It belongs to:", patient_info[['Effective Patient ID', 'Patient DOB']].to_string(index=False))


        print("\nScript finished.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    SERIAL_DATA = 'Guardant Project - Serial Genomic Data.csv'
    MAPPING_FILE = 'final_genomic_cohorts/complete_patient_mapping.csv'
    MASTER_LIST = 'patient_master_list.csv' # The original, most complete list from script 01
    
    analyze_mapping_completeness(SERIAL_DATA, MAPPING_FILE, MASTER_LIST) 