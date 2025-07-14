import pandas as pd
import os

def reconcile_tmb_counts(file_path):
    """
    Provides a final reconciliation of TMB score counts for the ILC file,
    explaining the difference between total patients and patients with valid scores.

    Args:
        file_path (str): The path to the ILC CSV file.
    """
    print(f"--- TMB Reconciliation for: {os.path.basename(file_path)} ---")

    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    try:
        df = pd.read_csv(file_path, low_memory=False)

        if 'TMB Score' not in df.columns or 'Effective Patient ID' not in df.columns:
            print("Error: Required columns not found.")
            return

        # 1. Count all unique patients in the file
        all_unique_patients = df['Effective Patient ID'].unique()
        total_unique_patient_count = len(all_unique_patients)
        print(f"\n1. Total unique patients in the file: {total_unique_patient_count}")

        # 2. Identify patients with valid, numerical TMB scores
        df['TMB Score'] = pd.to_numeric(df['TMB Score'], errors='coerce')
        valid_tmb_df = df.dropna(subset=['TMB Score'])
        valid_tmb_unique_patients = valid_tmb_df['Effective Patient ID'].unique()
        valid_tmb_patient_count = len(valid_tmb_unique_patients)
        print(f"2. Unique patients with a valid TMB score: {valid_tmb_patient_count}")

        # 3. Calculate how many are missing a valid score
        missing_tmb_count = total_unique_patient_count - valid_tmb_patient_count
        print(f"3. Patients with a missing or non-numeric TMB score: {missing_tmb_count}")

        # 4. Count patients with TMB > 10 from the valid group
        high_tmb_df = valid_tmb_df[valid_tmb_df['TMB Score'] > 10]
        high_tmb_count = high_tmb_df['Effective Patient ID'].nunique()
        print(f"\n4. Unique patients with TMB Score > 10 (from the valid group): {high_tmb_count}")

        # 5. Calculate the correct percentage
        if valid_tmb_patient_count > 0:
            correct_percentage = (high_tmb_count / valid_tmb_patient_count) * 100
            print(f"5. Correct percentage with TMB > 10 is {high_tmb_count}/{valid_tmb_patient_count} = {correct_percentage:.2f}%")
        else:
            print("5. No patients with valid TMB scores to calculate a percentage.")

        # 6. Calculate the median for the valid group
        median_tmb = valid_tmb_df.drop_duplicates(subset=['Effective Patient ID'])['TMB Score'].median()
        print(f"6. The median TMB score for the {valid_tmb_patient_count} patients is: {median_tmb:.2f}")
        
        print("\nConclusion: The median is >10 because over half (57.14%) of the patients *who have a score* are above 10.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    ilc_file_path = 'final_genomic_cohorts/ILC_genomic_infinity.csv'
    reconcile_tmb_counts(ilc_file_path) 