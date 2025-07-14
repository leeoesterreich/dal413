import pandas as pd
import sys

def find_intersection_overlap(genomic_file, ilc_file, idc_file):
    """
    Identifies and lists patients who are present in the genomic intersection
    AND in both the ILC and IDC histology cohorts.

    Args:
        genomic_file (str): Path to the first-true-alteration genomic data.
        ilc_file (str): Path to the ILC cohort data.
        idc_file (str): Path to the IDC cohort data.
    """
    try:
        df_genomic = pd.read_csv(genomic_file)
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # Define the patient ID sets
    genomic_ids = set(df_genomic['Effective Patient ID'].unique())
    ilc_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_ids = set(df_idc['Effective Patient ID'].unique())

    # Find the main intersection (histology patients with true alterations)
    histology_intersection = genomic_ids.intersection(ilc_ids.union(idc_ids))
    
    # Find the overlap between the two histology cohorts
    histology_overlap = ilc_ids.intersection(idc_ids)

    # Find which of the overlapping histology patients are ALSO in the main intersection
    final_overlap = histology_overlap.intersection(histology_intersection)

    print("--- Overlap Investigation Report ---")
    print(f"Total unique patients in the main intersection: {len(histology_intersection)}")
    print(f"Total unique patients in ILC intersection: {len(genomic_ids.intersection(ilc_ids))}")
    print(f"Total unique patients in IDC intersection: {len(genomic_ids.intersection(idc_ids))}")
    print("-" * 40)
    print(f"Number of patients in BOTH ILC and IDC cohorts (and in the intersection): {len(final_overlap)}")
    
    if final_overlap:
        print("Patient IDs present in both ILC and IDC cohorts:")
        for patient_id in sorted(list(final_overlap)):
            print(f"  - {patient_id}")

if __name__ == "__main__":
    GENOMIC_FILE = 'first_true_test_per_patient.csv'
    ILC_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_FILE = 'idc_cohort_with_patient_id.csv'
    
    find_intersection_overlap(GENOMIC_FILE, ILC_FILE, IDC_FILE) 