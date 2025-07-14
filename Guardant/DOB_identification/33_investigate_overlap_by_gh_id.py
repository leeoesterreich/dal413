import pandas as pd
import sys

def investigate_overlap_by_gh_id(ilc_file, idc_file):
    """
    Investigates the 7 overlapping patients by looking at their original GH_ID
    and cohort source.

    Args:
        ilc_file (str): Path to the ILC cohort data with patient IDs.
        idc_file (str): Path to the IDC cohort data with patient IDs.
    """
    try:
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # Find the overlapping Effective Patient IDs
    ilc_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_ids = set(df_idc['Effective Patient ID'].unique())
    overlap_ids = sorted(list(ilc_ids.intersection(idc_ids)))

    print("--- Overlapping Patient Investigation by GH_ID ---")
    print(f"Found {len(overlap_ids)} patients present in both ILC and IDC cohorts.")
    print("-" * 50)

    for patient_id in overlap_ids:
        print(f"\nInvestigating Patient ID: {patient_id}")
        
        # Get records from ILC cohort
        ilc_records = df_ilc[df_ilc['Effective Patient ID'] == patient_id][['GH_ID']]
        if not ilc_records.empty:
            for gh_id in ilc_records['GH_ID'].unique():
                print(f"  - Found in ILC cohort with GH_ID: {gh_id}")

        # Get records from IDC cohort
        idc_records = df_idc[df_idc['Effective Patient ID'] == patient_id][['GH_ID']]
        if not idc_records.empty:
            for gh_id in idc_records['GH_ID'].unique():
                print(f"  - Found in IDC cohort with GH_ID: {gh_id}")

if __name__ == "__main__":
    ILC_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_FILE = 'idc_cohort_with_patient_id.csv'
    
    investigate_overlap_by_gh_id(ILC_FILE, IDC_FILE) 