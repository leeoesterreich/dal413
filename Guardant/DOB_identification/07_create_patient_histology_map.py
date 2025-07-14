import pandas as pd

def create_patient_map(master_list_path, ilc_cohort_path, idc_cohort_path, output_path):
    """
    Creates a mapping file containing GH_ID, Effective Patient ID, Histological Subtype, and DOB.

    Args:
        master_list_path (str): Path to the patient master list.
        ilc_cohort_path (str): Path to the ILC cohort with IDs.
        idc_cohort_path (str): Path to the IDC cohort with IDs.
        output_path (str): Path to save the output CSV.
    """
    try:
        master_list = pd.read_csv(master_list_path)
        ilc_cohort = pd.read_csv(ilc_cohort_path)
        idc_cohort = pd.read_csv(idc_cohort_path)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found.")
        print(e)
        return

    # Add histology information
    ilc_cohort['Histological Subtype'] = 'ILC'
    idc_cohort['Histological Subtype'] = 'IDC'

    # Combine histology information
    ilc_histology = ilc_cohort[['Effective Patient ID', 'Histological Subtype']].drop_duplicates()
    idc_histology = idc_cohort[['Effective Patient ID', 'Histological Subtype']].drop_duplicates()
    
    combined_histology = pd.concat([ilc_histology, idc_histology])

    # For patients in multiple cohorts, aggregate the subtype into a single string
    aggregated_histology = combined_histology.groupby('Effective Patient ID')['Histological Subtype'].apply(lambda x: ', '.join(sorted(x))).reset_index()

    # Merge with the master list to get GH_ID and DOB
    final_map = pd.merge(master_list, aggregated_histology, on='Effective Patient ID', how='left')
    
    # Reorder and select final columns
    final_map = final_map[['GH_ID', 'Effective Patient ID', 'Patient DOB', 'Histological Subtype']]

    # Fill missing histology with a placeholder
    final_map['Histological Subtype'].fillna('No Histology Found', inplace=True)

    final_map.to_csv(output_path, index=False)
    print(f"Patient map saved to '{output_path}'")
    print(f"Total patients in map: {len(final_map)}")
    print("\\nValue counts for Histological Subtype:")
    print(final_map['Histological Subtype'].value_counts(dropna=False))


if __name__ == "__main__":
    MASTER_LIST_FILE = "patient_master_list.csv"
    ILC_COHORT_FILE = "ilc_cohort_with_ids.csv"
    IDC_COHORT_FILE = "idc_cohort_with_ids.csv"
    OUTPUT_FILE = "patient_histology_map.csv"
    create_patient_map(MASTER_LIST_FILE, ILC_COHORT_FILE, IDC_COHORT_FILE, OUTPUT_FILE) 