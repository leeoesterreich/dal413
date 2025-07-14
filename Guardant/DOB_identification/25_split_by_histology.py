import pandas as pd
import sys

def split_and_enrich_by_histology(main_file, ilc_file, idc_file, ilc_output, idc_output):
    """
    Splits genomic data by cancer histology (ILC vs. IDC) and enriches it
    with histology descriptions.

    Args:
        main_file (str): Path to the main genomic data CSV.
        ilc_file (str): Path to the ILC cohort CSV.
        idc_file (str): Path to the IDC cohort CSV.
        ilc_output (str): Path to save the ILC genomic data.
        idc_output (str): Path to save the IDC genomic data.
    """
    try:
        df_main = pd.read_csv(main_file)
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # Prepare histology mappings from cohort files
    ilc_histology_map = df_ilc.set_index('Effective Patient ID')['Histo/Behavior ICD-O-3-Desc'].to_dict()
    idc_histology_map = df_idc.set_index('Effective Patient ID')['Histo/Behavior ICD-O-3-Desc'].to_dict()

    # Create the histology description column in the main dataframe
    df_main['Histo/Behavior ICD-O-3-Desc'] = df_main['Effective Patient ID'].map(ilc_histology_map).fillna(
        df_main['Effective Patient ID'].map(idc_histology_map)
    )

    # Split the main dataframe into ILC and IDC based on the patient IDs in each cohort
    df_ilc_genomic = df_main[df_main['Effective Patient ID'].isin(ilc_histology_map.keys())].copy()
    df_idc_genomic = df_main[df_main['Effective Patient ID'].isin(idc_histology_map.keys())].copy()

    # Function to reorder columns
    def reorder_and_save(df, output_path, cancer_type_str):
        if 'Cancer Type' in df.columns:
            # Find the index of 'Cancer Type'
            cancer_type_idx = df.columns.get_loc('Cancer Type')
            # Create the new column order
            new_cols = list(df.columns)
            histology_col = new_cols.pop(new_cols.index('Histo/Behavior ICD-O-3-Desc'))
            new_cols.insert(cancer_type_idx + 1, histology_col)
            df = df[new_cols]
        
        df.to_csv(output_path, index=False)
        print(f"Successfully created {cancer_type_str} genomic file with {len(df)} records.")
        print(f"Data saved to {output_path}")

    # Process and save ILC data
    reorder_and_save(df_ilc_genomic, ilc_output, "ILC")

    # Process and save IDC data
    reorder_and_save(df_idc_genomic, idc_output, "IDC")


if __name__ == "__main__":
    MAIN_GENOMIC_FILE = 'first_true_test_per_patient.csv'
    ILC_COHORT_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_COHORT_FILE = 'idc_cohort_with_patient_id.csv'
    ILC_OUTPUT_FILE = 'ILC_genomic.csv'
    IDC_OUTPUT_FILE = 'IDC_genomic.csv'
    
    split_and_enrich_by_histology(MAIN_GENOMIC_FILE, ILC_COHORT_FILE, IDC_COHORT_FILE, ILC_OUTPUT_FILE, IDC_OUTPUT_FILE) 