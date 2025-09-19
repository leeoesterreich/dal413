import pandas as pd
import sys

def add_histology_and_reorder(report_file, mapping_file):
    """
    Adds histology data to an existing report and reorders the columns
    to place Histology and Effective Patient ID at the front. Overwrites
    the original report file.

    Args:
        report_file (str): Absolute path to the report to be modified.
        mapping_file (str): Absolute path to the patient mapping data.
    """
    print(f"--- Modifying Report: {report_file.split('/')[-1]} ---")
    try:
        report_df = pd.read_csv(report_file)
        mapping_df = pd.read_csv(mapping_file)
    except Exception as e:
        print(f"Error reading source files: {e}", file=sys.stderr)
        return

    # To add histology, we only need the Patient ID and Histology from the mapping file.
    # We drop duplicates to ensure we only have one histology per patient.
    histology_map_df = mapping_df[['Effective Patient ID', 'Histology']].drop_duplicates()

    # Merge the histology data into our report
    # We use a left merge to keep all rows from the original report.
    modified_df = pd.merge(report_df, histology_map_df, on='Effective Patient ID', how='left')

    # Check if 'Histology' was successfully merged before reordering
    if 'Histology' in modified_df.columns:
        # Reorder the columns
        all_cols = modified_df.columns.tolist()
        preferred_cols = ['Histology', 'Effective Patient ID']
        other_cols = [col for col in all_cols if col not in preferred_cols]
        final_col_order = preferred_cols + other_cols
        
        final_df = modified_df[final_col_order]
    else:
        print("Warning: 'Histology' column could not be added. The report will not be reordered.")
        final_df = modified_df

    try:
        # Overwrite the original file with the modified data
        final_df.to_csv(report_file, index=False)
        print(f"Successfully modified the report.")
        if 'Histology' in final_df.columns:
            print(f"Added 'Histology' and moved columns to the front.")
    except Exception as e:
        print(f"Error saving the modified report: {e}", file=sys.stderr)


if __name__ == '__main__':
    base = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification'
    report_to_modify = f'{base}/ERBB2_analysis_summary/multi_hit_non_vus_erbb2_patients.csv'
    mapping_data_file = f'{base}/final_genomic_cohorts/complete_patient_mapping.csv'

    add_histology_and_reorder(report_to_modify, mapping_data_file) 