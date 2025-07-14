import pandas as pd
import sys

def clean_and_reorder_report(report_file):
    """
    Cleans up a report by consolidating duplicate 'Histology' columns
    and reordering the final columns. Overwrites the original file.

    Args:
        report_file (str): Absolute path to the report to be cleaned.
    """
    print(f"--- Cleaning Report: {report_file.split('/')[-1]} ---")
    try:
        report_df = pd.read_csv(report_file)
    except Exception as e:
        print(f"Error reading source file: {e}", file=sys.stderr)
        return

    # Consolidate the two histology columns.
    # It seems one might be empty (NaN), so we use .fillna() to merge them.
    # This takes the value from _x, and if it's missing, takes the value from _y.
    if 'Histology_x' in report_df.columns and 'Histology_y' in report_df.columns:
        report_df['Histology'] = report_df['Histology_x'].fillna(report_df['Histology_y'])
        # Drop the old, messy columns
        report_df.drop(columns=['Histology_x', 'Histology_y'], inplace=True)
        print("Successfully consolidated 'Histology' columns.")
    else:
        print("Warning: Could not find duplicate histology columns to clean.")


    # Reorder the columns to place the clean Histology and Patient ID at the front
    all_cols = report_df.columns.tolist()
    preferred_cols = ['Histology', 'Effective Patient ID']
    other_cols = [col for col in all_cols if col not in preferred_cols]
    
    # Ensure the preferred columns actually exist before trying to order them
    final_col_order = [col for col in preferred_cols if col in report_df.columns] + other_cols
    
    final_df = report_df[final_col_order]

    try:
        # Overwrite the original file with the cleaned data
        final_df.to_csv(report_file, index=False)
        print(f"Successfully cleaned and reordered the report.")
    except Exception as e:
        print(f"Error saving the cleaned report: {e}", file=sys.stderr)


if __name__ == '__main__':
    base = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification'
    report_to_clean = f'{base}/ERBB2_analysis_summary/multi_hit_non_vus_erbb2_patients.csv'

    clean_and_reorder_report(report_to_clean) 