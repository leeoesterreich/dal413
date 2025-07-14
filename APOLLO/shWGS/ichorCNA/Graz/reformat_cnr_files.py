import os
import pandas as pd
import numpy as np
import glob

def reformat_cnr_files():
    """
    Reads .annotated.cnr files from an input directory, reformats them to a standard structure,
    and saves them to an output directory.
    """
    input_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/cna_seg_files"
    output_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/cna_seg_files_reformatted"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cnr_files = glob.glob(os.path.join(input_dir, '*.annotated.cnr'))
    print("Found {} files to reformat.".format(len(cnr_files)))

    for file_path in cnr_files:
        file_name = os.path.basename(file_path)
        print("Processing {}...".format(file_name))

        try:
            df = pd.read_csv(file_path, sep='\t')
            
            # Determine sample ID from a column name (e.g., 'B123_1.copy.number')
            sample_id_col = [c for c in df.columns if '.copy.number' in c]
            if not sample_id_col:
                print("  - ERROR: Could not determine sample ID from columns in {}. Skipping.".format(file_name))
                continue
            sample_id = sample_id_col[0].replace('.copy.number', '')

            # --- Column Renaming ---
            df.rename(columns={'chr': 'chromosome', 'genes': 'gene'}, inplace=True, errors='ignore')
            
            logR_col = "{}.logR".format(sample_id)
            if logR_col in df.columns:
                df.rename(columns={logR_col: 'log2'}, inplace=True)
            else:
                df['log2'] = np.nan
                print("  - WARNING: '{}' not found. Added 'log2' column with NaN values.".format(logR_col))

            # --- Data Manipulation ---
            if 'gene' in df.columns:
                df['gene'] = df['gene'].str.replace(';', ',')
            else:
                df['gene'] = '-'
                print("  - WARNING: 'genes' column not found. Added 'gene' column with '-' values.")

            # --- Column Reordering ---
            core_cols = ['chromosome', 'start', 'end', 'gene', 'log2']
            existing_cols = df.columns.tolist()
            
            new_order = [col for col in core_cols if col in existing_cols]
            other_cols = [col for col in existing_cols if col not in new_order]
            new_order.extend(other_cols)
            
            df = df[new_order]

            # --- Save file ---
            output_path = os.path.join(output_dir, file_name.replace(".annotated.cnr", ".cnr"))
            df.to_csv(output_path, sep='\t', index=False, na_rep='NA')
            print("  - Saved reformatted file to {}".format(output_path))

        except Exception as e:
            print("  - ERROR processing {}: {}".format(file_name, e))

if __name__ == "__main__":
    reformat_cnr_files() 