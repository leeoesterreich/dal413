import os
import pandas as pd
import numpy as np
import re

def main():
    """
    Collapses longitudinal segment data to the patient level by selecting the
    segments from the sample with the highest tumor fraction for each patient.
    The output is a GISTIC-compatible .seg file.
    """
    print("--- Collapsing Segments to Patient Level ---")

    # --- 1. Load Tumor Fraction and Patient-Sample Mapping ---
    # This logic is identical to our previous analysis scripts.
    file_path = 'results/TFx_hg38_1000kb.xlsx'
    df_tf = pd.read_excel(file_path, sheet_name=0, skiprows=1)
    df_tf.columns = df_tf.columns.map(str)
    df_tf = df_tf.loc[:, ~df_tf.columns.str.startswith('Tumor_fraction_1000kb')]
    df_tf.rename(columns=lambda x: 'Tumor_fraction' if x.startswith('Tumor_fraction_500kb') else x, inplace=True)
    
    patient_to_samples = {}
    for _, row in df_tf.iterrows():
        patient_code = str(row["Patient_Code"])
        tuples = []
        for i in range(1, len(df_tf.columns) // 2 + 1):
            if i == 1: sample_label_col = f"{i}st progression"
            elif i == 2: sample_label_col = f"{i}nd progression"
            elif i == 3: sample_label_col = f"{i}rd progression"
            else: sample_label_col = f"{i}th progression"
            if sample_label_col not in df_tf.columns: continue
            
            sample_label = row[sample_label_col]
            tumor_fraction = row.iloc[df_tf.columns.get_loc(sample_label_col) + 1]

            if pd.notna(sample_label) and pd.notna(tumor_fraction):
                # Use the full sample label from the excel sheet
                tuples.append((sample_label, tumor_fraction))
        patient_to_samples[patient_code] = tuples
    
    print(f"Loaded tumor fraction data for {len(patient_to_samples)} patients.")

    # --- Exclude specified patient ---
    patient_to_exclude = '4908'
    if patient_to_exclude in patient_to_samples:
        del patient_to_samples[patient_to_exclude]
        print(f"Patient {patient_to_exclude} has been excluded from the analysis.")
    else:
        print(f"Warning: Patient {patient_to_exclude} to exclude was not found.")

    # --- 2. Find Max TF Sample and Select Segments for Each Patient ---
    
    base_dir = "results/ichorCNA"
    collapsed_segments_list = []
    
    for patient, samples in patient_to_samples.items():
        if not samples:
            continue
            
        # Find the sample with the highest tumor fraction for this patient
        max_tf_sample_id, max_tf = max(samples, key=lambda item: item[1])

        # Construct the path to the corresponding .cna.seg file
        # The file is often named with a suffix, e.g., SAMPLE_ID_ichorCNA
        # We need to find the directory that starts with the sample ID.
        
        matching_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.startswith(max_tf_sample_id)]
        
        if not matching_dirs:
            print(f"Warning: No ichorCNA directory found for max TF sample {max_tf_sample_id} of patient {patient}. Skipping.")
            continue

        # Assuming the first match is the correct one
        sample_dir = matching_dirs[0]
        seg_file_path = os.path.join(base_dir, sample_dir, f"{sample_dir}.cna.seg")

        if not os.path.exists(seg_file_path):
            print(f"Warning: Seg file not found at {seg_file_path}. Skipping.")
            continue

        try:
            # Read the representative .seg file
            df_seg = pd.read_csv(seg_file_path, sep='\t')
            
            # Find the logR column (which is the segment mean)
            logr_col = [col for col in df_seg.columns if 'logR' in col][0]
            
            # Create the GISTIC-formatted DataFrame
            df_gistic = pd.DataFrame({
                'Sample': patient, # Use patient ID instead of sample ID
                'Chromosome': df_seg['chr'].str.replace('chr', ''),
                'Start': df_seg['start'],
                'End': df_seg['end'],
                'Num_Probes': df_seg['end'].sub(df_seg['start']).div(1000).round().astype(int), # Num markers is not always present, so we estimate
                'Segment_Mean': df_seg[logr_col]
            })
            
            collapsed_segments_list.append(df_gistic)
            print(f"Patient {patient}: Using segments from sample {max_tf_sample_id} (TF={max_tf:.3f})")

        except Exception as e:
            print(f"Error processing file {seg_file_path} for patient {patient}: {e}")

    # --- 3. Combine and Save Final Collapsed .seg File ---
    if not collapsed_segments_list:
        print("No segments could be processed. Exiting.")
        return

    final_collapsed_df = pd.concat(collapsed_segments_list, ignore_index=True)

    # GISTIC format requires specific chromosome notation (X=23, Y=24)
    final_collapsed_df['Chromosome'] = final_collapsed_df['Chromosome'].replace({'X': '23', 'Y': '24'})
    
    # Ensure chromosome column is integer for correct sorting
    final_collapsed_df['Chromosome'] = pd.to_numeric(final_collapsed_df['Chromosome'])
    final_collapsed_df = final_collapsed_df.sort_values(by=['Sample', 'Chromosome', 'Start'])

    output_file = "results/all_segments_collapsed.seg"
    final_collapsed_df.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"\nSuccessfully created collapsed segment file for GISTIC.")
    print(f"Total patients processed: {len(collapsed_segments_list)}")
    print(f"Total segments written: {len(final_collapsed_df)}")
    print(f"Output file: {output_file}")

if __name__ == "__main__":
    main() 