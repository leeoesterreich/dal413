import pandas as pd
import os
import glob

def create_gistic_file():
    # Paths
    clinical_data_path = 'CtDNAHRHER2Novartis-ClinicalNina_DATA_LABELS_2025-05-14_1651.csv'
    seg_files_dir = 'ichor_result/'
    output_dir = 'GISTIC/'
    output_file = os.path.join(output_dir, 'all_IDC_segments.txt')

    # Load clinical data and filter for IDC patients
    try:
        clinical_df = pd.read_csv(clinical_data_path)
    except FileNotFoundError:
        print("Error: Clinical data file not found at {}".format(clinical_data_path))
        return

    # Assuming 'Histology' and 'Patient ID' are the correct column names.
    # The user can modify this if the script fails.
    if 'Histology' not in clinical_df.columns or 'HG-ID' not in clinical_df.columns:
        print("Error: 'Histology' or 'HG-ID' column not found in clinical data file.")
        print("Available columns: {}".format(clinical_df.columns.tolist()))
        return
        
    idc_patients = clinical_df[clinical_df['Histology'] == 'NST']['HG-ID'].unique()
    print("Found {} NST patients.".format(len(idc_patients)))

    # Find all seg.txt files
    seg_files = glob.glob(os.path.join(seg_files_dir, '*.seg.txt'))
    print("Found {} .seg.txt files.".format(len(seg_files)))

    all_segments_df = pd.DataFrame()

    # Process each seg file
    for seg_file in seg_files:
        filename = os.path.basename(seg_file)
        
        # Extract patient ID from filename
        patient_id_from_file = filename.split('_')[0].split('-')[0]
        sample_id = filename.replace('.seg.txt', '')

        if patient_id_from_file in idc_patients:
            try:
                seg_df = pd.read_csv(seg_file, sep='\t', engine='python')
                
                # Check for expected columns
                # GISTIC format: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
                # ichorCNA seg file format: ID, chrom, start, end, num.mark, seg.median.logR
                required_cols = {'ID', 'chrom', 'start', 'end', 'num.mark', 'seg.median.logR'}
                if not required_cols.issubset(seg_df.columns):
                    print("Warning: Skipping {}. Missing one or more required columns: {}".format(filename, required_cols))
                    print("Available columns: {}".format(seg_df.columns.tolist()))
                    continue

                # Reformat for GISTIC
                gistic_df = seg_df[['chrom', 'start', 'end', 'num.mark', 'seg.median.logR']].copy()
                gistic_df.insert(0, 'Sample', sample_id)
                
                all_segments_df = pd.concat([all_segments_df, gistic_df], ignore_index=True)

            except Exception as e:
                print("Error processing file {}: {}".format(filename, e))

    if not all_segments_df.empty:
        # Rename columns to match GISTIC format example
        all_segments_df.columns = ['Sample', 'Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean']
        
        # Save to file
        all_segments_df.to_csv(output_file, sep='\t', index=False, header=True)
        print("Successfully created GISTIC input file at: {}".format(output_file))
        print("Total segments: {}".format(len(all_segments_df)))
        print("Number of unique patients in final file: {}".format(all_segments_df['Sample'].nunique()))
    else:
        print("No data was processed to generate the GISTIC file.")

if __name__ == '__main__':
    create_gistic_file() 