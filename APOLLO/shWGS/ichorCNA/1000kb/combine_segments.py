#!/usr/bin/env python

import os
import pandas as pd
import glob

def process_cna_seg_file(file_path):
    """Process a single .cna.seg file and return a DataFrame with required columns"""
    # Extract sample name from the path
    sample = os.path.basename(os.path.dirname(file_path))
    
    # Read the .seg file
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        # Get the logR column name (it contains the sample name)
        logr_col = [col for col in df.columns if 'logR' in col and not 'Copy_Number' in col][0]
        copy_number_col = [col for col in df.columns if 'copy.number' in col][0]
        
        # Select and rename required columns
        df_selected = pd.DataFrame({
            'sample': sample,
            'chr': df['chr'],
            'start': df['start'],
            'end': df['end'],
            'bins': df['end'].sub(df['start']).div(1e6).round().astype(int),  # Calculate bins as (end-start)/1e6
            'median': df[logr_col]  # Using logR as median
        })
        
        return df_selected
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

def process_seg_file(file_path):
    """Process a single .seg file and return a DataFrame with required columns"""
    try:
        # Read the .seg file
        df = pd.read_csv(file_path, sep='\t')
        
        # For .seg files, we expect these columns: sample, chr, start, end, event, copy.number, bins, median
        df_selected = pd.DataFrame({
            'Sample': df['sample'],
            'Chromosome': df['chr'].str.replace('chr', '').map(lambda x: '23' if x == 'X' else ('24' if x == 'Y' else x)),
            'Start': df['start'],
            'End': df['end'],
            'Num_Probes': df['bins'],
            'Seg.CN': df['median']  # Using median logR value
        })
        
        return df_selected
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

def process_files(file_pattern, process_func, output_file):
    """Process all files matching the pattern and save to output file"""
    # Find all matching files
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files found matching pattern: {file_pattern}")
        return
    
    print(f"Found {len(files)} files matching pattern: {file_pattern}")
    
    # Process all files
    dfs = []
    for file in files:
        df = process_func(file)
        if df is not None:
            dfs.append(df)
    
    if not dfs:
        print("No data could be processed!")
        return
    
    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save to output file
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined data saved to {output_file}")
    print(f"Total number of segments: {len(combined_df)}")

def main():
    # Base directory containing ichorCNA results
    base_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/ichorCNA"
    
    # Process .cna.seg files
    cna_pattern = os.path.join(base_dir, "*/*.cna.seg")
    cna_output = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/all_segments_cna.txt"
    process_files(cna_pattern, process_cna_seg_file, cna_output)
    
    # Process .seg files (excluding .cna.seg files)
    seg_files = []
    for f in glob.glob(os.path.join(base_dir, "*/*.seg")):
        if not f.endswith('.cna.seg'):
            seg_files.append(f)
    
    print(f"Found {len(seg_files)} .seg files (excluding .cna.seg files)")
    
    # Process all files
    dfs = []
    for file in seg_files:
        df = process_seg_file(file)
        if df is not None:
            dfs.append(df)
    
    if not dfs:
        print("No .seg data could be processed!")
        return
    
    # Combine all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save to output file - NO HEADER for GISTIC
    seg_output = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/all_segments_seg.txt"
    combined_df.to_csv(seg_output, sep='\t', index=False, header=False)
    print(f"Combined .seg data saved to {seg_output}")
    print(f"Total number of segments: {len(combined_df)}")

if __name__ == "__main__":
    main() 