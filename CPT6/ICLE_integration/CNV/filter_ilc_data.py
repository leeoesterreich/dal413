#!/usr/bin/env python3

import pandas as pd
import os
import sys

# File paths
tcga_sample_file = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/TCGA_sample.csv'
copy_number_file = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.txt'
output_file = '/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/CNV/ILC_data_CNA_gistic.txt'

def main():
    try:
        print("Step 1: Reading TCGA sample file...")
        # Read the TCGA sample file
        tcga_samples = pd.read_csv(tcga_sample_file)
        
        # Filter for Invasive lobular carcinoma
        print("Step 2: Filtering for Invasive lobular carcinoma samples...")
        ilc_samples = tcga_samples[tcga_samples.iloc[:, 1] == "Invasive lobular carcinoma"]
        
        if ilc_samples.empty:
            print("Error: No samples with 'Invasive lobular carcinoma' found!")
            sys.exit(1)
        
        # Extract IDs from the first column
        ilc_ids_raw = ilc_samples.iloc[:, 0].tolist()
        print(f"Found {len(ilc_ids_raw)} samples with Invasive lobular carcinoma")
        
        # Process IDs to match the format in the copy number file
        # The IDs in TCGA_sample.csv have "-01A" suffix, but the copy number file has "-01" suffix
        ilc_ids_processed = set()
        for raw_id in ilc_ids_raw:
            # Original ID
            ilc_ids_processed.add(raw_id)
            
            # Convert from TCGA-XX-XXXX-01A to TCGA-XX-XXXX-01
            if raw_id.endswith("-01A"):
                base_id = raw_id[:-1]  # Remove the 'A' at the end
                ilc_ids_processed.add(base_id)
            
            # Also add the base ID without any suffix
            parts = raw_id.split('-')
            if len(parts) >= 4:
                base_id = '-'.join(parts[:4])
                ilc_ids_processed.add(base_id)
        
        print(f"Processed {len(ilc_ids_processed)} variations of sample IDs for matching")
        
        # Read the copy number file header to get column names
        print("Step 3: Reading copy number data file header...")
        with open(copy_number_file, 'r') as f:
            header = f.readline().strip()
        
        cn_columns = header.split('\t')
        gene_symbol_col = cn_columns[0]  # First column is Gene Symbol
        
        # More efficient column matching
        print("Step 4: Identifying matching columns in copy number data...")
        # Create a list of columns to keep
        filtered_columns = [gene_symbol_col]
        matched_ids = []
        
        # For each column in the copy number file
        for col in cn_columns[1:]:
            # Check if this column matches any of our processed IDs
            for ilc_id in ilc_ids_processed:
                if ilc_id == col or col.startswith(ilc_id + '-') or ilc_id.startswith(col + '-'):
                    filtered_columns.append(col)
                    matched_ids.append(ilc_id)
                    break
        
        if len(filtered_columns) <= 1:
            print("Error: No matching columns found in copy number data!")
            print("Sample IDs in TCGA_sample.csv (first 5):", ilc_ids_raw[:5])
            print("Column headers in copy number file (first 5):", cn_columns[1:6])
            sys.exit(1)
            
        print(f"Found {len(filtered_columns) - 1} matching columns in copy number data")
        
        # Now read the actual data with only the columns we need
        print("Step 5: Reading and filtering copy number data...")
        # Use low_memory=False to avoid DtypeWarning for mixed types
        copy_number_data = pd.read_csv(copy_number_file, sep='\t', usecols=filtered_columns, low_memory=False)
        
        # Save the filtered data to the output file
        print(f"Step 6: Saving filtered data to {output_file}...")
        copy_number_data.to_csv(output_file, sep='\t', index=False)
        
        print(f"Done! Filtered data saved to {output_file}")
        print(f"Original samples: {len(ilc_ids_raw)}, Matched samples in output: {len(filtered_columns) - 1}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 