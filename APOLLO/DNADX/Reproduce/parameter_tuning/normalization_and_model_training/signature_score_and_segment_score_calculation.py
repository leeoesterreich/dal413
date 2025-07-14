'''
Calculates RNA and CNA signature scores from raw data using specified GMT files.
Allows skipping log2 transformation for RNA and assumes CNA is already appropriately scaled.
Accepts input file paths and an output directory via command-line arguments.
'''
import pandas as pd
import numpy as np
import os
import argparse # Added for command-line arguments
from helper import (
    exp_preprocess_for_rna, # This helper might need adjustment or careful use if we skip log2
    calculate_signature_scores,
    load_gmt, 
    convert_symbol_to_entrez_for_data_index
)

# Define a gene ID mapping function
def load_gene_id_mapping(mapping_file_path):
    return pd.read_csv(mapping_file_path, sep='\t', header=None, names=['GeneID_from_file', 'Symbol_from_file'])

def process_dataset_pair(
    rna_input_path, 
    cna_input_path, 
    mapping_df, 
    rna_gmt_file, 
    cna_gmt_file, 
    output_dir, 
    dataset_prefix,
    g_rna_skip_log2_transform, # Global flag from args
    g_cna_already_log2_assumed # Global flag from args
    ):
    """
    Processes a pair of RNA and CNA datasets (e.g., training or validation).
    g_rna_skip_log2_transform: If True, RNA data is NOT log2 transformed by this script.
    g_cna_already_log2_assumed: If True, CNA data is assumed to be already log2 (or not needing it), 
                                so this script does NOT log2 transform it.
    """
    print(f"\n--- Processing {dataset_prefix} Data ---M--")

    # --- RNA Signature Score Calculation ---
    print(f"--- {dataset_prefix} RNA Processing ---")
    if rna_input_path:
        print(f"Loading {dataset_prefix} RNA expression data from {rna_input_path}...")
        # Check if rna_input_path is a directory (like HiSeqV2) or a file
        if os.path.isdir(rna_input_path):
            # Assuming HiSeqV2 format: a directory containing a file, or the file itself needs to be found
            # This part might need specific logic if HiSeqV2 is a dir of multiple files.
            # For now, assuming it's a single file or a path that pandas can read directly.
            # If it's a directory, this pd.read_csv will likely fail.
            # The original script used: rna_file = os.path.join(training_data_dir, "HiSeqV2")
            # And then pd.read_csv(rna_file). This implies HiSeqV2 was treated as a file path.
            # We'll stick to that assumption for now. If it's a dir, user needs to point to the actual data file.
             print(f"Warning: {rna_input_path} is a directory. Attempting to read it as a file path directly. This may fail if it's not a single data file.")

        edata = pd.read_csv(rna_input_path, sep='\t', index_col=0)
        print(f"Raw {dataset_prefix} RNA data loaded. Shape: {edata.shape}")

        # Drop the 'Entrez_Gene_Id' column if it exists, as it's not a sample
        # This column typically arises when the original second column (Entrez IDs)
        # gets this name after the first column (Hugo Symbols) is set as index.
        if 'Entrez_Gene_Id' in edata.columns:
            edata.drop(columns=['Entrez_Gene_Id'], inplace=True)
            print(f"Dropped 'Entrez_Gene_Id' column. New shape for {dataset_prefix} RNA: {edata.shape}")
        elif edata.shape[1] > 0 and edata.columns[0] == 'Entrez_Gene_Id': # Fallback
            print(f"Warning: First column of RNA data is 'Entrez_Gene_Id'. Attempting to drop it.")
            edata.drop(columns=[edata.columns[0]], inplace=True)
            print(f"Dropped first data column (assumed Entrez_Gene_Id). New shape for {dataset_prefix} RNA: {edata.shape}")

        rna_data_for_scoring = edata.copy() # Use a copy

        # exp_preprocess_for_rna internally does:
        # 1. Gene ID conversion (Symbol to Entrez)
        # 2. Filters genes not in mapping or signature
        # 3. Log2 transformation (if not already log2)
        # 4. Median centering (genes)
        # 5. Z-scoring (genes)
        # We want to skip log2, median centering, and z-scoring for RNA scores here.
        # So, we should NOT use exp_preprocess_for_rna if we want raw scores.
        # We will rely on calculate_signature_scores to handle necessary gene mapping.
        
        print(f"Calculating {dataset_prefix} RNA signature scores...")
        # If g_rna_skip_log2_transform is True, we tell calculate_signature_scores the data is already log2
        # so it doesn't attempt to do it. This effectively uses the raw (non-log2) values.
        rna_data_is_effectively_log2_for_helper = g_rna_skip_log2_transform

        rna_signature_score = calculate_signature_scores(
            data_matrix=rna_data_for_scoring, # Raw data
            gmt_file_path=rna_gmt_file,
            aggregation_method='median',
            mapping_df=mapping_df,
            data_index_uses_entrez=False, # Input RNA data index is Symbol
            gmt_already_entrez=True,    # RNA GMT contains Entrez IDs
            data_already_log2=rna_data_is_effectively_log2_for_helper 
        )
        print(f"{dataset_prefix} RNA signature scores calculated. Shape: {rna_signature_score.shape}")

        if not rna_signature_score.empty and rna_signature_score.shape[1] > 0:
            initial_signature_count = rna_signature_score.shape[0]
            if rna_signature_score.columns[0] in rna_signature_score.columns:
                na_in_first_sample = rna_signature_score[rna_signature_score.columns[0]].isna()
                rna_signature_score = rna_signature_score[~na_in_first_sample]
                print(f"Removed {initial_signature_count - rna_signature_score.shape[0]} {dataset_prefix} RNA signatures with NA in the first sample.")
            else:
                print(f"First sample column not found for NA check, skipping {dataset_prefix} RNA signature NA removal.")
        else:
            print(f"{dataset_prefix} RNA signature score DataFrame is empty or has no columns, skipping NA removal.")
        
        print(f"Final {dataset_prefix} RNA signature scores. Shape: {rna_signature_score.shape}")
        rna_output_filename = f"{dataset_prefix}_rna_signature_scores_recalculated.pkl"
        rna_signature_output_file = os.path.join(output_dir, rna_output_filename)
        rna_signature_score.to_pickle(rna_signature_output_file)
        print(f"{dataset_prefix} RNA signature scores saved to {rna_signature_output_file}")
    else:
        print(f"No RNA input path provided for {dataset_prefix}. Skipping RNA processing.")

    # --- CNA Segment Score Calculation ---
    print(f"\n--- {dataset_prefix} CNA Processing ---")
    if cna_input_path:
        print(f"Loading {dataset_prefix} CNA data from {cna_input_path}...")
        cna_data = pd.read_csv(cna_input_path, sep='\t', index_col=0)
        cna_data.columns = cna_data.columns.str.replace('-', '.') # Sample name consistency
        print(f"Raw {dataset_prefix} CNA data loaded. Shape: {cna_data.shape}")

        # Drop the 'Entrez_Gene_Id' column if it exists, as it's not a sample
        if 'Entrez_Gene_Id' in cna_data.columns:
            cna_data.drop(columns=['Entrez_Gene_Id'], inplace=True)
            print(f"Dropped 'Entrez_Gene_Id' column. New shape for {dataset_prefix} CNA: {cna_data.shape}")
        elif cna_data.shape[1] > 0 and cna_data.columns[0] == 'Entrez_Gene_Id': # Fallback
            print(f"Warning: First column of CNA data is 'Entrez_Gene_Id'. Attempting to drop it.")
            cna_data.drop(columns=[cna_data.columns[0]], inplace=True)
            print(f"Dropped first data column (assumed Entrez_Gene_Id). New shape for {dataset_prefix} CNA: {cna_data.shape}")

        cna_data_for_scoring = cna_data.copy()

        # Convert CNA data index to Entrez IDs (always do this)
        print(f"Converting {dataset_prefix} CNA data index from gene symbols to Entrez IDs...")
        original_cna_index_symbols = cna_data_for_scoring.index
        cna_data_for_scoring.index = convert_symbol_to_entrez_for_data_index(original_cna_index_symbols, mapping_df)
        cna_data_for_scoring = cna_data_for_scoring[cna_data_for_scoring.index.notna()]
        cna_data_for_scoring = cna_data_for_scoring[~cna_data_for_scoring.index.duplicated(keep='first')]
        print(f"{dataset_prefix} CNA data index converted to Entrez IDs and filtered. Shape: {cna_data_for_scoring.shape}")
        
        # Regarding log2 for CNA:
        # If g_cna_already_log2_assumed is True, we tell calculate_signature_scores it's "already log2"
        # to prevent it from applying any transformation.
        # GISTIC data is usually already log2 ratios. Other CNA data might vary.
        # For "no normalization and further processing", we skip any additional log2 here.
        cna_data_is_effectively_log2_for_helper = g_cna_already_log2_assumed

        print(f"Calculating {dataset_prefix} CNA segment scores...")
        cna_segment_score = calculate_signature_scores(
            data_matrix=cna_data_for_scoring, # CNA data with Entrez IDs
            gmt_file_path=cna_gmt_file,
            aggregation_method='mean',
            mapping_df=mapping_df,
            data_index_uses_entrez=True, # CNA data index is now Entrez
            gmt_already_entrez=True,   # CNA GMT also contains Entrez IDs
            data_already_log2=cna_data_is_effectively_log2_for_helper
        )
        print(f"{dataset_prefix} CNA segment scores calculated. Shape: {cna_segment_score.shape}")

        if not cna_segment_score.empty and cna_segment_score.shape[1] > 0:
            initial_segment_count = cna_segment_score.shape[0]
            if cna_segment_score.columns[0] in cna_segment_score.columns:
                na_in_first_sample_cna = cna_segment_score[cna_segment_score.columns[0]].isna()
                cna_segment_score = cna_segment_score[~na_in_first_sample_cna]
                print(f"Removed {initial_segment_count - cna_segment_score.shape[0]} {dataset_prefix} CNA segments with NA in the first sample.")
            else:
                print(f"First sample column not found for NA check, skipping {dataset_prefix} CNA segment NA removal.")
        else:
            print(f"{dataset_prefix} CNA segment score DataFrame is empty or has no columns, skipping NA removal.")

        print(f"Final {dataset_prefix} CNA segment scores. Shape: {cna_segment_score.shape}")
        cna_output_filename = f"{dataset_prefix}_cna_segment_scores_recalculated.pkl"
        cna_segment_output_file = os.path.join(output_dir, cna_output_filename)
        cna_segment_score.to_pickle(cna_segment_output_file)
        print(f"{dataset_prefix} CNA segment scores saved to {cna_segment_output_file}")
    else:
        print(f"No CNA input path provided for {dataset_prefix}. Skipping CNA processing.")


def main():
    parser = argparse.ArgumentParser(description='Calculate RNA and CNA signature scores with options for log2 transformation.')
    parser.add_argument('--train-cna', type=str, required=True, help='Path to training CNA data file.')
    parser.add_argument('--train-rna', type=str, required=True, help='Path to training RNA data file or directory.')
    parser.add_argument('--valid-cna', type=str, required=True, help='Path to validation CNA data file.')
    parser.add_argument('--valid-rna', type=str, required=True, help='Path to validation RNA data file.')
    parser.add_argument('--output-dir', type=str, required=True, help='Directory to save the output PKL files.')
    parser.add_argument('--mapping-file', type=str, default='gene_id_mapping.txt', help='Path to gene ID mapping file (Symbol to Entrez). Default: gene_id_mapping.txt')
    parser.add_argument('--rna-gmt', type=str, help='Path to RNA GMT file. Default uses path from original script.')
    parser.add_argument('--cna-gmt', type=str, help='Path to CNA GMT file. Default uses path from original script.')
    
    # Flags to control log2 transformation.
    # For RNA: if True, this script will NOT apply log2. calculate_signature_scores is told data_already_log2=True.
    parser.add_argument('--rna-skip-log2-transform', action='store_true', 
                        help='If set, RNA data will NOT be log2 transformed by this script. Assumes raw input is desired or already log2.')
    # For CNA: if True, this script will NOT apply log2. calculate_signature_scores is told data_already_log2=True.
    parser.add_argument('--cna-assume-already-log2', action='store_true',
                        help='If set, CNA data is assumed to be appropriately scaled (e.g., already log2 transformed or not needing it). This script will NOT apply log2 transformation.')

    args = parser.parse_args()

    print("Starting signature and segment score calculation with specified inputs...")
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output will be saved to: {os.path.abspath(args.output_dir)}")

    # Define GMT file paths (use defaults from original script if not provided)
    original_data_base_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data"
    rna_gmt_file = args.rna_gmt if args.rna_gmt else os.path.join(original_data_base_path, "gene_signatures_20170111.gmt")
    # Use the UTF-8 encoded CNA GMT file from the current working directory or specified path
    cna_gmt_file = args.cna_gmt if args.cna_gmt else "CNA_segments.utf8.gmt" # Assuming it's in the run directory or base_dir equivalent

    print(f"Using RNA GMT: {rna_gmt_file}")
    print(f"Using CNA GMT: {cna_gmt_file}")

    # Load gene ID mapping file
    print(f"Loading gene ID mapping from {args.mapping_file}...")
    raw_mapping_df = load_gene_id_mapping(args.mapping_file)
    mapping_df = pd.DataFrame({
        'Symbol': raw_mapping_df['Symbol_from_file'],
        'GeneID': raw_mapping_df['GeneID_from_file']
    })
    print(f"Mapping file processed. Shape: {mapping_df.shape}. First 5 symbols: {list(mapping_df['Symbol'][:5])}, First 5 GeneIDs: {list(mapping_df['GeneID'][:5])}")

    # Process Training Data
    process_dataset_pair(
        rna_input_path=args.train_rna,
        cna_input_path=args.train_cna,
        mapping_df=mapping_df,
        rna_gmt_file=rna_gmt_file,
        cna_gmt_file=cna_gmt_file,
        output_dir=args.output_dir,
        dataset_prefix="training",
        g_rna_skip_log2_transform=args.rna_skip_log2_transform,
        g_cna_already_log2_assumed=args.cna_assume_already_log2
    )

    # Process Validation Data
    process_dataset_pair(
        rna_input_path=args.valid_rna,
        cna_input_path=args.valid_cna,
        mapping_df=mapping_df,
        rna_gmt_file=rna_gmt_file,
        cna_gmt_file=cna_gmt_file,
        output_dir=args.output_dir,
        dataset_prefix="validation",
        g_rna_skip_log2_transform=args.rna_skip_log2_transform,
        g_cna_already_log2_assumed=args.cna_assume_already_log2
    )

    print("\nRecalculation of signature and segment scores finished.")

if __name__ == '__main__':
    main() 