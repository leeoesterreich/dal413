import pandas as pd
import numpy as np
import os
from helper import calc_signatures, calc_segments

def load_gene_id_mapping(mapping_file):
    """Load gene ID mapping file"""
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    return mapping_df

def convert_symbol_to_entrez_for_data_index(index_symbols, mapping_df):
    """Convert gene symbols to Entrez IDs for data index"""
    entrez_ids = []
    for symbol in index_symbols:
        match = mapping_df[mapping_df['Symbol'] == symbol]
        if not match.empty:
            entrez_ids.append(str(match.iloc[0]['GeneID']))
        else:
            entrez_ids.append(np.nan)
    return entrez_ids

def main():
    print("Starting signature and segment score calculation...")
    # Define file paths
    base_dir = '.' # Assuming script is run from the project root
    training_data_dir = os.path.join(base_dir, 'training_data')
    os.makedirs(training_data_dir, exist_ok=True) # Ensure training_data_dir exists

    # Original paper's data paths
    original_data_base_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data"

    rna_file = os.path.join(training_data_dir, "HiSeqV2") # Keep using local training HiSeqV2
    cna_file = os.path.join(training_data_dir, "Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt") # Keep using local training Gistic
    mapping_file = os.path.join(base_dir, "gene_id_mapping.txt") # Local mapping file
    
    # Use ORIGINAL GMT files from the paper
    rna_gmt_file = os.path.join(original_data_base_path, "gene_signatures_20170111.gmt")
    # Use the UTF-8 encoded CNA GMT file from the current working directory
    cna_gmt_file = os.path.join(base_dir, "CNA_segments.utf8.gmt")

    print("Using RNA GMT: {}".format(rna_gmt_file))
    print("Using CNA GMT: {}".format(cna_gmt_file))

    # Load gene ID mapping file
    print("Loading gene ID mapping from {}...".format(mapping_file))
    raw_mapping_df = load_gene_id_mapping(mapping_file)
    # Create the mapping_df in the expected format for helper functions: 'Symbol' and 'GeneID' columns
    mapping_df = pd.DataFrame({
        'Symbol': raw_mapping_df['Symbol_from_file'],
        'GeneID': raw_mapping_df['GeneID_from_file']
    })
    print("Mapping file processed. Shape: {}. First 5 symbols: {}, First 5 GeneIDs: {}".format(
        mapping_df.shape, list(mapping_df['Symbol'][:5]), list(mapping_df['GeneID'][:5])))

    # --- RNA Signature Score Calculation ---
    print("\n--- Processing RNA Data ---")
    print("Loading RNA expression data from {}...".format(rna_file))
    edata = pd.read_csv(rna_file, sep='\t', index_col=0)
    print("Raw RNA data loaded. Shape: {}".format(edata.shape))

    # Convert RNA data index to Entrez IDs
    print("Converting RNA data index from gene symbols to Entrez IDs...")
    original_rna_index_symbols = edata.index
    edata.index = convert_symbol_to_entrez_for_data_index(original_rna_index_symbols, mapping_df)

    # Filter RNA data: keep rows with valid (non-NaN) Entrez IDs and remove duplicates
    edata = edata[edata.index.notna()]
    edata = edata[~edata.index.duplicated(keep='first')]
    print("RNA data index converted to Entrez IDs and filtered. Shape: {}".format(edata.shape))

    # Calculate RNA signature scores
    print("Calculating RNA signature scores...")
    rna_signature_score = calc_signatures(edata, rna_gmt_file, method="median")
    print("RNA signature scores calculated. Shape: {}".format(rna_signature_score.shape))
    
    print("Final RNA signature scores. Shape: {}".format(rna_signature_score.shape))
    rna_signature_output_file = os.path.join(training_data_dir, "rna_signature_score_median_no_norm.pkl")
    rna_signature_score.to_pickle(rna_signature_output_file)
    print("RNA signature scores (median aggregation, no normalization) saved to {}".format(rna_signature_output_file))

    # --- CNA Segment Score Calculation ---
    print("\n--- Processing CNA Data ---")
    print("Loading CNA data from {}...".format(cna_file))
    cna_data = pd.read_csv(cna_file, sep='\t', index_col=0)
    cna_data.columns = cna_data.columns.str.replace('-', '.')
    print("Raw CNA data loaded. Shape: {}".format(cna_data.shape))

    # Convert CNA data index to Entrez IDs
    print("Converting CNA data index from gene symbols to Entrez IDs...")
    original_cna_index_symbols = cna_data.index
    cna_data.index = convert_symbol_to_entrez_for_data_index(original_cna_index_symbols, mapping_df)

    # Filter CNA data: keep rows with valid (non-NaN) Entrez IDs and remove duplicates
    cna_data = cna_data[cna_data.index.notna()]
    cna_data = cna_data[~cna_data.index.duplicated(keep='first')]
    print("CNA data index converted to Entrez IDs and filtered. Shape: {}".format(cna_data.shape))

    # Calculate CNA segment scores
    print("Calculating CNA segment scores...")
    cna_segment_score = calc_segments(cna_data, cna_gmt_file, method="mean")
    print("CNA segment scores calculated. Shape: {}".format(cna_segment_score.shape))

    # Remove rows with NA in first sample for consistency
    if not cna_segment_score.empty and cna_segment_score.shape[1] > 0:
        initial_segment_count = cna_segment_score.shape[0]
        if cna_segment_score.columns[0] in cna_segment_score.columns:
            na_in_first_sample_cna = cna_segment_score[cna_segment_score.columns[0]].isna()
            cna_segment_score = cna_segment_score[~na_in_first_sample_cna]
            print("Removed {} CNA segments with NA in the first sample.".format(
                initial_segment_count - cna_segment_score.shape[0]))
        else:
            print("First sample column not found for NA check, skipping CNA segment NA removal.")
    else:
        print("CNA segment score DataFrame is empty or has no columns, skipping NA removal.")

    print("Final CNA segment scores (mean aggregation, no Z-score). Shape: {}".format(cna_segment_score.shape))
    cna_segment_output_file = os.path.join(training_data_dir, "cna_segment_score_mean_no_norm.pkl")
    cna_segment_score.to_pickle(cna_segment_output_file)
    print("CNA segment scores (mean aggregation, no normalization) saved to {}".format(cna_segment_output_file))

if __name__ == "__main__":
    main() 