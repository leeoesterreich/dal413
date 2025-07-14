'''
Calculates RNA and CNA signature scores for METABRIC data using specified GMT files.
Follows the same logic as the original score calculation script.
'''
import pandas as pd
import numpy as np
import os
from helper import (
    exp_preprocess_for_rna,
    calculate_signature_scores,
    load_gmt,
    convert_symbol_to_entrez_for_data_index
)

def load_gene_id_mapping(mapping_file_path):
    """Loads the gene ID mapping file, assuming no header and specific column order: EntrezID then Symbol."""
    return pd.read_csv(mapping_file_path, sep='\t', header=None, names=['GeneID_from_file', 'Symbol_from_file'])

def main():
    print("Starting METABRIC signature and segment score calculation...")
    
    # Define file paths
    base_dir = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training'
    validation_dir = os.path.join(base_dir, 'validation_data')
    output_dir = os.path.join(validation_dir, 'metabric_scores')
    os.makedirs(output_dir, exist_ok=True)

    # Input files
    rna_file = os.path.join(validation_dir, "data_mrna_agilent_microarray.txt")
    cna_file = os.path.join(validation_dir, "data_cna_METABRIC.txt")
    mapping_file = os.path.join(base_dir, "gene_id_mapping.txt")
    
    # Original GMT files from the paper's directory
    original_data_base_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data"
    rna_gmt_file = os.path.join(original_data_base_path, "gene_signatures_20170111.gmt")
    cna_gmt_file = os.path.join(base_dir, "CNA_segments.utf8.gmt")

    print("Using RNA GMT: {}".format(rna_gmt_file))
    print("Using CNA GMT: {}".format(cna_gmt_file))

    # Load gene ID mapping file
    print("Loading gene ID mapping from {}...".format(mapping_file))
    raw_mapping_df = load_gene_id_mapping(mapping_file)
    mapping_df = pd.DataFrame({
        'Symbol': raw_mapping_df['Symbol_from_file'],
        'GeneID': raw_mapping_df['GeneID_from_file']
    })
    print("Mapping file processed. Shape: {}".format(mapping_df.shape))

    # --- RNA Signature Score Calculation ---
    print("\n--- Processing METABRIC RNA Data ---")
    print("Loading RNA expression data from {}...".format(rna_file))
    edata = pd.read_csv(rna_file, sep='\t', index_col=0)
    print("Raw RNA data loaded. Shape: {}".format(edata.shape))

    # Preprocess RNA data
    processed_edata = exp_preprocess_for_rna(edata, mapping_df)
    print("RNA data preprocessed. Shape: {}".format(processed_edata.shape))

    # Calculate RNA signature scores
    rna_signature_score = calculate_signature_scores(
        data_matrix=edata,
        gmt_file_path=rna_gmt_file,
        aggregation_method='median',
        mapping_df=mapping_df,
        data_index_uses_entrez=False,
        gmt_already_entrez=True,
        data_already_log2=True
    )
    print("Raw RNA signature scores calculated. Shape: {}".format(rna_signature_score.shape))

    # Remove rows with NA in first sample
    if not rna_signature_score.empty and rna_signature_score.shape[1] > 0:
        initial_signature_count = rna_signature_score.shape[0]
        if rna_signature_score.columns[0] in rna_signature_score.columns:
            na_in_first_sample = rna_signature_score[rna_signature_score.columns[0]].isna()
            rna_signature_score = rna_signature_score[~na_in_first_sample]
            print("Removed {} RNA signatures with NA in the first sample.".format(
                initial_signature_count - rna_signature_score.shape[0]))

    print("Final RNA signature scores. Shape: {}".format(rna_signature_score.shape))
    rna_signature_output_file = os.path.join(output_dir, "rna_signature_score.pkl")
    rna_signature_score.to_pickle(rna_signature_output_file)
    print("RNA signature scores saved to {}".format(rna_signature_output_file))

    # --- CNA Segment Score Calculation ---
    print("\n--- Processing METABRIC CNA Data ---")
    print("Loading CNA data from {}...".format(cna_file))
    cna_data = pd.read_csv(cna_file, sep='\t', index_col=0)
    cna_data.columns = cna_data.columns.str.replace('-', '.')
    print("Raw CNA data loaded. Shape: {}".format(cna_data.shape))

    # Convert CNA data index to Entrez IDs
    print("Converting CNA data index from gene symbols to Entrez IDs...")
    original_cna_index_symbols = cna_data.index
    cna_data.index = convert_symbol_to_entrez_for_data_index(original_cna_index_symbols, mapping_df)

    # Filter CNA data
    cna_data = cna_data[cna_data.index.notna()]
    cna_data = cna_data[~cna_data.index.duplicated(keep='first')]
    print("CNA data index converted to Entrez IDs and filtered. Shape: {}".format(cna_data.shape))

    # Calculate CNA segment scores
    cna_segment_score = calculate_signature_scores(
        data_matrix=cna_data,
        gmt_file_path=cna_gmt_file,
        aggregation_method='mean',
        mapping_df=mapping_df,
        data_index_uses_entrez=True,
        gmt_already_entrez=True,
        data_already_log2=True
    )
    print("CNA segment scores calculated. Shape: {}".format(cna_segment_score.shape))

    # Remove rows with NA in first sample
    if not cna_segment_score.empty and cna_segment_score.shape[1] > 0:
        initial_segment_count = cna_segment_score.shape[0]
        if cna_segment_score.columns[0] in cna_segment_score.columns:
            na_in_first_sample_cna = cna_segment_score[cna_segment_score.columns[0]].isna()
            cna_segment_score = cna_segment_score[~na_in_first_sample_cna]
            print("Removed {} CNA segments with NA in the first sample.".format(
                initial_segment_count - cna_segment_score.shape[0]))

    print("Final CNA segment scores. Shape: {}".format(cna_segment_score.shape))
    cna_segment_output_file = os.path.join(output_dir, "cna_signature_score.pkl")
    cna_segment_score.to_pickle(cna_segment_output_file)
    print("CNA segment scores saved to {}".format(cna_segment_output_file))

    print("\nMETABRIC signature and segment score calculation finished.")

if __name__ == '__main__':
    main() 