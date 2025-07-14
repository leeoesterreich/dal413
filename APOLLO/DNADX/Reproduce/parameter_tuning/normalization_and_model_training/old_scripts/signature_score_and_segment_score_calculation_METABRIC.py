import pandas as pd
import numpy as np
import pickle
import os
from helper import calc_signatures, calc_segments
from sklearn.preprocessing import StandardScaler

def clean_column_names(df):
    """Clean column names by removing any special characters and converting to uppercase."""
    # Don't modify the column names as they are sample IDs
    return df

def convert_symbol_to_entrez(index_symbols, mapping_file):
    """Convert gene symbols to Entrez IDs using mapping file"""
    # Load mapping file
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Convert symbols to Entrez IDs
    entrez_ids = []
    for symbol in index_symbols:
        match = mapping_df[mapping_df['Symbol_from_file'] == symbol]
        if not match.empty:
            entrez_ids.append(str(match.iloc[0]['GeneID_from_file']))
        else:
            entrez_ids.append(np.nan)
    return entrez_ids

def zscore_normalize(df):
    """Normalize data using z-score standardization"""
    scaler = StandardScaler()
    # Transpose to normalize samples (columns) instead of genes (rows)
    normalized_data = scaler.fit_transform(df.T)
    # Transpose back to original orientation
    return pd.DataFrame(normalized_data.T, index=df.index, columns=df.columns)

def main():
    # Define input and output paths
    validation_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/validation_data"
    rna_file = os.path.join(validation_dir, "data_mrna_agilent_microarray.txt")
    cna_file = os.path.join(validation_dir, "data_cna_METABRIC.txt")
    mapping_file = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/gene_id_mapping.txt"
    output_dir = os.path.join(validation_dir, "metabric_scores")
    os.makedirs(output_dir, exist_ok=True)

    # Load RNA expression data
    print("Loading RNA data...")
    edata = pd.read_csv(rna_file, sep='\t', index_col=0)
    print(f"RNA data shape: {edata.shape}")
    print("First few RNA samples:", list(edata.columns)[:5])
    edata = clean_column_names(edata)

    # Convert gene symbols to Entrez IDs for expression data
    print("\nConverting gene symbols to Entrez IDs for expression data...")
    edata.index = convert_symbol_to_entrez(edata.index, mapping_file)

    # Load CNA data
    print("\nLoading CNA data...")
    cna_data = pd.read_csv(cna_file, sep='\t', index_col=0)
    print(f"CNA data shape: {cna_data.shape}")
    print("First few CNA samples:", list(cna_data.columns)[:5])
    cna_data = clean_column_names(cna_data)

    # Convert gene symbols to Entrez IDs for CNA data
    print("\nConverting gene symbols to Entrez IDs for CNA data...")
    cna_data.index = convert_symbol_to_entrez(cna_data.index, mapping_file)

    # Find common samples
    common_samples = list(set(edata.columns) & set(cna_data.columns))
    print(f"\nFound {len(common_samples)} common samples")
    if len(common_samples) > 0:
        print("First few common samples:", common_samples[:5])

    # Filter data to common samples
    edata = edata[common_samples]
    cna_data = cna_data[common_samples]

    # Print first 10 row labels
    print("\nFirst 10 row labels in edata:")
    print(list(edata.index)[:10])
    print("\nFirst 10 row labels in CNA data:")
    print(list(cna_data.index)[:10])

    # Normalize both datasets using z-score
    print("\nNormalizing RNA data...")
    edata_normalized = zscore_normalize(edata)
    print("\nNormalizing CNA data...")
    cna_data_normalized = zscore_normalize(cna_data)

    # Calculate signature scores for RNA data
    print("\nCalculating signature scores for RNA data...")
    rna_signature_score = calc_signatures(edata_normalized, "gene_signatures_20170111.gmt", method="median")

    # Calculate segment scores for CNA data
    print("\nCalculating segment scores for CNA data...")
    cna_signature_score = calc_segments(cna_data_normalized, "CNA_segments.utf8.gmt", method="median")

    print("\nFirst 10 row labels in RNA signature_score:")
    print(list(rna_signature_score.index)[:10])
    print("\nFirst 10 row labels in CNA signature_score:")
    print(list(cna_signature_score.index)[:10])

    # Save signature scores
    rna_signature_score.to_pickle(os.path.join(output_dir, "rna_signature_score.pkl"))
    cna_signature_score.to_pickle(os.path.join(output_dir, "cna_signature_score.pkl"))

    # Also save as CSV for easier inspection
    rna_signature_score.to_csv(os.path.join(output_dir, "rna_signature_score.csv"))
    cna_signature_score.to_csv(os.path.join(output_dir, "cna_signature_score.csv"))

    print("\nDone! Files saved in:", output_dir)
    print("RNA signature scores shape:", rna_signature_score.shape)
    print("CNA signature scores shape:", cna_signature_score.shape)

if __name__ == "__main__":
    main() 