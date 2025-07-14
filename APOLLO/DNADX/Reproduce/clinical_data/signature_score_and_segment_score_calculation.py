import pandas as pd
import numpy as np
from helper import (calc_signatures, assign_diff_score_dwd, GHI_RS,
                   calc_segments)
import os

def convert_symbol_to_entrez(symbols, id_mapping_file):
    """Convert gene symbols to Entrez IDs using the mapping file"""
    mapping = pd.read_csv(id_mapping_file, sep='\t', header=None, names=['EntrezID', 'Symbol'])
    # Create symbol to EntrezID mapping
    symbol_to_entrez = dict(zip(mapping['Symbol'], mapping['EntrezID'].astype(str)))
    # Convert symbols to EntrezIDs
    return [symbol_to_entrez.get(symbol, symbol) for symbol in symbols]

def convert_entrez_to_symbol(entrez_ids, id_mapping_file):
    """Convert Entrez IDs to gene symbols using the mapping file"""
    mapping = pd.read_csv(id_mapping_file, sep='\t', header=None, names=['EntrezID', 'Symbol'])
    # Create EntrezID to symbol mapping
    entrez_to_symbol = dict(zip(mapping['EntrezID'].astype(str), mapping['Symbol']))
    # Convert EntrezIDs to symbols
    return [entrez_to_symbol.get(entrez, entrez) for entrez in entrez_ids]

def main():
    # Define base paths
    base_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce"
    data_path = os.path.join(base_path, "DNA-based-predictors-of-non-genetic-cancer-phenotypes/data")
    
    # Load input data
    print("Loading input data...")
    edata = pd.read_csv(os.path.join(base_path, "HiSeqV2"), sep='\t', index_col=0)
    
    # Convert gene symbols to Entrez IDs for expression data
    print("Converting gene symbols to Entrez IDs for expression data...")
    edata.index = convert_symbol_to_entrez(edata.index, os.path.join(base_path, "gene_id_mapping.txt"))
    
    # Load and process CNA data
    print("\nLoading CNA data...")
    cna_data = pd.read_csv(os.path.join(base_path, "Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt"), sep='\t')
    
    # Convert gene symbols to Entrez IDs for CNA data
    print("Converting gene symbols to Entrez IDs for CNA data...")
    cna_data['Gene Symbol'] = convert_symbol_to_entrez(cna_data['Gene Symbol'], 
                                                      os.path.join(base_path, "gene_id_mapping.txt"))
    
    # Set gene symbols as index
    cna_data.set_index('Gene Symbol', inplace=True)
    
    # Print first 10 row labels
    print("\nFirst 10 row labels in edata:")
    print(edata.index[:10].tolist())
    print("\nFirst 10 row labels in CNA data:")
    print(cna_data.index[:10].tolist())
    
    # Calculate signature scores
    print("\nCalculating signature scores...")
    signature_score = calc_signatures(edata, 
                                    os.path.join(data_path, "gene_signatures_20170111.gmt"),
                                    method="median")
    
    print("\nFirst 10 row labels in signature_score:")
    print(signature_score.index[:10].tolist())
    
    # Calculate CD103_Ratio
    print("Calculating CD103_Ratio...")
    CD103_pos = signature_score.loc["CD103_Positive_Median_Cancer.Cell.2014_PMID.25446897"]
    CD103_neg = signature_score.loc["CD103_Negative_Median_Cancer.Cell.2014_PMID.25446897"]
    CD103_ratio = CD103_pos - CD103_neg
    signature_score.loc["CD103_Ratio_Cancer.Cell.2014_PMID.25446897"] = CD103_ratio
    
    # Calculate differentiation score
    print("Calculating differentiation score...")
    diff_centroid = pd.read_csv(os.path.join(data_path, "special_gene_signature_training_sets/UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035.txt"),
                               sep='\t', index_col=0)
    diff_score = assign_diff_score_dwd(diff_centroid, edata)
    signature_score.loc["UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035"] = diff_score
    
    # Calculate Oncotype DX score
    print("Calculating Oncotype DX score...")
    oncotype = GHI_RS(edata)
    signature_score.loc["GHI_RS_Model_NJEM.2004_PMID.15591335"] = oncotype
    
    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)
    
    # Save signature scores
    print("Saving signature scores...")
    signature_score.to_pickle("results/signature_score.pkl")
    
    # Calculate segment scores
    print("\nCalculating segment scores...")
    segment_score = calc_segments(cna_data, os.path.join(data_path, 'CNA_segments.utf8.gmt'), method='mean')
    
    # Save segment scores
    print("Saving segment scores...")
    segment_score.to_pickle("results/segment_score.pkl")
    
    print("\nAll calculations completed successfully!")

if __name__ == "__main__":
    main() 