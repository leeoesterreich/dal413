import pandas as pd
import numpy as np
from helper import GHI_RS

def convert_symbol_to_entrez(symbols, id_mapping_file):
    """Convert gene symbols to Entrez IDs using the mapping file"""
    mapping = pd.read_csv(id_mapping_file, sep='\t', header=None, names=['EntrezID', 'Symbol'])
    # Create symbol to EntrezID mapping
    symbol_to_entrez = dict(zip(mapping['Symbol'], mapping['EntrezID'].astype(str)))
    # Convert symbols to EntrezIDs
    return [symbol_to_entrez.get(symbol, symbol) for symbol in symbols]

def main():
    # Load expression data
    print("Loading expression data...")
    edata = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_mrna_agilent_microarray.txt", 
                        sep='\t', index_col=0)
    
    # Convert gene symbols to Entrez IDs
    print("Converting gene symbols to Entrez IDs...")
    edata.index = convert_symbol_to_entrez(edata.index, "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/gene_id_mapping.txt")
    
    # Calculate Oncotype DX score
    print("Calculating Oncotype DX score...")
    oncotype_score = GHI_RS(edata)
    
    # Convert to DataFrame and save
    print("Saving results...")
    result_df = pd.DataFrame({
        'GHI_RS_Model_NJEM.2004_PMID.15591335': oncotype_score
    })
    result_df.to_csv('oncotype_dx_score.csv')
    print("Done! Results saved to oncotype_dx_score.csv")

if __name__ == "__main__":
    main() 