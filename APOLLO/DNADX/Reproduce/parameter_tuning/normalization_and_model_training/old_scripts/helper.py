import pandas as pd
import numpy as np

def load_gene_id_mapping(mapping_file):
    """Load gene ID mapping file"""
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    return mapping_df

def convert_symbol_to_entrez(index_symbols, mapping_file):
    """Convert gene symbols to Entrez IDs using mapping file directly"""
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

def convert_symbol_to_entrez_for_data_index(index_symbols, mapping_df):
    """Convert gene symbols to Entrez IDs using pre-loaded mapping DataFrame"""
    entrez_ids = []
    for symbol in index_symbols:
        match = mapping_df[mapping_df['Symbol'] == symbol]
        if not match.empty:
            entrez_ids.append(str(match.iloc[0]['GeneID']))
        else:
            entrez_ids.append(np.nan)
    return entrez_ids

def calc_signatures(data_matrix, gmt_file, method="median"):
    """Calculate signature scores from expression data using gene sets."""
    signatures = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:  # Skip malformed lines
                continue
            signature_name = parts[0]
            genes = parts[2:]  # Skip name and description
            signatures[signature_name] = genes
    
    scores = {}
    for name, genes in signatures.items():
        # Get intersection of genes in signature and data
        common_genes = list(set(genes) & set(data_matrix.index))
        if common_genes:
            if method == "median":
                scores[name] = data_matrix.loc[common_genes].median()
            elif method == "mean":
                scores[name] = data_matrix.loc[common_genes].mean()
    
    return pd.DataFrame.from_dict(scores, orient='index')

def calc_segments(data_matrix, gmt_file, method="mean"):
    """Calculate segment scores from CNA data using gene sets."""
    segments = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:  # Skip malformed lines
                continue
            segment_name = parts[0]
            genes = parts[2:]  # Skip name and description
            segments[segment_name] = genes
    
    scores = {}
    for name, genes in segments.items():
        # Get intersection of genes in segment and data
        common_genes = list(set(genes) & set(data_matrix.index))
        if common_genes:
            if method == "median":
                scores[name] = data_matrix.loc[common_genes].median()
            elif method == "mean":
                scores[name] = data_matrix.loc[common_genes].mean()
    
    return pd.DataFrame.from_dict(scores, orient='index')

def clean_feature_names(feature_names):
    """Clean feature names by removing special characters and standardizing format.
    This function is used to ensure consistent feature naming between training and validation."""
    cleaned_names = []
    for name in feature_names:
        # Convert to string if not already
        name = str(name)
        # Remove any special characters that might cause issues
        name = name.replace('/', '_').replace('\\', '_').replace(' ', '_')
        cleaned_names.append(name)
    return cleaned_names

def calculate_signature_scores(data_matrix, gmt_file_path, aggregation_method='median', 
                             mapping_df=None, data_index_uses_entrez=False, 
                             gmt_already_entrez=False, data_already_log2=False):
    """Calculate signature scores from a data matrix using gene sets from a GMT file.
    
    Args:
        data_matrix: DataFrame with genes as index and samples as columns
        gmt_file_path: Path to GMT file containing gene sets
        aggregation_method: Method to aggregate gene scores ('mean' or 'median')
        mapping_df: DataFrame with gene symbol to Entrez ID mapping
        data_index_uses_entrez: Whether data_matrix index already uses Entrez IDs
        gmt_already_entrez: Whether GMT file already uses Entrez IDs
        data_already_log2: Whether data is already log2 transformed
    
    Returns:
        DataFrame with signature scores (signatures as index, samples as columns)
    """
    print(f"Loading gene sets from {gmt_file_path}...")
    gene_sets = {}
    with open(gmt_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:  # Skip malformed lines
                continue
            signature_name = parts[0]
            genes = parts[2:]  # Skip name and description
            gene_sets[signature_name] = genes
    
    print(f"Loaded {len(gene_sets)} gene sets")
    
    # Make a copy of the data matrix to avoid modifying the original
    working_matrix = data_matrix.copy()
    
    # Convert gene symbols to Entrez IDs if needed
    if not data_index_uses_entrez and mapping_df is not None:
        print("Converting data matrix index from symbols to Entrez IDs...")
        original_shape = working_matrix.shape
        working_matrix.index = convert_symbol_to_entrez_for_data_index(working_matrix.index, mapping_df)
        working_matrix = working_matrix[working_matrix.index.notna()]
        working_matrix = working_matrix[~working_matrix.index.duplicated(keep='first')]
        print(f"After conversion and filtering: {working_matrix.shape} (original: {original_shape})")
    
    # Log2 transform if needed
    if not data_already_log2:
        print("Applying log2(x+1) transformation...")
        working_matrix = np.log2(working_matrix + 1)
    
    # Initialize results matrix
    signature_scores = pd.DataFrame(index=gene_sets.keys(), columns=working_matrix.columns)
    
    # Calculate signature scores
    print("\nCalculating scores for each signature...")
    for signature_name, genes_in_signature in gene_sets.items():
        # Find overlapping genes
        genes_in_signature_and_data = list(set(genes_in_signature) & set(working_matrix.index))
        
        if not genes_in_signature_and_data:
            print(f"Warning: No overlapping genes found for signature '{signature_name}'")
            print(f"  Signature genes (first 5): {genes_in_signature[:5]}")
            print(f"  Data matrix genes (first 5): {list(working_matrix.index[:5])}")
            score_row = pd.Series([np.nan] * working_matrix.shape[1], 
                                index=working_matrix.columns, 
                                name=signature_name)
        else:
            # Calculate signature score
            signature_matrix = working_matrix.loc[genes_in_signature_and_data]
            
            if aggregation_method == 'mean':
                score_row = signature_matrix.mean(axis=0)
            elif aggregation_method == 'median':
                score_row = signature_matrix.median(axis=0)
            else:
                raise ValueError(f"Unknown aggregation method: {aggregation_method}")
            
            score_row.name = signature_name
            
            # Print overlap stats periodically
            if len(signature_scores) % 50 == 0:
                print(f"  {signature_name}: Found {len(genes_in_signature_and_data)} overlapping genes out of {len(genes_in_signature)} signature genes")
        
        signature_scores.loc[signature_name] = score_row
    
    # Print final statistics
    total_na = signature_scores.isna().sum().sum()
    total_elements = signature_scores.size
    print(f"\nFinal signature scores shape: {signature_scores.shape}")
    print(f"Total NaN values: {total_na} ({(total_na/total_elements)*100:.2f}%)")
    
    return signature_scores 