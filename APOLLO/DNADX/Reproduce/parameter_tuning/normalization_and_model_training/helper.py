import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import gseapy as gp
from scipy import stats
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegressionCV

def calc_signatures(x, gmt_file, method="mean", scale_data=False):
    """
    Calculate gene signatures from expression data
    
    Parameters:
    -----------
    x : pandas.DataFrame
        Gene expression data matrix (genes x samples)
    gmt_file : str
        Path to GMT file containing gene sets
    method : str
        Method to calculate signature scores ('mean', 'median', 'pca', 'gsa')
    scale_data : bool
        Whether to scale the data before calculation
        
    Returns:
    --------
    pandas.DataFrame
        Signature scores matrix (signatures x samples)
    """
    # Read GMT file (now always utf-8)
    gmt = gp.parser.read_gmt(gmt_file)
 
    # Scale data if requested
    if scale_data:
        xs = pd.DataFrame(scale(x.T).T, index=x.index, columns=x.columns)
    else:
        xs = x.copy()  # added .copy() for safety
    
    # Initialize results matrix
    val = pd.DataFrame(index=gmt.keys(), columns=x.columns, dtype=float)  # dtype() for safety
    
    # Calculate scores for each gene set
    for set_name, genes in gmt.items():
        gene_set = [g for g in genes if g in xs.index]
        
        if len(gene_set) > 1:
            if method == "mean":
                val.loc[set_name] = xs.loc[gene_set].mean(axis=0)
            elif method == "median":
                val.loc[set_name] = xs.loc[gene_set].median(axis=0)
            elif method == "pca":
                pca = PCA(n_components=1)
                pca.fit(xs.loc[gene_set].T)  # fit PCA on samples-by-genes (transpose of xs.loc[gene_set])
                loadings = pca.components_[0]  # extract first PC loadings (length = #genes in gene_set)
                val.loc[set_name] = xs.loc[gene_set].T.dot(loadings)
            elif method == "gsa":
                # Implement GSA method if needed
                pass
        elif len(gene_set) == 1:
            val.loc[set_name] = xs.loc[gene_set[0]]
    
    return val

def overlap_sets(x, y):
    """
    Find overlapping genes between two datasets
    """
    common_genes = list(set(x.index) & set(y.index))
    x_subset = x.loc[common_genes].sort_index()
    y_subset = y.loc[common_genes].sort_index()
    return x_subset, y_subset

def assign_diff_score_dwd(x, y):
    """
    Calculate differentiation score using DWD method
    """
    x_subset, y_subset = overlap_sets(x, y)
    
    # Normalize vectors
    x_norm = x_subset.apply(lambda col: np.sign(col) * np.sqrt(col**2 / np.sum(col**2)))
    y_norm = y_subset.apply(lambda col: np.sign(col) * np.sqrt(col**2 / np.sum(col**2)))
    
    # Project samples
    ms_proj = np.dot(y_norm.T, x_norm.iloc[:, 0])
    ml_proj = np.dot(y_norm.T, x_norm.iloc[:, 1])
    
    # Calculate differentiation score
    diff_score = ml_proj - ms_proj
    return pd.Series(diff_score, index=y.columns)

def GHI_RS(edata):
    """
    Calculate Oncotype DX score
    """
    def get_gene(edata, gene):
        return edata.loc[gene] if gene in edata.index else pd.Series(0, index=edata.columns)
    
    # Gene IDs
    genes = [2597, 2990, 60, 7037, 6175, 2886, 2064, 596, 5241, 57758, 2099, 
             6790, 4605, 891, 332, 4288, 4320, 1515, 968, 2944, 573]
    
    # Get gene expression values
    gene_expr = {f'X{gene}': get_gene(edata, str(gene)) for gene in genes}
    
    # Calculate reference average
    reference_avg = sum(gene_expr[f'X{gene}'] for gene in [2597, 2990, 60, 7037, 6175]) / 5
    
    # Normalize genes
    def ref_norm(x, ref):
        x = x - ref
        x = x - x.min()
        x = x * 15 / x.max()
        return x
    
    for gene in genes:
        gene_expr[f'X{gene}'] = ref_norm(gene_expr[f'X{gene}'], reference_avg)
    
    # Calculate group scores
    GRB7_Group = 0.9 * gene_expr['X2886'] + 0.1 * gene_expr['X2064']
    GRB7_Group = GRB7_Group.clip(lower=8)
    
    ER_Group = (gene_expr['X596'] + 1.2 * gene_expr['X5241'] + 
                gene_expr['X57758'] + 0.8 * gene_expr['X2099']) / 4
    
    Prolif_Group = (gene_expr['X6790'] + gene_expr['X4605'] + gene_expr['X891'] + 
                    gene_expr['X332'] + gene_expr['X4288']) / 5
    Prolif_Group = Prolif_Group.clip(lower=6.5)
    
    Invasion_Group = (gene_expr['X4320'] + gene_expr['X1515']) / 2
    
    CD68 = gene_expr['X968']
    GSTM1 = gene_expr['X2944']
    BAG1 = gene_expr['X573']
    
    # Calculate final score
    RSU = (0.47 * GRB7_Group - 0.34 * ER_Group + 1.04 * Prolif_Group + 
           0.10 * Invasion_Group + 0.05 * CD68 - 0.08 * GSTM1 - 0.07 * BAG1)
    RS = 20 * (RSU - 6.7)
    
    return RS

def fisher_test(score, CN):
    """
    Perform Fisher's exact test between signature and gene
    """
    top_q = np.quantile(score, 0.75)
    module = pd.Series('low', index=score.index)
    module[score >= top_q] = 'high'
    
    # Create contingency table
    table = pd.crosstab(module, CN)
    
    # Perform Fisher's exact test
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    
    return p_value, odds_ratio

def sig_CN_test(score, CN_score, CN_gain, CN_loss):
    """
    Test associations between signature and all genes
    """
    results = []
    
    for gene in CN_score.index:
        CN = CN_score.loc[gene]
        
        # Spearman correlation
        spearman_cor, spearman_pos_p = stats.spearmanr(score, CN, alternative='greater')
        _, spearman_neg_p = stats.spearmanr(score, CN, alternative='less')
        
        # Linear model
        X = pd.DataFrame({
            'CN': CN,
            'basal': subtype_variable['basal'],
            'her2': subtype_variable['her2'],
            'lumA': subtype_variable['lumA'],
            'lumB': subtype_variable['lumB']
        })
        y = score
        
        model = sm.OLS(y, sm.add_constant(X)).fit()
        beta = model.params['CN']
        p_value = model.pvalues['CN']
        
        if beta > 0:
            lm_pos = p_value / 2
            lm_neg = 1 - p_value / 2
        else:
            lm_pos = 1 - p_value / 2
            lm_neg = p_value / 2
        
        # Fisher's exact test for CN gain/loss
        gain_p, gain_or = fisher_test(score, CN_gain.loc[gene])
        loss_p, loss_or = fisher_test(score, CN_loss.loc[gene])
        
        results.append({
            'spearman_pos': spearman_pos_p,
            'spearman_neg': spearman_neg_p,
            'CN_gain': gain_p,
            'CN_loss': loss_p,
            'lm_pos': lm_pos,
            'lm_neg': lm_neg
        })
    
    return pd.DataFrame(results, index=CN_score.index)

def calc_segments(x, gmt_file, method="mean", scale_data=False):
    """
    Calculate segment scores from CNA data
    """
    return calc_signatures(x, gmt_file, method, scale_data)

def median_ctr(x):
    """
    Center data by median
    """
    medians = x.median(axis=1)
    x_centered = x.sub(medians, axis=0)
    return x_centered

def standardize(x):
    """
    Standardize data
    """
    return pd.DataFrame(scale(x), index=x.index, columns=x.columns)

def exp_wrap(edata):
    """
    Process expression data
    """
    keep = ['29126', '1493', '5133']  # PD1, PDL1, CTLA4
    
    # Filter data
    edata70 = edata[edata.le(2).sum(axis=1) < (0.3 * edata.shape[1])]
    
    # Add required genes if missing
    for gene in keep:
        if gene not in edata70.index:
            edata70 = pd.concat([edata70, edata.loc[[gene]]])
    
    # Process data
    edata70[edata70 <= 2] = 0
    edata70log2 = np.log2(edata70)
    edata70log2[edata70log2 == -np.inf] = 0
    
    # Normalize
    exp = median_ctr(edata70log2)
    exp = standardize(exp)
    
    return exp

def caret_wrap(trainX, trainY, testX, testY, bi=True,
               alphas=None, l1_ratio=0.5, cv=5):
    # Pre-scale to match glmnet's standardize=TRUE (ddof=1)
    X_centered = trainX - trainX.mean(axis=0)
    X_scaled = X_centered / trainX.std(axis=0, ddof=1)
    testX_centered = testX - trainX.mean(axis=0)
    testX_scaled = testX_centered / trainX.std(axis=0, ddof=1)
 
    # Try the Elastic Net with CV
    if not bi:
        # regression
        model = ElasticNetCV(alphas=alphas,
                            l1_ratio=l1_ratio,
                            cv=cv,
                            normalize=False)  # normalized by R instead of by python
        model.fit(X_scaled, trainY)
    else:
        # classification-logistic
        model = LogisticRegressionCV(penalty="elasticnet",
                                    solver="saga",
                                    l1_ratios=[l1_ratio],
                                    Cs=None if alphas is None else 1/np.array(alphas),
                                    cv=cv,
                                    scoring="accuracy",
                                    max_iter=5000)
        model.fit(X_scaled, trainY)
 
    return model

def plot_ROC(perf1, perf2, a1, a2, main):
    """
    Plot ROC curves
    """
    plt.figure(figsize=(1.5, 1.5))
    plt.plot(perf1[0], perf1[1], 'r-', label=f'AUC = {a1:.2f}')
    plt.plot(perf2[0], perf2[1], 'b-', label=f'AUC = {a2:.2f}')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(main)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{main}.tiff', dpi=300, bbox_inches='tight')
    plt.close() 

def convert_symbol_to_entrez_for_gmt(gene_list, mapping_df):
    """Converts a list of gene symbols from GMT to Entrez IDs."""
    print(f"    DEBUG (convert_symbol_to_entrez_for_gmt): Input gene_list (first 5): {gene_list[:5]}")
    symbol_to_entrez = pd.Series(mapping_df.GeneID.astype(str).values, index=mapping_df.Symbol.astype(str)).to_dict()
    # Print a small part of the created map for verification
    temp_map_debug_count = 0
    print("      DEBUG (convert_symbol_to_entrez_for_gmt): Sample of symbol_to_entrez map:")
    for k, v in symbol_to_entrez.items():
        if temp_map_debug_count < 3:
            print(f"        '{k}': '{v}'")
            temp_map_debug_count += 1
        else:
            break

    entrez_ids = []
    mapped_count = 0
    unmapped_count = 0
    for gene_symbol_from_gmt in gene_list:
        s_gene_symbol = str(gene_symbol_from_gmt).strip()
        entrez_id = symbol_to_entrez.get(s_gene_symbol) # Direct lookup

        if entrez_id is None: # Attempt case-insensitive if direct fails
            for map_sym, map_eid in symbol_to_entrez.items():
                if s_gene_symbol.lower() == str(map_sym).lower():
                    entrez_id = map_eid
                    break
        
        if entrez_id is not None and str(entrez_id).lower() != 'nan':
            entrez_ids.append(str(entrez_id).split('.')[0])
            mapped_count +=1
        else:
            unmapped_count +=1 # Not adding None, just counting unmapped for debug
            
    print(f"    DEBUG (convert_symbol_to_entrez_for_gmt): Mapped {mapped_count} GMT genes, Unmapped: {unmapped_count}")
    print(f"    DEBUG (convert_symbol_to_entrez_for_gmt): Resulting Entrez IDs (first 5): {entrez_ids[:5]}")
    return entrez_ids # Return only successfully mapped and cleaned Entrez IDs

def convert_symbol_to_entrez_for_data_index(index_series, mapping_df):
    """Converts a pandas Index of gene symbols to Entrez IDs."""
    print("\nDebugging convert_symbol_to_entrez_for_data_index:")
    print(f"  First 5 input gene symbols from data index: {list(index_series[:5])}")
    
    # Create both forward and reverse mappings
    symbol_to_entrez = pd.Series(mapping_df.GeneID.astype(str).values, index=mapping_df.Symbol.astype(str)).to_dict()
    entrez_to_symbol = pd.Series(mapping_df.Symbol.astype(str).values, index=mapping_df.GeneID.astype(str)).to_dict()
    
    print(f"  First 5 entries in symbol_to_entrez map (from mapping_df):")
    for i, (k, v) in enumerate(symbol_to_entrez.items()):
        if i < 5:
            print(f"    '{k}': '{v}'")
        else:
            break
            
    new_index = []
    mapped_count = 0
    unmapped_count = 0
    for gene_symbol in index_series:
        # Ensure gene_symbol is a string and strip whitespace
        s_gene_symbol = str(gene_symbol).strip()
        entrez_id = None
        
        # Try direct lookup first
        if s_gene_symbol in symbol_to_entrez:
            entrez_id = symbol_to_entrez[s_gene_symbol]
        
        # Try case-insensitive lookup if direct fails
        if entrez_id is None:
            for map_sym, map_eid in symbol_to_entrez.items():
                if s_gene_symbol.lower() == str(map_sym).lower():
                    entrez_id = map_eid
                    break
        
        # Try reverse lookup (in case input is already Entrez ID)
        if entrez_id is None and s_gene_symbol in entrez_to_symbol:
            entrez_id = s_gene_symbol
        
        # Try numeric check (in case input is numeric Entrez ID as string)
        if entrez_id is None and s_gene_symbol.isdigit():
            entrez_id = s_gene_symbol

        if entrez_id is not None and str(entrez_id).lower() != 'nan':
            new_index.append(str(entrez_id).split('.')[0]) # Take part before any decimal, ensure string
            mapped_count += 1
        else:
            new_index.append(None) # Keep None for unmapped to filter later
            unmapped_count +=1
            
    print(f"  Mapped {mapped_count} gene symbols to Entrez IDs.")
    print(f"  Failed to map {unmapped_count} gene symbols.")
    print(f"  First 5 resulting Entrez IDs (or None): {new_index[:5]}")
    return pd.Index(new_index)

def load_gmt(gmt_file_path, mapping_df=None, convert_to_entrez=False, gmt_contains_entrez_ids=False):
    """Load gene sets from a GMT file.
    
    Parameters:
    -----------
    gmt_file_path : str
        Path to GMT file
    mapping_df : pd.DataFrame, optional
        DataFrame with 'Symbol' and 'GeneID' columns for mapping
    convert_to_entrez : bool
        Whether to convert gene symbols to Entrez IDs
    gmt_contains_entrez_ids : bool
        Whether GMT file already contains Entrez IDs
        
    Returns:
    --------
    dict
        Dictionary mapping signature names to lists of genes
    """
    print(f"\nLoading gene sets from {gmt_file_path}")
    print(f"  convert_to_entrez: {convert_to_entrez}")
    print(f"  gmt_contains_entrez_ids: {gmt_contains_entrez_ids}")
    
    gene_sets = {}
    unmapped_genes = set()
    total_genes = 0
    mapped_genes = 0
    
    with open(gmt_file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:  # Skip malformed lines
                continue
                
            signature_name = fields[0]
            description = fields[1]  # Usually same as name or description
            genes = fields[2:]
            total_genes += len(genes)
            
            # Clean and process genes
            processed_genes = []
            for gene in genes:
                gene = str(gene).strip()
                if not gene:  # Skip empty strings
                    continue
                    
                if gmt_contains_entrez_ids:
                    # If GMT contains Entrez IDs, just clean and validate
                    if gene.isdigit():
                        processed_genes.append(gene)
                    else:
                        print(f"Warning: Expected Entrez ID but got '{gene}' in signature {signature_name}")
                elif convert_to_entrez and mapping_df is not None:
                    # Convert symbol to Entrez ID
                    entrez_id = None
                    # Try direct lookup
                    symbol_matches = mapping_df[mapping_df.Symbol.astype(str).str.strip() == gene]
                    if not symbol_matches.empty:
                        entrez_id = str(symbol_matches.iloc[0].GeneID)
                    else:
                        # Try case-insensitive lookup
                        symbol_matches = mapping_df[mapping_df.Symbol.astype(str).str.strip().str.lower() == gene.lower()]
                        if not symbol_matches.empty:
                            entrez_id = str(symbol_matches.iloc[0].GeneID)
                    
                    if entrez_id and str(entrez_id).lower() != 'nan':
                        processed_genes.append(entrez_id)
                        mapped_genes += 1
                    else:
                        unmapped_genes.add(gene)
                else:
                    # Keep original symbols
                    processed_genes.append(gene)
                    mapped_genes += 1
            
            if processed_genes:  # Only add signatures with at least one valid gene
                gene_sets[signature_name] = processed_genes
    
    # Print statistics
    print(f"\nLoaded {len(gene_sets)} gene sets")
    print(f"Total genes processed: {total_genes}")
    print(f"Successfully mapped/validated genes: {mapped_genes}")
    if convert_to_entrez:
        print(f"Unmapped genes: {len(unmapped_genes)}")
        if len(unmapped_genes) > 0:
            print(f"Sample of unmapped genes: {list(unmapped_genes)[:5]}")
    
    return gene_sets

def median_center_rows(df):
    """Median-centers each row of the DataFrame."""
    medians = df.median(axis=1)
    return df.subtract(medians, axis=0)

def standardize_columns(df):
    """Z-score normalizes each column of the DataFrame, handling NaN std."""
    # Ensure mean and std are calculated correctly even with NaNs
    # R's scale uses N-1 for std, pandas default is N for .std(), use ddof=1 for N-1
    return df.apply(lambda x: (x - np.nanmean(x)) / (np.nanstd(x, ddof=1) if pd.notna(np.nanstd(x, ddof=1)) and np.nanstd(x, ddof=1) != 0 else 1), axis=0)

def exp_preprocess_for_rna(edata, mapping_df, genes_to_keep_symbols=None):
    """Preprocesses RNA expression data similar to the R script's exp_wrap function.

    Handles:
    - Gene ID conversion (Symbol to Entrez, if mapping_df provided).
    - Filtering for genes expressed in a minimum percentage of samples (similar to R script logic if implemented).
    - Log2 transformation (log2(x+1)).
    - Median centering of rows (genes).
    - Standardization of columns (samples) to mean 0, std 1.

    Parameters:
    -----------
    edata : pd.DataFrame
        Raw RNA expression data (genes as index, samples as columns).
        Assumed to have gene symbols as index if mapping_df is used for conversion.
    mapping_df : pd.DataFrame
        DataFrame with 'Symbol' and 'GeneID' columns for mapping.
    genes_to_keep_symbols : list, optional
        A list of gene symbols that should be preferentially kept if any filtering logic
        would otherwise remove them (similar to 'keep' in R's exp_wrap).

    Returns:
    --------
    pd.DataFrame
        Processed RNA expression data.
    """
    print("    Entering exp_preprocess_for_rna...")
    processed_data = edata.copy()

    # --- Gene ID Conversion (Symbol to Entrez) ---
    if mapping_df is not None:
        print("        Converting RNA data index from Symbol to Entrez ID...")
        original_index_symbols = processed_data.index
        processed_data.index = convert_symbol_to_entrez_for_data_index(original_index_symbols, mapping_df)
        
        # Keep only rows where Entrez ID conversion was successful (non-NaN)
        # and remove duplicates that might arise from mapping
        initial_rows = len(processed_data)
        processed_data = processed_data[processed_data.index.notna()]
        processed_data = processed_data[~processed_data.index.duplicated(keep='first')]
        print(f"        RNA data index converted. Rows before: {initial_rows}, Rows after: {len(processed_data)}")
    else:
        print("        No mapping_df provided, skipping Entrez ID conversion for RNA data.")

    # --- Filtering (Placeholder for R's edata70 logic) ---
    # R script filters genes: edata70 <- edata[rowSums(edata<=2)<(0.3*ncol(edata)),]
    # This means keep genes where the count of samples with expression <= 2 is less than 30% of total samples.
    # And then adds back certain 'keep' genes: PD1, PDL1, CTLA4 (by Entrez ID)
    # This step is complex and requires careful implementation if fully replicated.
    # For now, we assume this pre-filtering might have happened upstream or is omitted for simplification.
    # If genes_to_keep_symbols is provided, ensure they are present (assuming index is Symbol before conversion or Entrez after)
    # print("        Filtering step (placeholder)...")

    # --- Log2 Transformation ---
    # R script does: edata70log2 <- log(edata70,base = 2); edata70log2[edata70log2=="-Inf"] <- 0
    # We'll do log2(x+1) which is common to handle zeros and avoid -Inf directly.
    print("        Applying log2(x+1) transformation...")
    processed_data = np.log2(processed_data + 1)

    # --- Median Centering (Rows/Genes) ---
    # R script: exp <- medianCtr(edata70log2)
    print("        Median centering rows (genes)...")
    processed_data = median_center_rows(processed_data)

    # --- Standardization (Columns/Samples) ---
    # R script: exp <- standardize(exp)
    print("        Standardizing columns (samples) to mean 0, std 1...")
    processed_data = standardize_columns(processed_data)
    
    print(f"    Finished exp_preprocess_for_rna. Processed data shape: {processed_data.shape}")
    return processed_data

def calculate_signature_scores(data_matrix, gmt_file_path, aggregation_method='median', mapping_df=None, data_index_uses_entrez=False, gmt_already_entrez=False, data_already_log2=False):
    """Calculate signature scores from a data matrix using gene sets from a GMT file.
    
    Parameters:
    -----------
    data_matrix : pd.DataFrame
        Gene expression or CNA data matrix (genes x samples)
    gmt_file_path : str
        Path to GMT file containing gene sets
    aggregation_method : str
        Method to aggregate gene-level scores ('mean' or 'median')
    mapping_df : pd.DataFrame, optional
        DataFrame with 'Symbol' and 'GeneID' columns for mapping
    data_index_uses_entrez : bool
        Whether data_matrix index already uses Entrez IDs
    gmt_already_entrez : bool
        Whether GMT file already contains Entrez IDs
    data_already_log2 : bool
        Whether data is already log2 transformed
        
    Returns:
    --------
    pd.DataFrame
        Signature scores matrix (signatures x samples)
    """
    print(f"\nCalculating signature scores using {aggregation_method} aggregation...")
    print(f"Data matrix shape: {data_matrix.shape}")
    print(f"First few data matrix index entries: {list(data_matrix.index[:5])}")
    
    # Load gene sets from GMT file
    gene_sets = load_gmt(gmt_file_path, mapping_df, 
                        convert_to_entrez=not gmt_already_entrez,
                        gmt_contains_entrez_ids=gmt_already_entrez)
    print(f"Loaded {len(gene_sets)} gene sets from GMT file")
    
    # Convert data matrix index to Entrez IDs if needed
    working_matrix = data_matrix.copy()
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