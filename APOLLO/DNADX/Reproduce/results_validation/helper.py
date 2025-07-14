import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
import gseapy as gp
from scipy import stats
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

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
    gene_sets = {k: v for k, v in gmt.items()}
    
    # Scale data if requested
    if scale_data:
        xs = pd.DataFrame(scale(x.T).T, index=x.index, columns=x.columns)
    else:
        xs = x
    
    # Initialize results matrix
    val = pd.DataFrame(index=gene_sets.keys(), columns=x.columns)
    
    # Calculate scores for each gene set
    for i, (set_name, genes) in enumerate(gene_sets.items()):
        # Find genes in the set
        gene_set = [g for g in genes if g in xs.index]
        
        if len(gene_set) > 1:
            if method == "mean":
                val.loc[set_name] = xs.loc[gene_set].mean()
            elif method == "median":
                val.loc[set_name] = xs.loc[gene_set].median()
            elif method == "pca":
                pca = PCA(n_components=1)
                val.loc[set_name] = pca.fit_transform(xs.loc[gene_set].T)[:, 0]
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

def caret_wrap(trainX, trainY, testX, testY, bi=True):
    """
    Wrapper for Elastic Net modeling
    """
    if not bi:
        # Regression
        model = ElasticNet(alpha=0.5, l1_ratio=0.5)
        model.fit(trainX, trainY)
        return model
    else:
        # Classification
        model = ElasticNet(alpha=0.5, l1_ratio=0.5)
        model.fit(trainX, trainY)
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