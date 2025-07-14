import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

def calc_signatures(x, gmt_file, method="mean", scale=False):
    """
    Calculate gene signatures from expression data using GMT file
    Args:
        x: Expression data matrix (genes x samples)
        gmt_file: Path to GMT file
        method: Method to calculate signature ("mean", "median", "pca")
        scale: Whether to scale the data
    Returns:
        Signature scores matrix
    """
    # Read GMT file
    genesets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            genesets[parts[0]] = parts[2:]
    
    genenames = x.index
    np = len(genesets)
    
    if scale:
        xs = pd.DataFrame(StandardScaler().fit_transform(x.T).T, 
                         index=x.index, columns=x.columns)
    else:
        xs = x
    
    val = np.zeros((np, x.shape[1]))
    geneset_names = list(genesets.keys())
    
    for i, (name, gene_set) in enumerate(genesets.items()):
        gene_set = [g for g in gene_set if g in genenames]
        if len(gene_set) > 1:
            if method == "mean":
                val[i] = xs.loc[gene_set].mean()
            elif method == "median":
                val[i] = xs.loc[gene_set].median()
            elif method == "pca":
                pca = PCA(n_components=1)
                val[i] = pca.fit_transform(xs.loc[gene_set].T)[:, 0]
        elif len(gene_set) == 1:
            val[i] = xs.loc[gene_set].values[0]
    
    return pd.DataFrame(val, index=geneset_names, columns=x.columns)

def overlap_sets(x, y):
    """Find overlapping genes between two datasets"""
    common_genes = list(set(x.index) & set(y.index))
    x = x.loc[common_genes].sort_index()
    y = y.loc[common_genes].sort_index()
    return x, y

def fisher_test(score, CN):
    """Perform Fisher's exact test between signature score and CN status"""
    top_q = np.percentile(score, 75)
    module = np.where(score >= top_q, "high", "low")
    
    table = np.zeros((2, 2))
    table[0, 0] = np.sum((module == "high") & (CN == "mut"))
    table[0, 1] = np.sum((module == "high") & (CN == "wt"))
    table[1, 0] = np.sum((module == "low") & (CN == "mut"))
    table[1, 1] = np.sum((module == "low") & (CN == "wt"))
    
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')
    return p_value, odds_ratio

def sig_CN_test(score, CN_score, CN_gain, CN_loss, subtype_variable=None):
    """
    Test associations between signature and CN data
    Args:
        score: Signature scores
        CN_score: Continuous CN scores
        CN_gain: Binary CN gain matrix
        CN_loss: Binary CN loss matrix
        subtype_variable: Optional subtype information
    Returns:
        DataFrame with p-values and statistics
    """
    spearman_pos = []
    spearman_neg = []
    spearman_cor = []
    lm_pos = []
    lm_neg = []
    beta_coeff = []
    r_squared = []
    gain = []
    loss = []
    
    for j in range(len(CN_score)):
        CN = CN_score.iloc[j]
        
        # Spearman correlation
        pos_cor, pos_p = stats.spearmanr(score, CN, alternative='greater')
        neg_cor, neg_p = stats.spearmanr(score, CN, alternative='less')
        cor = stats.spearmanr(score, CN)[0]
        
        spearman_pos.append(pos_p)
        spearman_neg.append(neg_p)
        spearman_cor.append(cor)
        
        # Linear model
        if subtype_variable is not None:
            X = pd.concat([CN, subtype_variable], axis=1)
        else:
            X = pd.DataFrame(CN)
        
        model = sm.OLS(score, sm.add_constant(X)).fit()
        beta = model.params[1]  # CN coefficient
        p = model.pvalues[1]
        
        if beta > 0:
            lm_pos.append(p/2)
            lm_neg.append(1-p/2)
        else:
            lm_pos.append(1-p/2)
            lm_neg.append(p/2)
            
        beta_coeff.append(beta)
        r_squared.append(model.rsquared)
        
        # Fisher's exact test
        gain_p, gain_or = fisher_test(score, CN_gain.iloc[j])
        loss_p, loss_or = fisher_test(score, CN_loss.iloc[j])
        gain.append(gain_p)
        loss.append(loss_p)
    
    results = pd.DataFrame({
        'spearman_pos': spearman_pos,
        'spearman_neg': spearman_neg,
        'CN_gain': gain,
        'CN_loss': loss,
        'lm_pos': lm_pos,
        'lm_neg': lm_neg,
        'beta_coeff': beta_coeff,
        'r_squared': r_squared
    })
    
    return results

def median_center(x):
    """Center data by median"""
    medians = x.median(axis=1)
    x_centered = x.sub(medians, axis=0)
    return x_centered

def standardize(x):
    """Standardize data"""
    return pd.DataFrame(
        StandardScaler().fit_transform(x.T).T,
        index=x.index,
        columns=x.columns
    )

def exp_wrap(edata):
    """Process expression data"""
    keep = ['29126', '1493', '5133']  # PD1, PDL1, CTLA4
    
    # Filter low expression
    edata70 = edata[edata.le(2).sum(axis=1) < (0.3 * edata.shape[1])]
    
    # Add back important genes if filtered out
    for gene in keep:
        if gene not in edata70.index:
            edata70 = pd.concat([edata70, edata.loc[[gene]]])
    
    # Handle missing data
    edata70[edata70 <= 2] = 0
    edata70log2 = np.log2(edata70)
    edata70log2[edata70log2 == -np.inf] = 0
    
    # Center and standardize
    exp = median_center(edata70log2)
    exp = standardize(exp)
    return exp 