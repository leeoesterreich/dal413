import pandas as pd
import numpy as np
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pickle
import re

# Parameters
CNA_FILE = 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt'
SIGNATURE_SCORE_FILE = 'signature_score.pkl'
GAIN_THRESHOLD = 0.3
LOSS_THRESHOLD = -0.3
QVAL_CUTOFF = 0.01

sig_name = 'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN'
sig_short = 'Estrogen_signaling'
driver_genes = ['MAPK3', 'CREBBP', 'ZFHX3']
out_csv = '{}_association_results.csv'.format(sig_short)
out_png = '{}_association_landscape_genomewide.png'.format(sig_short)

def sanitize_filename(name):
    return re.sub(r'[^\w\-_\. ]', '_', name)

print('Loading CNA data...')
df_cna = pd.read_csv(CNA_FILE, sep='\t', index_col=0)
print('CNA data loaded. Shape:', df_cna.shape)
print('Generating gain/loss matrices...')
df_gain = (df_cna > GAIN_THRESHOLD).astype(int)
df_loss = (df_cna < LOSS_THRESHOLD).astype(int)
print('Gain/loss matrices generated.')

print('Loading signature scores...')
with open(SIGNATURE_SCORE_FILE, 'rb') as f:
    sig_score_df = pickle.load(f)
if not isinstance(sig_score_df, pd.DataFrame):
    raise ValueError('Signature score file format not recognized.')

if sig_name not in sig_score_df.index:
    print('Signature {} not found in signature score file!'.format(sig_name))
    exit(1)
sig_score = sig_score_df.loc[sig_name]
sig_score.index = sig_score.index.astype(str)

# Intersect samples
overlap_samples = list(set(df_cna.columns) & set(sig_score.index))
df_cna_sub = df_cna[overlap_samples]
df_gain_sub = df_gain[overlap_samples]
df_loss_sub = df_loss[overlap_samples]
sig_score_sub = sig_score.loc[overlap_samples]
print('Number of overlapping samples:', len(overlap_samples))

# Association analysis
results = []
for i, gene in enumerate(df_cna_sub.index):
    if i % 500 == 0:
        print('  Processing gene {}/{}: {}'.format(i+1, len(df_cna_sub.index), gene))
    cna_vec = df_cna_sub.loc[gene]
    gain_vec = df_gain_sub.loc[gene]
    loss_vec = df_loss_sub.loc[gene]
    # Spearman correlation
    pos_corr, pos_p = spearmanr(sig_score_sub, cna_vec, alternative='greater')
    neg_corr, neg_p = spearmanr(sig_score_sub, cna_vec, alternative='less')
    # Fisher's exact test (gain)
    high = sig_score_sub >= sig_score_sub.quantile(0.75)
    gain_table = pd.crosstab(high, gain_vec)
    if gain_table.shape == (2,2):
        gain_p = fisher_exact(gain_table, alternative='greater')[1]
    else:
        gain_p = np.nan
    # Fisher's exact test (loss)
    loss_table = pd.crosstab(high, loss_vec)
    if loss_table.shape == (2,2):
        loss_p = fisher_exact(loss_table, alternative='greater')[1]
    else:
        loss_p = np.nan
    results.append([gene, pos_corr, pos_p, neg_corr, neg_p, gain_p, loss_p])
res_df = pd.DataFrame(results, columns=['gene', 'spearman_pos_corr', 'spearman_pos_p', 'spearman_neg_corr', 'spearman_neg_p', 'gain_p', 'loss_p'])

# Benjamini-Hochberg correction
for col in ['spearman_pos_p', 'spearman_neg_p', 'gain_p', 'loss_p']:
    mask = ~res_df[col].isna()
    res_df.loc[mask, col+'_adj'] = multipletests(res_df.loc[mask, col], method='fdr_bh')[1]
res_df.to_csv(out_csv, index=False)
print('Results saved as', out_csv)

# Print information about significant genes
sig_mask = ((res_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (res_df['gain_p_adj'] < QVAL_CUTOFF)) | \
           ((res_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (res_df['loss_p_adj'] < QVAL_CUTOFF))
sig_df = res_df[sig_mask].copy()

print('Number of significant genes:', len(sig_df))
print('Driver genes in significant set:', [g for g in driver_genes if g in sig_df['gene'].values])

# Plot all genes (genomewide) instead of just significant ones
plt.figure(figsize=(24, 6))

# Create masks for different categories across all genes
pos_mask = (res_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (res_df['gain_p_adj'] < QVAL_CUTOFF)
neg_mask = (res_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (res_df['loss_p_adj'] < QVAL_CUTOFF)
gain_mask = (res_df['gain_p_adj'] < QVAL_CUTOFF)
loss_mask = (res_df['loss_p_adj'] < QVAL_CUTOFF)

print('Number of genes with positive correlations:', sum(pos_mask))
print('Number of genes with negative correlations:', sum(neg_mask))
print('Number of genes with significant gains:', sum(gain_mask))
print('Number of genes with significant losses:', sum(loss_mask))

# Use full x range (all genes)
x_pos = np.arange(len(res_df))
bar_width = 0.8

# Plot all genes with significant values
if sum(pos_mask) > 0:
    print('Plotting positive correlations...')
    plt.bar(x_pos[pos_mask], -np.log10(res_df.loc[pos_mask, 'spearman_pos_p_adj']), 
            width=bar_width, color='red', alpha=0.7, label='Positive correlation')

if sum(neg_mask) > 0:
    print('Plotting negative correlations...')
    plt.bar(x_pos[neg_mask], np.log10(res_df.loc[neg_mask, 'spearman_neg_p_adj']), 
            width=bar_width, color='darkblue', alpha=0.7, label='Negative correlation')

if sum(gain_mask) > 0:
    print('Plotting gains...')
    plt.bar(x_pos[gain_mask], -np.log10(res_df.loc[gain_mask, 'gain_p_adj']), 
            width=bar_width, color='orange', alpha=0.5, label='Gain')

if sum(loss_mask) > 0:
    print('Plotting losses...')
    plt.bar(x_pos[loss_mask], np.log10(res_df.loc[loss_mask, 'loss_p_adj']), 
            width=bar_width, color='lightblue', alpha=0.5, label='Loss')

# Annotate driver genes with repelled labels and larger font
texts = []
for gene in driver_genes:
    if gene in res_df['gene'].values:
        idx = np.where(res_df['gene'] == gene)[0][0]
        texts.append(plt.text(idx, 0, gene, fontsize=16, fontweight='bold', ha='center', va='bottom'))

try:
    if texts:
        print('Adjusting text positions...')
        adjust_text(texts, arrowprops=dict(arrowstyle='-|>', color='black'), force_text=1.5, force_points=1.5)
except Exception as e:
    print(f"Warning: Error adjusting text positions: {str(e)}")
    print("Continuing with unadjusted text...")

# Add significance threshold lines
plt.axhline(-np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
plt.axhline(np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')

# Labels and formatting
plt.xlabel('Genes (genomic order)', fontsize=18)
plt.ylabel('-log10(q)', fontsize=18)
plt.title('{} association landscape (genome-wide)'.format(sig_short), fontsize=20)
plt.xticks([])  # Hide x-ticks since there are too many genes
plt.yticks(fontsize=14)
plt.legend(fontsize=14, loc='upper right')
plt.tight_layout()

print('Saving plot to', out_png)
plt.savefig(out_png, dpi=300, bbox_inches='tight')
print('Plot saved as', out_png) 