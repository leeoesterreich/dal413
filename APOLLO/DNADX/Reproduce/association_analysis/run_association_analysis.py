import pandas as pd
import numpy as np
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import pickle

# Parameters
CNA_FILE = 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt'
SIGNATURE_SCORE_FILE = 'signature_score.pkl'
GENE_SEGMENTS_FILE = 'gene_segments.bed'
DRIVER_GENES = ['SOS1', 'E2F3', 'CCND3', 'CDK6', 'E2F5', 'MYC', 'CCND1', 'CCND2', 'RB1', 'E2F1']
SIGNATURE_NAME = 'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450'
GAIN_THRESHOLD = 0.3
LOSS_THRESHOLD = -0.3
QVAL_CUTOFF = 0.01

print('Step 1: Loading CNA data...')
df_cna = pd.read_csv(CNA_FILE, sep='\t', index_col=0)
print('CNA data loaded. Shape:', df_cna.shape)

print('Step 2: Generating gain/loss matrices...')
df_gain = (df_cna > GAIN_THRESHOLD).astype(int)
df_loss = (df_cna < LOSS_THRESHOLD).astype(int)
print('Gain/loss matrices generated.')

print('Step 3: Loading signature score...')
with open(SIGNATURE_SCORE_FILE, 'rb') as f:
    sig_score_df = pickle.load(f)
if isinstance(sig_score_df, pd.DataFrame):
    sig_score = sig_score_df.loc[SIGNATURE_NAME]
    sig_score.index = sig_score.index.astype(str)
    print('Signature score loaded. Length:', len(sig_score))
else:
    raise ValueError('Signature score file format not recognized.')

print('Step 4: Intersecting samples...')
overlap_samples = list(set(df_cna.columns) & set(sig_score.index))
df_cna = df_cna[overlap_samples]
df_gain = df_gain[overlap_samples]
df_loss = df_loss[overlap_samples]
sig_score = sig_score.loc[overlap_samples]
print('Sample intersection complete. Number of overlapping samples:', len(overlap_samples))

print('Step 5: Running association analysis (Spearman and Fisher tests)...')
results = []
for i, gene in enumerate(df_cna.index):
    if i % 500 == 0:
        print(f'  Processing gene {i+1}/{len(df_cna.index)}: {gene}')
    cna_vec = df_cna.loc[gene]
    gain_vec = df_gain.loc[gene]
    loss_vec = df_loss.loc[gene]
    # Spearman correlation
    pos_corr, pos_p = spearmanr(sig_score, cna_vec, alternative='greater')
    neg_corr, neg_p = spearmanr(sig_score, cna_vec, alternative='less')
    # Fisher's exact test (gain)
    high = sig_score >= sig_score.quantile(0.75)
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
print('Association analysis complete.')

res_df = pd.DataFrame(results, columns=['gene', 'spearman_pos_corr', 'spearman_pos_p', 'spearman_neg_corr', 'spearman_neg_p', 'gain_p', 'loss_p'])

print('Step 6: Applying Benjamini-Hochberg correction...')
for col in ['spearman_pos_p', 'spearman_neg_p', 'gain_p', 'loss_p']:
    mask = ~res_df[col].isna()
    res_df.loc[mask, col+'_adj'] = multipletests(res_df.loc[mask, col], method='fdr_bh')[1]
print('Multiple testing correction complete.')

# Save results as CSV
res_df.to_csv('RB_LOH_association_results.csv', index=False)
print('Results DataFrame saved as RB_LOH_association_results.csv')

print('Step 7: Plotting results (only significant genes)...')
try:
    gene_segments = pd.read_csv(GENE_SEGMENTS_FILE, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])
    res_df = res_df.merge(gene_segments, on='gene', how='left')
    res_df['chrom_num'] = res_df['chrom'].str.replace('chr','').replace({'X':23,'Y':24}).astype('category')
    res_df = res_df.sort_values(['chrom_num','start'])
    print('Gene segment annotation loaded and merged.')
except Exception as e:
    print('Gene segment annotation not used in plot:', e)

# Only plot genes significant in both analyses (q < 0.01)
sig_mask = ((res_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (res_df['gain_p_adj'] < QVAL_CUTOFF)) | \
           ((res_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (res_df['loss_p_adj'] < QVAL_CUTOFF))
sig_df = res_df[sig_mask].copy()

plt.figure(figsize=(18,5))
# Plot positive Spearman (red) and negative (blue) for significant genes only
pos_mask = (sig_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (sig_df['gain_p_adj'] < QVAL_CUTOFF)
neg_mask = (sig_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (sig_df['loss_p_adj'] < QVAL_CUTOFF)
plt.bar(sig_df.index[pos_mask], -np.log10(sig_df.loc[pos_mask, 'spearman_pos_p_adj']), color='red', label='Positive correlation')
plt.bar(sig_df.index[neg_mask], np.log10(sig_df.loc[neg_mask, 'spearman_neg_p_adj']), color='darkblue', label='Negative correlation')
# Plot gain (orange) and loss (lightblue)
gain_mask = (sig_df['gain_p_adj'] < QVAL_CUTOFF)
loss_mask = (sig_df['loss_p_adj'] < QVAL_CUTOFF)
plt.bar(sig_df.index[gain_mask], -np.log10(sig_df.loc[gain_mask, 'gain_p_adj']), color='orange', alpha=0.5, label='Gain')
plt.bar(sig_df.index[loss_mask], np.log10(sig_df.loc[loss_mask, 'loss_p_adj']), color='lightblue', alpha=0.5, label='Loss')
# Annotate driver genes
for gene in DRIVER_GENES:
    if gene in sig_df['gene'].values:
        idx = sig_df.index[sig_df['gene']==gene][0]
        plt.annotate(gene, (idx, 0), xytext=(0,30), textcoords='offset points', arrowprops=dict(arrowstyle='-|>'))
plt.axhline(-np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
plt.axhline(np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
plt.xlabel('Significant genes (ordered by chromosome)')
plt.ylabel('-log10(q)')
plt.title('RB-LOH signature association landscape (significant genes only)')
plt.legend()
plt.tight_layout()
plt.savefig('RB_LOH_association_landscape.png', dpi=300)
print('Analysis complete. Plot saved as RB_LOH_association_landscape.png') 