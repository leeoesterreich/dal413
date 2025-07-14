import pandas as pd
import numpy as np
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests
# Set non-interactive backend first, before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pickle
import re
import sys
import os
import traceback

# Parameters
CNA_FILE = 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt'
SIGNATURE_SCORE_FILE = 'signature_score.pkl'
GENE_SEGMENTS_FILE = 'gene_segments.bed'
GAIN_THRESHOLD = 0.3
LOSS_THRESHOLD = -0.3
QVAL_CUTOFF = 0.01

SIGNATURES = [
    {
        'name': 'UNC_HER1_Cluster2_Median_BMC.Genomics.2007_PMID.17663798',
        'short': 'HER1_C2',
        'drivers': ['HER1', 'EGFR', 'HER2', 'HER4', 'TGFA', 'AREG', 'EGF', 'KRAS', 'PIK3ca', 'AKT1', 'MEK1', 'ERK1']
    },
    {
        'name': 'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP',
        'short': 'Basal_signaling',
        'drivers': ['RAD17', 'RAD50', 'PALB2', 'TP53', 'BRCA1']
    }
]

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

gene_segments = pd.read_csv(GENE_SEGMENTS_FILE, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])
gene_segments['chrom_num'] = gene_segments['chrom'].str.replace('chr','').replace({'X':23,'Y':24}).astype('category')
gene_segments['chrom_num'] = gene_segments['chrom_num'].cat.set_categories([str(i) for i in range(1,23)] + ['23','24'], ordered=True)
gene_segments = gene_segments.sort_values(['chrom_num','start'])

for sig in SIGNATURES:
    print('\n=== Processing signature: {} ==='.format(sig['name']))
    sig_name = sig['name']
    sig_short = sig['short']
    driver_genes = sig['drivers']
    out_csv = '{}_association_results.csv'.format(sanitize_filename(sig_short))
    out_png = '{}_association_landscape_genomewide.png'.format(sanitize_filename(sig_short))

    if sig_name not in sig_score_df.index:
        print('Signature {} not found in signature score file! Skipping.'.format(sig_name))
        continue
    sig_score = sig_score_df.loc[sig_name]
    sig_score.index = sig_score.index.astype(str)

    # Optimize overlap_samples subsetting
    overlap_samples = [s for s in df_cna.columns if s in sig_score.index]
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
    res_df = res_df.merge(gene_segments[['gene','chrom','start','chrom_num']], on='gene', how='left')
    res_df = res_df.sort_values(['chrom_num','start'])
    res_df.to_csv(out_csv, index=False)
    print('Results saved as', out_csv)

    try:
        # Only plot bars for genes where any q-value is < 0.01
        sig_mask = (
            (res_df['spearman_pos_p_adj'] < QVAL_CUTOFF) |
            (res_df['spearman_neg_p_adj'] < QVAL_CUTOFF) |
            (res_df['gain_p_adj'] < QVAL_CUTOFF) |
            (res_df['loss_p_adj'] < QVAL_CUTOFF)
        )
        sig_df = res_df[sig_mask].copy()
        
        if len(sig_df) == 0:
            print('No significant genes found for {}, skipping plot.'.format(sig_short))
            continue
            
        print('Plotting results for {} signature ({} significant genes)...'.format(sig_short, len(sig_df)))
        
        # Create a new figure
        plt.figure(figsize=(18, 5))
        
        # Create masks for different categories
        pos_mask = (sig_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (sig_df['gain_p_adj'] < QVAL_CUTOFF)
        neg_mask = (sig_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (sig_df['loss_p_adj'] < QVAL_CUTOFF)
        gain_mask = (sig_df['gain_p_adj'] < QVAL_CUTOFF)
        loss_mask = (sig_df['loss_p_adj'] < QVAL_CUTOFF)
        
        print('  Number of positive correlations:', sum(pos_mask))
        print('  Number of negative correlations:', sum(neg_mask))
        print('  Number of gains:', sum(gain_mask))
        print('  Number of losses:', sum(loss_mask))
        
        # Use scatter plot (faster than bars for large datasets)
        if sum(pos_mask) > 0:
            print('  Plotting positive correlations...')
            plt.scatter(
                sig_df.index[pos_mask], 
                -np.log10(sig_df.loc[pos_mask, 'spearman_pos_p_adj']), 
                color='red', s=30, label='Positive correlation'
            )
        
        if sum(neg_mask) > 0:
            print('  Plotting negative correlations...')
            plt.scatter(
                sig_df.index[neg_mask], 
                np.log10(sig_df.loc[neg_mask, 'spearman_neg_p_adj']), 
                color='darkblue', s=30, label='Negative correlation'
            )
        
        if sum(gain_mask) > 0:
            print('  Plotting gains...')
            plt.scatter(
                sig_df.index[gain_mask], 
                -np.log10(sig_df.loc[gain_mask, 'gain_p_adj']), 
                color='orange', s=20, alpha=0.7, label='Gain'
            )
        
        if sum(loss_mask) > 0:
            print('  Plotting losses...')
            plt.scatter(
                sig_df.index[loss_mask], 
                np.log10(sig_df.loc[loss_mask, 'loss_p_adj']), 
                color='lightblue', s=20, alpha=0.7, label='Loss'
            )
        
        # Annotate driver genes with repelled labels
        print('  Annotating driver genes...')
        texts = []
        for gene in driver_genes:
            if gene in sig_df['gene'].values:
                idx = sig_df.index[sig_df['gene']==gene][0]
                texts.append(plt.text(idx, 0, gene, fontsize=16, fontweight='bold', ha='center', va='bottom'))
        
        try:
            if texts:
                print('  Adjusting text positions...')
                adjust_text(texts, arrowprops=dict(arrowstyle='-|>', color='black'), force_text=1.5, force_points=1.5)
        except Exception as e:
            print('  Warning: Error adjusting text positions:', str(e))
            print('  Continuing with unadjusted text...')
        
        # Add significance threshold lines
        plt.axhline(-np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
        plt.axhline(np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
        
        # Labels and formatting
        plt.xlabel('Genes (ordered)', fontsize=18)
        plt.ylabel('-log10(q)', fontsize=18)
        plt.title('{} association landscape (significant genes)'.format(sig_short), fontsize=20)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        plt.tight_layout()
        
        # Save figure
        print('  Saving plot to {}...'.format(out_png))
        plt.savefig(out_png, dpi=300)
        plt.close()
        print('Plot saved as', out_png)
        
    except Exception as e:
        print('Error during plotting for {} signature:'.format(sig_short))
        print(traceback.format_exc())
        print('Continuing with next signature...')
        
print('\nAssociation analysis complete!') 