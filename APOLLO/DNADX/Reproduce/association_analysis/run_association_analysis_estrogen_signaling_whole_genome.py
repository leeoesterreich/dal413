import pandas as pd
import numpy as np
from scipy.stats import spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pickle

# Parameters
CNA_FILE = 'Gistic2_CopyNumber_Gistic2_all_data_by_genes.txt'
SIGNATURE_SCORE_FILE = 'signature_score.pkl'
GENE_SEGMENTS_FILE = 'gene_segments.bed'
DRIVER_GENES = ['MAPK3', 'CREBBP', 'ZFHX3']
SIGNATURE_NAME = 'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN'
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
res_df.to_csv('Estrogen_signaling_association_results.csv', index=False)
print('Results DataFrame saved as Estrogen_signaling_association_results.csv')

print('Step 7: Plotting whole genome results...')
try:
    # Load gene positions and merge with results
    gene_segments = pd.read_csv(GENE_SEGMENTS_FILE, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene'])
    merged_df = res_df.merge(gene_segments, on='gene', how='left')
    
    # Convert chromosome to numeric for sorting
    merged_df['chrom_num'] = merged_df['chrom'].str.replace('chr','').replace({'X':23,'Y':24}).astype('category')
    merged_df = merged_df.sort_values(['chrom_num','start'])
    
    # Create chromosome mapping for continuous coordinates
    chromosomes = merged_df['chrom_num'].unique()
    chromosomes = sorted(chromosomes)
    
    # Calculate chromosome lengths and offsets
    chrom_offsets = {}
    current_offset = 0
    chrom_centers = {}
    
    for chrom in chromosomes:
        chrom_data = merged_df[merged_df['chrom_num'] == chrom]
        if len(chrom_data) == 0:
            continue
        chrom_offsets[chrom] = current_offset
        
        # Get max position for chromosome
        max_pos = chrom_data['end'].max()
        chrom_length = max_pos
        
        # Store center position for labels
        chrom_centers[chrom] = current_offset + chrom_length / 2
        
        # Update offset for next chromosome (add buffer)
        current_offset += chrom_length + 25000000  # 25Mb buffer between chromosomes
    
    # Create continuous genomic coordinates
    merged_df['plot_pos'] = merged_df.apply(
        lambda row: row['start'] + chrom_offsets.get(row['chrom_num'], 0), 
        axis=1
    )
    
    print('Gene segment annotation loaded and merged.')
except Exception as e:
    print('Gene segment annotation not used in plot:', e)
    # If we can't use chromosomal positions, just use indices
    merged_df = res_df.copy()
    merged_df['plot_pos'] = np.arange(len(merged_df))
    chrom_offsets = {}
    chrom_centers = {}

# Create a wider figure for whole genome plot
plt.figure(figsize=(24, 8))

# Create significance masks for plotting (whole genome)
pos_mask = (merged_df['spearman_pos_p_adj'] < QVAL_CUTOFF) & (merged_df['gain_p_adj'] < QVAL_CUTOFF)
neg_mask = (merged_df['spearman_neg_p_adj'] < QVAL_CUTOFF) & (merged_df['loss_p_adj'] < QVAL_CUTOFF)
gain_mask = (merged_df['gain_p_adj'] < QVAL_CUTOFF)
loss_mask = (merged_df['loss_p_adj'] < QVAL_CUTOFF)

print(f"Number of genes with positive correlations: {sum(pos_mask)}")
print(f"Number of genes with negative correlations: {sum(neg_mask)}")
print(f"Number of genes with significant gains: {sum(gain_mask)}")
print(f"Number of genes with significant losses: {sum(loss_mask)}")

# Plot only the genes with significant values, but across the whole genome
if sum(pos_mask) > 0:
    print("Plotting positive correlations...")
    plt.bar(merged_df.loc[pos_mask, 'plot_pos'], 
            -np.log10(merged_df.loc[pos_mask, 'spearman_pos_p_adj']), 
            width=5000000,  # Fixed width for bars (5Mb)
            color='red', alpha=0.7, label='Positive correlation')

if sum(neg_mask) > 0:
    print("Plotting negative correlations...")
    plt.bar(merged_df.loc[neg_mask, 'plot_pos'], 
            np.log10(merged_df.loc[neg_mask, 'spearman_neg_p_adj']), 
            width=5000000,  # Fixed width for bars (5Mb)
            color='darkblue', alpha=0.7, label='Negative correlation')

if sum(gain_mask) > 0:
    print("Plotting gains...")
    plt.bar(merged_df.loc[gain_mask, 'plot_pos'], 
            -np.log10(merged_df.loc[gain_mask, 'gain_p_adj']), 
            width=5000000,  # Fixed width for bars (5Mb)
            color='orange', alpha=0.5, label='Gain')

if sum(loss_mask) > 0:
    print("Plotting losses...")
    plt.bar(merged_df.loc[loss_mask, 'plot_pos'], 
            np.log10(merged_df.loc[loss_mask, 'loss_p_adj']), 
            width=5000000,  # Fixed width for bars (5Mb)
            color='lightblue', alpha=0.5, label='Loss')

# Add chromosome boundary lines
if chrom_offsets:
    for chrom in chromosomes[1:]:  # Skip first chromosome
        if chrom in chrom_offsets:
            plt.axvline(chrom_offsets[chrom] - 12500000, color='gray', linestyle='--', alpha=0.5)
    
    # Add chromosome labels
    for chrom in chromosomes:
        if chrom in chrom_centers:
            chrom_label = str(chrom) if chrom <= 22 else ('X' if chrom == 23 else 'Y')
            plt.text(chrom_centers[chrom], -1.5, f'Chr {chrom_label}', 
                    ha='center', va='center', fontsize=10, fontweight='bold')

# Annotate driver genes with repelled labels
texts = []
for gene in DRIVER_GENES:
    if gene in merged_df['gene'].values:
        gene_idx = merged_df[merged_df['gene'] == gene].index[0]
        pos = merged_df.loc[gene_idx, 'plot_pos']
        texts.append(plt.text(pos, 0, gene, fontsize=16, fontweight='bold', ha='center', va='bottom'))

try:
    if texts:
        print("Adjusting text positions...")
        adjust_text(texts, arrowprops=dict(arrowstyle='-|>', color='black'), 
                    force_text=1.5, force_points=1.5)
except Exception as e:
    print(f"Warning: Error adjusting text positions: {str(e)}")
    print("Continuing with unadjusted text...")

# Add significance threshold lines
plt.axhline(-np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')
plt.axhline(np.log10(QVAL_CUTOFF), color='k', linestyle='dashed')

# Format plot
plt.xlabel('Genomic Position', fontsize=18)
plt.ylabel('-log10(q)', fontsize=18)
plt.title('Estrogen signaling signature association landscape (whole genome)', fontsize=20)
plt.xticks([])  # Hide x-ticks as we have chromosome labels
plt.yticks(fontsize=14)
plt.legend(fontsize=14, loc='upper right')
plt.tight_layout()

# Save the plot
plt.savefig('Estrogen_signaling_association_landscape_genomewide.png', dpi=300, bbox_inches='tight')
print('Analysis complete. Plot saved as Estrogen_signaling_association_landscape_genomewide.png') 