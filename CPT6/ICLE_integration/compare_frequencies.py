import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
import numpy as np

# Set style
plt.style.use('seaborn-whitegrid')
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# Read the gene list from agg_sorted.csv
agg_sorted = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/agg_sorted.csv')
cpt6_genes = agg_sorted['Gene Symbol'].unique().tolist()

# Read the data files
icle_data = pd.read_csv('results/ICLE_mutation_frequency_summary.txt', sep='\t')
tcga_data = pd.read_csv('maftools/TCGA2015_Metabirc2012_16_MSK2018_19_20_ILC.txt', sep='\t')

# Convert frequency strings to floats
icle_data['Freq'] = icle_data['Freq'].str.rstrip('%').astype(float)
tcga_data['Freq'] = tcga_data['Freq'].str.rstrip('%').astype(float)

# Get total sample counts
icle_total = 17  # ICLE sample count
tcga_total = 643  # TCGA sample count

# Filter for CPT6 genes and prepare data for plotting
genes = []
icle_freqs = []
tcga_freqs = []
pvalues = []
odds_ratios = []

for gene in cpt6_genes:
    if gene in icle_data['Gene'].values and gene in tcga_data['Gene'].values:
        genes.append(gene)
        
        # Get ICLE frequency and calculate counts
        icle_freq = icle_data[icle_data['Gene'] == gene]['Freq'].iloc[0]
        icle_count = int(np.round(icle_freq * icle_total / 100))  # Convert percentage to count
        icle_freqs.append(icle_freq)
        
        # Get TCGA frequency and calculate counts
        tcga_freq = tcga_data[tcga_data['Gene'] == gene]['Freq'].iloc[0]
        tcga_count = int(np.round(tcga_freq * tcga_total / 100))  # Convert percentage to count
        tcga_freqs.append(tcga_freq)
        
        # Fisher's exact test
        table = [[icle_count, icle_total - icle_count],
                [tcga_count, tcga_total - tcga_count]]
        odds_ratio, pvalue = fisher_exact(table)
        pvalues.append(pvalue)
        odds_ratios.append(odds_ratio)

# Sort genes by ICLE frequency (reversed order)
sorted_indices = np.argsort(icle_freqs)[::-1]  # Reverse to get descending order
genes = [genes[i] for i in sorted_indices]
icle_freqs = [icle_freqs[i] for i in sorted_indices]
tcga_freqs = [tcga_freqs[i] for i in sorted_indices]
pvalues = [pvalues[i] for i in sorted_indices]
odds_ratios = [odds_ratios[i] for i in sorted_indices]

# Create figure
fig, ax = plt.subplots(figsize=(10, 12))

# Plot bars
y_pos = np.arange(len(genes))
bar_height = 0.35

# ICLE (red)
cell_lines = ax.barh(y_pos + bar_height/2, icle_freqs, height=bar_height, 
                    color='red', alpha=0.6, label='ICLE (n={})'.format(icle_total))

# TCGA primary ILC tumor (blue)
primary = ax.barh(y_pos - bar_height/2, tcga_freqs, height=bar_height,
                 color='blue', alpha=0.6, label='TCGA primary ILC tumor (n={})'.format(tcga_total))

# Add significance stars for p < 0.05 and meaningful fold change
for i, (pval, freq, odds) in enumerate(zip(pvalues, icle_freqs, odds_ratios)):
    if pval < 0.05 and (odds > 2 or odds < 0.5):  # Significant and meaningful difference
        stars = '*' * sum([pval < cutoff for cutoff in [0.05, 0.01, 0.001]])
        ax.text(max(freq, tcga_freqs[i]) + 0.5, i, stars, va='center', ha='left', fontsize=12)

# Customize plot
ax.set_yticks(y_pos)
ax.set_yticklabels(genes, fontsize=10)
ax.set_xlabel('Alteration Frequency (%)', fontsize=12)
ax.set_title('ICLE vs TCGA primary ILC tumor', fontsize=14, pad=20)

# Set x-axis limits based on data
max_freq = max(max(icle_freqs), max(tcga_freqs)) * 1.2  # Add 20% margin
ax.set_xlim(0, max_freq)

# Add legend at the bottom right
ax.legend(loc='lower right', bbox_to_anchor=(1, 0), frameon=True)

# Add gridlines
ax.grid(True, linestyle='--', alpha=0.3, axis='x')
ax.set_axisbelow(True)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adjust layout with more space at the bottom for legend
plt.tight_layout()

# Save plot
plt.savefig('results/mutation_frequency_comparison.pdf', bbox_inches='tight', dpi=300)
plt.savefig('results/mutation_frequency_comparison.png', bbox_inches='tight', dpi=300)

# Print statistical analysis results
print("\nStatistical Analysis Results:")
print("ICLE samples: {}".format(icle_total))
print("TCGA primary ILC tumor samples: {}\n".format(tcga_total))

for gene, pval, icle_f, tcga_f, odds in zip(genes, pvalues, icle_freqs, tcga_freqs, odds_ratios):
    if pval < 0.05 and (odds > 2 or odds < 0.5):  # Significant and meaningful difference
        stars = '*' * sum([pval < cutoff for cutoff in [0.05, 0.01, 0.001]])
        effect = "higher" if odds > 1 else "lower"
        print("{}: ICLE={:.1f}%, TCGA={:.1f}%, {:.1f}-fold {}, p={:.2e} ({})".format(
            gene, icle_f, tcga_f, abs(odds if odds > 1 else 1/odds), effect, pval, stars))

# Print legend for p-value significance
print("\nSignificance levels:")
print("* p < 0.05")
print("** p < 0.01")
print("*** p < 0.001")
print("\nNote: Stars are only shown for significant differences (p < 0.05)")
print("with at least 2-fold change in frequency between groups.") 