import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

# Read the data
idc_data = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Tendo/tumor_fraction_cna_data.csv')
ilc_data = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv')

# Create dataframes with type and tumor fraction
idc_tf = pd.DataFrame({
    'cancer_type': 'IDC',
    'tumor_fraction': idc_data['Tumor_Fraction']
})

# Filter out ILC samples with 0 tumor fraction
ilc_tf = pd.DataFrame({
    'cancer_type': 'ILC',
    'tumor_fraction': ilc_data['Tumor Fraction']
})
ilc_tf = ilc_tf[ilc_tf['tumor_fraction'] > 0]  # Remove samples with 0 tumor fraction

# Combine the data
combined_data = pd.concat([idc_tf, ilc_tf])

# Multiply tumor fractions by 100 and apply log10 transformation
# Add a small constant (0.1) to handle zero values
combined_data['log10_tumor_fraction'] = np.log10((combined_data['tumor_fraction'] * 100) + 0.1)

# Perform Mann-Whitney U test on both original and log-transformed data
stat_orig, pval_orig = stats.mannwhitneyu(
    idc_tf['tumor_fraction'], 
    ilc_tf['tumor_fraction'],
    alternative='two-sided'
)

stat_log, pval_log = stats.mannwhitneyu(
    combined_data[combined_data['cancer_type'] == 'IDC']['log10_tumor_fraction'],
    combined_data[combined_data['cancer_type'] == 'ILC']['log10_tumor_fraction'],
    alternative='two-sided'
)

# Create the square boxplot with improved styling
plt.figure(figsize=(8, 8))
sns.set_style("whitegrid")

# Create boxplot with white boxes and black edges
ax = sns.boxplot(x='cancer_type', y='log10_tumor_fraction', data=combined_data, width=0.5,
                 color='white', linewidth=1.5)

# Add individual points as solid dark grey dots
sns.stripplot(x='cancer_type', y='log10_tumor_fraction', data=combined_data, 
              color='#303030', alpha=0.7, size=6, jitter=0.2)

# Customize the plot
plt.title('Log10-transformed Tumor Fraction Distribution\n(ILC vs IDC)', fontsize=14, pad=20)
plt.xlabel('Cancer Type', fontsize=12)
plt.ylabel('Log10(Tumor Fraction %)', fontsize=12)

# Get current y-axis limits
ymin, ymax = plt.ylim()

# Add statistics including original mean and median above the plot
for cancer_type in ['IDC', 'ILC']:
    data = combined_data[combined_data['cancer_type'] == cancer_type]
    orig_mean = data['tumor_fraction'].mean() * 100
    orig_median = data['tumor_fraction'].median() * 100
    log_mean = data['log10_tumor_fraction'].mean()
    log_median = data['log10_tumor_fraction'].median()
    
    stats_text = f'n = {len(data)}\n'
    stats_text += f'Original Mean = {orig_mean:.1f}%\n'
    stats_text += f'Original Median = {orig_median:.1f}%'
    
    x_pos = 0 if cancer_type == 'IDC' else 1
    plt.text(x_pos, ymax + 0.05 * (ymax - ymin), stats_text, 
             horizontalalignment='center', verticalalignment='bottom',
             fontsize=10, color='darkblue')

# Add Mann-Whitney U test results
mw_text = f'Mann-Whitney U test:\n'
mw_text += f'Original scale: p = {pval_orig:.2e}\n'
mw_text += f'Log10 scale: p = {pval_log:.2e}'
plt.text(0.5, ymax + 0.15 * (ymax - ymin), mw_text,
         horizontalalignment='center', verticalalignment='bottom',
         fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.8))

# Adjust y-axis limits to make room for the text
plt.ylim(ymin, ymax + 0.3 * (ymax - ymin))

# Make axis labels and ticks larger
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('tumor_fraction_log_comparison_filtered.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics and test results
print("\nMann-Whitney U Test Results:")
print(f"Original scale - U statistic: {stat_orig:.1f}, p-value: {pval_orig:.2e}")
print(f"Log10 scale - U statistic: {stat_log:.1f}, p-value: {pval_log:.2e}")

for cancer_type in ['IDC', 'ILC']:
    data = combined_data[combined_data['cancer_type'] == cancer_type]
    print(f"\n{cancer_type} Summary Statistics:")
    print(f"Number of samples: {len(data)}")
    print(f"Original Tumor Fraction:")
    print(f"  Mean: {data['tumor_fraction'].mean():.3f}")
    print(f"  Median: {data['tumor_fraction'].median():.3f}")
    print(f"  Range: {data['tumor_fraction'].min():.3f} - {data['tumor_fraction'].max():.3f}")
    print(f"Log10 Transformed (%):")
    print(f"  Mean: {data['log10_tumor_fraction'].mean():.3f}")
    print(f"  Median: {data['log10_tumor_fraction'].median():.3f}")
    print(f"  Range: {data['log10_tumor_fraction'].min():.3f} - {data['log10_tumor_fraction'].max():.3f}")

# Print number of excluded samples
n_excluded = len(ilc_data) - len(ilc_tf)
print(f"\nNumber of ILC samples excluded (tumor fraction = 0): {n_excluded}") 