import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Read the data
file_path = 'Concentration_summary_IDC.csv'
df = pd.read_csv(file_path)

# Clean up Spin_status (remove extra spaces, lowercase for consistency)
df['Spin_status'] = df['Spin_status'].str.strip().str.capitalize()

# Handle '<100' values and convert to float
concs = []
for val in df['Concentration']:
    if isinstance(val, str) and val.strip() == '<100':
        concs.append(100.0)
    else:
        concs.append(float(val))
df['Concentration_clean'] = concs

# Log10 transform for plotting
df['log10_conc'] = np.log10(df['Concentration_clean'])

# Print some debug information
print("\nDouble spin samples:")
double_samples = df[df['Spin_status'] == 'Double']
print(double_samples[['Sample_name', 'Concentration', 'Concentration_clean']].head())
print(f"\nDouble spin mean: {np.mean(double_samples['Concentration_clean']):.0f}")

print("\nSingle spin samples:")
single_samples = df[df['Spin_status'] == 'Single']
print(single_samples[['Sample_name', 'Concentration', 'Concentration_clean']].head())
print(f"\nSingle spin mean: {np.mean(single_samples['Concentration_clean']):.0f}")

# Filter groups
double_log = df[df['Spin_status'] == 'Double']['log10_conc']
single_log = df[df['Spin_status'] == 'Single']['log10_conc']

# Calculate means and medians in original scale
double_mean = np.mean(df[df['Spin_status'] == 'Double']['Concentration_clean'])
single_mean = np.mean(df[df['Spin_status'] == 'Single']['Concentration_clean'])
double_median = np.median(df[df['Spin_status'] == 'Double']['Concentration_clean'])
single_median = np.median(df[df['Spin_status'] == 'Single']['Concentration_clean'])

# Mann-Whitney U test
stat, p_val = mannwhitneyu(double_log, single_log, alternative='two-sided')

# Create violin plot with larger figure size
plt.figure(figsize=(10, 8))
violin = plt.violinplot([double_log, single_log], showmeans=True, showmedians=True)
plt.xticks([1, 2], ['Double', 'Single'])
plt.ylabel('Log10 Concentration (ng/mL)')
plt.title('Comparison of double/single spun plasma extraction')

# Customize violin plot colors
for pc in violin['bodies']:
    pc.set_facecolor('skyblue')
    pc.set_alpha(0.7)

# Add individual data points
x1 = 1 + np.random.normal(0, 0.05, len(double_log))
x2 = 2 + np.random.normal(0, 0.05, len(single_log))
plt.plot(x1, double_log, 'k.', alpha=0.3, markersize=8)
plt.plot(x2, single_log, 'k.', alpha=0.3, markersize=8)

# Add mean lines
plt.plot([1], [double_log.mean()], 'ko', label='Mean')
plt.plot([2], [single_log.mean()], 'ko')

# Calculate y-axis limits with more space for annotation
y_min = min(double_log.min(), single_log.min())
y_max = max(double_log.max(), single_log.max())
y_range = y_max - y_min

# Annotate p-value and means below x-axis with more space
plt.text(1.5, y_min - y_range * 0.2, 
         f'Mann-Whitney U test p = {p_val:.3g}\n'
         f'Double: mean = {double_mean:.0f} ng/mL, median = {double_median:.0f} ng/mL\n'
         f'Single: mean = {single_mean:.0f} ng/mL, median = {single_median:.0f} ng/mL', 
         ha='center', color='black', fontsize=12)

# Adjust y-axis limits to accommodate the annotation
plt.ylim(y_min - y_range * 0.25, y_max + y_range * 0.05)

# Adjust layout with more space at the bottom
plt.subplots_adjust(bottom=0.25)

# Save the plot with high DPI and ensure it's visible
plt.savefig('boxplot_spin_status.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f'Mann-Whitney U test result: U = {stat:.3f}, p = {p_val:.3g}')
print(f'\nMean concentration:')
print(f'Double spin: {double_mean:.0f} ng/mL')
print(f'Single spin: {single_mean:.0f} ng/mL')
print(f'Difference: {double_mean - single_mean:.0f} ng/mL')
print(f'\nMedian concentration:')
print(f'Double spin: {double_median:.0f} ng/mL')
print(f'Single spin: {single_median:.0f} ng/mL')
print(f'Difference: {double_median - single_median:.0f} ng/mL') 