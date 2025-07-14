import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the tumor fraction data
tumor_fraction = pd.read_csv('results/summary.csv')

# Create figure
plt.figure(figsize=(15, 8))

# Create scatter plot of tumor fraction vs ploidy
plt.scatter(tumor_fraction['Ploidy'], tumor_fraction['Tumor Fraction'], 
           alpha=0.6, s=100)

# Add horizontal line at 0.03
plt.axhline(y=0.03, color='red', linestyle='--', alpha=0.7, label='TF = 0.03')

# Add horizontal line at 0.05
plt.axhline(y=0.05, color='green', linestyle='--', alpha=0.7, label='TF = 0.05')

# Customize the plot
plt.xlabel('Ploidy', fontsize=12)
plt.ylabel('Tumor Fraction', fontsize=12)
plt.title('Tumor Fraction vs Ploidy', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend()

# Save the plot
plt.savefig('results/plots/tumor_fraction_vs_ploidy_thresholds.png', dpi=300, bbox_inches='tight')
plt.close()

# Calculate statistics
total_samples = len(tumor_fraction)
samples_above_003 = sum(tumor_fraction['Tumor Fraction'] > 0.03)
samples_above_005 = sum(tumor_fraction['Tumor Fraction'] > 0.05)
samples_zero = sum(tumor_fraction['Tumor Fraction'] == 0)

# Save statistics to file
with open('results/plots/tumor_fraction_threshold_stats.txt', 'w') as f:
    f.write('Tumor Fraction Statistics:\n\n')
    f.write(f'Total number of samples: {total_samples}\n')
    f.write(f'Number of samples with TF > 0.03: {samples_above_003} ({samples_above_003/total_samples*100:.1f}%)\n')
    f.write(f'Number of samples with TF > 0.05: {samples_above_005} ({samples_above_005/total_samples*100:.1f}%)\n')
    f.write(f'Number of samples with TF = 0: {samples_zero} ({samples_zero/total_samples*100:.1f}%)\n\n')
    
    f.write('Summary Statistics:\n')
    f.write(f'Mean tumor fraction: {tumor_fraction["Tumor Fraction"].mean():.4f}\n')
    f.write(f'Median tumor fraction: {tumor_fraction["Tumor Fraction"].median():.4f}\n')
    f.write(f'Standard deviation: {tumor_fraction["Tumor Fraction"].std():.4f}\n')
    f.write(f'Minimum: {tumor_fraction["Tumor Fraction"].min():.4f}\n')
    f.write(f'Maximum: {tumor_fraction["Tumor Fraction"].max():.4f}\n')

print("Analysis complete. Check results/plots/ for the new plot and statistics.") 