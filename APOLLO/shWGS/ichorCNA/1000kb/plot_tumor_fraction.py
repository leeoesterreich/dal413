import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import os

# Read the summary data
summary = pd.read_csv('results/summary.csv')

# Create a directory for plots if it doesn't exist
if not os.path.exists('results/plots'):
    os.makedirs('results/plots')

# Create a histogram of tumor fractions
plt.figure(figsize=(10, 6))
sns.histplot(data=summary, x='Tumor Fraction', bins=30)
plt.title('Distribution of Tumor Fractions')
plt.xlabel('Tumor Fraction')
plt.ylabel('Count')
plt.grid(True)
plt.tight_layout()
plt.savefig('results/plots/tumor_fraction_distribution.png')
plt.close()

# Create a box plot of tumor fractions
plt.figure(figsize=(10, 6))
sns.boxplot(data=summary, y='Tumor Fraction')
plt.title('Tumor Fraction Distribution')
plt.ylabel('Tumor Fraction')
plt.grid(True)
plt.tight_layout()
plt.savefig('results/plots/tumor_fraction_boxplot.png')
plt.close()

# Calculate and save statistics
with open('results/plots/tumor_fraction_stats.txt', 'w') as f:
    f.write('Tumor Fraction Statistics:\n\n')
    f.write(f'Number of samples: {len(summary)}\n')
    f.write(f'Mean tumor fraction: {summary["Tumor Fraction"].mean():.4f}\n')
    f.write(f'Median tumor fraction: {summary["Tumor Fraction"].median():.4f}\n')
    f.write(f'Standard deviation: {summary["Tumor Fraction"].std():.4f}\n')
    f.write(f'Minimum: {summary["Tumor Fraction"].min():.4f}\n')
    f.write(f'Maximum: {summary["Tumor Fraction"].max():.4f}\n')

# Create a scatter plot of tumor fraction vs ploidy if ploidy exists
if 'Ploidy' in summary.columns:
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=summary, x='Ploidy', y='Tumor Fraction')
    plt.title('Tumor Fraction vs Ploidy')
    plt.xlabel('Ploidy')
    plt.ylabel('Tumor Fraction')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('results/plots/tumor_fraction_vs_ploidy.png')
    plt.close()

print("Analysis complete. Check the results/plots directory for output files.") 