#!/usr/bin/env python3

import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# Read both summary files
df1 = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb_Broad/results/summary.csv')
df2 = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv')

# Clean up sample names if needed
df1['Sample'] = df1['Sample'].str.replace('.recal.bam', '')
df2['Sample'] = df2['Sample'].str.replace('.recal.bam', '')

# Merge the dataframes on Sample
merged = pd.merge(df1, df2, on='Sample', suffixes=('_broad', '_orig'))

print(f'Number of paired samples: {len(merged)}')

# Function to perform comparison for a parameter
def compare_parameter(param_name):
    param_broad = pd.to_numeric(merged[f'{param_name}_broad'], errors='coerce')
    param_orig = pd.to_numeric(merged[f'{param_name}_orig'], errors='coerce')
    
    t_stat, p_val = stats.ttest_rel(param_broad, param_orig)
    
    print(f'\n{param_name} Comparison:')
    print(f'Mean {param_name} (Broad): {param_broad.mean():.3f}')
    print(f'Mean {param_name} (Original): {param_orig.mean():.3f}')
    print(f'Std {param_name} (Broad): {param_broad.std():.3f}')
    print(f'Std {param_name} (Original): {param_orig.std():.3f}')
    print(f'Paired t-test p-value: {p_val:.6f}')
    
    if param_name == 'Tumor Fraction':
        print('\nZero Tumor Fraction Counts:')
        print(f'Number of zeros (Broad): {sum(param_broad == 0)}')
        print(f'Number of zeros (Original): {sum(param_orig == 0)}')
        
        # Create scatter plot
        plt.figure(figsize=(10, 10))
        plt.scatter(param_orig, param_broad, alpha=0.5)
        plt.plot([0, 1], [0, 1], 'r--')  # Add diagonal line
        plt.xlabel('Original Tumor Fraction')
        plt.ylabel('Broad Tumor Fraction')
        plt.title('Comparison of Tumor Fraction Estimates')
        plt.savefig('tumor_fraction_comparison.pdf')
        plt.close()
        
        # Create boxplot
        plot_data = pd.DataFrame({
            'Original': param_orig,
            'Broad': param_broad
        })
        plt.figure(figsize=(8, 6))
        sns.boxplot(data=plot_data)
        plt.title('Distribution of Tumor Fraction Estimates')
        plt.savefig('tumor_fraction_boxplot.pdf')
        plt.close()

# Compare all relevant parameters
parameters = ['Tumor Fraction', 'Ploidy', 'Subclone Fraction']
for param in parameters:
    compare_parameter(param)

# Save detailed results to a file
with open('comparison_results.txt', 'w') as f:
    f.write('Sample-by-sample comparison:\n\n')
    for idx, row in merged.iterrows():
        f.write(f"Sample: {row['Sample']}\n")
        f.write(f"Tumor Fraction (Broad): {row['Tumor Fraction_broad']}, (Original): {row['Tumor Fraction_orig']}\n")
        f.write(f"Ploidy (Broad): {row['Ploidy_broad']}, (Original): {row['Ploidy_orig']}\n")
        f.write(f"Subclone Fraction (Broad): {row['Subclone Fraction_broad']}, (Original): {row['Subclone Fraction_orig']}\n")
        f.write('-' * 50 + '\n') 