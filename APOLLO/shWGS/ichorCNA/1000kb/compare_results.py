import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Read the summary files
df_1000kb = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv')
df_500kb = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/500kb/results/summary.csv')

# Reset index to avoid duplicate index issues
df_1000kb = df_1000kb.reset_index(drop=True)
df_500kb = df_500kb.reset_index(drop=True)

# Add a column to identify the source
df_1000kb['Resolution'] = '1000kb'
df_500kb['Resolution'] = '500kb'

# Combine the dataframes
df_combined = pd.concat([df_1000kb, df_500kb], ignore_index=True)

# Create a figure with subplots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Comparison of Results between 1000kb and 500kb Analyses', fontsize=16)

# List of metrics to plot
metrics = [
    'Tumor Fraction',
    'Ploidy',
    'Subclone Fraction',
    'Fraction Genome Subclonal',
    'Fraction CNA Subclonal'
]

# Create box plots for each metric
for i, metric in enumerate(metrics):
    row = i // 3
    col = i % 3
    
    # Create box plot
    sns.boxplot(data=df_combined, x='Resolution', y=metric, ax=axes[row, col])
    axes[row, col].set_title(f'{metric} Comparison')
    axes[row, col].set_ylabel(metric)
    
    # Add individual points
    sns.stripplot(data=df_combined, x='Resolution', y=metric, 
                 color='black', alpha=0.3, ax=axes[row, col], jitter=True)

# Remove the empty subplot
axes[1, 2].remove()

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(top=0.9)

# Save the plot
plt.savefig('resolution_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSummary Statistics:")
for metric in metrics:
    print(f"\n{metric}:")
    print(df_combined.groupby('Resolution')[metric].describe()) 