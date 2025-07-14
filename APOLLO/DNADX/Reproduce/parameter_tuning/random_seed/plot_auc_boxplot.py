import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the CSV file
df = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/results/validation_summary/auc_comparison.csv')

# Create a figure with a larger size
plt.figure(figsize=(12, 6))

# Create box plot with blank filling
sns.boxplot(data=df[['train_auc', 'test_auc', 'validation_auc']], 
            palette=['white', 'white', 'white'],
            width=0.5,
            showfliers=False)  # Hide default outliers

# Add individual points in dark grey
sns.stripplot(data=df[['train_auc', 'test_auc', 'validation_auc']], 
              color='darkgrey',
              size=4,
              alpha=0.3)

# Add horizontal line at y=0.75
plt.axhline(y=0.75, color='gray', linestyle='--', alpha=0.5)

# Highlight specific signatures
highlight_signatures = {
    'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450': 'red',
    'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP': 'blue',
    'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN': 'green'
}

# Create legend handles
legend_handles = []
legend_labels = []

# Add highlighted points for each signature
for signature, color in highlight_signatures.items():
    # Find the row where the signature matches
    mask = df.apply(lambda x: x.astype(str).str.contains(signature, na=False)).any(axis=1)
    if any(mask):
        for i, col in enumerate(['train_auc', 'test_auc', 'validation_auc']):
            value = df.loc[mask, col].values[0]
            point = plt.plot(i, value, 'o', color=color, markersize=8, alpha=0.7, markeredgecolor='black')[0]
            if i == 0:  # Only add to legend once per signature
                legend_handles.append(point)
                # Create a shorter label for the legend
                if 'UNC_RB_LOH' in signature:
                    label = 'UNC_RB_LOH'
                elif 'Basal_signaling' in signature:
                    label = 'Basal signaling'
                else:
                    label = 'Estrogen signaling'
                legend_labels.append(label)

# Customize the plot
plt.title('AUC Comparison Across Different Sets', fontsize=14, pad=20)
plt.ylabel('AUC Score', fontsize=12)
plt.xticks([0, 1, 2], ['TCGA training', 'TCGA testing', 'METABRIC validation'], rotation=0)
plt.grid(True, linestyle='--', alpha=0.7)

# Add legend
plt.legend(legend_handles, legend_labels, loc='upper right', bbox_to_anchor=(1, 1))

# Set y-axis limits to ensure visibility of all points
plt.ylim(0.4, 1.0)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('auc_comparison_boxplot.png', dpi=300, bbox_inches='tight')
plt.close() 