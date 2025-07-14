import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set the style
plt.style.use('seaborn')

# Read the CSV files
df_rs100 = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/random_seed/comparison/auc_comparison_rs100.csv')
df_rs42 = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/random_seed/comparison/auc_comparison_rs_42.csv')
df_rs10 = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/random_seed/comparison/auc_comparison_rs_10.csv')

# Filter out rows where all three AUC columns are NaN in df_rs42
df_rs42 = df_rs42.dropna(subset=['train_auc', 'test_auc', 'validation_auc'], how='all')

# Create a figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# Function to create boxplot for each dataset
def create_boxplot(ax, df, title):
    # Create box plot with blank filling
    box = ax.boxplot([df['train_auc'], df['test_auc'], df['validation_auc']], 
                     patch_artist=True,
                     widths=0.5,
                     showfliers=False)
    
    # Set box colors to white
    for patch in box['boxes']:
        patch.set_facecolor('white')

    # Add individual points
    for i, col in enumerate(['train_auc', 'test_auc', 'validation_auc']):
        x = np.random.normal(i+1, 0.04, size=len(df))
        ax.scatter(x, df[col], color='darkgrey', alpha=0.3, s=20)

    # Add horizontal line at y=0.75
    ax.axhline(y=0.75, color='gray', linestyle='--', alpha=0.5)

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
                point = ax.plot(i+1, value, 'o', color=color, markersize=8, alpha=0.7, markeredgecolor='black')[0]
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
    ax.set_title(title, fontsize=12, pad=20)
    ax.set_ylabel('AUC Score', fontsize=10)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['TCGA training', 'TCGA testing', 'METABRIC validation'], rotation=0)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_ylim(0.4, 1.0)

    # Add legend only to the first subplot
    if ax == ax1:
        ax.legend(legend_handles, legend_labels, loc='upper right', bbox_to_anchor=(1, 1))

# Create boxplots for each dataset
create_boxplot(ax1, df_rs100, 'Random Seed 100')
create_boxplot(ax2, df_rs42, 'Random Seed 42')
create_boxplot(ax3, df_rs10, 'Random Seed 10')

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('auc_comparison_boxplot_all_seeds.png', dpi=300, bbox_inches='tight')
plt.close() 