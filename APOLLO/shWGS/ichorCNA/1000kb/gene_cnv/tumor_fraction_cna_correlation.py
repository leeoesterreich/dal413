import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

def read_tumor_fraction_data(file_path):
    """Read tumor fraction data from summary.csv file."""
    df = pd.read_csv(file_path)
    # Remove .bam extension from sample names
    df['Sample'] = df['Sample'].str.replace('.bam', '', regex=False)
    return df[['Sample', 'Tumor Fraction']]

def count_cna_genes(cnr_file):
    """Count the number of gained and deleted genes in a CNVkit annotated.cnr file."""
    try:
        df = pd.read_csv(cnr_file, sep='\t')
        call_column = [col for col in df.columns if col.endswith('.Corrected_Call')][0]
        genes = df[df['gene'] != '-']['gene'].str.split(',').explode()
        gain_genes = set(genes[df[call_column].isin(['GAIN', 'HLAMP'])].dropna())
        del_genes = set(genes[df[call_column].isin(['HETD', 'HOMD'])].dropna())
        return len(gain_genes), len(del_genes)
    except Exception as e:
        print(f"Error processing {cnr_file}: {str(e)}")
        return 0, 0

def main():
    # Read tumor fraction data
    tf_data = read_tumor_fraction_data('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv')
    print("\nTumor fraction data:")
    print(tf_data.head())
    print(f"Number of samples in tf_data: {len(tf_data)}")

    # Process CNVkit annotated.cnr files
    cnr_dir = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/CNVkit'
    cna_data = []
    for row in tf_data.itertuples():
        sample_name = row.Sample
        cnr_file = glob(os.path.join(cnr_dir, f"{sample_name}.annotated.cnr"))
        if cnr_file:
            print(f"\nProcessing {cnr_file[0]}")
            gain_count, del_count = count_cna_genes(cnr_file[0])
            print(f"Found {gain_count} gained genes and {del_count} deleted genes")
            cna_data.append({'Sample': sample_name, 'Gain_Genes': gain_count, 'Del_Genes': del_count})
        else:
            print(f"\nNo CNR file found for sample {sample_name}")
            cna_data.append({'Sample': sample_name, 'Gain_Genes': np.nan, 'Del_Genes': np.nan})
    
    cna_df = pd.DataFrame(cna_data)
    print("\nCNA data:")
    print(cna_df.head())
    print(f"Number of samples in cna_df: {len(cna_df)}")

    # Merge
    merged_df = pd.merge(tf_data, cna_df, on='Sample', how='inner')
    
    # Filter out samples with more than 15000 gained genes
    filtered_df = merged_df[merged_df['Gain_Genes'] <= 15000].copy()
    print(f"\nRemoved {len(merged_df) - len(filtered_df)} samples with >15000 gained genes")
    print("Samples removed:")
    print(merged_df[merged_df['Gain_Genes'] > 15000][['Sample', 'Tumor Fraction', 'Gain_Genes']])
    
    # Calculate Spearman correlations
    gain_corr, gain_pval = stats.spearmanr(filtered_df['Tumor Fraction'], filtered_df['Gain_Genes'], nan_policy='omit')
    del_corr, del_pval = stats.spearmanr(filtered_df['Tumor Fraction'], filtered_df['Del_Genes'], nan_policy='omit')
    
    # Plot
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    sns.regplot(data=filtered_df, x='Tumor Fraction', y='Gain_Genes', scatter_kws={'alpha':0.5}, line_kws={'color': 'blue'}, color='blue')
    plt.title(f'Tumor Fraction vs Gained Genes\nSpearman r={gain_corr:.3f}, p={gain_pval:.3e}')
    plt.xlabel('Tumor Fraction')
    plt.ylabel('Number of Gained Genes')
    plt.subplot(1, 2, 2)
    sns.regplot(data=filtered_df, x='Tumor Fraction', y='Del_Genes', scatter_kws={'alpha':0.5}, line_kws={'color': 'blue'}, color='blue')
    plt.title(f'Tumor Fraction vs Deleted Genes\nSpearman r={del_corr:.3f}, p={del_pval:.3e}')
    plt.xlabel('Tumor Fraction')
    plt.ylabel('Number of Deleted Genes')
    plt.tight_layout()
    plt.savefig('tumor_fraction_cna_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save merged data
    filtered_df.to_csv('tumor_fraction_cna_data.csv', index=False)
    # Print results
    print("\nSpearman Correlation Results:")
    print(f"Gains vs Tumor Fraction: r={gain_corr:.3f}, p={gain_pval:.3e}")
    print(f"Deletions vs Tumor Fraction: r={del_corr:.3f}, p={del_pval:.3e}")

if __name__ == "__main__":
    main() 