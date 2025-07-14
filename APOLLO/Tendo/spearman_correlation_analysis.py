import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import os
from glob import glob

def read_tumor_fraction_data(file_path):
    """Read tumor fraction data from CSV file."""
    df = pd.read_csv(file_path, skiprows=1)  # Skip the first row (Reference genome info)
    
    # Create a list to store processed data
    processed_data = []
    
    # Process each row
    for _, row in df.iterrows():
        patient = row['Patient_Code']
        
        # Process each progression
        for i in range(1, 7):  # 6 possible progressions
            sample_col = f'{i}st progression' if i == 1 else \
                        f'{i}nd progression' if i == 2 else \
                        f'{i}rd progression' if i == 3 else \
                        f'{i}th progression'
            tf_col = 'Tumor_fraction_500kb bin'  # Using 500kb bin data
            
            if pd.notna(row[sample_col]):
                sample_name = row[sample_col].split()[0]  # Get sample name without the number in parentheses
                tf_value = row[tf_col]
                
                if pd.notna(tf_value):
                    processed_data.append({
                        'Patient': patient,
                        'Sample': sample_name,
                        'Tumor_Fraction': float(tf_value)
                    })
    
    return pd.DataFrame(processed_data)

def count_cna_genes(cnr_file):
    """Count the number of gained and deleted genes in a CNR file."""
    try:
        df = pd.read_csv(cnr_file, sep='\t')
        
        # Get the column name that contains the call information (ends with .Corrected_Call)
        call_column = [col for col in df.columns if col.endswith('.Corrected_Call')][0]
        
        # Count unique genes for gains and deletions
        # Split the gene column by comma and create a list of all genes
        genes = df[df['gene'] != '-']['gene'].str.split(',').explode()
        
        # Count genes in gain regions
        gain_genes = set(genes[df[call_column].isin(['GAIN', 'HLAMP'])].dropna())
        
        # Count genes in deletion regions
        del_genes = set(genes[df[call_column].isin(['HETD', 'HOMD'])].dropna())
        
        return len(gain_genes), len(del_genes)
    except Exception as e:
        print(f"Error processing {cnr_file}: {str(e)}")
        return 0, 0

def main():
    # Read tumor fraction data
    tf_data = read_tumor_fraction_data('/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Clinical_information/TFx_hg38_1000_500kb.csv')
    
    # Process CNR files
    cnr_dir = '/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Tendo/Annotated_cnr_file_500Kb'
    cna_data = []
    
    for sample_data in tf_data.itertuples():
        sample_name = sample_data.Sample
        cnr_file = glob(os.path.join(cnr_dir, f"{sample_name}.annotated.cnr"))
        
        if cnr_file:
            gain_count, del_count = count_cna_genes(cnr_file[0])
            cna_data.append({
                'Sample': sample_name,
                'Gain_Genes': gain_count,
                'Del_Genes': del_count
            })
    
    # Create CNA dataframe
    cna_df = pd.DataFrame(cna_data)
    
    # Merge tumor fraction and CNA data
    merged_df = pd.merge(tf_data, cna_df, on='Sample', how='inner')
    
    # Calculate Spearman correlations
    gain_corr, gain_pval = stats.spearmanr(merged_df['Tumor_Fraction'], merged_df['Gain_Genes'])
    del_corr, del_pval = stats.spearmanr(merged_df['Tumor_Fraction'], merged_df['Del_Genes'])
    
    # Set style for better visualization
    sns.set_style("whitegrid")
    
    # Create scatter plots with regression lines
    plt.figure(figsize=(12, 5))
    
    # Plot for gains
    plt.subplot(1, 2, 1)
    sns.regplot(data=merged_df, x='Tumor_Fraction', y='Gain_Genes',
                scatter_kws={'alpha':0.5}, 
                line_kws={'color': 'blue'},
                color='blue')
    plt.title(f'Tumor Fraction vs Gained Genes\nSpearman r={gain_corr:.3f}, p={gain_pval:.3e}')
    plt.xlabel('Tumor Fraction')
    plt.ylabel('Number of Gained Genes')
    
    # Plot for deletions
    plt.subplot(1, 2, 2)
    sns.regplot(data=merged_df, x='Tumor_Fraction', y='Del_Genes',
                scatter_kws={'alpha':0.5}, 
                line_kws={'color': 'blue'},
                color='blue')
    plt.title(f'Tumor Fraction vs Deleted Genes\nSpearman r={del_corr:.3f}, p={del_pval:.3e}')
    plt.xlabel('Tumor Fraction')
    plt.ylabel('Number of Deleted Genes')
    
    plt.tight_layout()
    plt.savefig('tumor_fraction_cna_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print results
    print("\nSpearman Correlation Results:")
    print(f"Gains vs Tumor Fraction: r={gain_corr:.3f}, p={gain_pval:.3e}")
    print(f"Deletions vs Tumor Fraction: r={del_corr:.3f}, p={del_pval:.3e}")
    
    # Save detailed results to CSV
    merged_df.to_csv('tumor_fraction_cna_data.csv', index=False)

if __name__ == "__main__":
    main() 