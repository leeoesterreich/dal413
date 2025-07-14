import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import sys

def read_cnr_file(file_path):
    """Read and process a CNR file."""
    print("Reading file:", file_path)
    df = pd.read_csv(file_path, sep='\t')
    print("Columns in file:", df.columns.tolist())
    
    # Get the column name that contains the call information
    call_column = [col for col in df.columns if col.endswith('.Corrected_Call')]
    copy_number_column = [col for col in df.columns if col.endswith('.Corrected_Copy_Number')]
    
    if not call_column or not copy_number_column:
        print("Warning: Could not find required columns in", file_path)
        print("Available columns:", df.columns.tolist())
        return None
        
    call_column = call_column[0]
    copy_number_column = copy_number_column[0]
    
    # Split genes into individual rows and filter out '-'
    genes = df[df['gene'] != '-']['gene'].str.split(',').explode()
    copy_numbers = df[df['gene'] != '-'][copy_number_column]
    calls = df[df['gene'] != '-'][call_column]
    
    # Create DataFrame with gene-level information
    gene_df = pd.DataFrame({
        'gene': genes,
        'copy_number': copy_numbers,
        'call': calls
    }).dropna()
    
    print("Processed", len(gene_df), "gene entries")
    return gene_df

def calculate_odds_ratio(gene, ilc_samples, idc_samples, threshold=3, is_deletion=False):
    """Calculate odds ratio for a gene between ILC and IDC samples."""
    # For deletions, we want CN < threshold
    # For amplifications, we want CN >= threshold
    if is_deletion:
        ilc_pos = sum(1 for df in ilc_samples if gene in set(df[df['copy_number'] < threshold]['gene']))
        idc_pos = sum(1 for df in idc_samples if gene in set(df[df['copy_number'] < threshold]['gene']))
    else:
        ilc_pos = sum(1 for df in ilc_samples if gene in set(df[df['copy_number'] >= threshold]['gene']))
        idc_pos = sum(1 for df in idc_samples if gene in set(df[df['copy_number'] >= threshold]['gene']))
    
    ilc_neg = len(ilc_samples) - ilc_pos
    idc_neg = len(idc_samples) - idc_pos
    
    # Add 0.5 to all cells to handle zero counts
    contingency = np.array([[ilc_pos + 0.5, ilc_neg + 0.5],
                          [idc_pos + 0.5, idc_neg + 0.5]])
    
    # Calculate odds ratio (ILC vs IDC) and confidence interval
    odds_ratio = (contingency[0,0] * contingency[1,1]) / (contingency[0,1] * contingency[1,0])
    log_odds = np.log(odds_ratio)
    se = np.sqrt(sum(1/contingency.flatten()))
    ci_lower = np.exp(log_odds - 1.96*se)
    ci_upper = np.exp(log_odds + 1.96*se)
    
    # Fisher's exact test
    _, pvalue = stats.fisher_exact(contingency)
    
    return {
        'gene': gene,
        'odds_ratio': odds_ratio,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'pvalue': pvalue,
        'ilc_count': ilc_pos,
        'ilc_total': len(ilc_samples),
        'idc_count': idc_pos,
        'idc_total': len(idc_samples),
        'ilc_freq': ilc_pos/len(ilc_samples),
        'idc_freq': idc_pos/len(idc_samples)
    }

def plot_forest(results_df, title, output_file, is_deletion=False):
    """Create a forest plot for odds ratios."""
    plt.figure(figsize=(18, len(results_df)*0.6 + 4))  # Further increased figure size
    
    # Set font weight to bold for all text
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    
    # Convert to log scale for plotting
    y_pos = np.arange(len(results_df))
    log_odds = np.log10(results_df['odds_ratio'])
    log_ci_lower = np.log10(results_df['ci_lower'])
    log_ci_upper = np.log10(results_df['ci_upper'])
    
    # Plot points and lines
    highlight_color = 'blue' if is_deletion else 'red'
    for i, (lo, hi, odds, sig) in enumerate(zip(log_ci_lower, log_ci_upper, log_odds, results_df['reject'])):
        color = highlight_color if sig else 'grey'
        plt.plot([lo, hi], [i, i], color=color, linewidth=3, zorder=1)  # Even thicker lines
        plt.scatter(odds, i, color=color, zorder=2, s=200)  # Even larger points
    
    # Add vertical line at x=0 (odds ratio = 1)
    plt.axvline(x=0, color='black', linestyle='--', linewidth=2)
    
    # Customize plot with larger and bold fonts
    plt.yticks(y_pos, results_df['gene'], fontsize=24, fontweight='bold')  # Even larger gene names
    plt.xlabel('log10(Odds Ratio)', fontsize=22, fontweight='bold')  # Larger axis label
    plt.title(title, fontsize=26, pad=30, fontweight='bold')  # Larger title
    plt.grid(True, axis='x', alpha=0.3)  # Light grid
    
    # Increase tick label font size
    plt.xticks(fontsize=20, fontweight='bold')
    
    # Add more whitespace around the plot
    plt.margins(y=0.03)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Reset font weight to default for subsequent plots
    plt.rcParams['font.weight'] = 'normal'
    plt.rcParams['axes.labelweight'] = 'normal'
    plt.rcParams['axes.titleweight'] = 'normal'

def main():
    # Define paths
    ilc_dir = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/CNVkit'
    idc_dir = '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/tendo_ctDNA_paper/hg38/Annotated_cnr_file_500Kb'
    
    print("Checking directories:")
    print("ILC dir exists:", os.path.exists(ilc_dir))
    print("IDC dir exists:", os.path.exists(idc_dir))
    
    # Process ILC samples
    ilc_samples = []
    print("\nProcessing ILC samples from:", ilc_dir)
    for file in os.listdir(ilc_dir):
        if file.endswith('.annotated.cnr'):
            try:
                df = read_cnr_file(os.path.join(ilc_dir, file))
                if df is not None:
                    ilc_samples.append(df)
            except Exception as e:
                print("Error processing {0}: {1}".format(file, str(e)))
    
    # Process IDC samples
    idc_samples = []
    print("\nProcessing IDC samples from:", idc_dir)
    for file in os.listdir(idc_dir):
        if file.endswith('.cnr'):
            try:
                df = read_cnr_file(os.path.join(idc_dir, file))
                if df is not None:
                    idc_samples.append(df)
            except Exception as e:
                print("Error processing {0}: {1}".format(file, str(e)))
    
    print("\nProcessed {0} ILC samples and {1} IDC samples".format(len(ilc_samples), len(idc_samples)))
    
    if len(ilc_samples) == 0 or len(idc_samples) == 0:
        print("Error: No samples were processed successfully")
        sys.exit(1)
        
    # Define gene sets
    gain_geneset = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", 
                      "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", 
                      "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", 
                      "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", 
                      "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", 
                      "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", 
                      "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4",
                      "SDHC", "DDR2", "EGFR", "PTPRT", "NUF2", "TERT", "BRCA1", 
                      "HOXB13", "SMYD3", "CCND2"]
    
    loss_geneset = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", 
                      "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", 
                      "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]
    
    # Calculate odds ratios for amplifications
    amp_results = []
    for gene in gain_geneset:
        result = calculate_odds_ratio(gene, ilc_samples, idc_samples, threshold=3, is_deletion=False)
        amp_results.append(result)
    
    # Calculate odds ratios for deletions
    del_results = []
    for gene in loss_geneset:
        result = calculate_odds_ratio(gene, ilc_samples, idc_samples, threshold=2, is_deletion=True)
        del_results.append(result)
    
    # Convert to DataFrames and add multiple testing correction
    amp_df = pd.DataFrame(amp_results)
    del_df = pd.DataFrame(del_results)
    
    # FDR correction
    amp_df['reject'] = multipletests(amp_df['pvalue'], method='fdr_bh')[0]
    del_df['reject'] = multipletests(del_df['pvalue'], method='fdr_bh')[0]
    
    # Sort by odds ratio
    amp_df = amp_df.sort_values('odds_ratio', ascending=False)
    del_df = del_df.sort_values('odds_ratio', ascending=False)
    
    # Create forest plots
    plot_forest(amp_df, 'Amplification Odds Ratios (ILC vs IDC)', 'ilc_idc_amplification_forest.png', is_deletion=False)
    plot_forest(del_df, 'Deletion Odds Ratios (ILC vs IDC)', 'ilc_idc_deletion_forest.png', is_deletion=True)
    
    # Save results with frequencies
    amp_df.to_csv('ilc_idc_amplification_odds.csv', index=False)
    del_df.to_csv('ilc_idc_deletion_odds.csv', index=False)
    
    # Print PTEN deletion details
    if 'PTEN' in loss_geneset:
        pten_result = del_df[del_df['gene'] == 'PTEN'].iloc[0]
        print("\nPTEN deletion details:")
        print("ILC: {}/{} samples ({:.1f}%)".format(
            pten_result['ilc_count'], pten_result['ilc_total'], 
            100 * pten_result['ilc_freq']))
        print("IDC: {}/{} samples ({:.1f}%)".format(
            pten_result['idc_count'], pten_result['idc_total'], 
            100 * pten_result['idc_freq']))
        print("Odds ratio (ILC vs IDC): {:.2f}".format(pten_result['odds_ratio']))
        print("P-value: {:.2e}".format(pten_result['pvalue']))

if __name__ == "__main__":
    main() 