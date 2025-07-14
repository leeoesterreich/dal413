# analyze_patient_cnv.py

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import re
from collections import Counter
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats
import itertools
import seaborn as sns
import matplotlib.gridspec as gridspec

# Create output directory
output_dir = "."
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- Data Loading and Preprocessing ---

# Load and process .cnr files
folder_path = "../../CNVkit"  # Adjusted for script location
sample_dataframes = {}
for file_name in os.listdir(folder_path):
    if file_name.endswith(".cnr"):
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep="\t")
        sample_id_cols = [col for col in df.columns if col.startswith("TP") and "Corrected_Copy_Number" in col]
        if sample_id_cols:
            sample_id = sample_id_cols[0].split(".")[0]
            df = df[["gene", sample_id + ".Corrected_Copy_Number"]]
            df.columns = ["gene_symbol", "copy_number"]
            df = df.assign(gene_symbol=df["gene_symbol"].str.split(",")).explode("gene_symbol").reset_index(drop=True)
            df = df[df["gene_symbol"] != "-"]
            df["copy_number"] = pd.to_numeric(df["copy_number"], errors="coerce")
            df = df.groupby("gene_symbol", as_index=False)["copy_number"].max()
            sample_dataframes[sample_id] = df

print(f"Processed {len(sample_dataframes)} samples.")

# Reformat sample IDs to remove suffixes
reformatted_sample_dataframes = {key.split("_")[0]: value for key, value in sample_dataframes.items()}
print(f"Reformatted to {len(reformatted_sample_dataframes)} unique sample IDs.")

# Load tumor fraction Excel file and build patient_dict
file_path = '../../results/TFx_hg38_1000kb.xlsx'
df_tf = pd.read_excel(file_path, sheet_name=0, skiprows=1)
df_tf.columns = df_tf.columns.map(str)
df_tf = df_tf.loc[:, ~df_tf.columns.str.startswith('Tumor_fraction_1000kb')]
df_tf.rename(columns=lambda x: 'Tumor_fraction' if x.startswith('Tumor_fraction_500kb') else x, inplace=True)
patient_dict = {}
for _, row in df_tf.iterrows():
    patient_code = row["Patient_Code"]
    tuples = []
    for i in range(1, len(df_tf.columns) // 2 + 1):
        if i == 1:
            sample_label_col = f"{i}st progression"
        elif i == 2:
            sample_label_col = f"{i}nd progression"
        elif i == 3:
            sample_label_col = f"{i}rd progression"
        else:
            sample_label_col = f"{i}th progression"
        
        if sample_label_col not in df_tf.columns:
            continue
            
        sample_label_col_index = df_tf.columns.get_loc(sample_label_col)
        tumor_fraction_col_index = sample_label_col_index + 1
        sample_label = row[sample_label_col]
        tumor_fraction = row.iloc[tumor_fraction_col_index]
        if pd.notna(sample_label) and pd.notna(tumor_fraction):
            tuples.append((sample_label, tumor_fraction))
    patient_dict[patient_code] = tuples

# Reformat patient_dict sample IDs
reformatted_data = {
    key: [(item[0].split(" ")[0].split("_")[0], item[1]) for item in value]
    for key, value in patient_dict.items()
}

# --- Exclude specified patient ---
patient_to_exclude = 4908 # Note: patient_dict uses integer keys
if patient_to_exclude in reformatted_data:
    del reformatted_data[patient_to_exclude]
    print(f"Patient {patient_to_exclude} has been excluded from the analysis.")

print("Data loading and preprocessing complete.")

# --- Phase 1: Patient-Level Collapse ---

CNA_AMP_THRESHOLD = 3
CNA_DEL_THRESHOLD = 2

def is_amplified(copy_number):
    return copy_number >= CNA_AMP_THRESHOLD

def is_deleted(copy_number):
    return copy_number < CNA_DEL_THRESHOLD

patient_amplifications = {}
patient_deletions = {}
patient_to_max_tf = {}

for patient_id, samples in reformatted_data.items():
    patient_amp_genes = set()
    patient_del_genes = set()
    max_tf = 0
    
    valid_samples = [s for s in samples if s[0] in reformatted_sample_dataframes]
    if not valid_samples:
        continue

    for sample_id, tumor_fraction in valid_samples:
        # Update max tumor fraction for the patient
        if tumor_fraction > max_tf:
            max_tf = tumor_fraction
        
        # Get the sample's dataframe
        sample_df = reformatted_sample_dataframes[sample_id]
        
        # Find amplified and deleted genes for the current sample
        amp_genes = set(sample_df[sample_df['copy_number'].apply(is_amplified)]['gene_symbol'])
        del_genes = set(sample_df[sample_df['copy_number'].apply(is_deleted)]['gene_symbol'])
        
        # Add to the patient's overall set of altered genes
        patient_amp_genes.update(amp_genes)
        patient_del_genes.update(del_genes)
        
    patient_amplifications[patient_id] = patient_amp_genes
    patient_deletions[patient_id] = patient_del_genes
    patient_to_max_tf[patient_id] = max_tf

print(f"\nCollapsed data for {len(patient_amplifications)} patients.")
print(f"Example: Patient {list(patient_amplifications.keys())[0]} has {len(patient_amplifications[list(patient_amplifications.keys())[0]])} unique amplified genes.")

# --- Phase 2: Patient-Level Odds Ratio Analysis ---

def fisher_confint(table, alpha=0.05):
    """Calculate confidence intervals for Fisher's exact test."""
    oddsratio, pvalue = fisher_exact(table, alternative='two-sided')
    a, b = table[0]
    c, d = table[1]
    # Add 0.5 to avoid division by zero errors in case of zero counts
    if 0 in [a, b, c, d]:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    oratio = (a * d) / (b * c)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    z = stats.norm.ppf(1 - alpha/2)
    log_or = np.log(oratio)
    ci_low = np.exp(log_or - z * se)
    ci_high = np.exp(log_or + z * se)
    return ci_low, ci_high

gain_geneset = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4"]
loss_geneset = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]

# Segregate patients by median tumor fraction
median_tf = np.median(list(patient_to_max_tf.values()))
high_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf >= median_tf}
low_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf < median_tf}

# --- Amplification Analysis ---
amp_results = []
for gene in gain_geneset:
    # Count patients in high TF group with the amplification
    high_tf_amp = sum(1 for p in high_tf_patients if gene in patient_amplifications.get(p, set()))
    high_tf_no_amp = len(high_tf_patients) - high_tf_amp
    
    # Count patients in low TF group with the amplification
    low_tf_amp = sum(1 for p in low_tf_patients if gene in patient_amplifications.get(p, set()))
    low_tf_no_amp = len(low_tf_patients) - low_tf_amp
    
    table = [[high_tf_amp, high_tf_no_amp], [low_tf_amp, low_tf_no_amp]]
    oddsratio, pvalue = fisher_exact(table)
    ci_low, ci_high = fisher_confint(table)
    
    amp_results.append({
        'gene': gene,
        'oddsratio': oddsratio,
        'pvalue': pvalue,
        'ci_low': ci_low,
        'ci_high': ci_high
    })

amp_results_df = pd.DataFrame(amp_results)
if not amp_results_df.empty:
    reject, pvals_corrected, _, _ = multipletests(amp_results_df['pvalue'], method='fdr_bh')
    amp_results_df['pvalue_adj'] = pvals_corrected
    amp_results_df['reject'] = reject
    amp_results_df.to_csv(os.path.join(output_dir, "patient_level_amp_propensity.csv"), index=False)
    print("\nSaved patient-level amplification propensity results.")

    # Generate Forest Plot for Amplifications
    plt.figure(figsize=(10, 12))
    # Define colors based on significance
    colors = ['red' if reject else 'grey' for reject in amp_results_df['reject']]
    plt.errorbar(np.log2(amp_results_df['oddsratio']), np.arange(len(amp_results_df)), 
                 xerr=[np.log2(amp_results_df['oddsratio']) - np.log2(amp_results_df['ci_low']), 
                       np.log2(amp_results_df['ci_high']) - np.log2(amp_results_df['oddsratio'])],
                 fmt='o', color='black', ecolor='grey', elinewidth=1, capsize=3)
    # Add colored points for significance
    plt.scatter(np.log2(amp_results_df['oddsratio']), np.arange(len(amp_results_df)), 
                color=colors, zorder=10)
    plt.yticks(np.arange(len(amp_results_df)), amp_results_df['gene'])
    plt.axvline(0, linestyle='--', color='black')
    plt.xlabel("Log2(Odds Ratio)")
    plt.title("Patient-Level Propensity for Gene Amplification (High vs. Low TF)")
    plt.grid(axis='x')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "patient_level_amp_forest_plot.png"), dpi=300)
    plt.close()
    print("Generated patient-level amplification forest plot.")

# --- Deletion Analysis ---
del_results = []
for gene in loss_geneset:
    high_tf_del = sum(1 for p in high_tf_patients if gene in patient_deletions.get(p, set()))
    high_tf_no_del = len(high_tf_patients) - high_tf_del
    
    low_tf_del = sum(1 for p in low_tf_patients if gene in patient_deletions.get(p, set()))
    low_tf_no_del = len(low_tf_patients) - low_tf_del
    
    table = [[high_tf_del, high_tf_no_del], [low_tf_del, low_tf_no_del]]
    oddsratio, pvalue = fisher_exact(table)
    ci_low, ci_high = fisher_confint(table)
    
    del_results.append({
        'gene': gene,
        'oddsratio': oddsratio,
        'pvalue': pvalue,
        'ci_low': ci_low,
        'ci_high': ci_high
    })

del_results_df = pd.DataFrame(del_results)
if not del_results_df.empty:
    reject, pvals_corrected, _, _ = multipletests(del_results_df['pvalue'], method='fdr_bh')
    del_results_df['pvalue_adj'] = pvals_corrected
    del_results_df['reject'] = reject
    del_results_df.to_csv(os.path.join(output_dir, "patient_level_del_propensity.csv"), index=False)
    print("\nSaved patient-level deletion propensity results.")

    # Generate Forest Plot for Deletions
    plt.figure(figsize=(10, 8))
    # Define colors based on significance
    colors = ['blue' if reject else 'grey' for reject in del_results_df['reject']]
    plt.errorbar(np.log2(del_results_df['oddsratio']), np.arange(len(del_results_df)), 
                 xerr=[np.log2(del_results_df['oddsratio']) - np.log2(del_results_df['ci_low']), 
                       np.log2(del_results_df['ci_high']) - np.log2(del_results_df['oddsratio'])],
                 fmt='o', color='black', ecolor='grey', elinewidth=1, capsize=3)
    # Add colored points for significance
    plt.scatter(np.log2(del_results_df['oddsratio']), np.arange(len(del_results_df)),
                color=colors, zorder=10)
    plt.yticks(np.arange(len(del_results_df)), del_results_df['gene'])
    plt.axvline(0, linestyle='--', color='black')
    plt.xlabel("Log2(Odds Ratio)")
    plt.title("Patient-Level Propensity for Gene Deletion (High vs. Low TF)")
    plt.grid(axis='x')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "patient_level_del_forest_plot.png"), dpi=300)
    plt.close()
    print("Generated patient-level deletion forest plot.")

# --- Phase 3: Co-occurrence and Mutual Exclusivity Analysis ---
print("\nStarting co-occurrence and mutual exclusivity analysis...")

# Focus on the gain_geneset for this analysis
genes_of_interest = sorted(list(set(gain_geneset)))
co_occurrence_results = []

# Get all unique pairs of genes
gene_pairs = list(itertools.combinations(genes_of_interest, 2))

for gene1, gene2 in gene_pairs:
    # Count how many patients have both, one, or neither gene amplified
    g1_and_g2 = sum(1 for pid, amps in patient_amplifications.items() if gene1 in amps and gene2 in amps)
    g1_not_g2 = sum(1 for pid, amps in patient_amplifications.items() if gene1 in amps and gene2 not in amps)
    g2_not_g1 = sum(1 for pid, amps in patient_amplifications.items() if gene1 not in amps and gene2 in amps)
    neither = len(patient_amplifications) - (g1_and_g2 + g1_not_g2 + g2_not_g1)
    
    table = [[g1_and_g2, g1_not_g2], [g2_not_g1, neither]]
    
    # An odds ratio > 1 suggests co-occurrence, < 1 suggests mutual exclusivity
    oddsratio, pvalue = fisher_exact(table, alternative='two-sided')
    
    co_occurrence_results.append({
        'gene1': gene1,
        'gene2': gene2,
        'oddsratio': oddsratio,
        'pvalue': pvalue
    })

if co_occurrence_results:
    co_occurrence_df = pd.DataFrame(co_occurrence_results)
    
    # Adjust p-values for multiple comparisons
    reject, pvals_corrected, _, _ = multipletests(co_occurrence_df['pvalue'], method='fdr_bh')
    co_occurrence_df['pvalue_adj'] = pvals_corrected
    
    # Save the raw results
    co_occurrence_df.sort_values('pvalue').to_csv(os.path.join(output_dir, "co_occurrence_analysis.csv"), index=False)
    print("Saved co-occurrence analysis results.")

    # --- Phase 4: Visualization ---
    # We will visualize the log2 of the odds ratio for significant pairs
    # Create a matrix for the heatmap
    co_occurrence_df['log2_odds'] = np.log2(co_occurrence_df['oddsratio'].replace(0, 1e-10)) # replace 0 to avoid log errors
    
    # Filter for significant results to make the plot readable
    significant_df = co_occurrence_df[co_occurrence_df['pvalue_adj'] < 0.25] # Using a less strict cutoff for visualization
    
    if not significant_df.empty:
        heatmap_data = significant_df.pivot(index='gene1', columns='gene2', values='log2_odds')
        
        # Make the matrix symmetric
        all_genes = sorted(list(set(significant_df['gene1']) | set(significant_df['gene2'])))
        heatmap_data = heatmap_data.reindex(index=all_genes, columns=all_genes)
        
        # Fill the other half of the matrix
        heatmap_data = heatmap_data.add(heatmap_data.T, fill_value=0)

        # Plotting the heatmap
        plt.figure(figsize=(20, 18))
        cmap = sns.diverging_palette(240, 10, as_cmap=True) # Blue for negative (exclusivity), Red for positive (co-occurrence)
        sns.heatmap(heatmap_data, cmap=cmap, center=0, annot=False, fmt=".2f", linewidths=.5)
        plt.title('Co-occurrence and Mutual Exclusivity of Gene Amplifications (Log2 Odds Ratio)', fontsize=16)
        plt.xlabel('')
        plt.ylabel('')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "co_occurrence_heatmap.png"), dpi=300)
        plt.close()
        print("Generated co-occurrence heatmap.")
    else:
        print("No significant co-occurrence or mutual exclusivity found to plot.") 

# --- Phase 5: Patient-Level Oncoplot Visualization ---

def generate_patient_oncoplot(patient_alterations, patient_tf, geneset, title, filename, cmap, freq_bar_color):
    """Generates an oncoplot at the patient level with frequency bars."""
    
    # Create a patient-by-gene matrix
    all_patients = sorted([str(p) for p in patient_alterations.keys()])
    df = pd.DataFrame(0, index=geneset, columns=all_patients)
    for patient, genes in patient_alterations.items():
        altered_in_set = genes.intersection(geneset)
        if altered_in_set:
            df.loc[list(altered_in_set), str(patient)] = 1

    # Filter out genes with no alterations
    df = df.loc[(df.sum(axis=1) > 0)]
    if df.empty:
        print(f"No alterations found for {filename}. Skipping plot.")
        return

    # Sort genes by frequency
    gene_freq = df.sum(axis=1).sort_values(ascending=False)
    df = df.loc[gene_freq.index]

    # Sort patients by number of collections and then by alterations
    patient_sample_counts = {str(p): len(s) for p, s in reformatted_data.items()}
    patient_alteration_counts = df.sum(axis=0)
    sorted_patients = sorted(
        df.columns,
        key=lambda p: (patient_sample_counts.get(p, 0), patient_alteration_counts.get(p, 0)),
        reverse=True
    )
    df = df[sorted_patients]
    
    # Get corresponding TF values and gene frequencies
    patient_tf_str_keys = {str(k): v for k, v in patient_tf.items()}
    tf_values = [patient_tf_str_keys.get(p, 0) for p in sorted_patients]
    gene_freq_percent = (df.sum(axis=1) / df.shape[1]) * 100

    # --- Plotting with GridSpec ---
    fig = plt.figure(figsize=(24, 20))
    gs = gridspec.GridSpec(2, 2, width_ratios=[20, 3], height_ratios=[2, 10],
                           wspace=0.05, hspace=0.05)

    ax_heatmap = fig.add_subplot(gs[1, 0])
    # Restore sharex/sharey for the most robust alignment
    ax_tf = fig.add_subplot(gs[0, 0], sharex=ax_heatmap)
    ax_freq = fig.add_subplot(gs[1, 1], sharey=ax_heatmap)

    # Main heatmap
    sns.heatmap(df, ax=ax_heatmap, cbar=False, cmap=cmap, linewidths=.5)
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=90, fontsize=12)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), fontsize=16, rotation=0)
    ax_heatmap.set_xlabel("Patients", fontsize=14)
    ax_heatmap.set_ylabel("")
    ax_heatmap.tick_params(axis='y', length=0)

    # Top bar plot for tumor fraction
    # Shift bars by 0.5 to align centers with heatmap cells
    ax_tf.bar(np.arange(len(sorted_patients)) + 0.5, tf_values, width=0.8, color='grey', edgecolor='black', linewidth=0.5)
    ax_tf.set_ylabel("Max Tumor Fraction", fontsize=12)
    ax_tf.set_ylim(0, 0.8)
    plt.setp(ax_tf.get_xticklabels(), visible=False)
    ax_tf.tick_params(axis='x', bottom=False, length=0)
    ax_tf.set_title(title, fontsize=22)
    for spine in ['top', 'right', 'left']:
        ax_tf.spines[spine].set_visible(False)

    # Right bar plot for gene frequency
    # Shift bars by 0.5 and set height to 1.0 for perfect alignment
    bars = ax_freq.barh(np.arange(len(gene_freq_percent)) + 0.5, gene_freq_percent.values, height=0.8, color=freq_bar_color)
    ax_freq.set_xlabel("% Patients", fontsize=12)
    ax_freq.set_xlim(0, 110)
    plt.setp(ax_freq.get_yticklabels(), visible=False)
    ax_freq.tick_params(axis='y', left=False, length=0)
    for spine in ['top', 'right', 'bottom']:
        ax_freq.spines[spine].set_visible(False)
    
    # Add text annotations
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax_freq.text(width + 2, bar.get_y() + bar.get_height()/2.,
                         f'{int(round(width))}%', ha='left', va='center', fontsize=11)
    
    # Adjust bottom margin to prevent patient labels from being cut off
    fig.subplots_adjust(bottom=0.15)

    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()
    print(f"Generated enhanced patient-level oncoplot: {filename}")

# Generate plots
print("\nGenerating patient-level oncoplots...")
try:
    generate_patient_oncoplot(patient_amplifications, patient_to_max_tf, gain_geneset, 
                              "Patient-Level Gene Amplifications (CN >= 3)", 
                              "patient_level_amplification_oncoplot.png", 
                              "Reds", "darkred")
except Exception as e:
    print(f"Could not generate amplification oncoplot. Error: {e}")

try:
    generate_patient_oncoplot(patient_deletions, patient_to_max_tf, loss_geneset,
                              "Patient-Level Gene Deletions (CN < 2)", 
                              "patient_level_deletion_oncoplot.png",
                              "Blues", "darkblue")
except Exception as e:
    print(f"Could not generate deletion oncoplot. Error: {e}") 