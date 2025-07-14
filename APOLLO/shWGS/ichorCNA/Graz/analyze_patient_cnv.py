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
output_dir = "./cnv_analysis_results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- Function Definitions ---
def is_amplified(copy_number):
    return copy_number >= 3

def is_deleted(copy_number):
    return copy_number < 2

def generate_patient_oncoplot(patient_alterations, patient_tf, patient_dict, geneset, title, filename, cmap, freq_bar_color):
    """
    Generates and saves a patient-level oncoplot.
    """
    print(f"Generating oncoplot: {title}")

    # Create a binary matrix for the oncoplot
    onco_matrix = pd.DataFrame(index=geneset, columns=list(patient_tf.keys()))
    for gene in geneset:
        for patient in onco_matrix.columns:
            onco_matrix.loc[gene, patient] = 1 if gene in patient_alterations.get(patient, set()) else 0
    onco_matrix = onco_matrix.astype(int)

    # Calculate alteration frequency for sorting
    alteration_freq = onco_matrix.sum(axis=1) / len(onco_matrix.columns)
    onco_matrix = onco_matrix.loc[alteration_freq.sort_values(ascending=False).index]
    
    # Sort patients by max tumor fraction
    sorted_patients = sorted(patient_tf.keys(), key=lambda p: patient_tf[p], reverse=True)
    onco_matrix = onco_matrix[sorted_patients]

    # Create figure with GridSpec for layout
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 8], width_ratios=[18, 2], 
                           wspace=0.05, hspace=0.05)

    ax_onco = plt.subplot(gs[1, 0])
    ax_bar = plt.subplot(gs[1, 1], sharey=ax_onco)
    ax_tf = plt.subplot(gs[0, 0], sharex=ax_onco)

    # Oncoprint
    sns.heatmap(onco_matrix, ax=ax_onco, cbar=False, cmap=cmap, linewidths=.5, linecolor='lightgrey')
    ax_onco.set_yticks(np.arange(len(onco_matrix.index)) + 0.5)
    ax_onco.set_yticklabels(onco_matrix.index, rotation=0)
    ax_onco.tick_params(left=True, bottom=False)
    
    # Hide patient labels on the main heatmap
    ax_onco.set_xticklabels([])

    # Bar plot for alteration frequency
    freq = (onco_matrix.sum(axis=1) / len(onco_matrix.columns)) * 100
    bars = ax_bar.barh(np.arange(len(freq)) + 0.5, freq, color=freq_bar_color, height=0.7)
    ax_bar.set_xlabel('% Altered')
    ax_bar.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.spines['left'].set_visible(False)
    ax_bar.set_xlim(0, 110) # Make space for labels

    # Add percentage annotation to bars
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax_bar.text(width + 2, bar.get_y() + bar.get_height() / 2,
                        f'{width:.0f}%', ha='left', va='center', fontsize=8)

    # Top plot for tumor fraction
    tf_values = [patient_tf[p] for p in onco_matrix.columns]
    ax_tf.bar(np.arange(len(tf_values)) + 0.5, tf_values, color='grey', width=1.0)
    ax_tf.set_ylabel("log10(TFx z-score)")
    ax_tf.spines['top'].set_visible(False)
    ax_tf.spines['right'].set_visible(False)
    ax_tf.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    # Add dividing line for TF cutoff
    TF_CUTOFF = np.log10(3)
    ax_tf.axhline(y=TF_CUTOFF, color='red', linestyle='--', linewidth=1)

    fig.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make room for suptitle
    
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved oncoplot to {output_path}")


# --- Data Loading and Preprocessing ---

# Load clinical data to filter for NST patients
clinical_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/CtDNAHRHER2Novartis-ClinicalNina_DATA_LABELS_2025-05-14_1651.csv'
df_clinical = pd.read_csv(clinical_file_path)
nst_patients = set(df_clinical[df_clinical['Histology'] == 'NST']['HG-ID'])
print(f"Loaded {len(nst_patients)} patients with NST histology.")

# Load and process .cnr files from the reformatted directory
folder_path = "cna_seg_files_reformatted"
sample_dataframes = {}

for file_name in os.listdir(folder_path):
    if file_name.endswith(".cnr"):
        long_sample_id = file_name.replace(".cnr", "")
        patient_id_from_file = long_sample_id.split('_')[0]

        if patient_id_from_file in nst_patients:
            file_path = os.path.join(folder_path, file_name)
            try:
                df = pd.read_csv(file_path, sep="\t")
                copy_number_col = f"{long_sample_id}.Corrected_Copy_Number"

                if 'gene' in df.columns and copy_number_col in df.columns:
                    df_proc = df[["gene", copy_number_col]]
                    df_proc.columns = ["gene_symbol", "copy_number"]
                    df_proc = df_proc[df_proc['gene_symbol'].notna()]
                    df_proc = df_proc.assign(gene_symbol=df_proc["gene_symbol"].str.split(",")).explode("gene_symbol").reset_index(drop=True)
                    df_proc = df_proc[df_proc["gene_symbol"] != "-"]
                    df_proc["copy_number"] = pd.to_numeric(df_proc["copy_number"], errors="coerce")
                    df_proc.dropna(subset=['copy_number'], inplace=True)
                    df_proc = df_proc.groupby("gene_symbol", as_index=False)["copy_number"].max()
                    sample_dataframes[long_sample_id] = df_proc
                else:
                    print(f"Could not find required columns in file {file_name}")
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")

print(f"Processed {len(sample_dataframes)} total samples from NST patients.")

# Load tumor fraction Excel file and build patient_dict
tf_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/Novartis TFx.csv'
df_tf = pd.read_csv(tf_file_path)
df_tf.columns = df_tf.columns.str.replace('\\n', ' ', regex=True)
tf_col_name = 'GW z-score (mFAST-SeqS)'
df_tf.rename(columns={col: tf_col_name for col in df_tf.columns if 'GW' in col and 'z-score' in col}, inplace=True)

# Use log10 of the z-score
# Filter for positive z-scores to avoid log10(x<=0) errors and apply transformation
df_tf = df_tf[df_tf[tf_col_name] > 0].copy()
df_tf[tf_col_name] = np.log10(df_tf[tf_col_name])

patient_dict = {}
for _, row in df_tf.iterrows():
    db_id = row["DB-ID"]
    patient_id = db_id.split('_')[0]
    
    if patient_id in nst_patients:
        if tf_col_name in df_tf.columns and pd.notna(row[tf_col_name]):
            tumor_fraction = row[tf_col_name]
            if patient_id not in patient_dict:
                patient_dict[patient_id] = []
            patient_dict[patient_id].append((db_id, tumor_fraction))

print(f"Built patient dictionary for {len(patient_dict)} NST patients from tumor fraction file.")

print("Data loading and preprocessing complete.")

# --- Phase 1: Patient-Level Collapse ---

CNA_AMP_THRESHOLD = 3
CNA_DEL_THRESHOLD = 2

patient_amplifications = {}
patient_deletions = {}
patient_to_max_tf = {}

# First, populate patient_to_max_tf from the Novartis TFx file
for patient_id in nst_patients:
    patient_tfx = df_tf[df_tf['DB-ID'].str.startswith(f"{patient_id}_")]
    if not patient_tfx.empty:
        patient_to_max_tf[patient_id] = patient_tfx[tf_col_name].max()

# Now, collapse CNV data for each patient
for patient_id in nst_patients:
    patient_amp_genes = set()
    patient_del_genes = set()

    # Find all samples belonging to this patient
    patient_sample_ids = [sid for sid in sample_dataframes.keys() if sid.startswith(f"{patient_id}_")]

    if not patient_sample_ids:
        continue # Skip patients with no CNV data

    for sample_id in patient_sample_ids:
        sample_df = sample_dataframes[sample_id]
        
        amp_genes = set(sample_df[sample_df['copy_number'] >= CNA_AMP_THRESHOLD]['gene_symbol'])
        del_genes = set(sample_df[sample_df['copy_number'] < CNA_DEL_THRESHOLD]['gene_symbol'])
        
        patient_amp_genes.update(amp_genes)
        patient_del_genes.update(del_genes)
        
    patient_amplifications[patient_id] = patient_amp_genes
    patient_deletions[patient_id] = patient_del_genes

# Remove patients from patient_to_max_tf if they had no CNV data to be collapsed
valid_patients = set(patient_amplifications.keys())
patient_to_max_tf = {p: tf for p, tf in patient_to_max_tf.items() if p in valid_patients}

print(f"\nCollapsed data for {len(patient_amplifications)} patients.")
if patient_amplifications:
    first_patient_key = list(patient_amplifications.keys())[0]
    print(f"Example: Patient {first_patient_key} has {len(patient_amplifications[first_patient_key])} unique amplified genes.")

# --- Phase 2: Patient-Level Odds Ratio Analysis ---

gain_geneset = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4"]
loss_geneset = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]

# Segregate patients by tumor fraction cutoff of log10(3)
TF_CUTOFF = np.log10(3)
high_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf >= TF_CUTOFF}
low_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf < TF_CUTOFF}

print(f"\nSegregated patients into High TF (>= log10(3)) and Low TF (< log10(3)) groups.")
print(f"High TF patients: {len(high_tf_patients)}")
print(f"Low TF patients: {len(low_tf_patients)}")

# --- Generate Oncoplots First ---
print("\nGenerating oncoplots...")

# Oncoplot for Amplifications
generate_patient_oncoplot(
    patient_alterations=patient_amplifications,
    patient_tf=patient_to_max_tf,
    patient_dict=patient_dict,
    geneset=gain_geneset,
    title="Frequency of Gene Amplification (CN >= 3) in NST Patients (Collapsed)",
    filename="patient_level_gains_oncoplot.png",
    cmap=plt.cm.get_cmap("Reds"),
    freq_bar_color='lightcoral'
)

# Oncoplot for Deletions
generate_patient_oncoplot(
    patient_alterations=patient_deletions,
    patient_tf=patient_to_max_tf,
    patient_dict=patient_dict,
    geneset=loss_geneset,
    title="Frequency of Gene Deletion (CN < 2) in NST Patients (Collapsed)",
    filename="patient_level_losses_oncoplot.png",
    cmap=plt.cm.get_cmap("Blues"),
    freq_bar_color='lightblue'
)

print("Oncoplot generation complete.")

# --- Amplification Analysis ---
amp_results = []
if high_tf_patients and low_tf_patients:
    for gene in gain_geneset:
        high_tf_amp = sum(1 for p in high_tf_patients if gene in patient_amplifications.get(p, set()))
        high_tf_no_amp = len(high_tf_patients) - high_tf_amp
        
        low_tf_amp = sum(1 for p in low_tf_patients if gene in patient_amplifications.get(p, set()))
        low_tf_no_amp = len(low_tf_patients) - low_tf_amp
        
        table = [[int(high_tf_amp), int(high_tf_no_amp)], [int(low_tf_amp), int(low_tf_no_amp)]]
        
        if table[0][0] + table[0][1] == 0 or table[1][0] + table[1][1] == 0:
            continue

        try:
            oddsratio, pvalue = fisher_exact(table)
            
            # Calculate confidence interval for odds ratio
            se = np.sqrt(1/(table[0][0]+0.5) + 1/(table[0][1]+0.5) + 1/(table[1][0]+0.5) + 1/(table[1][1]+0.5))
            ci_low = np.exp(np.log(oddsratio) - 1.96 * se)
            ci_high = np.exp(np.log(oddsratio) + 1.96 * se)
            
            amp_results.append({
                'gene': gene,
                'oddsratio': oddsratio,
                'pvalue': pvalue,
                'ci_low': ci_low,
                'ci_high': ci_high
            })
        except ValueError:
            # This can happen if a row or column sums to zero, even with the check above for some edge cases.
            print(f"Skipping Fisher's exact test for gene {gene} due to zero counts in a row/column.")
            continue

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
else:
    print("Skipping amplification analysis and plotting due to no data or invalid group separation.")

# --- Deletion Analysis ---
del_results = []
if high_tf_patients and low_tf_patients:
    for gene in loss_geneset:
        high_tf_del = sum(1 for p in high_tf_patients if gene in patient_deletions.get(p, set()))
        high_tf_no_del = len(high_tf_patients) - high_tf_del
        
        low_tf_del = sum(1 for p in low_tf_patients if gene in patient_deletions.get(p, set()))
        low_tf_no_del = len(low_tf_patients) - low_tf_del
        
        table = [[int(high_tf_del), int(high_tf_no_del)], [int(low_tf_del), int(low_tf_no_del)]]

        if table[0][0] + table[0][1] == 0 or table[1][0] + table[1][1] == 0:
            continue

        try:
            oddsratio, pvalue = fisher_exact(table)
            
            # Calculate confidence interval for odds ratio
            se = np.sqrt(1/(table[0][0]+0.5) + 1/(table[0][1]+0.5) + 1/(table[1][0]+0.5) + 1/(table[1][1]+0.5))
            ci_low = np.exp(np.log(oddsratio) - 1.96 * se)
            ci_high = np.exp(np.log(oddsratio) + 1.96 * se)
            
            del_results.append({
                'gene': gene,
                'oddsratio': oddsratio,
                'pvalue': pvalue,
                'ci_low': ci_low,
                'ci_high': ci_high
            })
        except ValueError:
            print(f"Skipping Fisher's exact test for gene {gene} due to zero counts in a row/column.")
            continue

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
else:
    print("Skipping deletion propensity analysis and plotting due to lack of data.")

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

# --- Phase 4: Oncoplot Visualization ---

print("\nGenerating patient-level oncoplots...")

# --- Generate Deletion Plot First ---
generate_patient_oncoplot(
    patient_deletions, patient_to_max_tf, patient_dict, loss_geneset,
    "Patient-Level Gene Deletions (CN < 2)",
    "patient_level_deletion_oncoplot.png",
    matplotlib.colors.ListedColormap(['white', 'blue']), "darkblue"
)

# --- Generate Amplification Plot Second ---
generate_patient_oncoplot(
    patient_amplifications, patient_to_max_tf, patient_dict, gain_geneset,
    "Patient-Level Gene Amplifications (CN >= 3)",
    "patient_level_amplification_oncoplot.png",
    matplotlib.colors.ListedColormap(['white', 'red']), "darkred"
)

print("\nScript finished.") 