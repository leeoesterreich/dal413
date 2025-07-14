# Oncoplot_propensity.py
# Converted from Oncoplot_propensity.ipynb

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from collections import Counter
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
import scipy.stats as stats

# --- Cell 0: Load and process .cnr files ---
folder_path = "../CNVkit/annotated"  # Updated path to annotated CNR files
sample_dataframes = {}
for file_name in os.listdir(folder_path):
    if file_name.endswith(".annotated.cnr"):  # Updated file extension
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep="\t")
        sample_id = file_name.split(".")[0]  # Get sample ID from filename
        # Extract copy number and gene information
        df = df[["gene", f"{sample_id}.Corrected_Copy_Number"]]
        df.columns = ["gene_symbol", "copy_number"]
        df = df.assign(gene_symbol=df["gene_symbol"].str.split(",")).explode("gene_symbol").reset_index(drop=True)
        df = df[df["gene_symbol"] != "-"]
        df["copy_number"] = pd.to_numeric(df["copy_number"], errors="coerce")
        df = df.groupby("gene_symbol", as_index=False)["copy_number"].max()
        sample_dataframes[sample_id] = df
print("Processed sample IDs:", sample_dataframes.keys())

# --- Cell 1: Reformat sample IDs ---
reformatted_sample_dataframes = {key: value for key, value in sample_dataframes.items()}
print("Reformatted sample IDs:", reformatted_sample_dataframes.keys())

# --- Cell 2: Load tumor fraction Excel file and build patient_dict ---
file_path = '../results/TFx_hg38_1000kb.xlsx'
df = pd.read_excel(file_path, sheet_name=0, skiprows=1)
df.columns = df.columns.map(str)
df = df.loc[:, ~df.columns.str.startswith('Tumor_fraction_1000kb')]
df.rename(columns=lambda x: 'Tumor_fraction' if x.startswith('Tumor_fraction_500kb') else x, inplace=True)
patient_dict = {}
for _, row in df.iterrows():
    patient_code = row["Patient_Code"]
    tuples = []
    for i in range(1, len(df.columns) // 2 + 1):
        if i == 1:
            sample_label_col = f"{i}st progression"
        elif i == 2:
            sample_label_col = f"{i}nd progression"
        elif i == 3:
            sample_label_col = f"{i}rd progression"
        else:
            sample_label_col = f"{i}th progression"
        if sample_label_col not in df.columns:
            continue
        sample_label_col_index = df.columns.get_loc(sample_label_col)
        tumor_fraction_col_index = sample_label_col_index + 1
        sample_label = row[sample_label_col]
        tumor_fraction = row.iloc[tumor_fraction_col_index]
        if pd.notna(sample_label) and pd.notna(tumor_fraction):
            tuples.append((sample_label, tumor_fraction))
    patient_dict[patient_code] = tuples

# --- Cell 3: Reformat patient_dict ---
reformatted_data = {
    key: [(item[0].split(" ")[0], item[1]) for item in value]
    for key, value in patient_dict.items()
}

# --- Cell 4: Top amplified genes (copy number >= 3) ---
gene_counter = Counter()
for sample_id, df in reformatted_sample_dataframes.items():
    amplified_genes = set(df[df['copy_number'] >= 3]['gene_symbol'])
    gene_counter.update(amplified_genes)
top_amplified_genes_list = gene_counter.most_common()
print("Top amplified genes (copy number >= 3):")
for gene, count in top_amplified_genes_list:
    print(f"{gene}: {count}")

# --- Cell 5: Top deleted genes (copy number < 2) ---
gene_counter = Counter()
for sample_id, df in reformatted_sample_dataframes.items():
    deleted_genes = set(df[df['copy_number'] < 2]['gene_symbol'])
    gene_counter.update(deleted_genes)
top_deleted_genes_list = gene_counter.most_common()
print("Top deleted genes (copy number < 2):")
for gene, count in top_deleted_genes_list:
    print(f"{gene}: {count}")

# --- Cell 6: Top amplified genes across patients ---
gene_counter = Counter()
for patient, samples in reformatted_data.items():
    patient_amplified_genes = set()
    for sample, _ in samples:
        if sample in reformatted_sample_dataframes:
            df = reformatted_sample_dataframes[sample]
            amplified_genes = df[df['copy_number'] >= 3]['gene_symbol']
            patient_amplified_genes.update(amplified_genes)
    gene_counter.update(patient_amplified_genes)
top_amplified_genes_list = gene_counter.most_common()
print("Top amplified genes (copy number >= 3 across patients):")
for gene, count in top_amplified_genes_list:
    if gene == "MYC":
        print(f"{gene}: {count}")

# --- Cell 7: Gain geneset ---
gain_geneset = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4"]
print(gain_geneset)

# --- Cell 8: Oncoplot for amplifications ---
def generate_oncoplot_with_sample_and_patient_labels(reformatted_sample_dataframes, gain_geneset, reformatted_data):
    def parse_sample_id(sample_id):
        match = re.match(r"TP(\d+)-M(\d+)", sample_id)
        if match:
            year, sample_num = map(int, match.groups())
            return (year, sample_num)
        return (0, 0)
    def extract_patient_num(patient_id):
        match = re.search(r'\d+', str(patient_id))
        return int(match.group()) if match else float('inf')
    top_genes = gain_geneset
    # Sort patients by number of samples (descending), then by patient number
    patient_sample_counts = {pid: sum(1 for sample, _ in reformatted_data[pid] if sample in reformatted_sample_dataframes) for pid in reformatted_data}
    unique_patients = sorted(
        reformatted_data.keys(),
        key=lambda pid: (-patient_sample_counts[pid], extract_patient_num(pid))
    )
    # Prepare sample order and all data containers
    sample_order = []
    sample_to_patient = {}
    sample_to_fraction = {}
    for patient in unique_patients:
        patient_samples = [
            (sample, fraction) for sample, fraction in reformatted_data[patient]
            if sample in reformatted_sample_dataframes
        ]
        patient_samples.sort(key=lambda x: parse_sample_id(x[0]))
        for sample, fraction in patient_samples:
            sample_order.append(sample)
            sample_to_patient[sample] = patient
            sample_to_fraction[sample] = fraction
    # Build oncoplot_data in this sample order
    oncoplot_data = pd.DataFrame(0, index=top_genes, columns=sample_order)
    for sample in sample_order:
        df = reformatted_sample_dataframes[sample]
        df_amplified = df[(df['copy_number'] >= 3) & (df['gene_symbol'].isin(gain_geneset))]
        amplified_genes = set(df_amplified['gene_symbol'])
        oncoplot_data.loc[list(amplified_genes), sample] = 1
    tumor_fraction = np.array([sample_to_fraction[sample] for sample in sample_order])
    sample_positions = [sample_to_patient[sample] for sample in sample_order]
    sample_labels = sample_order
    gene_frequencies = oncoplot_data.sum(axis=1).sort_values(ascending=False)
    oncoplot_data = oncoplot_data.loc[gene_frequencies.index]
    if oncoplot_data.empty:
        print("⚠️ No genes to plot. Exiting.")
        return
    # Calculate figure size for square boxes
    n_samples = len(sample_order)
    n_genes = len(oncoplot_data.index)
    cell_size = 0.6
    width = max(10, n_samples * cell_size)
    height = max(8, n_genes * cell_size)
    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        sharex=True,
        gridspec_kw={'height_ratios': [0.5, 4]},
        figsize=(width, height)
    )
    x = np.arange(len(sample_order))
    ax1.bar(x, tumor_fraction, color='grey', edgecolor="black", linewidth=0.5)
    ax1.set_ylabel("Tumor Fraction", fontsize=16)
    ax1.set_xlim(-0.5, len(sample_order) - 0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(sample_labels, rotation=90, fontsize=12)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax1.set_yticks([0.2, 0.4, 0.6])
    ax1.tick_params(axis='y', labelsize=12)
    cax = ax2.matshow(oncoplot_data, cmap="Reds", aspect="auto", zorder=1, alpha=0.7)
    ax2.set_xlim(-0.5, len(sample_order) - 0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(sample_labels, rotation=90, fontsize=12)
    ax2.set_yticks(np.arange(len(oncoplot_data.index)))
    ax2.set_yticklabels(oncoplot_data.index, fontsize=24, fontweight='bold')
    ax1.yaxis.grid(True, which='major', linestyle=':', color='gray', alpha=0.7)
    ax2.set_xticks(x - 0.5, minor=True)
    ax2.set_yticks(np.arange(-0.5, len(oncoplot_data.index), 1), minor=True)
    ax2.grid(which="minor", color="gray", linestyle=":", linewidth=0.5, zorder=0)
    patient_ticks = []
    unique_patient_labels = []
    current_patient = sample_positions[0]
    group_start = 0
    for i in range(1, len(sample_positions)):
        if sample_positions[i] != current_patient:
            ax2.axvline(i - 0.5, color="black", linestyle="-", linewidth=1.5, zorder=3)
            midpoint = (group_start + i - 1) / 2
            patient_ticks.append(midpoint)
            unique_patient_labels.append(f"Patient {current_patient}")
            group_start = i
            current_patient = sample_positions[i]
    # Add last patient
    midpoint = (group_start + len(sample_positions) - 1) / 2
    patient_ticks.append(midpoint)
    unique_patient_labels.append(f"Patient {current_patient}")
    ax2.set_xticks(patient_ticks)
    ax2.set_xticklabels(unique_patient_labels, rotation=90, fontsize=12)
    plt.tight_layout()
    plt.savefig("oncoplot_amplifications.pdf", bbox_inches='tight')
    plt.close()

# --- Cell 9: Loss geneset ---
loss_geneset = ["PTEN", "RB1", "CDKN2A", "CDKN2B", "ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", "PBRM1", "BAP1", "SETD2", "KMT2C", "KMT2D", "KDM6A", "CREBBP", "EP300", "ASXL1", "ASXL2", "EZH2", "NF1", "TSC1", "TSC2", "STK11", "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54L", "XRCC2", "XRCC3"]
print(loss_geneset)

# --- Cell 10: Oncoplot for deletions ---
def generate_oncoplot_with_sample_and_patient_labels_loss(reformatted_sample_dataframes, loss_geneset, reformatted_data):
    def parse_sample_id(sample_id):
        match = re.match(r"TP(\d+)-M(\d+)", sample_id)
        if match:
            year, sample_num = map(int, match.groups())
            return (year, sample_num)
        return (0, 0)
    def extract_patient_num(patient_id):
        match = re.search(r'\d+', str(patient_id))
        return int(match.group()) if match else float('inf')
    top_genes = loss_geneset
    # Sort patients by number of samples (descending), then by patient number
    patient_sample_counts = {pid: sum(1 for sample, _ in reformatted_data[pid] if sample in reformatted_sample_dataframes) for pid in reformatted_data}
    unique_patients = sorted(
        reformatted_data.keys(),
        key=lambda pid: (-patient_sample_counts[pid], extract_patient_num(pid))
    )
    # Prepare sample order and all data containers
    sample_order = []
    sample_to_patient = {}
    sample_to_fraction = {}
    for patient in unique_patients:
        patient_samples = [
            (sample, fraction) for sample, fraction in reformatted_data[patient]
            if sample in reformatted_sample_dataframes
        ]
        patient_samples.sort(key=lambda x: parse_sample_id(x[0]))
        for sample, fraction in patient_samples:
            sample_order.append(sample)
            sample_to_patient[sample] = patient
            sample_to_fraction[sample] = fraction
    # Build oncoplot_data in this sample order
    oncoplot_data = pd.DataFrame(0, index=top_genes, columns=sample_order)
    for sample in sample_order:
        df = reformatted_sample_dataframes[sample]
        df_deleted = df[(df['copy_number'] < 2) & (df['gene_symbol'].isin(loss_geneset))]
        deleted_genes = set(df_deleted['gene_symbol'])
        oncoplot_data.loc[list(deleted_genes), sample] = 1
    tumor_fraction = np.array([sample_to_fraction[sample] for sample in sample_order])
    sample_positions = [sample_to_patient[sample] for sample in sample_order]
    sample_labels = sample_order
    gene_frequencies = oncoplot_data.sum(axis=1).sort_values(ascending=False)
    oncoplot_data = oncoplot_data.loc[gene_frequencies.index]
    if oncoplot_data.empty:
        print("⚠️ No genes to plot. Exiting.")
        return
    # Calculate figure size for square boxes
    n_samples = len(sample_order)
    n_genes = len(oncoplot_data.index)
    cell_size = 0.6
    width = max(10, n_samples * cell_size)
    height = max(8, n_genes * cell_size)
    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        sharex=True,
        gridspec_kw={'height_ratios': [0.5, 4]},
        figsize=(width, height)
    )
    x = np.arange(len(sample_order))
    ax1.bar(x, tumor_fraction, color='grey', edgecolor="black", linewidth=0.5)
    ax1.set_ylabel("Tumor Fraction", fontsize=16)
    ax1.set_xlim(-0.5, len(sample_order) - 0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(sample_labels, rotation=90, fontsize=12)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax1.set_yticks([0.2, 0.4, 0.6])
    ax1.tick_params(axis='y', labelsize=12)
    cax = ax2.matshow(oncoplot_data, cmap="Blues", aspect="auto", zorder=1, alpha=0.7)
    ax2.set_xlim(-0.5, len(sample_order) - 0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(sample_labels, rotation=90, fontsize=12)
    ax2.set_yticks(np.arange(len(oncoplot_data.index)))
    ax2.set_yticklabels(oncoplot_data.index, fontsize=24, fontweight='bold')
    ax1.yaxis.grid(True, which='major', linestyle=':', color='gray', alpha=0.7)
    ax2.set_xticks(x - 0.5, minor=True)
    ax2.set_yticks(np.arange(-0.5, len(oncoplot_data.index), 1), minor=True)
    ax2.grid(which="minor", color="gray", linestyle=":", linewidth=0.5, zorder=0)
    patient_ticks = []
    unique_patient_labels = []
    current_patient = sample_positions[0]
    group_start = 0
    for i in range(1, len(sample_positions)):
        if sample_positions[i] != current_patient:
            ax2.axvline(i - 0.5, color="black", linestyle="-", linewidth=1.5, zorder=3)
            midpoint = (group_start + i - 1) / 2
            patient_ticks.append(midpoint)
            unique_patient_labels.append(f"Patient {current_patient}")
            group_start = i
            current_patient = sample_positions[i]
    # Add last patient
    midpoint = (group_start + len(sample_positions) - 1) / 2
    patient_ticks.append(midpoint)
    unique_patient_labels.append(f"Patient {current_patient}")
    ax2.set_xticks(patient_ticks)
    ax2.set_xticklabels(unique_patient_labels, rotation=90, fontsize=12)
    plt.tight_layout()
    plt.savefig("oncoplot_deletions.pdf", bbox_inches='tight')
    plt.close()

# --- Cell 11: Generate plots ---
generate_oncoplot_with_sample_and_patient_labels(reformatted_sample_dataframes, gain_geneset, reformatted_data)
generate_oncoplot_with_sample_and_patient_labels_loss(reformatted_sample_dataframes, loss_geneset, reformatted_data)

# --- Cell 12: Statistical analysis ---
def get_amplified_genes(df):
    return set(df[df['copy_number'] >= 3]['gene_symbol'])

def get_deleted_genes(df):
    return set(df[df['copy_number'] < 2]['gene_symbol'])

def fisher_confint(table, alpha=0.05):
    """Calculate confidence interval for odds ratio using Fisher's exact test"""
    oddsratio, pvalue = fisher_exact(table)
    # Calculate standard error of log odds ratio
    if table[0][0] * table[1][1] * table[0][1] * table[1][0] == 0:
        return oddsratio, (float('nan'), float('nan'))
    se = np.sqrt(sum(1/x for x in table.flatten() if x != 0))
    # Calculate confidence interval
    ci_lower = np.exp(np.log(oddsratio) - stats.norm.ppf(1-alpha/2) * se)
    ci_upper = np.exp(np.log(oddsratio) + stats.norm.ppf(1-alpha/2) * se)
    return oddsratio, (ci_lower, ci_upper)

# Analyze amplifications
all_genes = set()
for df in reformatted_sample_dataframes.values():
    all_genes.update(df['gene_symbol'].unique())

results = []
for gene in all_genes:
    group1_amp = sum(1 for df in reformatted_sample_dataframes.values() if gene in get_amplified_genes(df))
    group1_not_amp = len(reformatted_sample_dataframes) - group1_amp
    
    contingency_table = np.array([[group1_amp, group1_not_amp]])
    if contingency_table[0][0] + contingency_table[0][1] == 0:
        continue
        
    oddsratio, (ci_lower, ci_upper) = fisher_confint(contingency_table)
    results.append({
        'Gene': gene,
        'Amplified_Count': group1_amp,
        'Not_Amplified_Count': group1_not_amp,
        'Odds_Ratio': oddsratio,
        'CI_Lower': ci_lower,
        'CI_Upper': ci_upper
    })

# Create DataFrame and sort by amplification frequency
results_df = pd.DataFrame(results)
results_df['Frequency'] = results_df['Amplified_Count'] / len(reformatted_sample_dataframes)
results_df = results_df.sort_values('Frequency', ascending=False)

# Save results
results_df.to_csv('gene_amplification_statistics.csv', index=False)

# Analyze deletions
results = []
for gene in all_genes:
    group1_del = sum(1 for df in reformatted_sample_dataframes.values() if gene in get_deleted_genes(df))
    group1_not_del = len(reformatted_sample_dataframes) - group1_del
    
    contingency_table = np.array([[group1_del, group1_not_del]])
    if contingency_table[0][0] + contingency_table[0][1] == 0:
        continue
        
    oddsratio, (ci_lower, ci_upper) = fisher_confint(contingency_table)
    results.append({
        'Gene': gene,
        'Deleted_Count': group1_del,
        'Not_Deleted_Count': group1_not_del,
        'Odds_Ratio': oddsratio,
        'CI_Lower': ci_lower,
        'CI_Upper': ci_upper
    })

# Create DataFrame and sort by deletion frequency
results_df = pd.DataFrame(results)
results_df['Frequency'] = results_df['Deleted_Count'] / len(reformatted_sample_dataframes)
results_df = results_df.sort_values('Frequency', ascending=False)

# Save results
results_df.to_csv('gene_deletion_statistics.csv', index=False) 