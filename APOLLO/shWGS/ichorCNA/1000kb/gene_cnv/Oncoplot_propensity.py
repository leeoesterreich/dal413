# Oncoplot_propensity.py
# Converted from Oncoplot_propensity.ipynb

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

# Create output directory for CN >= 4 results
output_dir = "CN_GE_4"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def fisher_confint(table, alpha=0.05):
    """Calculate confidence intervals for Fisher's exact test."""
    oddsratio, pvalue = fisher_exact(table, alternative='two-sided')
    a, b = table[0]
    c, d = table[1]
    if 0 in [a, b, c, d]:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    oratio = (a * d) / (b * c)
    se = (1/a + 1/b + 1/c + 1/d) ** 0.5
    z = stats.norm.ppf(1 - alpha/2)
    log_or = np.log(oratio)
    ci_low = np.exp(log_or - z * se)
    ci_high = np.exp(log_or + z * se)
    return ci_low, ci_high

# --- Cell 0: Load and process .cnr files ---
folder_path = "../CNVkit"  # Adjusted for script location
sample_dataframes = {}
for file_name in os.listdir(folder_path):
    if file_name.endswith(".cnr"):
        file_path = os.path.join(folder_path, file_name)
        df = pd.read_csv(file_path, sep="\t")
        sample_id = [col for col in df.columns if col.startswith("TP") and "Corrected_Copy_Number" in col]
        if sample_id:
            sample_id = sample_id[0].split(".")[0]
            df = df[["gene", sample_id + ".Corrected_Copy_Number"]]
            df.columns = ["gene_symbol", "copy_number"]
            df = df.assign(gene_symbol=df["gene_symbol"].str.split(",")).explode("gene_symbol").reset_index(drop=True)
            df = df[df["gene_symbol"] != "-"]
            df["copy_number"] = pd.to_numeric(df["copy_number"], errors="coerce")
            df = df.groupby("gene_symbol", as_index=False)["copy_number"].max()
            sample_dataframes[sample_id] = df
print("Processed sample IDs:", sample_dataframes.keys())

# --- Cell 1: Reformat sample IDs ---
reformatted_sample_dataframes = {key.split("_")[0]: value for key, value in sample_dataframes.items()}
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
    key: [(item[0].split(" ")[0].split("_")[0], item[1]) for item in value]
    for key, value in patient_dict.items()
}

# --- Cell 4: Top amplified genes (copy number >= 4) ---
gene_counter = Counter()
for sample_id, df in reformatted_sample_dataframes.items():
    amplified_genes = set(df[df['copy_number'] >= 4]['gene_symbol'])
    gene_counter.update(amplified_genes)
top_amplified_genes_list = gene_counter.most_common()
print("Top amplified genes (copy number >= 4):")
with open(os.path.join(output_dir, "top_amplified_genes_CN_GE_4.txt"), "w") as f:
    f.write("Gene\tCount\n")
    for gene, count in top_amplified_genes_list:
        print(f"{gene}: {count}")
        f.write(f"{gene}\t{count}\n")

# --- Cell 5: Top deleted genes (copy number < 2) ---
gene_counter = Counter()
for sample_id, df in reformatted_sample_dataframes.items():
    deleted_genes = set(df[df['copy_number'] < 2]['gene_symbol'])
    gene_counter.update(deleted_genes)
top_deleted_genes_list = gene_counter.most_common()
print("\nTop deleted genes (copy number < 2):")
with open(os.path.join(output_dir, "top_deleted_genes.txt"), "w") as f:
    f.write("Gene\tCount\n")
    for gene, count in top_deleted_genes_list:
        print(f"{gene}: {count}")
        f.write(f"{gene}\t{count}\n")

# --- Cell 6: Top amplified genes across patients ---
gene_counter = Counter()
for patient, samples in reformatted_data.items():
    patient_amplified_genes = set()
    for sample, _ in samples:
        if sample in reformatted_sample_dataframes:
            df = reformatted_sample_dataframes[sample]
            amplified_genes = df[df['copy_number'] >= 4]['gene_symbol']
            patient_amplified_genes.update(amplified_genes)
    gene_counter.update(patient_amplified_genes)
top_amplified_genes_list = gene_counter.most_common()
print("\nTop amplified genes (copy number >= 4 across patients):")
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
        df_amplified = df[(df['copy_number'] >= 4) & (df['gene_symbol'].isin(gain_geneset))]
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
    midpoint = (group_start + len(sample_positions) - 1) / 2
    patient_ticks.append(midpoint)
    unique_patient_labels.append(f"Patient {current_patient}")
    ax_patient = ax2.secondary_xaxis('bottom')
    ax_patient.set_xticks(patient_ticks)
    ax_patient.set_xticklabels(unique_patient_labels, rotation=45, fontsize=16, ha="right")
    ax_patient.spines['bottom'].set_position(('outward', 80))
    plt.subplots_adjust(bottom=0.35)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(output_dir, "oncoplot_amplification.png"), dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to free memory

generate_oncoplot_with_sample_and_patient_labels(reformatted_sample_dataframes, gain_geneset, reformatted_data)

# --- Cell 9: Loss geneset ---
loss_geneset = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]

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
    midpoint = (group_start + len(sample_positions) - 1) / 2
    patient_ticks.append(midpoint)
    unique_patient_labels.append(f"Patient {current_patient}")
    ax_patient = ax2.secondary_xaxis('bottom')
    ax_patient.set_xticks(patient_ticks)
    ax_patient.set_xticklabels(unique_patient_labels, rotation=45, fontsize=16, ha="right")
    ax_patient.spines['bottom'].set_position(('outward', 80))
    plt.subplots_adjust(bottom=0.35)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(output_dir, "oncoplot_deletion.png"), dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to free memory

generate_oncoplot_with_sample_and_patient_labels_loss(reformatted_sample_dataframes, loss_geneset, reformatted_data)

# --- Cell 11: Statistical analysis ---
def get_amplified_genes(df):
    return set(df.loc[df['copy_number'] >= 4, 'gene_symbol'])

# Calculate median tumor fraction
sample_tumor_fraction = {}
for patient, samples in reformatted_data.items():
    for sample, fraction in samples:
        if sample in reformatted_sample_dataframes:
            sample_tumor_fraction[sample] = fraction

median_tumor_fraction = np.median(list(sample_tumor_fraction.values()))
high_tumor_fraction_samples = [s for s, tf in sample_tumor_fraction.items() if tf >= median_tumor_fraction]
low_tumor_fraction_samples = [s for s, tf in sample_tumor_fraction.items() if tf < median_tumor_fraction]

# Analyze amplifications
genes_of_interest = gain_geneset
final_results = []

for gene in genes_of_interest:
    # High tumor fraction group
    high_tf_amplified = sum(1 for sample in high_tumor_fraction_samples 
                           if gene in get_amplified_genes(reformatted_sample_dataframes[sample]))
    high_tf_not_amplified = len(high_tumor_fraction_samples) - high_tf_amplified
    
    # Low tumor fraction group
    low_tf_amplified = sum(1 for sample in low_tumor_fraction_samples 
                          if gene in get_amplified_genes(reformatted_sample_dataframes[sample]))
    low_tf_not_amplified = len(low_tumor_fraction_samples) - low_tf_amplified
    
    contingency_table = [[high_tf_amplified, high_tf_not_amplified],
                        [low_tf_amplified, low_tf_not_amplified]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    ci_lower, ci_upper = fisher_confint(contingency_table)
    
    final_results.append({
        'gene': gene,
        'a_amplified_gene': high_tf_amplified,
        'b_not_amplified_gene': high_tf_not_amplified,
        'c_amplified_others': low_tf_amplified,
        'd_not_amplified_others': low_tf_not_amplified,
        'oddsratio': oddsratio,
        'pvalue': pvalue,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper
    })

if final_results:
    results_df = pd.DataFrame(final_results)
    results_df['pvalue_adj'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]
    results_df['reject'] = results_df['pvalue_adj'] < 0.05
    
    print("\nFinal Results:")
    print(results_df[['gene', 'a_amplified_gene', 'b_not_amplified_gene', 'c_amplified_others', 'd_not_amplified_others', 'oddsratio', 'pvalue', 'pvalue_adj', 'reject']])
    results_df.to_csv(os.path.join(output_dir, "propensity_results.csv"), index=False)
    
    # Plot Odds Ratios
    plt.figure(figsize=(12, 8))
    plt.scatter(range(len(results_df)), np.log10(results_df['oddsratio']), 
               c=['red' if x else 'blue' for x in results_df['reject']])
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    plt.xticks(range(len(results_df)), results_df['gene'], rotation=90)
    plt.ylabel('log10(Odds Ratio)', fontsize=12)
    plt.title("log10(Odds Ratios) for Gene Amplifications (High vs Low Tumor Fraction)\nSignificant Genes Highlighted", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "amplification_oddsratio.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Forest plot
    plt.figure(figsize=(12, len(results_df) * 0.3))
    y_pos = np.arange(len(results_df))
    plt.errorbar(np.log10(results_df['oddsratio']), y_pos,
                xerr=[np.log10(results_df['oddsratio']) - np.log10(results_df['ci_lower']),
                      np.log10(results_df['ci_upper']) - np.log10(results_df['oddsratio'])],
                fmt='o', capsize=5,
                c=['red' if x else 'blue' for x in results_df['reject']])
    plt.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    plt.yticks(y_pos, results_df['gene'])
    plt.xlabel('log10(Odds Ratio)', fontsize=12)
    plt.title('Forest Plot: log10(Odds Ratios) for Gene Amplifications (High vs Low Tumor Fraction)\n(95% CI, Red = FDR Significant)', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'amplification_oddsratio_forest.png'), dpi=300, bbox_inches='tight')
    plt.close()

# Analyze deletions similarly
genes_of_interest = loss_geneset
final_results = []

for gene in genes_of_interest:
    # High tumor fraction group
    high_tf_deleted = sum(1 for sample in high_tumor_fraction_samples 
                         if gene in get_deleted_genes(reformatted_sample_dataframes[sample]))
    high_tf_not_deleted = len(high_tumor_fraction_samples) - high_tf_deleted
    
    # Low tumor fraction group
    low_tf_deleted = sum(1 for sample in low_tumor_fraction_samples 
                        if gene in get_deleted_genes(reformatted_sample_dataframes[sample]))
    low_tf_not_deleted = len(low_tumor_fraction_samples) - low_tf_deleted
    
    contingency_table = [[high_tf_deleted, high_tf_not_deleted],
                        [low_tf_deleted, low_tf_not_deleted]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    ci_lower, ci_upper = fisher_confint(contingency_table)
    
    final_results.append({
        'gene': gene,
        'a_deleted_gene': high_tf_deleted,
        'b_not_deleted_gene': high_tf_not_deleted,
        'c_deleted_others': low_tf_deleted,
        'd_not_deleted_others': low_tf_not_deleted,
        'oddsratio': oddsratio,
        'pvalue': pvalue,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper
    })

if final_results:
    results_df_del = pd.DataFrame(final_results)
    results_df_del['pvalue_adj'] = multipletests(results_df_del['pvalue'], method='fdr_bh')[1]
    results_df_del['reject'] = results_df_del['pvalue_adj'] < 0.05
    
    print("\nFinal Deletion Results:")
    print(results_df_del[['gene', 'a_deleted_gene', 'b_not_deleted_gene', 'c_deleted_others', 'd_not_deleted_others', 'oddsratio', 'pvalue', 'pvalue_adj', 'reject']])
    results_df_del.to_csv(os.path.join(output_dir, "propensity_results_deletion.csv"), index=False)
    
    # Plot Odds Ratios
    plt.figure(figsize=(12, 8))
    plt.scatter(range(len(results_df_del)), np.log10(results_df_del['oddsratio']), 
               c=['blue' if x else 'gray' for x in results_df_del['reject']])
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    plt.xticks(range(len(results_df_del)), results_df_del['gene'], rotation=90)
    plt.ylabel('log10(Odds Ratio)', fontsize=12)
    plt.title("log10(Odds Ratios) for Gene Deletions (High vs Low Tumor Fraction)\nSignificant Genes Highlighted", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "deletion_oddsratio.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Forest plot
    plt.figure(figsize=(12, len(results_df_del) * 0.3))
    y_pos = np.arange(len(results_df_del))
    plt.errorbar(np.log10(results_df_del['oddsratio']), y_pos,
                xerr=[np.log10(results_df_del['oddsratio']) - np.log10(results_df_del['ci_lower']),
                      np.log10(results_df_del['ci_upper']) - np.log10(results_df_del['oddsratio'])],
                fmt='o', capsize=5,
                c=['blue' if x else 'gray' for x in results_df_del['reject']])
    plt.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    plt.yticks(y_pos, results_df_del['gene'])
    plt.xlabel('log10(Odds Ratio)', fontsize=12)
    plt.title('Forest Plot: log10(Odds Ratios) for Gene Deletions (High vs Low Tumor Fraction)\n(95% CI, Blue = FDR Significant)', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'deletion_oddsratio_forest.png'), dpi=300, bbox_inches='tight')
    plt.close() 