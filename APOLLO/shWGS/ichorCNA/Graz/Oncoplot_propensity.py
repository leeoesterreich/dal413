# Oncoplot_propensity_longitudinal.py

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
from collections import Counter

# --- Configuration & Setup ---
output_dir = "longitudinal_analysis"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

GAIN_THRESHOLD = 3
LOSS_THRESHOLD = 2
GAIN_GENESET = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4"]
LOSS_GENESET = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]


# --- Data Loading ---
print("Loading and preprocessing data...")

# 1. Load tumor fraction data to define the patient cohort
tf_file_path = 'Novartis TFx.csv'
df_tf = pd.read_csv(tf_file_path)
df_tf.columns = df_tf.columns.str.replace('\\n', ' ', regex=True)
tf_col_name = 'GW z-score (mFAST-SeqS)'
df_tf.rename(columns={col: tf_col_name for col in df_tf.columns if 'GW' in col and 'z-score' in col}, inplace=True)
df_tf = df_tf[df_tf[tf_col_name] > 0].copy()
df_tf[tf_col_name] = np.log10(df_tf[tf_col_name])

# Create a dictionary to hold patient longitudinal data
patient_dict = {}
for _, row in df_tf.iterrows():
    db_id = row["DB-ID"]
    patient_id = db_id.split('_')[0]
    if patient_id not in patient_dict:
        patient_dict[patient_id] = []
    patient_dict[patient_id].append({
        'sample_id': db_id,
        'timepoint': int(db_id.split('_')[1]),
        'tf': row[tf_col_name]
    })
    
cohort_patients = set(patient_dict.keys())
print(f"Defined cohort of {len(cohort_patients)} patients from the tumor fraction file.")

# 2. Load reformatted .cnr files for the cohort
folder_path = "cna_seg_files_reformatted"
sample_dataframes = {}
for file_name in os.listdir(folder_path):
    if file_name.endswith(".cnr"):
        long_sample_id = file_name.replace(".cnr", "")
        patient_id = long_sample_id.split('_')[0]
        short_sample_id = "_".join(long_sample_id.split("_")[:2])

        if patient_id in cohort_patients:
            try:
                df = pd.read_csv(os.path.join(folder_path, file_name), sep="\t")
                copy_number_col = f"{long_sample_id}.Corrected_Copy_Number"
                if 'gene' in df.columns and copy_number_col in df.columns:
                    df_filtered = df[["gene", copy_number_col]].copy()
                    df_filtered.columns = ["gene_symbol", "copy_number"]
                    sample_dataframes[short_sample_id] = df_filtered
            except Exception as e:
                print(f"Could not process file {file_name}: {e}")

print(f"Loaded data for {len(sample_dataframes)} samples.")


# --- Plotting Function ---
def generate_longitudinal_oncoplot(patient_data, sample_dfs, geneset, threshold, alteration_type, output_filename):
    print(f"Generating oncoplot for {alteration_type}...")

    # Sort patients by the number of longitudinal samples
    sorted_patients = sorted(patient_data.keys(), key=lambda p: len(patient_data[p]), reverse=True)
    
    # Prepare sample order and metadata
    sample_order = []
    sample_to_patient = {}
    sample_to_tf = {}
    
    for patient in sorted_patients:
        # Sort samples for each patient by timepoint
        sorted_samples = sorted(patient_data[patient], key=lambda x: x['timepoint'])
        for s in sorted_samples:
            if s['sample_id'] in sample_dfs:
                sample_order.append(s['sample_id'])
                sample_to_patient[s['sample_id']] = patient
                sample_to_tf[s['sample_id']] = s['tf']

    if not sample_order:
        print(f"No valid samples found for {alteration_type} plot. Skipping.")
        return

    # Build the oncoplot matrix
    oncoplot_data = pd.DataFrame(0, index=geneset, columns=sample_order)
    for sample_id in sample_order:
        df_raw = sample_dfs[sample_id]
        # Explode the gene symbols string into one gene per row
        df = df_raw.assign(gene_symbol=df_raw['gene_symbol'].str.split(',')).explode('gene_symbol')
        
        if alteration_type == 'gain':
            altered_genes = set(df[df['copy_number'] >= threshold]['gene_symbol'])
        else: # loss
            altered_genes = set(df[df['copy_number'] < threshold]['gene_symbol'])
        
        genes_in_set = altered_genes.intersection(geneset)
        oncoplot_data.loc[list(genes_in_set), sample_id] = 1

    # Filter and sort data for plotting
    oncoplot_data = oncoplot_data.loc[oncoplot_data.sum(axis=1) > 0]
    gene_frequencies = oncoplot_data.sum(axis=1).sort_values(ascending=False)
    oncoplot_data = oncoplot_data.loc[gene_frequencies.index]

    if oncoplot_data.empty:
        print(f"No alterations found in the given geneset for {alteration_type}. Skipping plot.")
        return

    tumor_fraction = np.array([sample_to_tf[s] for s in oncoplot_data.columns])
    
    # --- Plotting ---
    fig = plt.figure(figsize=(max(20, len(oncoplot_data.columns) * 0.4), max(10, len(oncoplot_data.index) * 0.4)))
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 5], hspace=0.05)
    
    # Tumor fraction bar plot
    ax1 = plt.subplot(gs[0])
    ax1.bar(np.arange(len(oncoplot_data.columns)), tumor_fraction, color='grey', edgecolor="black", linewidth=0.5)
    ax1.set_ylabel("log10(TFx)")
    ax1.set_xlim(-0.5, len(oncoplot_data.columns) - 0.5)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    # Main oncoplot
    ax2 = plt.subplot(gs[1], sharex=ax1)
    cmap = "Reds" if alteration_type == 'gain' else "Blues"
    ax2.matshow(oncoplot_data, cmap=cmap, aspect="auto")
    ax2.set_yticks(np.arange(len(oncoplot_data.index)))
    ax2.set_yticklabels(oncoplot_data.index, fontsize=12)
    ax2.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    ax2.set_xticks(np.arange(len(oncoplot_data.columns)))
    ax2.set_xticklabels(oncoplot_data.columns, rotation=90, fontsize=8)
    
    # Add patient separators and labels at the bottom
    patient_ticks = []
    unique_patient_labels = []
    sample_to_patient_list = [sample_to_patient[s] for s in oncoplot_data.columns]
    
    if not sample_to_patient_list:
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, output_filename), dpi=300)
        plt.close()
        return

    current_patient = sample_to_patient_list[0]
    group_start = 0
    for i in range(1, len(sample_to_patient_list)):
        if sample_to_patient_list[i] != current_patient:
            # Draw line between different patients
            ax2.axvline(x=i - 0.5, color='black', linewidth=1.5)
            # Calculate position for the previous patient's label
            midpoint = (group_start + i - 1) / 2
            patient_ticks.append(midpoint)
            unique_patient_labels.append(current_patient)
            # Start new group
            group_start = i
            current_patient = sample_to_patient_list[i]
    
    # Add the last patient group
    midpoint = (group_start + len(sample_to_patient_list) - 1) / 2
    patient_ticks.append(midpoint)
    unique_patient_labels.append(current_patient)

    # Create and configure the secondary axis for patient labels
    ax_patient = ax2.secondary_xaxis('bottom')
    ax_patient.set_xticks(patient_ticks)
    ax_patient.set_xticklabels(unique_patient_labels, rotation=45, ha="right", fontsize=10)
    ax_patient.spines['bottom'].set_position(('outward', 40))
    ax_patient.tick_params(axis='x', which='major', bottom=False)

    plt.subplots_adjust(bottom=0.2)
    plt.savefig(os.path.join(output_dir, output_filename), dpi=300)
    plt.close()
    print(f"Successfully generated {output_filename}")


# --- Main Execution ---
if __name__ == "__main__":
    generate_longitudinal_oncoplot(
        patient_data=patient_dict,
        sample_dfs=sample_dataframes,
        geneset=GAIN_GENESET,
        threshold=GAIN_THRESHOLD,
        alteration_type='gain',
        output_filename='longitudinal_gains_oncoplot.png'
    )

    generate_longitudinal_oncoplot(
        patient_data=patient_dict,
        sample_dfs=sample_dataframes,
        geneset=LOSS_GENESET,
        threshold=LOSS_THRESHOLD,
        alteration_type='loss',
        output_filename='longitudinal_losses_oncoplot.png'
    ) 