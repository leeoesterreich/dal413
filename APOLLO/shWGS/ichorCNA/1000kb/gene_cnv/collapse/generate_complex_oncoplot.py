import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns

def get_alteration_code(cn):
    """Converts a copy number value to a categorical code."""
    if cn >= 5:
        return 4  # High-Level Amplification
    if 4 <= cn < 5:
        return 3  # Amplification
    if 3 <= cn < 4:
        return 2  # Gain
    if cn < 2:
        return 1  # Deletion
    return 0  # No alteration

def main():
    """
    Generates a complex, multi-layered oncoplot for the top 50 most frequently
    altered genes, based on collapsed patient-level copy number data.
    """
    print("--- Generating Complex Oncoplot ---")

    # --- 1. Load and Process Data (Reusing verified logic) ---
    # Load .cnr files
    folder_path = "CNVkit"
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
    
    reformatted_sample_dataframes = {key.split("_")[0]: value for key, value in sample_dataframes.items()}

    # Load tumor fraction and patient mapping
    file_path_tf = 'results/TFx_hg38_1000kb.xlsx'
    df_tf = pd.read_excel(file_path_tf, sheet_name=0, skiprows=1)
    df_tf.columns = df_tf.columns.map(str)
    df_tf = df_tf.loc[:, ~df_tf.columns.str.startswith('Tumor_fraction_1000kb')]
    df_tf.rename(columns=lambda x: 'Tumor_fraction' if x.startswith('Tumor_fraction_500kb') else x, inplace=True)
    patient_dict = {}
    for _, row in df_tf.iterrows():
        patient_code = row["Patient_Code"]
        tuples = []
        for i in range(1, len(df_tf.columns) // 2 + 1):
            if i == 1: sample_label_col = f"{i}st progression"
            elif i == 2: sample_label_col = f"{i}nd progression"
            elif i == 3: sample_label_col = f"{i}rd progression"
            else: sample_label_col = f"{i}th progression"
            if sample_label_col not in df_tf.columns: continue
            sample_label_col_index = df_tf.columns.get_loc(sample_label_col)
            tumor_fraction_col_index = sample_label_col_index + 1
            sample_label = row[sample_label_col]
            if pd.notna(sample_label):
                tuples.append(sample_label)
        patient_dict[patient_code] = tuples
    
    # Invert patient_dict to map sample to patient
    sample_to_patient = {}
    for patient, samples in patient_dict.items():
        for sample in samples:
            sample_base = sample.split(" ")[0].split("_")[0]
            sample_to_patient[sample_base] = patient

    # Exclude patient 4908
    patient_to_exclude = 4908
    samples_to_exclude = {s for s, p in sample_to_patient.items() if p == patient_to_exclude}
    reformatted_sample_dataframes = {s: df for s, df in reformatted_sample_dataframes.items() if s not in samples_to_exclude}
    print(f"Patient {patient_to_exclude} and their samples have been excluded.")
    
    # Create a unified DataFrame: [patient, gene_symbol, copy_number]
    all_cn_data = []
    for sample_id, df in reformatted_sample_dataframes.items():
        patient_id = sample_to_patient.get(sample_id)
        if patient_id:
            df['patient_id'] = patient_id
            all_cn_data.append(df)
    
    if not all_cn_data:
        print("No data to process after filtering. Exiting.")
        return
        
    master_df = pd.concat(all_cn_data, ignore_index=True)

    # --- 2. Collapse to Patient-Level Max Copy Number ---
    patient_level_cn = master_df.groupby(['patient_id', 'gene_symbol'])['copy_number'].max().unstack()

    # --- 3. Use Predefined Gene List ---
    # Apply the alteration coding
    patient_level_coded = patient_level_cn.apply(np.vectorize(get_alteration_code))
    
    # Define the gene sets from the propensity analysis
    gain_geneset = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4"]
    loss_geneset = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]
    genes_of_interest = sorted(list(set(gain_geneset) | set(loss_geneset)))

    # Filter the coded matrix to only the genes of interest that are present in our data
    available_genes = [gene for gene in genes_of_interest if gene in patient_level_coded.columns]
    filtered_coded_df = patient_level_coded[available_genes]

    # Calculate frequency of any alteration for sorting purposes
    alteration_freq = (filtered_coded_df > 0).sum(axis=0) / filtered_coded_df.shape[0]
    sorted_genes = alteration_freq.sort_values(ascending=False).index
    
    # Final sorted DataFrame for plotting
    final_plotting_df = filtered_coded_df[sorted_genes]
    
    # Sort patients based on a sort key (e.g., number of alterations) for better visualization
    patient_sort_key = (final_plotting_df > 0).sum(axis=1)
    sorted_patients = patient_sort_key.sort_values(ascending=False).index
    final_plotting_df = final_plotting_df.loc[sorted_patients]

    print(f"Processed data for {len(sorted_patients)} patients and {len(sorted_genes)} genes.")

    # --- 4. Generate the Oncoplot ---
    fig = plt.figure(figsize=(20, 24))
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 10], width_ratios=[18, 5], 
                           wspace=0.05, hspace=0.05)
    
    ax_heatmap = fig.add_subplot(gs[1, 0])
    ax_top_bar = fig.add_subplot(gs[0, 0], sharex=ax_heatmap)
    ax_right_bar = fig.add_subplot(gs[1, 1], sharey=ax_heatmap)

    # Define colors
    colors = ['#c0c0c0', '#0072b2', '#f0e442', '#e69f00', '#d55e00']
    cmap = plt.cm.colors.ListedColormap(colors[1:] + [colors[0]]) # Reorder for heatmap: No Alteration is last

    # Main Heatmap
    sns.heatmap(final_plotting_df.T, ax=ax_heatmap, cmap=cmap, cbar=False, linewidths=0.5, linecolor='white', vmin=0, vmax=4)
    ax_heatmap.set_xlabel("Patients", fontsize=16)
    ax_heatmap.set_ylabel("")
    plt.setp(ax_heatmap.get_xticklabels(), rotation=90)
    
    # Y-tick labels with gene name and frequency
    gene_freq_percent = alteration_freq[sorted_genes] * 100
    ax_heatmap.set_yticks(np.arange(len(sorted_genes)) + 0.5)
    ax_heatmap.set_yticklabels([f"{gene}   {freq:.0f}%" for gene, freq in gene_freq_percent.items()], fontsize=14)

    # Top Stacked Bar Plot (Alterations per patient)
    patient_counts = pd.DataFrame({
        'Deletion': (final_plotting_df == 1).sum(axis=1),
        'Gain': (final_plotting_df == 2).sum(axis=1),
        'Amplification': (final_plotting_df == 3).sum(axis=1),
        'High-Level Amp': (final_plotting_df == 4).sum(axis=1)
    }).loc[sorted_patients]
    
    bottom = np.zeros(len(sorted_patients))
    bar_positions = np.arange(len(sorted_patients)) + 0.5 # Center the bars
    for i, col in enumerate(['Deletion', 'Gain', 'Amplification', 'High-Level Amp']):
        if col in patient_counts.columns:
            ax_top_bar.bar(bar_positions, patient_counts[col], bottom=bottom, width=0.8,
                           color=colors[i+1], label=col)
            bottom += patient_counts[col]

    ax_top_bar.set_title("CNA calling of cfDNA ULP-WGS", fontsize=24)
    plt.setp(ax_top_bar.get_xticklabels(), visible=False)
    ax_top_bar.tick_params(axis='x', length=0)
    ax_top_bar.set_ylabel("Altered Genes", fontsize=14)

    # Right Stacked Bar Plot (Frequency per gene)
    gene_counts_percent_stacked = pd.DataFrame({
        'Deletion': (final_plotting_df == 1).sum(axis=0),
        'Gain': (final_plotting_df == 2).sum(axis=0),
        'Amplification': (final_plotting_df == 3).sum(axis=0),
        'High-Level Amp': (final_plotting_df == 4).sum(axis=0)
    }).T
    gene_counts_percent_stacked = (gene_counts_percent_stacked / len(sorted_patients)) * 100
    gene_counts_percent_stacked = gene_counts_percent_stacked.reindex(sorted_genes) # Ensure order matches heatmap
    
    left = np.zeros(len(sorted_genes))
    bar_positions_h = np.arange(len(sorted_genes)) + 0.5 # Center the bars
    for i, col in enumerate(['Deletion', 'Gain', 'Amplification', 'High-Level Amp']):
        if col in gene_counts_percent_stacked.columns:
            ax_right_bar.barh(bar_positions_h, gene_counts_percent_stacked[col], left=left, height=0.8,
                              color=colors[i+1], label=col)
            left += gene_counts_percent_stacked[col]
    
    plt.setp(ax_right_bar.get_yticklabels(), visible=False)
    ax_right_bar.tick_params(axis='y', length=0)
    ax_right_bar.set_xlabel("% Patients", fontsize=14)

    # Legend
    legend_elements = [
        Patch(facecolor=colors[1], label='Heterozygous Deletion'),
        Patch(facecolor=colors[2], label='Gain'),
        Patch(facecolor=colors[3], label='Amplification'),
        Patch(facecolor=colors[4], label='High Level Amplification')
    ]
    ax_right_bar.legend(handles=legend_elements, loc='upper right', title='Alterations', fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    output_filename = "gene_cnv/collapse/complex_oncoplot.png"
    plt.savefig(output_filename, dpi=300)
    plt.close()
    
    print(f"\nComplex oncoplot saved to {output_filename}")

if __name__ == "__main__":
    main() 