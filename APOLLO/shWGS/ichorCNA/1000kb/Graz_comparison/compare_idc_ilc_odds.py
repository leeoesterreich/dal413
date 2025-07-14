import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import sys

# --- Gene Sets ---
# Using the gene sets found in both analysis scripts
GAIN_GENESET = ["CCND1", "FGF19", "FGF4", "FGF3", "FGFR1", "NSD3", "PAK1", "ERBB2", "MYC", "CDK12", "RAD21", "PPM1D", "RECQL4", "NBN", "GNAS", "BRIP1", "RPS6KB2", "RAD51C", "MDM2", "MCL1", "CD79B", "AURKA", "ELOC", "SPOP", "PRKAR1A", "MDM4", "AGO2", "INPPL1", "ELF3", "FOXA1", "RNF43", "AXIN2", "RARA", "RTEL1", "IKBKE", "IL10", "IGF1R", "PREX2", "MSI2", "SOX17", "PRDM14", "PRDM1", "GATA3", "CDKN1B", "LYN", "FGFR2", "NCOA3", "ESR1", "SOX9", "CDK4", "SDHC", "DDR2", "EGFR", "PTPRT", "NUF2", "TERT", "BRCA1", "HOXB13", "SMYD3", "CCND2"]
LOSS_GENESET = ["CDKN2B", "CDKN2A", "PTEN", "MAP2K4", "RB1", "DUSP4", "RAC2", "NCOR1", "NF1", "TEK", "BIRC3", "ZFHX3", "EPHA7", "PRDM1", "TP53", "CRLF2", "FYN", "FAT1", "PTPRD", "MAP3K1"]

CNA_AMP_THRESHOLD = 3
CNA_DEL_THRESHOLD = 2

# --- ILC Data Processing ---
def get_ilc_alterations():
    """
    Loads and processes ILC data to get per-patient gene alterations.
    Logic adapted from /gene_cnv/collapse/analyze_patient_cnv.py
    """
    print("--- Processing ILC Data ---")
    
    # Load and process .cnr files
    folder_path = "CNVkit"
    sample_dataframes = {}
    print(f"Reading ILC .cnr files from: {folder_path}")
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
    print(f"Processed {len(reformatted_sample_dataframes)} unique ILC samples.")

    # Load tumor fraction Excel file to get patient-sample mapping
    file_path = 'results/TFx_hg38_1000kb.xlsx'
    print(f"Reading ILC patient info from: {file_path}")
    df_tf = pd.read_excel(file_path, sheet_name=0, skiprows=1)
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
            tumor_fraction = row.iloc[tumor_fraction_col_index]
            if pd.notna(sample_label) and pd.notna(tumor_fraction):
                tuples.append((sample_label, tumor_fraction))
        patient_dict[patient_code] = tuples

    reformatted_data = {
        key: [(item[0].split(" ")[0].split("_")[0], item[1]) for item in value]
        for key, value in patient_dict.items()
    }
    patient_to_exclude = 4908
    if patient_to_exclude in reformatted_data:
        del reformatted_data[patient_to_exclude]

    # Collapse to patient-level alterations
    patient_amplifications = {}
    patient_deletions = {}
    for patient_id, samples in reformatted_data.items():
        patient_amp_genes = set()
        patient_del_genes = set()
        valid_samples = [s for s in samples if s[0] in reformatted_sample_dataframes]
        if not valid_samples: continue
        for sample_id, _ in valid_samples:
            sample_df = reformatted_sample_dataframes[sample_id]
            amp_genes = set(sample_df[sample_df['copy_number'] >= CNA_AMP_THRESHOLD]['gene_symbol'])
            del_genes = set(sample_df[sample_df['copy_number'] < CNA_DEL_THRESHOLD]['gene_symbol'])
            patient_amp_genes.update(amp_genes)
            patient_del_genes.update(del_genes)
        patient_amplifications[patient_id] = patient_amp_genes
        patient_deletions[patient_id] = patient_del_genes
        
    print(f"Collapsed to {len(patient_amplifications)} ILC patients.")
    return patient_amplifications, patient_deletions

# --- IDC Data Processing ---
def get_idc_alterations():
    """
    Loads and processes IDC data to get per-patient gene alterations.
    Logic adapted from /Graz/analyze_patient_cnv.py
    """
    print("\n--- Processing IDC Data ---")
    
    # Load clinical data to filter for NST patients
    clinical_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/CtDNAHRHER2Novartis-ClinicalNina_DATA_LABELS_2025-05-14_1651.csv'
    print(f"Reading IDC clinical data from: {clinical_file_path}")
    df_clinical = pd.read_csv(clinical_file_path)
    nst_patients = set(df_clinical[df_clinical['Histology'] == 'NST']['HG-ID'])
    print(f"Found {len(nst_patients)} patients with NST histology (IDC).")

    # Load and process .cnr files
    folder_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/Graz/cna_seg_files_reformatted/"
    print(f"Reading IDC .cnr files from: {folder_path}")
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
                except Exception as e:
                    print(f"Could not process file {file_name}: {e}", file=sys.stderr)

    print(f"Processed {len(sample_dataframes)} total IDC samples.")

    # Collapse to patient-level alterations
    patient_amplifications = {}
    patient_deletions = {}
    for patient_id in nst_patients:
        patient_amp_genes = set()
        patient_del_genes = set()
        patient_sample_ids = [sid for sid in sample_dataframes.keys() if sid.startswith(f"{patient_id}_")]
        if not patient_sample_ids: continue
        for sample_id in patient_sample_ids:
            sample_df = sample_dataframes[sample_id]
            amp_genes = set(sample_df[sample_df['copy_number'] >= CNA_AMP_THRESHOLD]['gene_symbol'])
            del_genes = set(sample_df[sample_df['copy_number'] < CNA_DEL_THRESHOLD]['gene_symbol'])
            patient_amp_genes.update(amp_genes)
            patient_del_genes.update(del_genes)
        patient_amplifications[patient_id] = patient_amp_genes
        patient_deletions[patient_id] = patient_del_genes

    print(f"Collapsed to {len(patient_amplifications)} IDC patients.")
    return patient_amplifications, patient_deletions

# --- Odds Ratio Calculation and Plotting ---
def calculate_odds_ratio(gene, ilc_alts, idc_alts):
    """Calculate odds ratio for a gene between ILC and IDC."""
    ilc_pos = sum(1 for patient_id in ilc_alts if gene in ilc_alts[patient_id])
    idc_pos = sum(1 for patient_id in idc_alts if gene in idc_alts[patient_id])
    
    ilc_total = len(ilc_alts)
    idc_total = len(idc_alts)
    
    ilc_neg = ilc_total - ilc_pos
    idc_neg = idc_total - idc_pos
    
    contingency = np.array([[ilc_pos, ilc_neg], [idc_pos, idc_neg]])
    
    # Add 0.5 to all cells for stability if any cell is 0
    if np.any(contingency == 0):
        contingency = contingency + 0.5

    odds_ratio, pvalue = fisher_exact(contingency)
    
    log_odds = np.log(odds_ratio)
    se = np.sqrt(np.sum(1.0 / contingency))
    ci_lower = np.exp(log_odds - 1.96 * se)
    ci_upper = np.exp(log_odds + 1.96 * se)
    
    return {
        'gene': gene,
        'odds_ratio': odds_ratio,
        'pvalue': pvalue,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'ilc_count': ilc_pos,
        'ilc_total': ilc_total,
        'idc_count': idc_pos,
        'idc_total': idc_total,
        'ilc_freq': ilc_pos / ilc_total if ilc_total > 0 else 0,
        'idc_freq': idc_pos / idc_total if idc_total > 0 else 0,
    }

def plot_forest(results_df, title, output_file, is_deletion=False):
    """Create a forest plot for odds ratios."""
    plt.figure(figsize=(12, max(8, len(results_df) * 0.4)))
    
    results_df = results_df.sort_values('odds_ratio', ascending=True).reset_index(drop=True)
    
    y_pos = np.arange(len(results_df))
    log_odds = np.log2(results_df['odds_ratio'].astype(float))
    log_ci_lower = np.log2(results_df['ci_lower'].astype(float))
    log_ci_upper = np.log2(results_df['ci_upper'].astype(float))
    
    highlight_color = 'blue' if is_deletion else 'red'
    
    for i in range(len(results_df)):
        # Skip non-finite values that can cause plotting errors
        if not all(np.isfinite([log_odds[i], log_ci_lower[i], log_ci_upper[i]])):
            print(f"Skipping gene {results_df['gene'][i]} due to non-finite values for plotting.")
            continue

        color = highlight_color if results_df['reject'][i] else 'grey'
        plt.errorbar(log_odds[i], y_pos[i], 
                     xerr=[[log_odds[i] - log_ci_lower[i]], [log_ci_upper[i] - log_odds[i]]],
                     fmt='o', color='black', ecolor=color, elinewidth=1, capsize=3, markersize=5, zorder=2)
        plt.scatter(log_odds[i], y_pos[i], color=color, s=50, zorder=3)

    plt.axvline(x=0, color='black', linestyle='--', linewidth=1)
    
    plt.yticks(y_pos, results_df['gene'], fontsize=12)
    plt.xlabel('Log2(Odds Ratio) - ILC vs. IDC', fontsize=12)
    plt.title(title, fontsize=16, pad=20)
    plt.grid(True, axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved forest plot to {output_file}")


# --- Main Execution ---
def main():
    output_dir = "Graz_comparison"
    os.makedirs(output_dir, exist_ok=True)
    # Get patient-level alteration data
    ilc_amps_full, ilc_dels_full = get_ilc_alterations()
    idc_amps_full, idc_dels_full = get_idc_alterations()

    if not ilc_amps_full or not idc_amps_full:
        print("Error: Could not load data for one or both cohorts. Exiting.", file=sys.stderr)
        sys.exit(1)

    # --- Amplification Analysis ---
    print("\n--- Analyzing Amplifications (ILC vs. IDC) ---")
    
    # Filter for patients with at least one amplification
    ilc_amps_filtered = {p: g for p, g in ilc_amps_full.items() if g}
    idc_amps_filtered = {p: g for p, g in idc_amps_full.items() if g}
    print(f"Including {len(ilc_amps_filtered)} ILC patients with >=1 amplification.")
    print(f"Including {len(idc_amps_filtered)} IDC patients with >=1 amplification.")

    amp_results = [calculate_odds_ratio(gene, ilc_amps_filtered, idc_amps_filtered) for gene in GAIN_GENESET]
    amp_df = pd.DataFrame(amp_results)
    
    if not amp_df.empty:
        amp_df['reject'] = multipletests(amp_df['pvalue'], method='fdr_bh')[0]
        amp_df = amp_df.sort_values('odds_ratio', ascending=False)
        amp_output_csv = os.path.join(output_dir, 'ilc_vs_idc_amplification_odds.csv')
        amp_df.to_csv(amp_output_csv, index=False)
        print(f"Saved amplification results to {amp_output_csv}")
        
        plot_forest(amp_df, 'Amplification Odds Ratios (ILC vs. IDC)', os.path.join(output_dir, 'ilc_vs_idc_amplification_forest.png'), is_deletion=False)

    # --- Deletion Analysis ---
    print("\n--- Analyzing Deletions (ILC vs. IDC) ---")

    # Filter for patients with at least one deletion
    ilc_dels_filtered = {p: g for p, g in ilc_dels_full.items() if g}
    idc_dels_filtered = {p: g for p, g in idc_dels_full.items() if g}
    print(f"Including {len(ilc_dels_filtered)} ILC patients with >=1 deletion.")
    print(f"Including {len(idc_dels_filtered)} IDC patients with >=1 deletion.")

    del_results = [calculate_odds_ratio(gene, ilc_dels_filtered, idc_dels_filtered) for gene in LOSS_GENESET]
    del_df = pd.DataFrame(del_results)

    if not del_df.empty:
        del_df['reject'] = multipletests(del_df['pvalue'], method='fdr_bh')[0]
        del_df = del_df.sort_values('odds_ratio', ascending=False)
        del_output_csv = os.path.join(output_dir, 'ilc_vs_idc_deletion_odds.csv')
        del_df.to_csv(del_output_csv, index=False)
        print(f"Saved deletion results to {del_output_csv}")
        
        plot_forest(del_df, 'Deletion Odds Ratios (ILC vs. IDC)', os.path.join(output_dir, 'ilc_vs_idc_deletion_forest.png'), is_deletion=True)

    print("\nAnalysis complete.")

if __name__ == "__main__":
    main() 