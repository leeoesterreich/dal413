import os
import pandas as pd
import numpy as np
import re

def main():
    """
    Correctly examines the contingency tables for all genes with infinite odds ratios
    in both amplification and deletion results.
    """
    print(f"--- Examining All Genes with Infinite Odds Ratios (Corrected Logic) ---")

    # --- Data Loading and Preprocessing (Verified from main script) ---
    folder_path = "../../CNVkit"
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
    print("Data loading and preprocessing complete.")

    # --- Patient-Level Collapse (Verified from main script) ---
    CNA_AMP_THRESHOLD = 3
    CNA_DEL_THRESHOLD = 2
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

        for _, tumor_fraction in valid_samples:
            if tumor_fraction > max_tf:
                max_tf = tumor_fraction
        patient_to_max_tf[patient_id] = max_tf

        for sample_id, _ in valid_samples:
            sample_df = reformatted_sample_dataframes[sample_id]
            amp_genes = set(sample_df[sample_df['copy_number'] >= CNA_AMP_THRESHOLD]['gene_symbol'])
            del_genes = set(sample_df[sample_df['copy_number'] < CNA_DEL_THRESHOLD]['gene_symbol'])
            patient_amp_genes.update(amp_genes)
            patient_del_genes.update(del_genes)
            
        patient_amplifications[patient_id] = patient_amp_genes
        patient_deletions[patient_id] = patient_del_genes

    print(f"Collapsed data for {len(patient_amplifications)} patients.")

    # --- Identify Genes and Calculate Tables (Verified from main script) ---
    median_tf = np.median(list(patient_to_max_tf.values()))
    high_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf >= median_tf}
    low_tf_patients = {p for p, tf in patient_to_max_tf.items() if tf < median_tf}

    print(f"\nMedian Max Tumor Fraction: {median_tf:.3f}")
    print(f"Number of High TF Patients: {len(high_tf_patients)}")
    print(f"Number of Low TF Patients: {len(low_tf_patients)}")

    # Check Amplifications
    try:
        amp_results_df = pd.read_csv("patient_level_amp_propensity.csv")
        inf_amp_genes = amp_results_df[np.isinf(amp_results_df['oddsratio'])]['gene'].tolist()
        
        print("\n\n--- Infinite Odds Ratios in AMPLIFICATION Results ---")
        if not inf_amp_genes:
            print("None found.")
        else:
            for gene in inf_amp_genes:
                high_tf_amp = sum(1 for p in high_tf_patients if gene in patient_amplifications.get(p, set()))
                high_tf_no_amp = len(high_tf_patients) - high_tf_amp
                low_tf_amp = sum(1 for p in low_tf_patients if gene in patient_amplifications.get(p, set()))
                low_tf_no_amp = len(low_tf_patients) - low_tf_amp
                
                print(f"\n--- Contingency Table for Gene: {gene} (Amplification) ---")
                print(f"                     | High TF | Low TF |")
                print("---------------------------------------")
                print(f"  Amplified          |    {high_tf_amp:<2}    |    {low_tf_amp:<2}   |")
                print(f"  Not Amplified      |    {high_tf_no_amp:<2}    |   {low_tf_no_amp:<2}   |")
                print("---------------------------------------")
    except FileNotFoundError:
        print("\nCould not find patient_level_amp_propensity.csv")

    # Check Deletions
    try:
        del_results_df = pd.read_csv("patient_level_del_propensity.csv")
        inf_del_genes = del_results_df[np.isinf(del_results_df['oddsratio'])]['gene'].tolist()

        print("\n\n--- Infinite Odds Ratios in DELETION Results ---")
        if not inf_del_genes:
            print("None found.")
        else:
            for gene in inf_del_genes:
                high_tf_del = sum(1 for p in high_tf_patients if gene in patient_deletions.get(p, set()))
                high_tf_no_del = len(high_tf_patients) - high_tf_del
                low_tf_del = sum(1 for p in low_tf_patients if gene in patient_deletions.get(p, set()))
                low_tf_no_del = len(low_tf_patients) - low_tf_del

                print(f"\n--- Contingency Table for Gene: {gene} (Deletion) ---")
                print(f"                     | High TF | Low TF |")
                print("---------------------------------------")
                print(f"  Deleted            |    {high_tf_del:<2}    |    {low_tf_del:<2}   |")
                print(f"  Not Deleted        |    {high_tf_no_del:<2}    |   {low_tf_no_del:<2}   |")
                print("---------------------------------------")
    except FileNotFoundError:
        print("\nCould not find patient_level_del_propensity.csv")

if __name__ == "__main__":
    main() 