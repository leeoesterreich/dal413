import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import re

def parse_sample_id(sample_id):
    """Parse TP sample ID to get year and number for temporal ordering."""
    match = re.match(r"TP(\d+)-M(\d+)", sample_id)
    if match:
        year, sample_num = map(int, match.groups())
        return (year, sample_num)
    return (0, 0)

def load_sample_data(folder_path="../CNVkit"):
    """Load CNV data from .cnr files."""
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
    return sample_dataframes

def load_patient_data(excel_file='../results/TFx_hg38_1000kb.xlsx'):
    """Load patient timeline data."""
    df = pd.read_excel(excel_file, sheet_name=0, skiprows=1)
    df.columns = df.columns.map(str)
    df = df.loc[:, ~df.columns.str.startswith('Tumor_fraction_1000kb')]
    df.rename(columns=lambda x: 'Tumor_fraction' if x.startswith('Tumor_fraction_500kb') else x, inplace=True)
    
    patient_dict = {}
    for _, row in df.iterrows():
        patient_code = row["Patient_Code"]
        tuples = []
        for i in range(1, len(df.columns) // 2 + 1):
            if i == 1:
                sample_label_col = "1st progression"
            elif i == 2:
                sample_label_col = "2nd progression"
            elif i == 3:
                sample_label_col = "3rd progression"
            else:
                sample_label_col = "{}th progression".format(i)
            
            if sample_label_col not in df.columns:
                continue
            
            sample_label_col_index = df.columns.get_loc(sample_label_col)
            tumor_fraction_col_index = sample_label_col_index + 1
            sample_label = row[sample_label_col]
            tumor_fraction = row.iloc[tumor_fraction_col_index]
            
            if pd.notna(sample_label) and pd.notna(tumor_fraction):
                tuples.append((sample_label.split(" ")[0].split("_")[0], tumor_fraction))
        patient_dict[patient_code] = sorted(tuples, key=lambda x: parse_sample_id(x[0]))
    return patient_dict

def get_top_genes(sample_dataframes, n=100, gain_threshold=3, loss_threshold=2):
    """Get the top n most frequently gained and lost genes."""
    gain_counter = Counter()
    loss_counter = Counter()
    
    for df in sample_dataframes.values():
        gained_genes = set(df[df['copy_number'] >= gain_threshold]['gene_symbol'])
        lost_genes = set(df[df['copy_number'] < loss_threshold]['gene_symbol'])
        gain_counter.update(gained_genes)
        loss_counter.update(lost_genes)
    
    top_gained = [gene for gene, _ in gain_counter.most_common(n)]
    top_lost = [gene for gene, _ in loss_counter.most_common(n)]
    
    return top_gained, top_lost

def calculate_continuity_index(states):
    """Calculate continuity index for a sequence of states."""
    if len(states) < 2:
        return 1.0
    
    transitions = sum(1 for i in range(1, len(states)) if states[i] != states[i-1])
    max_possible_transitions = len(states) - 1
    
    continuity = 1 - (transitions / max_possible_transitions) if max_possible_transitions > 0 else 1.0
    return continuity

def analyze_gain_loss_continuity(sample_dataframes, patient_dict, top_gained, top_lost, gain_threshold=3, loss_threshold=2):
    """Calculate continuity indices for top gained and lost genes."""
    patient_continuity = {}
    
    for patient, samples in patient_dict.items():
        if len(samples) < 2:  # Skip patients with only one sample
            continue
        
        sample_dfs = []
        for sample, _ in samples:
            if sample in sample_dataframes:
                sample_dfs.append(sample_dataframes[sample])
        
        if not sample_dfs:
            continue
        
        # Calculate metrics for gains
        gain_transitions = 0
        gain_genes = 0
        for gene in top_gained:
            states = []
            for df in sample_dfs:
                if gene in df['gene_symbol'].values:
                    copy_number = df[df['gene_symbol'] == gene]['copy_number'].iloc[0]
                    states.append(1 if copy_number >= gain_threshold else 0)
                else:
                    states.append(0)
            
            if any(states):  # Only count genes that show gains at least once
                gain_genes += 1
                gain_transitions += sum(1 for i in range(1, len(states)) 
                                     if states[i] != states[i-1])
        
        # Calculate metrics for losses
        loss_transitions = 0
        loss_genes = 0
        for gene in top_lost:
            states = []
            for df in sample_dfs:
                if gene in df['gene_symbol'].values:
                    copy_number = df[df['gene_symbol'] == gene]['copy_number'].iloc[0]
                    states.append(1 if copy_number < loss_threshold else 0)
                else:
                    states.append(0)
            
            if any(states):  # Only count genes that show losses at least once
                loss_genes += 1
                loss_transitions += sum(1 for i in range(1, len(states)) 
                                     if states[i] != states[i-1])
        
        # Calculate continuity indices
        gain_max_transitions = (len(samples) - 1) * gain_genes if gain_genes > 0 else 1
        loss_max_transitions = (len(samples) - 1) * loss_genes if loss_genes > 0 else 1
        
        gain_continuity = 1 - (gain_transitions / gain_max_transitions)
        loss_continuity = 1 - (loss_transitions / loss_max_transitions)
        
        patient_continuity[patient] = {
            'gain_continuity': gain_continuity,
            'loss_continuity': loss_continuity,
            'gain_gene_count': gain_genes,
            'loss_gene_count': loss_genes,
            'sample_count': len(samples),
            'samples': [s[0] for s in samples]
        }
    
    return patient_continuity

def plot_gain_loss_continuity(patient_continuity, output_prefix=''):
    """Create visualizations for gain and loss continuity metrics."""
    # Prepare data for plotting
    patients = []
    gain_continuity = []
    loss_continuity = []
    sample_counts = []
    gain_genes = []
    loss_genes = []
    
    for patient, metrics in patient_continuity.items():
        patients.append(patient)
        gain_continuity.append(metrics['gain_continuity'])
        loss_continuity.append(metrics['loss_continuity'])
        sample_counts.append(metrics['sample_count'])
        gain_genes.append(metrics['gain_gene_count'])
        loss_genes.append(metrics['loss_gene_count'])
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12))
    
    # Sort indices for gain plot
    gain_sorted = np.argsort(gain_continuity)
    gain_patients = [patients[i] for i in gain_sorted]
    gain_values = [gain_continuity[i] for i in gain_sorted]
    gain_samples = [sample_counts[i] for i in gain_sorted]
    gain_gene_counts = [gain_genes[i] for i in gain_sorted]
    
    # Plot 1: Gain Continuity
    bars1 = ax1.bar(range(len(gain_patients)), gain_values)
    ax1.set_xlabel('Patients')
    ax1.set_ylabel('Gain Continuity Index')
    ax1.set_title('Patient Sample Continuity Index for Top 100 Gains\n(1.0 = perfectly continuous, 0.0 = maximum variation)')
    ax1.set_xticks(range(len(gain_patients)))
    ax1.set_xticklabels(gain_patients, rotation=45, ha='right')
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels for gains
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                '{:.2f}\n(n={}, g={})'.format(height, gain_samples[i], gain_gene_counts[i]),
                ha='center', va='bottom')
    
    # Sort indices for loss plot
    loss_sorted = np.argsort(loss_continuity)
    loss_patients = [patients[i] for i in loss_sorted]
    loss_values = [loss_continuity[i] for i in loss_sorted]
    loss_samples = [sample_counts[i] for i in loss_sorted]
    loss_gene_counts = [loss_genes[i] for i in loss_sorted]
    
    # Plot 2: Loss Continuity
    bars2 = ax2.bar(range(len(loss_patients)), loss_values)
    ax2.set_xlabel('Patients')
    ax2.set_ylabel('Loss Continuity Index')
    ax2.set_title('Patient Sample Continuity Index for Top 100 Losses\n(1.0 = perfectly continuous, 0.0 = maximum variation)')
    ax2.set_xticks(range(len(loss_patients)))
    ax2.set_xticklabels(loss_patients, rotation=45, ha='right')
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels for losses
    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                '{:.2f}\n(n={}, l={})'.format(height, loss_samples[i], loss_gene_counts[i]),
                ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('{}top100_gain_loss_continuity.png'.format(output_prefix), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Load data
    print("Loading sample data...")
    sample_dataframes = load_sample_data()
    
    print("Loading patient data...")
    patient_dict = load_patient_data()
    
    # Get top 100 gained and lost genes
    print("Finding top 100 gained and lost genes...")
    top_gained, top_lost = get_top_genes(sample_dataframes, n=100)
    
    # Calculate continuity
    print("Analyzing gain/loss continuity...")
    patient_continuity = analyze_gain_loss_continuity(sample_dataframes, patient_dict, top_gained, top_lost)
    
    # Save detailed results
    results = []
    for patient, metrics in patient_continuity.items():
        results.append({
            'Patient': patient,
            'Gain_Continuity': metrics['gain_continuity'],
            'Loss_Continuity': metrics['loss_continuity'],
            'Gain_Gene_Count': metrics['gain_gene_count'],
            'Loss_Gene_Count': metrics['loss_gene_count'],
            'Sample_Count': metrics['sample_count'],
            'Samples': ', '.join(metrics['samples'])
        })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv('top100_gain_loss_continuity_metrics.csv', index=False)
    
    # Plot results
    plot_gain_loss_continuity(patient_continuity)
    
    # Print summary statistics
    print("\nContinuity Summary (Top 100 genes):")
    print("Average Gain Continuity: {:.3f}".format(results_df['Gain_Continuity'].mean()))
    print("Average Loss Continuity: {:.3f}".format(results_df['Loss_Continuity'].mean()))
    print("Average Gained Genes: {:.1f}".format(results_df['Gain_Gene_Count'].mean()))
    print("Average Lost Genes: {:.1f}".format(results_df['Loss_Gene_Count'].mean()))
    
    print("\nTop 5 most continuous patients (Gains):")
    print(results_df.nlargest(5, 'Gain_Continuity')[['Patient', 'Gain_Continuity', 'Gain_Gene_Count', 'Sample_Count']])
    
    print("\nTop 5 most continuous patients (Losses):")
    print(results_df.nlargest(5, 'Loss_Continuity')[['Patient', 'Loss_Continuity', 'Loss_Gene_Count', 'Sample_Count']])

if __name__ == "__main__":
    main() 