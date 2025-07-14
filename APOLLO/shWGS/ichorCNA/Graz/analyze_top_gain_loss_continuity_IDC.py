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

def load_sample_data(folder_path="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Tendo/Annotated_cnr_file_500Kb"):
    """Load CNV data from .cnr files."""
    sample_dataframes = {}
    print("\nLoading .cnr files:")
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".cnr"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing {file_name}...")
            df = pd.read_csv(file_path, sep="\t")
            sample_id = [col for col in df.columns if col.startswith("TP") and "Corrected_Copy_Number" in col]
            if sample_id:
                sample_id = sample_id[0].split(".")[0]
                # Remove any suffix after underscore
                base_sample_id = sample_id.split("_")[0]
                df = df[["gene", sample_id + ".Corrected_Copy_Number"]]
                df.columns = ["gene_symbol", "copy_number"]
                df = df.assign(gene_symbol=df["gene_symbol"].str.split(",")).explode("gene_symbol").reset_index(drop=True)
                df = df[df["gene_symbol"] != "-"]
                df["copy_number"] = pd.to_numeric(df["copy_number"], errors="coerce")
                df = df.groupby("gene_symbol", as_index=False)["copy_number"].max()
                sample_dataframes[base_sample_id] = df
                print(f"  Loaded data for sample {base_sample_id}")
            else:
                print(f"  No sample ID found in {file_name}")
    return sample_dataframes

def load_patient_data(csv_file='/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Tendo/Clinical_information/TFx_hg38_1000_500kb.csv'):
    """Load patient timeline data from CSV."""
    df = pd.read_csv(csv_file, skiprows=1)  # Skip the first row with genome version
    
    patient_dict = {}
    for _, row in df.iterrows():
        patient_code = row["Patient_Code"]
        tuples = []
        
        # Process each progression point (1st through 6th)
        for i in range(1, 7):
            prog_suffix = {1: "st", 2: "nd", 3: "rd"}.get(i, "th")
            prog_col = f"{i}{prog_suffix} progression"
            tf_col = f"Tumor_fraction_1000kb bin.{i}"  # Column for the corresponding tumor fraction
            
            if prog_col in df.columns:
                sample_label = row[prog_col]
                if pd.notna(sample_label):
                    # Extract just the sample ID without the parenthetical number
                    sample_id = sample_label.split(" ")[0].strip()
                    # Get the tumor fraction from the third column after the progression column
                    tf_value = row.iloc[row.index.get_loc(prog_col) + 1]
                    if pd.notna(tf_value):
                        tuples.append((sample_id, tf_value))
        
        if tuples:  # Only add patients with valid samples
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
    
    print("\nAnalyzing patients:")
    for patient, samples in patient_dict.items():
        print(f"\nPatient {patient}:")
        print(f"  Samples: {samples}")
        
        if len(samples) < 2:  # Skip patients with only one sample
            print("  Skipped: Only one sample")
            continue
        
        sample_dfs = []
        for sample, _ in samples:
            if sample in sample_dataframes:
                sample_dfs.append(sample_dataframes[sample])
            else:
                print(f"  Warning: No data found for sample {sample}")
        
        if len(sample_dfs) < 2:  # Need at least 2 samples with data
            print("  Skipped: Not enough samples with data")
            continue
        
        print(f"  Processing {len(sample_dfs)} samples with data")
        
        # Calculate metrics for gains
        gain_transitions = 0
        gain_genes = 0
        gain_gene_details = []
        for gene in top_gained:
            states = []
            for df in sample_dfs:
                if gene in df['gene_symbol'].values:
                    copy_number = df[df['gene_symbol'] == gene]['copy_number'].iloc[0]
                    states.append(1 if copy_number >= gain_threshold else 0)
                else:
                    states.append(0)
            
            # Only count genes that show gains in at least 25% of samples
            if sum(states) >= len(states) * 0.25:
                gain_genes += 1
                # Count transitions between different states
                transitions = sum(1 for i in range(1, len(states)) if states[i] != states[i-1])
                gain_transitions += transitions
                # Store details for this gene
                continuity = 1 - (transitions / (len(states) - 1)) if len(states) > 1 else 1
                gain_gene_details.append((gene, continuity, sum(states)))
        
        # Calculate metrics for losses
        loss_transitions = 0
        loss_genes = 0
        loss_gene_details = []
        for gene in top_lost:
            states = []
            for df in sample_dfs:
                if gene in df['gene_symbol'].values:
                    copy_number = df[df['gene_symbol'] == gene]['copy_number'].iloc[0]
                    states.append(1 if copy_number < loss_threshold else 0)
                else:
                    states.append(0)
            
            # Only count genes that show losses in at least 25% of samples
            if sum(states) >= len(states) * 0.25:
                loss_genes += 1
                # Count transitions between different states
                transitions = sum(1 for i in range(1, len(states)) if states[i] != states[i-1])
                loss_transitions += transitions
                # Store details for this gene
                continuity = 1 - (transitions / (len(states) - 1)) if len(states) > 1 else 1
                loss_gene_details.append((gene, continuity, sum(states)))
        
        # Calculate continuity indices
        gain_max_transitions = (len(sample_dfs) - 1) * gain_genes if gain_genes > 0 else 1
        loss_max_transitions = (len(sample_dfs) - 1) * loss_genes if loss_genes > 0 else 1
        
        gain_continuity = 1 - (gain_transitions / gain_max_transitions) if gain_genes > 0 else 0
        loss_continuity = 1 - (loss_transitions / loss_max_transitions) if loss_genes > 0 else 0
        
        # Sort gene details by continuity
        gain_gene_details.sort(key=lambda x: (-x[1], -x[2]))  # Sort by continuity (desc) then occurrence (desc)
        loss_gene_details.sort(key=lambda x: (-x[1], -x[2]))
        
        print(f"  Results: {gain_genes} gained genes, {loss_genes} lost genes")
        print(f"  Continuity: Gain={gain_continuity:.3f}, Loss={loss_continuity:.3f}")
        
        if gain_genes > 0:
            print("  Top 5 most continuous gained genes:")
            for gene, cont, occur in gain_gene_details[:5]:
                print(f"    {gene}: continuity={cont:.3f}, occurrences={occur}/{len(sample_dfs)}")
        
        if loss_genes > 0:
            print("  Top 5 most continuous lost genes:")
            for gene, cont, occur in loss_gene_details[:5]:
                print(f"    {gene}: continuity={cont:.3f}, occurrences={occur}/{len(sample_dfs)}")
        
        patient_continuity[patient] = {
            'gain_continuity': gain_continuity,
            'loss_continuity': loss_continuity,
            'gain_gene_count': gain_genes,
            'loss_gene_count': loss_genes,
            'sample_count': len(samples),
            'samples': [s[0] for s in samples],
            'gain_gene_details': gain_gene_details[:10],  # Store top 10 genes
            'loss_gene_details': loss_gene_details[:10]
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
    ax1.set_title('IDC Patient Sample Continuity Index for Top 100 Gains\n(1.0 = perfectly continuous, 0.0 = maximum variation)')
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
    ax2.set_title('IDC Patient Sample Continuity Index for Top 100 Losses\n(1.0 = perfectly continuous, 0.0 = maximum variation)')
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
    plt.savefig('{}IDC_top100_gain_loss_continuity.png'.format(output_prefix), dpi=300, bbox_inches='tight')
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
    gene_results = []
    for patient, metrics in patient_continuity.items():
        results.append({
            'Patient': patient,
            'gain_continuity': metrics['gain_continuity'],
            'loss_continuity': metrics['loss_continuity'],
            'gain_gene_count': metrics['gain_gene_count'],
            'loss_gene_count': metrics['loss_gene_count'],
            'sample_count': metrics['sample_count'],
            'samples': ', '.join(metrics['samples'])
        })
        
        # Add gene-level details
        for gene, cont, occur in metrics['gain_gene_details']:
            gene_results.append({
                'Patient': patient,
                'Gene': gene,
                'Type': 'Gain',
                'Continuity': cont,
                'Occurrences': occur,
                'Sample_Count': metrics['sample_count']
            })
        
        for gene, cont, occur in metrics['loss_gene_details']:
            gene_results.append({
                'Patient': patient,
                'Gene': gene,
                'Type': 'Loss',
                'Continuity': cont,
                'Occurrences': occur,
                'Sample_Count': metrics['sample_count']
            })
    
    results_df = pd.DataFrame(results)
    gene_results_df = pd.DataFrame(gene_results)
    
    # Save results
    results_df.to_csv('IDC_top100_gain_loss_continuity_metrics.csv', index=False)
    gene_results_df.to_csv('IDC_top100_gain_loss_gene_details.csv', index=False)
    
    # Plot results
    plot_gain_loss_continuity(patient_continuity)
    
    # Print summary statistics
    print("\nContinuity Summary for IDC Samples (Top 100 genes):")
    print("Average Gain Continuity: {:.3f}".format(results_df['gain_continuity'].mean()))
    print("Average Loss Continuity: {:.3f}".format(results_df['loss_continuity'].mean()))
    print("Average Gained Genes: {:.1f}".format(results_df['gain_gene_count'].mean()))
    print("Average Lost Genes: {:.1f}".format(results_df['loss_gene_count'].mean()))
    
    print("\nTop 5 most continuous patients (Gains):")
    print(results_df.nlargest(5, 'gain_continuity')[['Patient', 'gain_continuity', 'gain_gene_count', 'sample_count']])
    
    print("\nTop 5 most continuous patients (Losses):")
    print(results_df.nlargest(5, 'loss_continuity')[['Patient', 'loss_continuity', 'loss_gene_count', 'sample_count']])
    
    # Print gene-level summary
    print("\nMost continuous genes across all patients:")
    for gene_type in ['Gain', 'Loss']:
        type_genes = gene_results_df[gene_results_df['Type'] == gene_type]
        avg_continuity = type_genes.groupby('Gene')['Continuity'].agg(['mean', 'count']).reset_index()
        avg_continuity = avg_continuity[avg_continuity['count'] >= 3]  # Only show genes found in at least 3 patients
        avg_continuity = avg_continuity.sort_values('mean', ascending=False)
        
        print(f"\nTop 10 most continuous {gene_type.lower()}ed genes (present in ≥3 patients):")
        print(avg_continuity.head(10)[['Gene', 'mean', 'count']].rename(
            columns={'mean': 'Avg_Continuity', 'count': 'Patient_Count'}))

if __name__ == "__main__":
    main() 