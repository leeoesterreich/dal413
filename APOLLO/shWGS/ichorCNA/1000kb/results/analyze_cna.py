import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Create output directory
output_dir = 'tumor_fraction/cna_analysis'
os.makedirs(output_dir, exist_ok=True)

# Read the data
cnv_data = pd.read_csv('gene_cnv_summary.txt', sep='\t')
summary_data = pd.read_csv('summary.csv')
blood_draw = pd.read_csv('Blood_draw.csv')

# Clean up sample names
summary_data['Sample'] = summary_data['Sample'].str.replace('.bam', '')

# Function to get patient ID from sample name
def get_patient_id(sample):
    for _, row in blood_draw.iterrows():
        if sample in row.values:
            return row['PT ID']
    return None

# Add patient ID to summary data
summary_data['Patient_ID'] = summary_data['Sample'].apply(get_patient_id)

# Process CNV data
def analyze_cna(cnv_data, threshold_amp=4, threshold_del=2):
    # Get amplified and deleted genes per sample
    amplified = cnv_data[cnv_data['CN'] >= threshold_amp].groupby(['Sample', 'Gene']).first().reset_index()
    deleted = cnv_data[cnv_data['CN'] < threshold_del].groupby(['Sample', 'Gene']).first().reset_index()
    
    # Count frequency of amplification/deletion per gene
    amp_freq = amplified['Gene'].value_counts()
    del_freq = deleted['Gene'].value_counts()
    
    return amp_freq, del_freq, amplified, deleted

# Analyze CNAs
amp_freq, del_freq, amplified_data, deleted_data = analyze_cna(cnv_data)

# Save top amplified and deleted genes
top_n = 50  # Number of top genes to save

with open(os.path.join(output_dir, 'top_amplified_genes.txt'), 'w') as f:
    f.write('Gene\tFrequency\n')
    for gene, freq in amp_freq.head(top_n).items():
        f.write(f'{gene}\t{freq}\n')

with open(os.path.join(output_dir, 'top_deleted_genes.txt'), 'w') as f:
    f.write('Gene\tFrequency\n')
    for gene, freq in del_freq.head(top_n).items():
        f.write(f'{gene}\t{freq}\n')

# Create oncoplot data
def prepare_oncoplot_data(amp_data, del_data, top_n=20):
    # Get top genes
    top_amp_genes = amp_freq.head(top_n).index.tolist()
    top_del_genes = del_freq.head(top_n).index.tolist()
    
    # Prepare matrix for amplifications
    amp_matrix = pd.DataFrame(0, 
                            index=cnv_data['Sample'].unique(),
                            columns=top_amp_genes)
    for _, row in amp_data.iterrows():
        if row['Gene'] in top_amp_genes:
            amp_matrix.loc[row['Sample'], row['Gene']] = 1
            
    # Prepare matrix for deletions
    del_matrix = pd.DataFrame(0,
                            index=cnv_data['Sample'].unique(),
                            columns=top_del_genes)
    for _, row in del_data.iterrows():
        if row['Gene'] in top_del_genes:
            del_matrix.loc[row['Sample'], row['Gene']] = -1
            
    return amp_matrix, del_matrix

# Create oncoplot
amp_matrix, del_matrix = prepare_oncoplot_data(amplified_data, deleted_data)

# Plot function
def plot_oncoplot(amp_matrix, del_matrix, summary_data, output_dir):
    # Combine amplification and deletion matrices
    combined_matrix = pd.concat([amp_matrix, del_matrix], axis=1)
    
    # Sort samples by tumor fraction
    tf_dict = dict(zip(summary_data['Sample'], summary_data['Tumor Fraction']))
    combined_matrix['TF'] = combined_matrix.index.map(tf_dict)
    combined_matrix = combined_matrix.sort_values('TF', ascending=False)
    combined_matrix = combined_matrix.drop('TF', axis=1)
    
    # Create the plot
    plt.figure(figsize=(20, 10))
    
    # Plot heatmap
    sns.heatmap(combined_matrix, 
                cmap='RdBu_r',
                center=0,
                vmin=-1,
                vmax=1,
                yticklabels=False)
    
    plt.title('CNA Profile Across Samples')
    plt.xlabel('Genes')
    plt.xticks(rotation=45, ha='right')
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'oncoplot.png'), dpi=300, bbox_inches='tight')
    plt.close()

# Create plots
plot_oncoplot(amp_matrix, del_matrix, summary_data, output_dir)

# Perform statistical analysis
def calculate_odds_ratios(amp_matrix, del_matrix, summary_data, threshold_tf=0.1):
    # Combine matrices
    combined_matrix = pd.concat([amp_matrix, del_matrix], axis=1)
    
    # Add tumor fraction information
    combined_matrix['TF'] = combined_matrix.index.map(
        dict(zip(summary_data['Sample'], summary_data['Tumor Fraction']))
    )
    
    # Calculate odds ratios
    results = []
    for gene in combined_matrix.columns[:-1]:  # Exclude TF column
        # Create contingency table
        table = pd.crosstab(
            combined_matrix['TF'] > threshold_tf,
            combined_matrix[gene] != 0
        )
        
        # Calculate odds ratio and p-value
        odds_ratio, p_value = stats.fisher_exact(table)
        
        results.append({
            'Gene': gene,
            'Odds_Ratio': odds_ratio,
            'P_Value': p_value,
            'Type': 'Amplification' if gene in amp_matrix.columns else 'Deletion'
        })
    
    return pd.DataFrame(results)

# Calculate and save odds ratios
odds_ratios = calculate_odds_ratios(amp_matrix, del_matrix, summary_data)
odds_ratios.to_csv(os.path.join(output_dir, 'odds_ratios.csv'), index=False)

# Add -log10(P_Value) column for plotting
odds_ratios['log10_p'] = -np.log10(odds_ratios['P_Value'])

# Plot odds ratios
plt.figure(figsize=(15, 8))
sns.scatterplot(data=odds_ratios, 
                x='Odds_Ratio',
                y='log10_p',
                hue='Type',
                alpha=0.6)

plt.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
plt.title('Odds Ratios vs Statistical Significance')
plt.ylabel('-log10(P-Value)')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'odds_ratios_plot.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Analysis complete. Results saved in:", output_dir)

# Read the gene CNV summary file
gene_cnv = pd.read_csv('results/gene_cnv_summary.txt', sep='\t')

# Read tumor fraction data
tumor_fraction = pd.read_csv('results/tumor_fraction.txt', sep='\t')

# Process data to ensure chronological order and group by patient
def get_patient_id(sample):
    return sample.split('_')[0]

def get_timepoint(sample):
    return int(''.join(filter(str.isdigit, sample.split('_')[1])))

# Group samples by patient and sort chronologically within each patient
samples_by_patient = {}
for sample in gene_cnv.columns[1:]:
    patient = get_patient_id(sample)
    if patient not in samples_by_patient:
        samples_by_patient[patient] = []
    samples_by_patient[patient].append(sample)

# Sort samples within each patient
for patient in samples_by_patient:
    samples_by_patient[patient].sort(key=get_timepoint)

# Create ordered list of all samples
all_samples = []
patient_boundaries = []
current_pos = 0

# Sort patients by their ID number
sorted_patients = sorted(samples_by_patient.keys(), 
                       key=lambda x: int(x.split('_')[1]))

for patient in sorted_patients:
    all_samples.extend(samples_by_patient[patient])
    current_pos += len(samples_by_patient[patient])
    if current_pos < len(gene_cnv.columns[1:]):
        patient_boundaries.append(current_pos - 0.5)

# Create the figure with specific size and layout
fig = plt.figure(figsize=(20, 12))
gs = plt.GridSpec(2, 1, height_ratios=[1, 4])

# Top subplot for tumor fraction
ax1 = fig.add_subplot(gs[0])
tumor_fraction_sorted = tumor_fraction.set_index('Sample').loc[all_samples]
sns.barplot(data=tumor_fraction_sorted.reset_index(), 
            x='Sample', y='TumorFraction',
            color='gray', ax=ax1)

# Add patient boundary lines to tumor fraction plot
for boundary in patient_boundaries:
    ax1.axvline(x=boundary, color='black', linestyle='-', linewidth=1)

ax1.set_ylim(0, 0.6)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
ax1.set_title('Tumor Fraction')
ax1.set_xlabel('')

# Bottom subplot for gene CNAs
ax2 = fig.add_subplot(gs[1])
gene_cnv_sorted = gene_cnv.set_index('Gene')[all_samples]

# Create heatmap with binary color scheme
sns.heatmap(gene_cnv_sorted, 
            cmap=['white', '#8B4513'],  # Using brown color for CNAs
            center=1.5,
            vmin=1,
            vmax=2,
            yticklabels=True,
            xticklabels=True,
            ax=ax2)

# Add patient boundary lines to heatmap
for boundary in patient_boundaries:
    ax2.axvline(x=boundary, color='black', linestyle='-', linewidth=1)

ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90)
ax2.set_title('Copy Number Alterations')

# Add patient labels at the bottom
ax2.set_xlabel('Patients')

# Adjust layout
plt.tight_layout()
plt.savefig('results/oncoplot.pdf', bbox_inches='tight', dpi=300)
plt.close() 