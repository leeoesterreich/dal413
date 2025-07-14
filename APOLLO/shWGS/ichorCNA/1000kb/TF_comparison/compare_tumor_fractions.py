import pandas as pd
import os
import glob
import re

# Paths
daisong_summary = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv"
rahul_base_dir = "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/ILC/ichorCNA/hg38_1000kb/results/ichorCNA"

def get_rahul_tumor_fraction(sample_id):
    params_file = os.path.join(rahul_base_dir, sample_id, f"{sample_id}.params.txt")
    if not os.path.exists(params_file):
        return None
    
    with open(params_file, 'r') as f:
        for line in f:
            if line.startswith("Tumor Fraction:"):
                return float(line.strip().split("\t")[1])
    return None

# Read Daisong's data
daisong_df = pd.read_csv(daisong_summary)
daisong_df['Sample'] = daisong_df['Sample'].str.replace('.bam', '')

# Create comparison dataframe
comparison_data = []
for sample in daisong_df['Sample']:
    daisong_tf = daisong_df[daisong_df['Sample'] == sample]['Tumor Fraction'].iloc[0]
    rahul_tf = get_rahul_tumor_fraction(sample)
    
    if rahul_tf is not None:
        comparison_data.append({
            'Sample': sample,
            'Daisong_TF': daisong_tf,
            'Rahul_TF': rahul_tf,
            'Absolute_Difference': abs(daisong_tf - rahul_tf),
            'Relative_Difference_Percent': abs(daisong_tf - rahul_tf) / ((daisong_tf + rahul_tf) / 2) * 100 if (daisong_tf + rahul_tf) > 0 else 0
        })

# Create comparison DataFrame and sort by absolute difference
comparison_df = pd.DataFrame(comparison_data)
comparison_df = comparison_df.sort_values('Absolute_Difference', ascending=False)

# Save results
output_file = "tumor_fraction_comparison.csv"
comparison_df.to_csv(output_file, index=False)

# Print summary statistics
print("\nComparison Summary:")
print(f"Total samples compared: {len(comparison_df)}")
print("\nSamples with largest differences:")
print(comparison_df.head().to_string(index=False))

print("\nOverall Statistics:")
print(f"Mean absolute difference: {comparison_df['Absolute_Difference'].mean():.4f}")
print(f"Median absolute difference: {comparison_df['Absolute_Difference'].median():.4f}")
print(f"Mean relative difference: {comparison_df['Relative_Difference_Percent'].mean():.2f}%") 