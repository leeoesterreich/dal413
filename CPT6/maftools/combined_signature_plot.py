#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import matplotlib
# Set the backend to a non-interactive one
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Path to the output folders
cpt6_output = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/CPT6_SigProfilerAssignmentOutput/filtered"
icle_output = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/CPT6_SigProfilerAssignmentOutput/ICLE"

# File paths for the activity files
cpt6_file_path = os.path.join(cpt6_output, 'Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
icle_file_path = os.path.join(icle_output, 'Assignment_Solution/Activities/Assignment_Solution_Activities.txt')

# Read the files into DataFrames
cpt6_df = pd.read_csv(cpt6_file_path, sep='\t')
icle_df = pd.read_csv(icle_file_path, sep='\t')

# Set 'Samples' column as the index
cpt6_df.set_index('Samples', inplace=True)
icle_df.set_index('Samples', inplace=True)

# Rename the CPT6 sample for clarity
cpt6_df.index = ['CPT6 (mm10)']

# Calculate the mean of ICLE samples
icle_mean = pd.DataFrame(icle_df.mean(axis=0)).T
icle_mean.index = ['ICLE (GRCh38)']

# Combine the two DataFrames
combined_df = pd.concat([cpt6_df, icle_mean])

# Drop columns where all values are zero
combined_df = combined_df.loc[:, (combined_df != 0).any(axis=0)]

# Sort columns by total contribution
column_sums = combined_df.sum()
sorted_columns = column_sums.sort_values(ascending=False).index.tolist()
combined_df = combined_df[sorted_columns]

# Keep only top 10 signatures for clarity
if len(combined_df.columns) > 10:
    top_signatures = sorted_columns[:10]
    other_signatures = sorted_columns[10:]
    
    # Sum the remaining signatures into an "Other" category
    combined_df['Other'] = combined_df[other_signatures].sum(axis=1)
    combined_df = combined_df[top_signatures + ['Other']]

# Calculate percentages for each row
combined_df_percent = combined_df.div(combined_df.sum(axis=1), axis=0) * 100

# Create a figure with two subplots - one for counts and one for percentages
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 16))

# Plot raw counts
combined_df.plot(kind='bar', stacked=True, ax=ax1, colormap='viridis')
ax1.set_ylabel('Mutation Counts', fontsize=18)
ax1.set_title('SBS Signature Counts Comparison', fontsize=20)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

# Plot percentages
combined_df_percent.plot(kind='bar', stacked=True, ax=ax2, colormap='viridis')
ax2.set_ylabel('Percentage (%)', fontsize=18)
ax2.set_title('SBS Signature Percentage Comparison', fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

# Rotate x-axis labels for better readability
plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')

# Add a note about different genome builds
fig.text(0.5, 0.01, 'Note: CPT6 uses mm10 (mouse) genome build, while ICLE uses GRCh38 (human) genome build.\n'
         'Direct comparison should be interpreted with caution.', 
         ha='center', fontsize=12, bbox=dict(facecolor='yellow', alpha=0.2))

plt.tight_layout(rect=[0, 0.03, 0.85, 0.98])  # Adjust layout to make room for the legend and note

# Create output directory if it doesn't exist
output_dir = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/ICLE/maftools/output_plots"
os.makedirs(output_dir, exist_ok=True)

# Save the plot in multiple formats
plt.savefig(os.path.join(output_dir, "combined_signatures.pdf"), format='pdf')
plt.savefig(os.path.join(output_dir, "combined_signatures.svg"), format='svg')
plt.savefig(os.path.join(output_dir, "combined_signatures.png"), format='png', dpi=300)

print(f"Combined signature plots saved to {output_dir}")

# Create a more detailed version with individual ICLE samples
# This will show CPT6 and all ICLE samples

# Create a new DataFrame with CPT6 and all ICLE samples
detailed_df = pd.concat([cpt6_df, icle_df])

# Rename the CPT6 sample
detailed_df = detailed_df.rename(index={'filtered_tumor_vs_liver': 'CPT6 (mm10)'})

# Add genome build information to ICLE samples
detailed_df.index = [f"{idx} (GRCh38)" if idx != "CPT6 (mm10)" else idx for idx in detailed_df.index]

# Drop columns where all values are zero
detailed_df = detailed_df.loc[:, (detailed_df != 0).any(axis=0)]

# Sort columns by total contribution
column_sums = detailed_df.sum()
sorted_columns = column_sums.sort_values(ascending=False).index.tolist()
detailed_df = detailed_df[sorted_columns]

# Keep only top 10 signatures for clarity
if len(detailed_df.columns) > 10:
    top_signatures = sorted_columns[:10]
    other_signatures = sorted_columns[10:]
    
    # Sum the remaining signatures into an "Other" category
    detailed_df['Other'] = detailed_df[other_signatures].sum(axis=1)
    detailed_df = detailed_df[top_signatures + ['Other']]

# Calculate percentages for each row
detailed_df_percent = detailed_df.div(detailed_df.sum(axis=1), axis=0) * 100

# Create a figure with two subplots for the detailed view
fig2, (ax3, ax4) = plt.subplots(2, 1, figsize=(15, 20))

# Plot raw counts for detailed view
detailed_df.plot(kind='bar', stacked=True, ax=ax3, colormap='viridis')
ax3.set_ylabel('Mutation Counts', fontsize=18)
ax3.set_title('SBS Signature Counts - All Samples', fontsize=20)
ax3.tick_params(axis='both', which='major', labelsize=12)
ax3.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

# Plot percentages for detailed view
detailed_df_percent.plot(kind='bar', stacked=True, ax=ax4, colormap='viridis')
ax4.set_ylabel('Percentage (%)', fontsize=18)
ax4.set_title('SBS Signature Percentage - All Samples', fontsize=20)
ax4.tick_params(axis='both', which='major', labelsize=12)
ax4.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

# Rotate x-axis labels for better readability
plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
plt.setp(ax4.get_xticklabels(), rotation=45, ha='right')

# Add a note about different genome builds
fig2.text(0.5, 0.01, 'Note: CPT6 uses mm10 (mouse) genome build, while ICLE samples use GRCh38 (human) genome build.\n'
         'Direct comparison should be interpreted with caution.', 
         ha='center', fontsize=12, bbox=dict(facecolor='yellow', alpha=0.2))

plt.tight_layout(rect=[0, 0.03, 0.85, 0.98])  # Adjust layout

# Save the detailed plot
plt.savefig(os.path.join(output_dir, "detailed_combined_signatures.pdf"), format='pdf')
plt.savefig(os.path.join(output_dir, "detailed_combined_signatures.svg"), format='svg')
plt.savefig(os.path.join(output_dir, "detailed_combined_signatures.png"), format='png', dpi=300)

print(f"Detailed signature plots saved to {output_dir}")

# Create a heatmap visualization for better comparison
plt.figure(figsize=(14, 10))
sns.heatmap(detailed_df_percent, cmap="viridis", annot=True, fmt=".1f", linewidths=.5)
plt.title('SBS Signature Percentage Heatmap', fontsize=20)
plt.ylabel('Samples', fontsize=18)
plt.xlabel('SBS Signatures', fontsize=18)
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(fontsize=12)

# Save the heatmap
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "signature_heatmap.pdf"), format='pdf')
plt.savefig(os.path.join(output_dir, "signature_heatmap.svg"), format='svg')
plt.savefig(os.path.join(output_dir, "signature_heatmap.png"), format='png', dpi=300)

print(f"Signature heatmap saved to {output_dir}") 