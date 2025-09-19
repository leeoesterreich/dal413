# -*- coding: utf-8 -*-
import pandas as pd

# Read the data files
histology_file = 'Guardant Project - Histology Data(IDC + ILC ONLY).csv'
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'

# Read histology data, skipping the first two rows
hist_df = pd.read_csv(histology_file, encoding='latin-1', skiprows=2)
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# Get the column names
gh_id_col = 'GH_ID'  # Column name in histology after skipping headers
genomic_gh_id_col = 'ï»¿GH_ID'  # Column name in genomic data

# Get sets of patient IDs
histology_ids = set(hist_df[gh_id_col])
genomic_ids = set(genomic_df[genomic_gh_id_col])

# Find patients with only genomic data
genomic_only_ids = genomic_ids - histology_ids

# Filter the genomic dataframe
genomic_only_df = genomic_df[genomic_df[genomic_gh_id_col].isin(genomic_only_ids)]

# Sort by GH_ID
genomic_only_df = genomic_only_df.sort_values(by=genomic_gh_id_col)

# Save to CSV
output_file = 'genomic_without_histology.csv'
genomic_only_df.to_csv(output_file, index=False)

# Print summary
print("\nSummary of patients with genomic data but no histology:")
print("---------------------------------------------------")
print(f"Total patients in genomic data: {len(genomic_ids)}")
print(f"Patients with only genomic data: {len(genomic_only_ids)} ({len(genomic_only_ids)/len(genomic_ids)*100:.1f}%)")
print(f"Total rows in output file: {len(genomic_only_df)}")
print(f"\nOutput saved to: {output_file}")

# Print sample of the data
print("\nSample of first few patients:")
print("---------------------------")
print(genomic_only_df[genomic_gh_id_col].head(10)) 