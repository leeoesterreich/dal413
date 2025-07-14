# -*- coding: utf-8 -*-

import pandas as pd

# Read the cohort files
idc_df = pd.read_csv('idc_cohort.csv')
ilc_df = pd.read_csv('ilc_cohort.csv')
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# Get column names and IDs
gh_id_col = 'Can delete duplicates of GH_ID'
genomic_gh_id_col = 'ï»¿GH_ID'

# Get sets of patient IDs
idc_ids = set(idc_df[gh_id_col])
ilc_ids = set(ilc_df[gh_id_col])
genomic_ids = set(genomic_df[genomic_gh_id_col])

# Find patients without genomic data
idc_without_genomic = idc_ids - genomic_ids
ilc_without_genomic = ilc_ids - genomic_ids

# Filter the original dataframes
idc_without_genomic_df = idc_df[idc_df[gh_id_col].isin(idc_without_genomic)]
ilc_without_genomic_df = ilc_df[ilc_df[gh_id_col].isin(ilc_without_genomic)]

# Add a column to indicate histology type
idc_without_genomic_df['Histology_Type'] = 'IDC'
ilc_without_genomic_df['Histology_Type'] = 'ILC'

# Combine the dataframes
combined_df = pd.concat([idc_without_genomic_df, ilc_without_genomic_df], ignore_index=True)

# Sort by GH_ID
combined_df = combined_df.sort_values(by=gh_id_col)

# Save to CSV
output_file = 'histology_without_genomic.csv'
combined_df.to_csv(output_file, index=False)

# Print summary
print("\nSummary of patients with histology but no genomic data:")
print("-----------------------------------------------------")
print(f"Total patients: {len(combined_df)}")
print(f"IDC patients: {len(idc_without_genomic)} ({len(idc_without_genomic)/len(combined_df)*100:.1f}%)")
print(f"ILC patients: {len(ilc_without_genomic)} ({len(ilc_without_genomic)/len(combined_df)*100:.1f}%)")
print(f"\nOutput saved to: {output_file}") 