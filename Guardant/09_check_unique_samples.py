# -*- coding: utf-8 -*-
import pandas as pd

# Read the histology data with proper encoding, skipping the first two rows
histology_file = 'Guardant Project - Histology Data(IDC + ILC ONLY).csv'
hist_df = pd.read_csv(histology_file, encoding='latin-1', skiprows=2)

# Focus on the GH_ID column
gh_id_col = 'GH_ID'  # The actual column name after skipping header rows

# Count unique GH_IDs
unique_ids = hist_df[gh_id_col].nunique()
total_rows = len(hist_df)

print("\nSummary of Histology Data:")
print("------------------------")
print(f"Total rows in file (excluding headers): {total_rows}")
print(f"Number of unique GH_IDs: {unique_ids}")
print(f"\nSample of GH_IDs:")
print("---------------")
print(hist_df[gh_id_col].head(10)) 