import pandas as pd

# Read the histology cohort files
idc_cohort = pd.read_csv('idc_cohort.csv', encoding='latin-1')
ilc_cohort = pd.read_csv('ilc_cohort.csv', encoding='latin-1')

# Read the genomic data files
idc_genomic = pd.read_csv('IDC_genomic.csv', encoding='latin-1')
ilc_genomic = pd.read_csv('ILC_genomic.csv', encoding='latin-1')

# Extract relevant columns from cohort files
idc_histology = idc_cohort[['GH_ID', 'Histo/Behavior ICD-O-3-Desc']]
ilc_histology = ilc_cohort[['GH_ID', 'Histo/Behavior ICD-O-3-Desc']]

# Combine histology data and remove any duplicates if they exist
histology_combined = pd.concat([idc_histology, ilc_histology])
histology_combined = histology_combined.rename(columns={'Histo/Behavior ICD-O-3-Desc': 'Histology Type'})
histology_combined = histology_combined.drop_duplicates(subset=['GH_ID'])

# Process IDC genomic data
idc_columns = idc_genomic.columns.tolist()
cancer_type_idx = idc_columns.index('Cancer Type')

# Count duplicates in genomic data before merge
idc_dup_counts = idc_genomic['GH_ID'].value_counts()
ilc_dup_counts = ilc_genomic['GH_ID'].value_counts()

# Merge IDC data
idc_merged = pd.merge(
    idc_genomic,
    histology_combined,
    on='GH_ID',
    how='left'
)

# Merge ILC data
ilc_merged = pd.merge(
    ilc_genomic,
    histology_combined,
    on='GH_ID',
    how='left'
)

# Reorder columns
columns_before = idc_columns[:cancer_type_idx + 1]
columns_after = idc_columns[cancer_type_idx + 1:]
new_column_order = columns_before + ['Histology Type'] + columns_after

idc_merged = idc_merged[new_column_order]
ilc_merged = ilc_merged[new_column_order]

# Save the merged data
idc_output = 'IDC_genomic_with_histology.csv'
ilc_output = 'ILC_genomic_with_histology.csv'

idc_merged.to_csv(idc_output, index=False, encoding='latin-1')
ilc_merged.to_csv(ilc_output, index=False, encoding='latin-1')

# Print summary statistics
print(f"IDC genomic data:")
print(f"Total records: {len(idc_genomic)}")
print(f"Unique GH_IDs: {idc_genomic['GH_ID'].nunique()}")
print(f"Records with histology: {idc_merged['Histology Type'].notna().sum()}")
print("\nIDC GH_IDs with multiple records:")
print(idc_dup_counts[idc_dup_counts > 1].to_string())

print(f"\nILC genomic data:")
print(f"Total records: {len(ilc_genomic)}")
print(f"Unique GH_IDs: {ilc_genomic['GH_ID'].nunique()}")
print(f"Records with histology: {ilc_merged['Histology Type'].notna().sum()}")
print("\nILC GH_IDs with multiple records:")
print(ilc_dup_counts[ilc_dup_counts > 1].to_string()) 