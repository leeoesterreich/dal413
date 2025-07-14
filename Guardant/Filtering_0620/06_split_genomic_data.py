import pandas as pd

# Read the cohort files
idc_df = pd.read_csv('idc_cohort.csv')
ilc_df = pd.read_csv('ilc_cohort.csv')

# Get the GH_IDs for each cohort and create histology mapping
gh_id_col = 'GH_ID'
histology_map = {}
for _, row in idc_df.iterrows():
    histology_map[row[gh_id_col]] = 'IDC'
for _, row in ilc_df.iterrows():
    histology_map[row[gh_id_col]] = 'ILC'

# Read the genomic data
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# The GH_ID column in genomic data has a special character prefix
genomic_gh_id_col = 'ï»¿GH_ID'

# Rename the problematic column
genomic_df = genomic_df.rename(columns={genomic_gh_id_col: 'GH_ID'})

# Add histology information
genomic_df['Histology'] = genomic_df['GH_ID'].map(lambda x: histology_map.get(x, 'Unknown'))

# Reorder columns to put Histology after Cancer Type
cols = list(genomic_df.columns)
cancer_type_idx = cols.index('Cancer Type')
cols.remove('Histology')
cols.insert(cancer_type_idx + 1, 'Histology')
genomic_df = genomic_df[cols]

# Filter genomic data for each cohort
idc_genomic = genomic_df[genomic_df['Histology'] == 'IDC']
ilc_genomic = genomic_df[genomic_df['Histology'] == 'ILC']

# Calculate excluded samples
all_genomic_ids = set(genomic_df['GH_ID'])
included_ids = set(histology_map.keys())
excluded_ids = all_genomic_ids - included_ids

# Save the filtered genomic data
idc_output = 'IDC_genomic.csv'
ilc_output = 'ILC_genomic.csv'
idc_genomic.to_csv(idc_output, index=False)
ilc_genomic.to_csv(ilc_output, index=False)

# Print summary
print("\nGenomic Data Summary:")
print(f"Total records in genomic data: {len(genomic_df)}")
print(f"Total unique GH_IDs in genomic data: {len(all_genomic_ids)}")
print(f"\nIDC genomic records: {len(idc_genomic)}")
print(f"Unique GH_IDs in IDC genomic data: {len(set(idc_genomic['GH_ID']))}")
print(f"\nILC genomic records: {len(ilc_genomic)}")
print(f"Unique GH_IDs in ILC genomic data: {len(set(ilc_genomic['GH_ID']))}")
print(f"\nNumber of GH_IDs excluded: {len(excluded_ids)}")

if len(excluded_ids) > 0:
    print("\nFirst 10 excluded GH_IDs (if any):")
    print(sorted(list(excluded_ids))[:10])

print(f"\nSaved IDC genomic data to {idc_output}")
print(f"Saved ILC genomic data to {ilc_output}")

# Verify that all records are properly classified
total_classified = len(idc_genomic) + len(ilc_genomic)
print(f"\nVerification:")
print(f"Total classified records: {total_classified}")
print(f"Records not in either cohort: {len(genomic_df) - total_classified}")

# Print sample of data distribution
print("\nSample of data distribution:")
print("\nIDC - First 3 unique GH_IDs and their alteration counts:")
for gh_id in list(set(idc_genomic['GH_ID']))[:3]:
    count = len(idc_genomic[idc_genomic['GH_ID'] == gh_id])
    print(f"GH_ID: {gh_id}, Number of alterations: {count}")

print("\nILC - First 3 unique GH_IDs and their alteration counts:")
for gh_id in list(set(ilc_genomic['GH_ID']))[:3]:
    count = len(ilc_genomic[ilc_genomic['GH_ID'] == gh_id])
    print(f"GH_ID: {gh_id}, Number of alterations: {count}") 