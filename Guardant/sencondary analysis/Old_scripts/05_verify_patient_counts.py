import pandas as pd

# Read the cohort files
idc_df = pd.read_csv('idc_cohort.csv')
ilc_df = pd.read_csv('ilc_cohort.csv')

# Read the genomic data
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# Get column names
gh_id_col = 'Can delete duplicates of GH_ID'
genomic_gh_id_col = 'ï»¿GH_ID'

# Get unique GH_IDs from each dataset
idc_ids = set(idc_df[gh_id_col])
ilc_ids = set(ilc_df[gh_id_col])
genomic_ids = set(genomic_df[genomic_gh_id_col])

# Calculate matches
idc_with_genomic = idc_ids.intersection(genomic_ids)
ilc_with_genomic = ilc_ids.intersection(genomic_ids)

print("\nPatient Count Summary:")
print("----------------------")
print(f"Total IDC patients: {len(idc_ids)}")
print(f"Total ILC patients: {len(ilc_ids)}")
print(f"Total patients in histology cohorts: {len(idc_ids) + len(ilc_ids)}")

print("\nGenomic Data Matching:")
print("----------------------")
print(f"Total unique patients in genomic data: {len(genomic_ids)}")
print(f"IDC patients with genomic data: {len(idc_with_genomic)} ({len(idc_with_genomic)/len(idc_ids)*100:.1f}%)")
print(f"ILC patients with genomic data: {len(ilc_with_genomic)} ({len(ilc_with_genomic)/len(ilc_ids)*100:.1f}%)")
print(f"Total patients with both histology and genomic data: {len(idc_with_genomic) + len(ilc_with_genomic)}")

# Calculate patients without genomic data
idc_without_genomic = idc_ids - genomic_ids
ilc_without_genomic = ilc_ids - genomic_ids

print("\nPatients Without Genomic Data:")
print("-----------------------------")
print(f"IDC patients without genomic data: {len(idc_without_genomic)}")
print(f"ILC patients without genomic data: {len(ilc_without_genomic)}")

# Calculate patients in genomic data but not in histology cohorts
genomic_only = genomic_ids - (idc_ids.union(ilc_ids))
print(f"\nPatients in genomic data but not in histology cohorts: {len(genomic_only)}") 