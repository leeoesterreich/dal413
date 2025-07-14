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

# Check for overlapping between IDC and ILC
overlap_idc_ilc = idc_ids.intersection(ilc_ids)

print("\nChecking for overlapping between cohorts:")
print("----------------------------------------")
print(f"Number of patients in IDC cohort: {len(idc_ids)}")
print(f"Number of patients in ILC cohort: {len(ilc_ids)}")
if overlap_idc_ilc:
    print(f"\nWARNING: Found {len(overlap_idc_ilc)} patients that appear in both IDC and ILC cohorts!")
    print("\nOverlapping GH_IDs:")
    for gh_id in sorted(overlap_idc_ilc):
        print(f"GH_ID: {gh_id}")
else:
    print("\nNo overlapping found between IDC and ILC cohorts - this is good!")

# Calculate various intersections with genomic data
idc_with_genomic = idc_ids.intersection(genomic_ids)
ilc_with_genomic = ilc_ids.intersection(genomic_ids)

# Check if any patient appears in both IDC and ILC genomic sets
overlap_genomic = idc_with_genomic.intersection(ilc_with_genomic)

print("\nGenomic Data Matching:")
print("----------------------")
print(f"Total unique patients in genomic data: {len(genomic_ids)}")
print(f"IDC patients with genomic data: {len(idc_with_genomic)}")
print(f"ILC patients with genomic data: {len(ilc_with_genomic)}")
if overlap_genomic:
    print(f"\nWARNING: Found {len(overlap_genomic)} patients that appear in both IDC and ILC genomic sets!")
    print("\nOverlapping GH_IDs in genomic data:")
    for gh_id in sorted(overlap_genomic):
        print(f"GH_ID: {gh_id}")
else:
    print("\nNo overlapping found in genomic data matches - this is good!")

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

# Final summary of unique patients
all_patients = idc_ids.union(ilc_ids).union(genomic_ids)
print(f"\nTotal number of unique patients across all datasets: {len(all_patients)}")

# Verify our counts add up
total_in_cohorts = len(idc_ids) + len(ilc_ids) - len(overlap_idc_ilc)
print("\nVerification of counts:")
print(f"Sum of IDC and ILC cohorts (accounting for overlap): {total_in_cohorts}")
print(f"Patients with genomic data only: {len(genomic_only)}")
print(f"Total unique patients (should match above): {total_in_cohorts + len(genomic_only)}") 