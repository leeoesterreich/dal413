import pandas as pd

# Read all files
print("Reading files...")
master_df = pd.read_csv('final_genomic_cohorts/complete_patient_mapping.csv')
idc_df = pd.read_csv('final_genomic_cohorts/IDC_genomic.csv')
ilc_df = pd.read_csv('final_genomic_cohorts/ILC_genomic.csv')

print("\nMaster File (complete_patient_mapping.csv):")
print("-" * 50)
master_counts = master_df['Histology'].value_counts()
print(f"Total patients: {len(master_df)}")
print(f"IDC patients: {master_counts.get('IDC', 0)}")
print(f"ILC patients: {master_counts.get('ILC', 0)}")

print("\nIDC Genomic File:")
print("-" * 50)
idc_unique = idc_df['Effective Patient ID'].nunique()
idc_records = len(idc_df)
print(f"Total records: {idc_records}")
print(f"Unique patients: {idc_unique}")
if 'Histology' in idc_df.columns:
    idc_hist_counts = idc_df['Histology'].value_counts()
    print("\nHistology breakdown:")
    print(idc_hist_counts)

print("\nILC Genomic File:")
print("-" * 50)
ilc_unique = ilc_df['Effective Patient ID'].nunique()
ilc_records = len(ilc_df)
print(f"Total records: {ilc_records}")
print(f"Unique patients: {ilc_unique}")
if 'Histology' in ilc_df.columns:
    ilc_hist_counts = ilc_df['Histology'].value_counts()
    print("\nHistology breakdown:")
    print(ilc_hist_counts)

# Check for any overlaps between IDC and ILC files
idc_patients = set(idc_df['Effective Patient ID'].unique())
ilc_patients = set(ilc_df['Effective Patient ID'].unique())
overlap = idc_patients & ilc_patients

print("\nOverlap Analysis:")
print("-" * 50)
print(f"Patients appearing in both IDC and ILC files: {len(overlap)}")
if len(overlap) > 0:
    print("\nOverlapping Patient IDs:")
    print(sorted(overlap)) 