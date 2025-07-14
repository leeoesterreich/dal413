import pandas as pd
import shutil
import os

# Read the master file
print("Reading master file...")
master_df = pd.read_csv('complete_patient_mapping.csv')

# Create histology info dictionary
histology_info = master_df.set_index('GH_ID')[['Histology', 'Histo/Behavior ICD-O-3-Desc']].to_dict('index')

# Process IDC file
print("\nProcessing IDC genomic file...")
idc_file = 'final_genomic_cohorts/IDC_genomic.csv'
idc_df = pd.read_csv(idc_file)

# Add histology columns to IDC
idc_df['Histology'] = idc_df['GH_ID'].map(lambda x: histology_info[x]['Histology'] if x in histology_info else None)
idc_df['Histo/Behavior ICD-O-3-Desc'] = idc_df['GH_ID'].map(lambda x: histology_info[x]['Histo/Behavior ICD-O-3-Desc'] if x in histology_info else None)

# Reorder columns for IDC
cols = idc_df.columns.tolist()
effective_pid_idx = cols.index('Effective Patient ID')
# Create new column order
new_cols = (
    cols[:effective_pid_idx + 1] +  # Keep columns up to Effective Patient ID
    ['Histology', 'Histo/Behavior ICD-O-3-Desc', 'GH_ID'] +  # Add histology columns and GH_ID
    [col for col in cols if col not in ['Effective Patient ID', 'Histology', 'Histo/Behavior ICD-O-3-Desc', 'GH_ID']]  # Add remaining columns
)
idc_df = idc_df[new_cols]

# Save updated IDC file
print("Saving updated IDC file...")
idc_df.to_csv(idc_file, index=False)

# Process ILC file
print("\nProcessing ILC genomic file...")
ilc_file = 'final_genomic_cohorts/ILC_genomic.csv'
ilc_df = pd.read_csv(ilc_file)

# Add histology columns to ILC
ilc_df['Histology'] = ilc_df['GH_ID'].map(lambda x: histology_info[x]['Histology'] if x in histology_info else None)
ilc_df['Histo/Behavior ICD-O-3-Desc'] = ilc_df['GH_ID'].map(lambda x: histology_info[x]['Histo/Behavior ICD-O-3-Desc'] if x in histology_info else None)

# Reorder columns for ILC
cols = ilc_df.columns.tolist()
effective_pid_idx = cols.index('Effective Patient ID')
# Create new column order
new_cols = (
    cols[:effective_pid_idx + 1] +  # Keep columns up to Effective Patient ID
    ['Histology', 'Histo/Behavior ICD-O-3-Desc', 'GH_ID'] +  # Add histology columns and GH_ID
    [col for col in cols if col not in ['Effective Patient ID', 'Histology', 'Histo/Behavior ICD-O-3-Desc', 'GH_ID']]  # Add remaining columns
)
ilc_df = ilc_df[new_cols]

# Save updated ILC file
print("Saving updated ILC file...")
ilc_df.to_csv(ilc_file, index=False)

# Copy master file to final_genomic_cohorts
print("\nCopying master file to final_genomic_cohorts...")
shutil.copy('complete_patient_mapping.csv', 'final_genomic_cohorts/complete_patient_mapping.csv')

print("Done!")

# Print statistics
print("\nStatistics:")
print(f"IDC file: {len(idc_df)} records, {idc_df['GH_ID'].nunique()} unique patients")
print(f"ILC file: {len(ilc_df)} records, {ilc_df['GH_ID'].nunique()} unique patients")
print(f"Records with missing histology in IDC: {idc_df['Histology'].isna().sum()}")
print(f"Records with missing histology in ILC: {ilc_df['Histology'].isna().sum()}")

# Print first few column names to verify order
print("\nColumn order in IDC file:")
print(", ".join(idc_df.columns[:5]))
print("\nColumn order in ILC file:")
print(", ".join(ilc_df.columns[:5])) 