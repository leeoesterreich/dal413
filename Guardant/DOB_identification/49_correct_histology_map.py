import pandas as pd

# Read the files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])

# Create initial mapping dataframe from genomic data
print("\nProcessing genomic data...")

# Get relevant columns from both genomic files
idc_map = idc_data[['Effective Patient ID', 'GH_ID', 'Patient DOB', 'Histo/Behavior ICD-O-3-Desc']].copy()
idc_map['Histological Subtype'] = 'IDC'
ilc_map = ilc_data[['Effective Patient ID', 'GH_ID', 'Patient DOB', 'Histo/Behavior ICD-O-3-Desc']].copy()
ilc_map['Histological Subtype'] = 'ILC'

# For overlapping patients, use the ILC records
combined_map = pd.concat([
    idc_map[~idc_map['Effective Patient ID'].isin(overlapping_ids)],  # Non-overlapping IDC patients
    ilc_map  # All ILC patients (including overlapping ones)
])

# Sort by Effective Patient ID
combined_map = combined_map.sort_values('Effective Patient ID')

# Save the corrected histology map
print("\nSaving corrected histology map...")
output_file = 'corrected_histology_map.csv'
combined_map.to_csv(output_file, index=False)

# Print statistics
print("\nStatistics:")
print("-" * 50)
print(f"Total unique patients: {len(combined_map)}")
print("\nHistological Subtype distribution:")
print(combined_map['Histological Subtype'].value_counts())

# Verify the corrections
print("\nVerifying corrections for overlapping patients:")
print("-" * 50)
corrected_patients = combined_map[combined_map['Effective Patient ID'].isin(overlapping_ids)]
print("\nCorrected patient records:")
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
print(corrected_patients)

print(f"\nThe corrected histology map has been saved to: {output_file}")
print("Changes made:")
print(f"- Reclassified {len(overlapping_ids)} patients from IDC to ILC")
print(f"- Final count: {len(combined_map)} total patients") 