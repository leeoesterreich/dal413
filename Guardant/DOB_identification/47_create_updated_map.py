import pandas as pd

# Read the files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Create a mapping dataframe from genomic data
print("\nProcessing genomic data...")

# Get relevant columns from both genomic files
idc_map = idc_data[['Effective Patient ID', 'GH_ID', 'Patient DOB', 'Histo/Behavior ICD-O-3-Desc']].copy()
idc_map['Histological Subtype'] = 'IDC'
ilc_map = ilc_data[['Effective Patient ID', 'GH_ID', 'Patient DOB', 'Histo/Behavior ICD-O-3-Desc']].copy()
ilc_map['Histological Subtype'] = 'ILC'

# Combine the data
combined_map = pd.concat([idc_map, ilc_map])

# Remove duplicates, keeping first occurrence (which maintains the original histological subtype)
combined_map = combined_map.drop_duplicates(subset=['Effective Patient ID'], keep='first')

# Sort by Effective Patient ID
combined_map = combined_map.sort_values('Effective Patient ID')

# Save the updated histology map
print("\nSaving updated histology map...")
combined_map.to_csv('new_histology_map.csv', index=False)

# Print statistics
print("\nStatistics:")
print("-" * 50)
print(f"Total unique patients: {len(combined_map)}")
print("\nHistological Subtype distribution:")
print(combined_map['Histological Subtype'].value_counts())
print("\nSample of the new mapping (first 5 rows):")
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
print(combined_map.head())

# Check for any patients with multiple records
duplicates = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])
if duplicates:
    print("\nWarning: Found patients present in both IDC and ILC files:")
    print(f"Number of overlapping patients: {len(duplicates)}")
    print("\nSample of overlapping patients:")
    overlap_sample = combined_map[combined_map['Effective Patient ID'].isin(duplicates)].head()
    print(overlap_sample[['Effective Patient ID', 'Histological Subtype', 'Histo/Behavior ICD-O-3-Desc']]) 