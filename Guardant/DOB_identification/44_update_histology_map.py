import pandas as pd

# Read the input files
histology_map = pd.read_csv('patient_histology_map.csv')
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Create dictionaries mapping Effective Patient ID to ICD-O-3 description
idc_desc_map = idc_data[['Effective Patient ID', 'Histo/Behavior ICD-O-3-Desc']].drop_duplicates()
ilc_desc_map = ilc_data[['Effective Patient ID', 'Histo/Behavior ICD-O-3-Desc']].drop_duplicates()

# Combine the description maps
desc_map = pd.concat([idc_desc_map, ilc_desc_map])
desc_map = desc_map.drop_duplicates('Effective Patient ID')

# Merge the descriptions with the histology map
updated_map = histology_map.merge(
    desc_map,
    on='Effective Patient ID',
    how='left'
)

# Save the updated map
updated_map.to_csv('updated_histology_map.csv', index=False)

# Print some statistics
total_patients = len(histology_map)
matched_patients = updated_map['Histo/Behavior ICD-O-3-Desc'].notna().sum()
print("\nStatistics:")
print("Total patients in histology map: {}".format(total_patients))
print("Patients matched with ICD-O-3 descriptions: {}".format(matched_patients))
print("Match rate: {:.1f}%".format(matched_patients/float(total_patients)*100)) 