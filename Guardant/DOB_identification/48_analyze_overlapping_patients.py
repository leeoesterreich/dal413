import pandas as pd
import numpy as np

# Read the files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])
print(f"\nFound {len(overlapping_ids)} overlapping patients")

# Get records for overlapping patients
idc_overlaps = idc_data[idc_data['Effective Patient ID'].isin(overlapping_ids)]
ilc_overlaps = ilc_data[ilc_data['Effective Patient ID'].isin(overlapping_ids)]

# Create a detailed comparison
print("\nDetailed Comparison:")
print("-" * 100)

for patient_id in sorted(overlapping_ids):
    print(f"\nPatient: {patient_id}")
    print("-" * 50)
    
    # Get patient records
    idc_record = idc_data[idc_data['Effective Patient ID'] == patient_id]
    ilc_record = ilc_data[ilc_data['Effective Patient ID'] == patient_id]
    
    # Compare key fields
    print("GH_ID:")
    print(f"  IDC file: {idc_record['GH_ID'].iloc[0]}")
    print(f"  ILC file: {ilc_record['GH_ID'].iloc[0]}")
    
    print("\nDOB:")
    print(f"  IDC file: {idc_record['Patient DOB'].iloc[0]}")
    print(f"  ILC file: {ilc_record['Patient DOB'].iloc[0]}")
    
    print("\nICD-O-3 Description:")
    print(f"  IDC file: {idc_record['Histo/Behavior ICD-O-3-Desc'].iloc[0]}")
    print(f"  ILC file: {ilc_record['Histo/Behavior ICD-O-3-Desc'].iloc[0]}")
    
    # Check if multiple records exist for this patient in either file
    idc_count = len(idc_record)
    ilc_count = len(ilc_record)
    if idc_count > 1 or ilc_count > 1:
        print(f"\nNumber of records:")
        print(f"  IDC file: {idc_count} records")
        print(f"  ILC file: {ilc_count} records")

# Create a summary dataframe
print("\nSummary Table:")
print("-" * 100)

summary_data = []
for patient_id in sorted(overlapping_ids):
    idc_record = idc_data[idc_data['Effective Patient ID'] == patient_id].iloc[0]
    ilc_record = ilc_data[ilc_data['Effective Patient ID'] == patient_id].iloc[0]
    
    summary_data.append({
        'Effective Patient ID': patient_id,
        'GH_ID (IDC)': idc_record['GH_ID'],
        'GH_ID (ILC)': ilc_record['GH_ID'],
        'DOB (IDC)': idc_record['Patient DOB'],
        'DOB (ILC)': ilc_record['Patient DOB'],
        'ICD-O-3 Desc (IDC)': idc_record['Histo/Behavior ICD-O-3-Desc'],
        'ICD-O-3 Desc (ILC)': ilc_record['Histo/Behavior ICD-O-3-Desc']
    })

summary_df = pd.DataFrame(summary_data)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
print(summary_df)

# Save the detailed analysis
print("\nSaving detailed analysis to 'overlapping_patients_analysis.csv'...")
summary_df.to_csv('overlapping_patients_analysis.csv', index=False) 