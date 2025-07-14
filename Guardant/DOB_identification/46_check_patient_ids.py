import pandas as pd

# Read all files
print("Reading files...")
histology_map = pd.read_csv('patient_histology_map.csv')
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Analyze patient IDs in each file
print("\nUnique Patient Counts:")
print("-" * 50)
print(f"Histology map: {len(histology_map['Effective Patient ID'].unique())} patients")
print(f"IDC genomic: {len(idc_data['Effective Patient ID'].unique())} patients")
print(f"ILC genomic: {len(ilc_data['Effective Patient ID'].unique())} patients")

# Check format of IDs in each file
print("\nSample Effective Patient IDs from each file:")
print("-" * 50)
print("Histology map:")
print(histology_map['Effective Patient ID'].head())
print("\nIDC genomic:")
print(idc_data['Effective Patient ID'].head())
print("\nILC genomic:")
print(ilc_data['Effective Patient ID'].head())

# Check DOB consistency
print("\nChecking DOB consistency:")
print("-" * 50)

# Function to get sample of DOB matches/mismatches
def check_dob_consistency(df1, df2, name1, name2):
    merged = df1.merge(df2, on='Effective Patient ID', suffixes=('_1', '_2'))
    merged['DOB_match'] = merged['Patient DOB_1'] == merged['Patient DOB_2']
    matches = merged[merged['DOB_match']]
    mismatches = merged[~merged['DOB_match']]
    
    print(f"\nComparison between {name1} and {name2}:")
    print(f"Total matches: {len(matches)}")
    print(f"Total mismatches: {len(mismatches)}")
    
    if len(mismatches) > 0:
        print("\nSample of mismatches:")
        sample = mismatches[['Effective Patient ID', 'Patient DOB_1', 'Patient DOB_2']].head()
        print(sample)

# Check DOB consistency between files
check_dob_consistency(histology_map, idc_data, "Histology map", "IDC genomic")
check_dob_consistency(histology_map, ilc_data, "Histology map", "ILC genomic") 